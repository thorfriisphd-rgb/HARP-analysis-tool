from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re

import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import capped_distance


@dataclass
class ContactTrajectoryResult:
    frame_times_ps: np.ndarray
    myht_frame_signal: np.ndarray
    pair_occupancy: np.ndarray
    myht_resids: np.ndarray
    myht_resnames: np.ndarray
    mg_resids: np.ndarray
    mg_resnames: np.ndarray


def _heavy_atom_group(universe, selection: str):
    ag = universe.select_atoms(selection)
    if len(ag) == 0:
        return ag
    names = np.array([str(name).upper() for name in ag.names], dtype=object)
    # Covers H, HA, HB2 and digit-prefixed GROMACS names such as 1HG1.
    keep = np.array([re.match(r"^\d*H", name) is None for name in names], dtype=bool)
    try:
        masses = np.asarray(ag.masses, dtype=float)
        if np.all(np.isfinite(masses)) and np.any(masses > 0):
            keep &= masses > 1.5
    except (AttributeError, mda.exceptions.NoDataError):
        pass
    return ag[keep]


def _local_residue_map(atom_group) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    residues = atom_group.residues
    global_to_local = {int(r.ix): i for i, r in enumerate(residues)}
    atom_to_local = np.array([global_to_local[int(ix)] for ix in atom_group.resindices], dtype=int)
    resids = np.array([int(r.resid) for r in residues], dtype=int)
    resnames = np.array([str(r.resname) for r in residues], dtype=object)
    return atom_to_local, resids, resnames


def compute_contact_trajectory(
    *,
    topology: str | Path,
    trajectory: str | Path,
    mg_selection: str,
    myht_selection: str,
    cutoff_angstrom: float = 4.5,
    contact_mode: str = "hard",
    smooth_power: int = 6,
    start: int | None = None,
    stop: int | None = None,
    step: int = 1,
) -> ContactTrajectoryResult:
    """Extract MG↔MyhT residue contacts from an MD trajectory.

    hard: a residue-pair score is 1 when any heavy-atom pair is within cutoff.
    smooth: score = 1/(1+(d_min/cutoff)^power) for pairs inside 2*cutoff.
    """
    if contact_mode not in {"hard", "smooth"}:
        raise ValueError("contact_mode must be 'hard' or 'smooth'")
    if step < 1:
        raise ValueError("step must be >= 1")

    u = mda.Universe(str(topology), str(trajectory))
    mg = _heavy_atom_group(u, mg_selection)
    myht = _heavy_atom_group(u, myht_selection)
    if len(mg) == 0 or len(myht) == 0:
        raise ValueError(
            f"Empty selection: MG atoms={len(mg)}, MyhT atoms={len(myht)}. "
            "Check topology chain/segid/resid naming."
        )

    mg_atom_res, mg_resids, mg_resnames = _local_residue_map(mg)
    my_atom_res, my_resids, my_resnames = _local_residue_map(myht)
    n_mg, n_my = len(mg_resids), len(my_resids)

    frame_signals: list[np.ndarray] = []
    pair_sum = np.zeros((n_mg, n_my), dtype=np.float64)
    times: list[float] = []

    traj_slice = u.trajectory[slice(start, stop, step)]
    search_cutoff = cutoff_angstrom if contact_mode == "hard" else 2.0 * cutoff_angstrom
    for ts in traj_slice:
        atom_pairs, distances = capped_distance(
            mg.positions,
            myht.positions,
            max_cutoff=search_cutoff,
            box=ts.dimensions,
            return_distances=True,
        )
        pair_weights = np.zeros(n_mg * n_my, dtype=np.float64)
        if len(atom_pairs):
            mg_local = mg_atom_res[atom_pairs[:, 0]]
            my_local = my_atom_res[atom_pairs[:, 1]]
            codes = mg_local * n_my + my_local
            if contact_mode == "hard":
                pair_weights[np.unique(codes)] = 1.0
            else:
                weights = 1.0 / (1.0 + (distances / cutoff_angstrom) ** smooth_power)
                np.maximum.at(pair_weights, codes, weights)

        pair_matrix = pair_weights.reshape(n_mg, n_my)
        pair_sum += pair_matrix
        # A MyhT residue signal is whether/degree to which any MG residue engages it.
        frame_signals.append(pair_matrix.max(axis=0).astype(np.float32))
        times.append(float(ts.time))

    if not frame_signals:
        raise ValueError("Trajectory slice contains no frames")
    frame_signal = np.vstack(frame_signals)
    pair_occupancy = pair_sum / frame_signal.shape[0]
    return ContactTrajectoryResult(
        frame_times_ps=np.asarray(times, dtype=float),
        myht_frame_signal=frame_signal,
        pair_occupancy=pair_occupancy,
        myht_resids=my_resids,
        myht_resnames=my_resnames,
        mg_resids=mg_resids,
        mg_resnames=mg_resnames,
    )
