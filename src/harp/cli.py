from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import MDAnalysis as mda

from .analysis import analyze_profile, run_from_config, run_panel
from .deepcoil import parse_deepcoil2
from .validation import validate_config, format_validation_report

AA3_TO_1 = {
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G",
    "HIS":"H","HSD":"H","HSE":"H","HSP":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M",
    "PHE":"F","PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V",
}


def _sequence(residues) -> str:
    return "".join(AA3_TO_1.get(str(r.resname).upper(), "X") for r in residues)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="harp",
        description="DeepCoil-anchored heptad analysis of IBAM–MyhT MD contact patterns",
    )
    sub = p.add_subparsers(dest="command", required=True)

    reg = sub.add_parser("register", help="Convert DeepCoil2 output into an explicit a–g register")
    reg.add_argument("--deepcoil", required=True)
    reg.add_argument("--out", required=True)
    reg.add_argument("--cc-threshold", type=float, default=0.5)
    reg.add_argument("--anchor-threshold", type=float, default=0.5)

    val = sub.add_parser("validate", help="Strictly validate HARP inputs without modifying them")
    val.add_argument("--config", required=True)
    val.add_argument("--report", help="Optional path for validation_report.json")

    ana = sub.add_parser("analyze", help="Run full trajectory analysis from YAML config")
    ana.add_argument("--config", required=True)

    prof = sub.add_parser("profile", help="Analyze a precomputed per-residue contact profile")
    prof.add_argument("--deepcoil", required=True)
    prof.add_argument("--profile", required=True, help="TSV/CSV with contact_occupancy column")
    prof.add_argument("--outdir", required=True)
    prof.add_argument("--n-permutations", type=int, default=9999)
    prof.add_argument("--block-size", type=int, default=4)
    prof.add_argument("--seed", type=int, default=20260801)

    inspect = sub.add_parser("inspect", help="Inspect topology segments or validate an MDAnalysis selection")
    inspect.add_argument("--topology", required=True)
    inspect.add_argument("--trajectory")
    inspect.add_argument("--selection")

    panel = sub.add_parser("panel", help="Test shared phase alignment across completed taxa")
    panel.add_argument("--manifest", required=True, help="CSV: taxon,summary_json")
    panel.add_argument("--outdir", default="results", help="Results root; HARP creates a timestamped n-taxa run directory")
    panel.add_argument("--n-permutations", type=int, default=9999)
    panel.add_argument("--seed", type=int, default=20260801)
    panel.add_argument("--direct-outdir", action="store_true")
    panel.add_argument("--run-name")
    return p


def main(argv=None) -> None:
    args = build_parser().parse_args(argv)
    if args.command == "validate":
        report = validate_config(args.config)
        print(format_validation_report(report))
        if args.report:
            report.write_json(args.report)
        if not report.passed:
            raise SystemExit(2)
    elif args.command == "register":
        reg = parse_deepcoil2(
            args.deepcoil,
            cc_threshold=args.cc_threshold,
            anchor_threshold=args.anchor_threshold,
        )
        reg.write_tsv(args.out)
        print(json.dumps({
            "sequence_length": reg.length,
            "cc_segment": [reg.cc_start, reg.cc_end],
            "a_origin": reg.a_origin,
            "d_offset": reg.d_offset,
            "register_quality": {
                "a_purity": reg.a_purity,
                "d_purity": reg.d_purity,
                "combined_purity": reg.register_purity,
                "n_a_calls": reg.n_a_calls,
                "n_d_calls": reg.n_d_calls,
                "a_modal_tie": reg.a_tie,
                "d_modal_tie": reg.d_tie,
            },
            "output": str(Path(args.out).resolve()),
        }, indent=2))
    elif args.command == "analyze":
        print(json.dumps(run_from_config(args.config), indent=2))
    elif args.command == "profile":
        sep = "\t" if Path(args.profile).suffix.lower() in {".tsv", ".txt"} else ","
        df = pd.read_csv(args.profile, sep=sep)
        if "contact_occupancy" not in df.columns:
            raise SystemExit("Profile file must have a contact_occupancy column")
        summary = analyze_profile(
            deepcoil_path=args.deepcoil,
            profile=df["contact_occupancy"].to_numpy(float),
            outdir=args.outdir,
            n_permutations=args.n_permutations,
            block_size=args.block_size,
            seed=args.seed,
        )
        print(json.dumps(summary, indent=2))
    elif args.command == "inspect":
        u = mda.Universe(args.topology, args.trajectory) if args.trajectory else mda.Universe(args.topology)
        if args.selection:
            ag = u.select_atoms(args.selection)
            print(json.dumps({
                "selection": args.selection,
                "n_atoms": len(ag),
                "n_residues": len(ag.residues),
                "resids": [int(r.resid) for r in ag.residues],
                "resnames": [str(r.resname) for r in ag.residues],
                "sequence": _sequence(ag.residues),
                "segids": sorted(set(str(x) for x in ag.segids)),
            }, indent=2))
        else:
            rows = []
            for seg in u.segments:
                resids = [int(r.resid) for r in seg.residues]
                rows.append({
                    "segid": str(seg.segid),
                    "n_atoms": len(seg.atoms),
                    "n_residues": len(seg.residues),
                    "resid_min": min(resids) if resids else None,
                    "resid_max": max(resids) if resids else None,
                    "sequence": _sequence(seg.residues),
                })
            print(json.dumps(rows, indent=2))
    elif args.command == "panel":
        print(json.dumps(run_panel(
            args.manifest,
            args.outdir,
            args.n_permutations,
            args.seed,
            direct_outdir=args.direct_outdir,
            run_name_override=args.run_name,
        ), indent=2))


if __name__ == "__main__":
    main()
