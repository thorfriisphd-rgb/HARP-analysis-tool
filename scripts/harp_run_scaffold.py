#!/usr/bin/env python3
"""
Operation HARP v2 scaffold

Adds normalized Socket2 integration on top of the original HARP scaffold.

New capabilities vs v1:
- load normalized Socket2 parser outputs
- compute conservative key a/d overlap against Socket2
- compute conservative KIH-support proxy against Socket2
- annotate residue-level Socket2 heptads in outputs
- preserve the rule that STRONG_SUPPORT is impossible without external confirmation

This remains intentionally cautious. It does not attempt exact structural re-derivation
of KIH geometry from coordinates yet. Instead it compares HARP assignments against the
normalized external evidence emitted by parse_socket2.py.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any, Set

try:
    import yaml  # type: ignore
except ImportError:
    yaml = None


HEPTAD = ("a", "b", "c", "d", "e", "f", "g")
HYDROPHOBIC = {"A", "V", "I", "L", "M", "F", "W", "Y", "C"}
CHARGED = {"D", "E", "K", "R", "H"}
POLAR = {"S", "T", "N", "Q", "Y", "C", "H"}


@dataclass
class RunConfig:
    run_id: str
    structure: str
    ibam_chain: str
    myht_chain: str
    ibam_range: List[int]
    myht_range: List[int]
    prco_table: str
    socket2_output: Optional[str]
    key_residue_mode: str
    key_residue_count: int
    weights: Dict[str, float]
    thresholds: Dict[str, float]


@dataclass
class ResidueRecord:
    chain: str
    resid: int
    resname: str
    prco_score: float = 0.0
    occupancy: float = 0.0
    is_key: bool = False
    assigned_heptad: Optional[str] = None
    socket2_heptad: Optional[str] = None
    ad_concordant: Optional[bool] = None
    kih_concordant: Optional[bool] = None


@dataclass
class CandidateModel:
    candidate_id: str
    phase_offset: int
    ibam_assignments: Dict[int, str]
    myht_assignments: Dict[int, str]
    contact_alignment_score: float
    hydrophobic_core_score: float
    charge_pattern_score: float
    periodicity_score: float
    total_internal_score: float
    socket2_ad_overlap: float = 0.0
    socket2_kih_overlap: float = 0.0
    socket2_status: str = "not_run"
    final_rank: Optional[int] = None


@dataclass
class FinalVerdict:
    verdict: str
    internal_margin: float
    ad_overlap: float
    kih_overlap: float
    rationale: str


def load_config(path: Path) -> RunConfig:
    if not path.exists():
        raise FileNotFoundError(f"Config not found: {path}")

    if path.suffix.lower() in {".yaml", ".yml"}:
        if yaml is None:
            raise RuntimeError("PyYAML is required for YAML configs. Install pyyaml or use JSON.")
        raw = yaml.safe_load(path.read_text())
    elif path.suffix.lower() == ".json":
        raw = json.loads(path.read_text())
    else:
        raise ValueError("Config must be .yaml, .yml, or .json")

    required = [
        "run_id", "structure", "ibam_chain", "myht_chain", "ibam_range", "myht_range",
        "prco_table", "key_residue_mode", "key_residue_count", "weights", "thresholds"
    ]
    missing = [k for k in required if k not in raw]
    if missing:
        raise ValueError(f"Missing config field(s): {', '.join(missing)}")

    return RunConfig(
        run_id=raw["run_id"],
        structure=raw["structure"],
        ibam_chain=raw["ibam_chain"],
        myht_chain=raw["myht_chain"],
        ibam_range=list(raw["ibam_range"]),
        myht_range=list(raw["myht_range"]),
        prco_table=raw["prco_table"],
        socket2_output=raw.get("socket2_output"),
        key_residue_mode=raw["key_residue_mode"],
        key_residue_count=int(raw["key_residue_count"]),
        weights={str(k): float(v) for k, v in raw["weights"].items()},
        thresholds={str(k): float(v) for k, v in raw["thresholds"].items()},
    )


def ensure_range_ok(name: str, span: List[int]) -> Tuple[int, int]:
    if len(span) != 2:
        raise ValueError(f"{name} must contain exactly two integers: [start, end]")
    start, end = int(span[0]), int(span[1])
    if end < start:
        raise ValueError(f"{name} end < start: {span}")
    return start, end

import csv
import re
from pathlib import Path

PARTNER_RE = re.compile(r"([A-Z]{3})(\d+)$")

def _parse_float(x):
    if x is None:
        return None
    s = str(x).strip()
    if not s:
        return None
    try:
        return float(s)
    except ValueError:
        return None


def _parse_partner_token(token):
    if token is None:
        return None

    s = str(token).strip()
    if not s:
        return None

    m = PARTNER_RE.match(s)

    if not m:
        return {
            "raw": s,
            "resname": None,
            "resid": None,
            "status": "unparsed",
        }

    return {
        "raw": s,
        "resname": m.group(1),
        "resid": int(m.group(2)),
        "status": "parsed",
    }


def load_prco_table(path, ibam_chain, myht_chain):

    path = Path(path)

    if not path.exists():
        raise FileNotFoundError(f"PRCO table not found: {path}")

    rows = []

    with path.open("r", newline="") as fh:

        reader = csv.DictReader(fh)

        headers = reader.fieldnames or []

        if not {"C12_resid", "C12_resname", "occupancy"}.issubset(headers):
            raise ValueError("PRCO file does not match expected Ovis schema")

        for i, rec in enumerate(reader, start=2):

            try:
                ibam_resid = int(str(rec["C12_resid"]).strip())
            except Exception:
                ibam_resid = None

            ibam_resname = str(rec.get("C12_resname", "")).strip() or None
            occ = _parse_float(rec.get("occupancy"))

            p1 = _parse_partner_token(rec.get("top_partner_1"))
            if p1:
                p1["occ"] = _parse_float(rec.get("partner1_occ"))
                p1["chain"] = myht_chain

            p2 = _parse_partner_token(rec.get("top_partner_2"))
            if p2:
                p2["occ"] = _parse_float(rec.get("partner2_occ"))
                p2["chain"] = myht_chain

            rows.append({
                "source_line": i,
                "ibam_chain": ibam_chain,
                "ibam_resid": ibam_resid,
                "ibam_resname": ibam_resname,
                "occupancy": occ,
                "partners": [p1, p2],
            })

    if not rows:
        raise ValueError("No PRCO rows parsed")

    return rows



def filter_span(records: List[ResidueRecord], start: int, end: int) -> List[ResidueRecord]:
    return [r for r in records if start <= r.resid <= end]


def mark_key_residues(records: List[ResidueRecord], mode: str, count: int) -> None:
    if mode != "top_prco":
        raise NotImplementedError(f"key_residue_mode '{mode}' is not implemented in scaffold")
    ranked = sorted(records, key=lambda r: (r.prco_score, r.occupancy), reverse=True)
    for rec in ranked[:count]:
        rec.is_key = True


def enumerate_phase_models(
    ibam_records: List[ResidueRecord],
    myht_records: List[ResidueRecord],
) -> List[Tuple[int, Dict[int, str], Dict[int, str]]]:
    results = []
    ibam_sorted = sorted(ibam_records, key=lambda r: r.resid)
    myht_sorted = sorted(myht_records, key=lambda r: r.resid)
    ibam_start = ibam_sorted[0].resid
    myht_start = myht_sorted[0].resid

    for offset in range(7):
        ibam_assignments = {}
        myht_assignments = {}
        for rec in ibam_sorted:
            idx = (rec.resid - ibam_start) % 7
            ibam_assignments[rec.resid] = HEPTAD[idx]
        for rec in myht_sorted:
            idx = ((rec.resid - myht_start) + offset) % 7
            myht_assignments[rec.resid] = HEPTAD[idx]
        results.append((offset, ibam_assignments, myht_assignments))
    return results


def score_contact_alignment(records: List[ResidueRecord], assignments: Dict[int, str]) -> float:
    interface_like = {"a", "d", "e", "g"}
    total_weight = sum(max(r.prco_score, 0.0) for r in records) or 1.0
    score = 0.0
    for r in records:
        pos = assignments.get(r.resid)
        if pos in interface_like:
            score += max(r.prco_score, 0.0)
    return score / total_weight


def score_hydrophobic_core(records: List[ResidueRecord], assignments: Dict[int, str]) -> float:
    total = 0.0
    max_total = 0.0
    for r in records:
        pos = assignments.get(r.resid)
        if pos not in {"a", "d"}:
            continue
        max_total += 1.0
        if r.resname in HYDROPHOBIC:
            total += 1.0
        elif r.resname in CHARGED:
            total -= 0.75
        else:
            total += 0.1
    if max_total == 0:
        return 0.0
    return (total / max_total + 1.0) / 2.0


def score_charge_pattern(records: List[ResidueRecord], assignments: Dict[int, str]) -> float:
    supportive = 0.0
    counted = 0.0
    for r in records:
        pos = assignments.get(r.resid)
        if pos is None:
            continue
        counted += 1.0
        if pos in {"e", "g", "b", "c", "f"} and (r.resname in CHARGED or r.resname in POLAR):
            supportive += 1.0
        elif pos in {"a", "d"} and r.resname in CHARGED:
            supportive -= 0.5
    if counted == 0:
        return 0.0
    return max(0.0, min(1.0, (supportive / counted + 0.5)))


def score_periodicity(records: List[ResidueRecord], assignments: Dict[int, str]) -> float:
    keyed = [r for r in records if r.prco_score > 0]
    if len(keyed) < 2:
        return 0.0

    weighted_hits = 0.0
    weighted_total = 0.0
    keyed_sorted = sorted(keyed, key=lambda r: r.resid)
    for i, r1 in enumerate(keyed_sorted):
        for r2 in keyed_sorted[i + 1:]:
            delta = abs(r2.resid - r1.resid)
            if delta == 0:
                continue
            pair_weight = max(r1.prco_score, 0.0) + max(r2.prco_score, 0.0)
            weighted_total += pair_weight
            same_class = assignments.get(r1.resid) == assignments.get(r2.resid)
            near_heptad = delta % 7 == 0
            if same_class and near_heptad:
                weighted_hits += pair_weight
    if weighted_total == 0:
        return 0.0
    return weighted_hits / weighted_total


def aggregate_score(contact_alignment: float, hydrophobic_core: float, charge_pattern: float, periodicity: float,
                    weights: Dict[str, float]) -> float:
    return (
        weights.get("contact_alignment", 1.0) * contact_alignment +
        weights.get("hydrophobic_core", 1.0) * hydrophobic_core +
        weights.get("charge_pattern", 0.7) * charge_pattern +
        weights.get("periodicity", 0.8) * periodicity
    )


def build_candidates(ibam_records: List[ResidueRecord], myht_records: List[ResidueRecord],
                     weights: Dict[str, float]) -> List[CandidateModel]:
    candidates: List[CandidateModel] = []
    for offset, ibam_assignments, myht_assignments in enumerate_phase_models(ibam_records, myht_records):
        ca = 0.5 * (
            score_contact_alignment(ibam_records, ibam_assignments) +
            score_contact_alignment(myht_records, myht_assignments)
        )
        hc = 0.5 * (
            score_hydrophobic_core(ibam_records, ibam_assignments) +
            score_hydrophobic_core(myht_records, myht_assignments)
        )
        cp = 0.5 * (
            score_charge_pattern(ibam_records, ibam_assignments) +
            score_charge_pattern(myht_records, myht_assignments)
        )
        ps = 0.5 * (
            score_periodicity(ibam_records, ibam_assignments) +
            score_periodicity(myht_records, myht_assignments)
        )
        total = aggregate_score(ca, hc, cp, ps, weights)
        candidates.append(CandidateModel(
            candidate_id=f"phase_{offset}",
            phase_offset=offset,
            ibam_assignments=ibam_assignments,
            myht_assignments=myht_assignments,
            contact_alignment_score=ca,
            hydrophobic_core_score=hc,
            charge_pattern_score=cp,
            periodicity_score=ps,
            total_internal_score=total,
        ))

    candidates.sort(key=lambda c: c.total_internal_score, reverse=True)
    for i, c in enumerate(candidates, start=1):
        c.final_rank = i
    return candidates


def load_socket2_normalized(socket2_path: Optional[Path]) -> Dict[str, Any]:
    if socket2_path is None:
        return {"status": "not_provided", "ad_map": {}, "kih_pairs": set(), "notes": []}
    if not socket2_path.exists():
        return {"status": "missing", "ad_map": {}, "kih_pairs": set(), "notes": [f"Missing path: {socket2_path}"]}

    if socket2_path.is_dir():
        base = socket2_path
        summary_path = base / "socket2_summary.json"
    else:
        if socket2_path.name == "socket2_summary.json":
            summary_path = socket2_path
            base = socket2_path.parent
        else:
            return {
                "status": "unsupported_path",
                "ad_map": {},
                "kih_pairs": set(),
                "notes": ["socket2_output must be a normalized parser directory or socket2_summary.json"],
            }

    if not summary_path.exists():
        return {
            "status": "missing_summary",
            "ad_map": {},
            "kih_pairs": set(),
            "notes": [f"Missing Socket2 summary: {summary_path}"],
        }

    try:
        summary = json.loads(summary_path.read_text())
    except Exception as exc:
        return {
            "status": "invalid_summary",
            "ad_map": {},
            "kih_pairs": set(),
            "notes": [f"Failed to parse Socket2 summary JSON: {exc}"],
        }

    ad_map: Dict[Tuple[str, int], str] = {}
    kih_pairs: Set[Tuple[Tuple[str, int], Tuple[str, int]]] = set()
    notes = list(summary.get("notes", []))

    assignments_path = base / "residue_assignments.tsv"
    if assignments_path.exists():
        with assignments_path.open("r", newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                try:
                    chain = str(row["chain"]).strip()
                    resid = int(str(row["resid"]).strip())
                    hep = str(row["socket2_heptad"]).strip().lower()
                except Exception:
                    continue
                if hep in HEPTAD:
                    ad_map[(chain, resid)] = hep

    kih_path = base / "kih_pairs.tsv"
    if kih_path.exists():
        with kih_path.open("r", newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                try:
                    kc = str(row["knob_chain"]).strip()
                    kr = int(str(row["knob_resid"]).strip())
                    hc = str(row["hole_chain"]).strip()
                    hr = int(str(row["hole_resid"]).strip())
                except Exception:
                    continue
                kih_pairs.add(((kc, kr), (hc, hr)))

    return {
        "status": str(summary.get("status", "unknown")),
        "ad_map": ad_map,
        "kih_pairs": kih_pairs,
        "notes": notes,
    }


def compute_socket2_ad_overlap(candidate: CandidateModel, ibam_records: List[ResidueRecord],
                               myht_records: List[ResidueRecord], socket2_blob: Dict[str, Any]
                               ) -> Tuple[float, Dict[Tuple[str, int], bool]]:
    ad_map = socket2_blob.get("ad_map", {})
    if not ad_map:
        return 0.0, {}

    key_records = [r for r in list(ibam_records) + list(myht_records) if r.is_key]
    if not key_records:
        return 0.0, {}

    matches: Dict[Tuple[str, int], bool] = {}
    comparable = 0
    agreed = 0
    ibam_chain = ibam_records[0].chain if ibam_records else ""

    for rec in key_records:
        harp_hep = candidate.ibam_assignments.get(rec.resid) if rec.chain == ibam_chain else candidate.myht_assignments.get(rec.resid)
        socket_hep = ad_map.get((rec.chain, rec.resid))
        if socket_hep is None or harp_hep is None:
            continue
        comparable += 1
        agree = (harp_hep in {"a", "d"}) and (socket_hep in {"a", "d"})
        matches[(rec.chain, rec.resid)] = agree
        if agree:
            agreed += 1

    return (agreed / comparable, matches) if comparable else (0.0, matches)


def compute_socket2_kih_overlap(candidate: CandidateModel, ibam_records: List[ResidueRecord],
                                myht_records: List[ResidueRecord], socket2_blob: Dict[str, Any]) -> float:
    kih_pairs = socket2_blob.get("kih_pairs", set())
    if not kih_pairs:
        return 0.0

    ibam_chain = ibam_records[0].chain if ibam_records else ""
    myht_chain = myht_records[0].chain if myht_records else ""
    supportive = 0
    total = 0

    for (knob, hole) in kih_pairs:
        (kc, kr), (hc, hr) = knob, hole
        if {kc, hc} != {ibam_chain, myht_chain}:
            continue

        if kc == ibam_chain:
            knob_hep = candidate.ibam_assignments.get(kr)
            hole_hep = candidate.myht_assignments.get(hr)
        else:
            knob_hep = candidate.myht_assignments.get(kr)
            hole_hep = candidate.ibam_assignments.get(hr)

        if knob_hep is None or hole_hep is None:
            continue

        total += 1
        if knob_hep in {"a", "d", "e", "g"} and hole_hep in {"a", "d", "e", "g"}:
            supportive += 1

    return supportive / total if total else 0.0


def compare_with_socket2(candidates: List[CandidateModel], ibam_records: List[ResidueRecord],
                         myht_records: List[ResidueRecord], socket2_blob: Dict[str, Any]) -> Dict[Tuple[str, int], bool]:
    residue_matches: Dict[Tuple[str, int], bool] = {}
    status = socket2_blob.get("status", "not_run")
    for idx, c in enumerate(candidates):
        ad_overlap, matches = compute_socket2_ad_overlap(c, ibam_records, myht_records, socket2_blob)
        kih_overlap = compute_socket2_kih_overlap(c, ibam_records, myht_records, socket2_blob)
        c.socket2_status = status
        c.socket2_ad_overlap = ad_overlap
        c.socket2_kih_overlap = kih_overlap
        if idx == 0:
            residue_matches = matches
    return residue_matches


def classify_result(candidates: List[CandidateModel], thresholds: Dict[str, float]) -> FinalVerdict:
    if not candidates:
        return FinalVerdict("UNSUPPORTED", 0.0, 0.0, 0.0, "No candidate models were generated.")

    top = candidates[0]
    second = candidates[1] if len(candidates) > 1 else None
    margin = top.total_internal_score - (second.total_internal_score if second else 0.0)
    strong_margin = thresholds.get("strong_margin", 0.20)
    ad_min = thresholds.get("ad_overlap_min", 0.60)
    kih_min = thresholds.get("kih_overlap_min", 0.60)

    if top.socket2_status in {"not_provided", "missing", "missing_summary", "unsupported_path", "invalid_summary", "partial", "unsupported_format", "empty", "unknown"}:
        if margin >= strong_margin:
            return FinalVerdict(
                "AMBIGUOUS", margin, top.socket2_ad_overlap, top.socket2_kih_overlap,
                "A preferred internal HARP model was identified, but external Socket2 confirmation is absent or incomplete."
            )
        return FinalVerdict(
            "UNSUPPORTED", margin, top.socket2_ad_overlap, top.socket2_kih_overlap,
            "Internal separation is weak and Socket2 confirmation is absent or incomplete."
        )

    if margin >= strong_margin and top.socket2_ad_overlap >= ad_min and top.socket2_kih_overlap >= kih_min:
        return FinalVerdict(
            "STRONG_SUPPORT", margin, top.socket2_ad_overlap, top.socket2_kih_overlap,
            "Top HARP model is well-separated and externally concordant with Socket2."
        )

    if margin > 0 and (top.socket2_ad_overlap > 0 or top.socket2_kih_overlap > 0):
        return FinalVerdict(
            "PARTIAL_NONCANONICAL", margin, top.socket2_ad_overlap, top.socket2_kih_overlap,
            "Some internal and external support exists, but not enough for strong support."
        )

    if margin > 0:
        return FinalVerdict(
            "AMBIGUOUS", margin, top.socket2_ad_overlap, top.socket2_kih_overlap,
            "A preferred model exists internally, but external concordance is weak."
        )

    return FinalVerdict(
        "UNSUPPORTED", margin, top.socket2_ad_overlap, top.socket2_kih_overlap,
        "No persuasive support for a stable register assignment."
    )


def write_candidate_table(path: Path, candidates: List[CandidateModel]) -> None:
    with path.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow([
            "candidate_id", "phase_offset", "contact_alignment_score", "hydrophobic_core_score",
            "charge_pattern_score", "periodicity_score", "total_internal_score",
            "socket2_ad_overlap", "socket2_kih_overlap", "socket2_status", "rank"
        ])
        for c in candidates:
            writer.writerow([
                c.candidate_id, c.phase_offset,
                f"{c.contact_alignment_score:.4f}",
                f"{c.hydrophobic_core_score:.4f}",
                f"{c.charge_pattern_score:.4f}",
                f"{c.periodicity_score:.4f}",
                f"{c.total_internal_score:.4f}",
                f"{c.socket2_ad_overlap:.4f}",
                f"{c.socket2_kih_overlap:.4f}",
                c.socket2_status,
                c.final_rank,
            ])


def write_residue_assignments(path: Path, ibam_records: List[ResidueRecord], myht_records: List[ResidueRecord],
                              top_candidate: CandidateModel) -> None:
    ibam_chain = ibam_records[0].chain if ibam_records else ""
    with path.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow([
            "chain", "resid", "resname", "prco_score", "occupancy", "is_key",
            "assigned_heptad", "socket2_heptad", "ad_concordant", "kih_concordant"
        ])
        for rec in list(ibam_records) + list(myht_records):
            assigned = top_candidate.ibam_assignments.get(rec.resid) if rec.chain == ibam_chain else top_candidate.myht_assignments.get(rec.resid)
            writer.writerow([
                rec.chain, rec.resid, rec.resname, f"{rec.prco_score:.4f}", f"{rec.occupancy:.4f}",
                int(rec.is_key), assigned, rec.socket2_heptad or "", rec.ad_concordant, rec.kih_concordant
            ])


def write_socket2_table(path: Path, candidates: List[CandidateModel]) -> None:
    top = candidates[0]
    with path.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["top_candidate", "socket2_status", "ad_overlap", "kih_overlap"])
        writer.writerow([top.candidate_id, top.socket2_status, f"{top.socket2_ad_overlap:.4f}", f"{top.socket2_kih_overlap:.4f}"])


def write_qc_log(path: Path, cfg: RunConfig, ibam_span: Tuple[int, int], myht_span: Tuple[int, int],
                 ibam_records: List[ResidueRecord], myht_records: List[ResidueRecord], socket2_blob: Dict[str, Any]) -> None:
    lines = [
        f"run_id: {cfg.run_id}",
        f"structure: {cfg.structure}",
        f"ibam_chain: {cfg.ibam_chain}",
        f"myht_chain: {cfg.myht_chain}",
        f"ibam_range: {ibam_span[0]}-{ibam_span[1]}",
        f"myht_range: {myht_span[0]}-{myht_span[1]}",
        f"prco_table: {cfg.prco_table}",
        f"socket2_output: {cfg.socket2_output}",
        f"socket2_status: {socket2_blob.get('status', 'unknown')}",
        f"socket2_assignment_count: {len(socket2_blob.get('ad_map', {}))}",
        f"socket2_kih_pair_count: {len(socket2_blob.get('kih_pairs', set()))}",
        f"ibam_records_in_span: {len(ibam_records)}",
        f"myht_records_in_span: {len(myht_records)}",
        f"key_residue_mode: {cfg.key_residue_mode}",
        f"key_residue_count: {cfg.key_residue_count}",
        "weights:",
    ]
    for k, v in cfg.weights.items():
        lines.append(f"  {k}: {v}")
    lines.append("thresholds:")
    for k, v in cfg.thresholds.items():
        lines.append(f"  {k}: {v}")
    for note in socket2_blob.get("notes", []):
        lines.append(f"socket2_note: {note}")
    path.write_text("\n".join(lines) + "\n")


def write_markdown_report(path: Path, cfg: RunConfig, candidates: List[CandidateModel], verdict: FinalVerdict,
                          ibam_records: List[ResidueRecord], myht_records: List[ResidueRecord], socket2_blob: Dict[str, Any]) -> None:
    top = candidates[0] if candidates else None
    key_ibam = [f"{r.resname}{r.resid}" for r in ibam_records if r.is_key]
    key_myht = [f"{r.resname}{r.resid}" for r in myht_records if r.is_key]

    lines = []
    lines.append(f"# HARP Report — {cfg.run_id}")
    lines.append("")
    lines.append("## Run summary")
    lines.append("")
    lines.append(f"- Structure: `{cfg.structure}`")
    lines.append(f"- IBAM chain/range: `{cfg.ibam_chain}` `{cfg.ibam_range[0]}-{cfg.ibam_range[1]}`")
    lines.append(f"- MyhT chain/range: `{cfg.myht_chain}` `{cfg.myht_range[0]}-{cfg.myht_range[1]}`")
    lines.append(f"- PRCO table: `{cfg.prco_table}`")
    lines.append(f"- Socket2 output: `{cfg.socket2_output}`")
    lines.append(f"- Socket2 status: `{socket2_blob.get('status', 'unknown')}`")
    lines.append("")
    lines.append("## Key residue definition")
    lines.append("")
    lines.append(f"- Mode: `{cfg.key_residue_mode}`")
    lines.append(f"- Count per chain (best effort within span): `{cfg.key_residue_count}`")
    lines.append(f"- IBAM key residues: {', '.join(key_ibam) if key_ibam else 'None'}")
    lines.append(f"- MyhT key residues: {', '.join(key_myht) if key_myht else 'None'}")
    lines.append("")
    lines.append("## Candidate ranking")
    lines.append("")
    lines.append("| Rank | Candidate | Phase | Contact | Hydro Core | Charge | Periodicity | Total | Socket2 a/d | Socket2 KIH | Status |")
    lines.append("|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---|")
    for c in candidates:
        lines.append(
            f"| {c.final_rank} | {c.candidate_id} | {c.phase_offset} | "
            f"{c.contact_alignment_score:.3f} | {c.hydrophobic_core_score:.3f} | "
            f"{c.charge_pattern_score:.3f} | {c.periodicity_score:.3f} | {c.total_internal_score:.3f} | "
            f"{c.socket2_ad_overlap:.3f} | {c.socket2_kih_overlap:.3f} | {c.socket2_status} |"
        )
    lines.append("")
    lines.append("## Preferred model")
    lines.append("")
    if top is not None:
        lines.append(f"- Top candidate: `{top.candidate_id}`")
        lines.append(f"- Phase offset: `{top.phase_offset}`")
        lines.append(f"- Total internal score: `{top.total_internal_score:.4f}`")
        lines.append(f"- Socket2 a/d overlap: `{top.socket2_ad_overlap:.4f}`")
        lines.append(f"- Socket2 KIH overlap: `{top.socket2_kih_overlap:.4f}`")
        lines.append(f"- Socket2 status: `{top.socket2_status}`")
    else:
        lines.append("No candidate model available.")
    lines.append("")
    lines.append("## Final verdict")
    lines.append("")
    lines.append(f"**{verdict.verdict}**")
    lines.append("")
    lines.append(f"- Internal margin: `{verdict.internal_margin:.4f}`")
    lines.append(f"- a/d overlap: `{verdict.ad_overlap:.4f}`")
    lines.append(f"- KIH overlap: `{verdict.kih_overlap:.4f}`")
    lines.append(f"- Rationale: {verdict.rationale}")
    if socket2_blob.get("notes"):
        lines.append("")
        lines.append("## Socket2 parser notes")
        lines.append("")
        for note in socket2_blob.get("notes", []):
            lines.append(f"- {note}")
    lines.append("")
    lines.append("## Notes")
    lines.append("")
    lines.append("This scaffold uses normalized Socket2 outputs from parse_socket2.py.")
    lines.append("KIH concordance is still conservative and should be treated as an external-support proxy rather than exact geometric equivalence.")
    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description="Run Operation HARP v2 scaffold")
    parser.add_argument("--config", required=True, help="Path to YAML or JSON config")
    parser.add_argument("--outdir", default="results", help="Base results directory")
    args = parser.parse_args()

    cfg = load_config(Path(args.config))
    ibam_start, ibam_end = ensure_range_ok("ibam_range", cfg.ibam_range)
    myht_start, myht_end = ensure_range_ok("myht_range", cfg.myht_range)

    structure_path = Path(cfg.structure)
    if not structure_path.exists():
        raise FileNotFoundError(f"Structure file not found: {structure_path}")

    ibam_all, myht_all = load_prco_table(Path(cfg.prco_table), cfg.ibam_chain, cfg.myht_chain)
    ibam_records = filter_span(ibam_all, ibam_start, ibam_end)
    myht_records = filter_span(myht_all, myht_start, myht_end)

    if not ibam_records:
        raise ValueError("No IBAM PRCO records fall within the declared IBAM span")
    if not myht_records:
        raise ValueError("No MyhT PRCO records fall within the declared MyhT span")

    mark_key_residues(ibam_records, cfg.key_residue_mode, cfg.key_residue_count)
    mark_key_residues(myht_records, cfg.key_residue_mode, cfg.key_residue_count)

    candidates = build_candidates(ibam_records, myht_records, cfg.weights)
    socket2_blob = load_socket2_normalized(Path(cfg.socket2_output) if cfg.socket2_output else None)
    residue_matches = compare_with_socket2(candidates, ibam_records, myht_records, socket2_blob)

    for rec in list(ibam_records) + list(myht_records):
        rec.socket2_heptad = socket2_blob.get("ad_map", {}).get((rec.chain, rec.resid))
        rec.ad_concordant = residue_matches.get((rec.chain, rec.resid))
        rec.kih_concordant = None

    verdict = classify_result(candidates, cfg.thresholds)

    result_dir = Path(args.outdir) / cfg.run_id
    result_dir.mkdir(parents=True, exist_ok=True)
    write_candidate_table(result_dir / "candidate_scores.tsv", candidates)
    write_residue_assignments(result_dir / "residue_assignments.tsv", ibam_records, myht_records, candidates[0])
    write_socket2_table(result_dir / "socket2_comparison.tsv", candidates)
    write_qc_log(result_dir / "qc.log", cfg, (ibam_start, ibam_end), (myht_start, myht_end), ibam_records, myht_records, socket2_blob)
    write_markdown_report(result_dir / "report.md", cfg, candidates, verdict, ibam_records, myht_records, socket2_blob)

    print(f"HARP v2 scaffold completed: {result_dir}")
    print(f"Top verdict: {verdict.verdict}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
