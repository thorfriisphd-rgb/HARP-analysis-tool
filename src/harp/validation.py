from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any
import json

import MDAnalysis as mda
from MDAnalysis.exceptions import NoDataError
import yaml

from .deepcoil import parse_deepcoil2


AA3_TO_1 = {
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G",
    "HIS":"H","HSD":"H","HSE":"H","HSP":"H","HID":"H","HIE":"H","HIP":"H",
    "ILE":"I","LEU":"L","LYS":"K","MET":"M","MSE":"M","PHE":"F","PRO":"P",
    "SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V","ASH":"D","GLH":"E",
    "LYN":"K","CYM":"C","CYX":"C",
}


@dataclass(frozen=True)
class ValidationIssue:
    code: str
    message: str
    details: dict[str, Any]


@dataclass(frozen=True)
class ValidationReport:
    status: str
    config: str
    deepcoil_file: str | None
    topology_file: str | None
    fasta_file: str | None
    sequence_length: int | None
    cc_segment: list[int] | None
    a_origin: int | None
    d_offset: int | None
    issues: list[ValidationIssue]

    @property
    def passed(self) -> bool:
        return self.status == "PASS"

    def to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        d["passed"] = self.passed
        return d

    def write_json(self, path: str | Path) -> None:
        Path(path).write_text(json.dumps(self.to_dict(), indent=2) + "\n", encoding="utf-8")


def _read_fasta(path: Path) -> str:
    if not path.is_file():
        raise FileNotFoundError(path)
    seq = "".join(
        line.strip().upper()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.startswith(">")
    )
    if not seq:
        raise ValueError("FASTA contains no sequence")
    return seq


def _residue_sequence(residues) -> tuple[str, list[dict[str, Any]]]:
    letters: list[str] = []
    unknown: list[dict[str, Any]] = []
    for residue in residues:
        name = str(residue.resname).strip().upper()
        aa = AA3_TO_1.get(name)
        if aa is None:
            unknown.append({"resid": int(residue.resid), "resname": name})
            letters.append("X")
        else:
            letters.append(aa)
    return "".join(letters), unknown


def _first_mismatch(a: str, b: str) -> dict[str, Any] | None:
    if len(a) != len(b):
        return {"length_a": len(a), "length_b": len(b)}
    for i, (aa, bb) in enumerate(zip(a, b), start=1):
        if aa != bb:
            return {"position": i, "a": aa, "b": bb}
    return None


def _candidate_molecules(universe):
    protein = universe.select_atoms("protein")
    seen: set[tuple[int, ...]] = set()
    if len(protein) == 0:
        return

    try:
        for molnum in sorted(set(int(x) for x in protein.molnums)):
            atoms = protein[protein.molnums == molnum]
            residues = atoms.residues
            key = tuple(int(r.ix) for r in residues)
            if key and key not in seen:
                seen.add(key)
                yield f"molnum={molnum}", residues
    except NoDataError:
        pass

    try:
        for number, fragment in enumerate(protein.fragments, start=1):
            residues = fragment.residues
            key = tuple(int(r.ix) for r in residues)
            if key and key not in seen:
                seen.add(key)
                yield f"fragment={number}", residues
    except NoDataError:
        pass

    for segment in universe.segments:
        atoms = segment.atoms.select_atoms("protein")
        if len(atoms) == 0:
            continue
        residues = atoms.residues
        key = tuple(int(r.ix) for r in residues)
        if key and key not in seen:
            seen.add(key)
            yield f"segid={segment.segid!r}", residues


def validate_config(config_path: str | Path) -> ValidationReport:
    """Validate a HARP v4.1 dataset without changing any input file.

    Validation is intentionally diagnostic-only: no reversal, trimming,
    threshold relaxation, sequence repair, or config rewriting is performed.
    """
    config_path = Path(config_path).resolve()
    cfg = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    base = config_path.parent

    issues: list[ValidationIssue] = []

    def resolve(value: str | None) -> Path | None:
        if not value:
            return None
        p = Path(value).expanduser()
        return p if p.is_absolute() else (base / p).resolve()

    deepcoil_cfg = cfg.get("deepcoil", {})
    traj_cfg = cfg.get("trajectory", {})

    deepcoil_path = resolve(deepcoil_cfg.get("path") or deepcoil_cfg.get("file"))
    topology_path = resolve(traj_cfg.get("topology"))
    fasta_path = resolve(
        cfg.get("myht_fasta")
        or cfg.get("fasta")
        or deepcoil_cfg.get("fasta")
    )

    for label, path in (("DEEPCOIL", deepcoil_path), ("TOPOLOGY", topology_path)):
        if path is None:
            issues.append(ValidationIssue(
                f"MISSING_{label}_PATH", f"{label.title()} path is not defined in config.", {}
            ))
        elif not path.is_file():
            issues.append(ValidationIssue(
                f"MISSING_{label}_FILE", f"{label.title()} file does not exist.",
                {"path": str(path)}
            ))

    reg = None
    if deepcoil_path is not None and deepcoil_path.is_file():
        try:
            reg = parse_deepcoil2(
                deepcoil_path,
                cc_threshold=float(deepcoil_cfg.get("cc_threshold", 0.5)),
                anchor_threshold=float(deepcoil_cfg.get("anchor_threshold", 0.5)),
            )
        except Exception as exc:
            text = str(exc)
            if "No residues exceed cc threshold" in text:
                code = "NO_VALID_CC_SEGMENT"
            elif "Too few confident DeepCoil" in text:
                code = "TOO_FEW_AD_ANCHORS"
            elif "not internally canonical" in text:
                code = "NONCANONICAL_AD_OFFSET"
            else:
                code = "DEEPCOIL_REGISTER_ERROR"
            issues.append(ValidationIssue(code, text, {"path": str(deepcoil_path)}))

    fasta_seq = None
    if fasta_path is not None:
        try:
            fasta_seq = _read_fasta(fasta_path)
        except Exception as exc:
            issues.append(ValidationIssue(
                "FASTA_READ_ERROR", str(exc), {"path": str(fasta_path)}
            ))

    if reg is not None and fasta_seq is not None and fasta_seq != reg.sequence:
        if fasta_seq == reg.sequence[::-1]:
            code = "DEEPCOIL_REVERSED"
            message = "DeepCoil sequence is the exact reverse of the FASTA."
        else:
            code = "FASTA_DEEPCOIL_MISMATCH"
            message = "FASTA and DeepCoil sequences are not identical."
        issues.append(ValidationIssue(
            code, message,
            {"fasta_length": len(fasta_seq),
             "deepcoil_length": len(reg.sequence),
             "first_difference": _first_mismatch(fasta_seq, reg.sequence)}
        ))

    if reg is not None and topology_path is not None and topology_path.is_file():
        try:
            u = mda.Universe(str(topology_path))
            exact: list[str] = []
            reverse: list[str] = []
            containing: list[dict[str, Any]] = []
            unknowns: list[dict[str, Any]] = []

            for label, residues in _candidate_molecules(u):
                seq, unk = _residue_sequence(residues)
                unknowns.extend({"group": label, **x} for x in unk)
                if "X" in seq:
                    continue
                if seq == reg.sequence:
                    exact.append(label)
                elif seq == reg.sequence[::-1]:
                    reverse.append(label)
                elif reg.sequence in seq:
                    start = seq.index(reg.sequence) + 1
                    containing.append({
                        "group": label, "start": start,
                        "stop": start + len(reg.sequence) - 1,
                        "group_length": len(seq),
                    })

            if unknowns:
                issues.append(ValidationIssue(
                    "UNKNOWN_TOPOLOGY_RESIDUE",
                    "One or more protein residues cannot be translated unambiguously.",
                    {"residues": unknowns},
                ))

            if len(exact) == 0:
                issues.append(ValidationIssue(
                    "NO_EXACT_TPR_MATCH",
                    "No unique complete protein molecule in the topology exactly matches the DeepCoil sequence.",
                    {"reverse_matches": reverse, "containing_matches": containing},
                ))
            elif len(exact) > 1:
                issues.append(ValidationIssue(
                    "MULTIPLE_TPR_MATCHES",
                    "More than one topology molecule exactly matches the DeepCoil sequence.",
                    {"exact_matches": exact},
                ))
        except Exception as exc:
            issues.append(ValidationIssue(
                "TOPOLOGY_VALIDATION_ERROR", str(exc), {"path": str(topology_path)}
            ))

    status = "PASS" if not issues else "FAIL"
    return ValidationReport(
        status=status,
        config=str(config_path),
        deepcoil_file=str(deepcoil_path) if deepcoil_path else None,
        topology_file=str(topology_path) if topology_path else None,
        fasta_file=str(fasta_path) if fasta_path else None,
        sequence_length=reg.length if reg is not None else None,
        cc_segment=[reg.cc_start, reg.cc_end] if reg is not None else None,
        a_origin=reg.a_origin if reg is not None else None,
        d_offset=reg.d_offset if reg is not None else None,
        issues=issues,
    )


def format_validation_report(report: ValidationReport) -> str:
    lines = [
        f"HARP v4.1 validation: {report.status}",
        f"Config: {report.config}",
    ]
    if report.sequence_length is not None:
        lines += [
            f"Sequence length: {report.sequence_length}",
            f"CC segment: {report.cc_segment[0]}-{report.cc_segment[1]}",
            f"a origin: {report.a_origin}",
            f"d offset: a+{report.d_offset}",
        ]
    if report.issues:
        lines.append("Issues:")
        for issue in report.issues:
            lines.append(f"  [{issue.code}] {issue.message}")
            if issue.details:
                lines.append(f"    {json.dumps(issue.details, sort_keys=True)}")
    else:
        lines.append("No validation failures detected.")
    return "\n".join(lines)
