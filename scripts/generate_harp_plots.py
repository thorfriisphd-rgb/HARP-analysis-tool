#!/usr/bin/env python3
"""
Generate the three manuscript HARP plots from a panel summary TSV/CSV.

Canonical output is SVG, produced with only the Python standard library.
Optional PNG/PDF conversion is supported via CairoSVG when installed:

    pip install cairosvg

Example:
    python scripts/generate_harp_plots.py \
        --summary results/panel_26/harp_panel_summary.tsv \
        --out-dir results/panel_26/plots \
        --formats svg,png,pdf

Notes:
    - SVG is always generated first and is treated as the canonical output.
    - PNG/PDF are conversions of the SVG files, not independently redrawn plots.
    - No HARP scoring/statistical logic is implemented here; this is visualization only.
    - SVG and converted PNG outputs use a white background for readability in dark-mode viewers.
"""

from __future__ import annotations

import argparse
import csv
import html
import math
from pathlib import Path
from typing import Dict, List, Tuple

try:
    import cairosvg  # type: ignore
    HAS_CAIROSVG = True
except Exception:
    cairosvg = None  # type: ignore
    HAS_CAIROSVG = False


REQUIRED = [
    "taxon",
    "full_best_phase",
    "full_best_score",
    "full_margin",
    "full_best_minus_flat",
    "full_phase_ranking",
    "full_null_p95_best",
    "full_p_best",
    "full_q_best",
    "full_null_stability_rate",
]

DEFAULT_FORMATS = "svg,png,pdf"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Generate standard SVG/PNG/PDF plots from HARP summary output."
    )
    p.add_argument("--summary", required=True, help="Path to harp_panel_summary.tsv or .csv")
    p.add_argument("--out-dir", required=True, help="Output directory for plots")
    p.add_argument("--top-n", type=int, default=26, help="Number of taxa to show in ranking plots")
    p.add_argument(
        "--formats",
        default=DEFAULT_FORMATS,
        help="Comma-separated output formats. Supported: svg,png,pdf. Default: svg,png,pdf",
    )
    p.add_argument(
        "--png-dpi",
        type=int,
        default=300,
        help="DPI for PNG conversion when PNG output is requested. Default: 300",
    )
    p.add_argument(
        "--margin-threshold",
        type=float,
        default=0.05,
        help="Visual guide threshold for phase margin in score-vs-margin plot. Default: 0.05",
    )
    p.add_argument(
        "--no-scatter-guides",
        action="store_true",
        help="Suppress quadrant labels and guide line in score-vs-margin plot.",
    )
    return p.parse_args()


def parse_formats(text: str) -> List[str]:
    allowed = {"svg", "png", "pdf"}
    formats = []
    for item in text.split(","):
        fmt = item.strip().lower()
        if not fmt:
            continue
        if fmt not in allowed:
            raise ValueError(f"Unsupported output format '{fmt}'. Supported formats: svg,png,pdf")
        if fmt not in formats:
            formats.append(fmt)
    if not formats:
        raise ValueError("At least one output format must be requested.")
    if "svg" not in formats:
        # SVG is required internally because PNG/PDF are generated from SVG.
        formats.insert(0, "svg")
    return formats


def detect_delimiter(path: Path) -> str:
    lines = path.read_text(errors="replace").splitlines()
    if not lines:
        raise ValueError(f"Summary file is empty: {path}")
    first = lines[0]
    return "\t" if "\t" in first else ","


def to_float(value: str, default: float = math.nan) -> float:
    try:
        return float(value)
    except Exception:
        return default


def to_int(value: str, default: int = -1) -> int:
    try:
        return int(float(value))
    except Exception:
        return default


def load_rows(path: Path) -> List[Dict[str, str]]:
    delim = detect_delimiter(path)
    with path.open("r", newline="") as fh:
        reader = csv.DictReader(fh, delimiter=delim)
        rows = list(reader)
        fieldnames = reader.fieldnames or []
    missing = [c for c in REQUIRED if c not in fieldnames]
    if missing:
        raise ValueError(f"Missing required column(s): {missing}")
    if not rows:
        raise ValueError(f"No data rows found in summary file: {path}")
    return rows


def esc(s: object) -> str:
    return html.escape(str(s), quote=True)


def svg_header(width: int, height: int) -> str:
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">\n'
        f'<rect width="{width}" height="{height}" fill="white"/>\n'
        '<style>\n'
        'text { font-family: Arial, Helvetica, sans-serif; font-size: 12px; }\n'
        '.title { font-size: 18px; font-weight: bold; }\n'
        '.subtitle { font-size: 12px; fill: #555; }\n'
        '.axis { stroke: #222; stroke-width: 1; }\n'
        '.grid { stroke: #ddd; stroke-width: 1; }\n'
        '.guide { stroke: #777; stroke-width: 1; stroke-dasharray: 4 4; }\n'
        '.baseline { stroke: #aa3333; stroke-width: 1; stroke-dasharray: 5 4; }\n'
        '.bar { fill: #6b8dbb; }\n'
        '.dot { fill: #333; }\n'
        '.null { fill: #111; font-size: 16px; }\n'
        '.label { font-size: 11px; }\n'
        '.annotation { font-size: 12px; fill: #555; }\n'
        '</style>\n'
    )


def svg_footer() -> str:
    return "</svg>\n"


def write(path: Path, text: str) -> None:
    path.write_text(text, encoding="utf-8")


def convert_svg(svg_path: Path, formats: List[str], png_dpi: int) -> List[Path]:
    """Convert SVG to PNG/PDF if requested and CairoSVG is available."""
    made: List[Path] = [svg_path]
    requested_conversion = any(fmt in formats for fmt in ("png", "pdf"))
    if not requested_conversion:
        return made

    if not HAS_CAIROSVG:
        print(
            f"[warn] CairoSVG is not installed; keeping SVG only for {svg_path.name}. "
            "Install with: pip install cairosvg"
        )
        return made

    if "png" in formats:
        png_path = svg_path.with_suffix(".png")
        try:
            cairosvg.svg2png(
                url=str(svg_path),
                write_to=str(png_path),
                dpi=png_dpi,
                background_color="white",
            )  # type: ignore[union-attr]
            made.append(png_path)
        except Exception as exc:
            print(f"[warn] PNG conversion failed for {svg_path.name}: {exc}")

    if "pdf" in formats:
        pdf_path = svg_path.with_suffix(".pdf")
        try:
            cairosvg.svg2pdf(url=str(svg_path), write_to=str(pdf_path))  # type: ignore[union-attr]
            made.append(pdf_path)
        except Exception as exc:
            print(f"[warn] PDF conversion failed for {svg_path.name}: {exc}")

    return made


def barh_plot(
    rows: List[Dict[str, str]],
    value_col: str,
    title: str,
    xlabel: str,
    out_path: Path,
    top_n: int,
    baseline: float | None = None,
) -> None:
    data: List[Tuple[str, float]] = []
    for r in rows:
        v = to_float(r.get(value_col, ""))
        if not math.isnan(v):
            data.append((r["taxon"], v))
    data = sorted(data, key=lambda x: x[1], reverse=True)[:top_n]
    data = list(reversed(data))

    left, right, top, bottom = 220, 45, 55, 65
    row_h = 24
    width = 940
    height = top + bottom + row_h * len(data)
    plot_w = width - left - right
    values = [v for _, v in data]
    xmin = min(0.0, min(values), baseline if baseline is not None else 0.0)
    xmax = max(values + ([baseline] if baseline is not None else [0.0]))
    if xmax == xmin:
        xmax = xmin + 1.0

    def xscale(v: float) -> float:
        return left + (v - xmin) / (xmax - xmin) * plot_w

    parts = [svg_header(width, height)]
    parts.append(f'<text x="{left}" y="28" class="title">{esc(title)}</text>\n')
    parts.append(f'<line x1="{left}" y1="{top-10}" x2="{left}" y2="{height-bottom+8}" class="axis"/>\n')
    parts.append(f'<line x1="{left}" y1="{height-bottom+8}" x2="{width-right}" y2="{height-bottom+8}" class="axis"/>\n')

    if baseline is not None:
        bx = xscale(baseline)
        parts.append(f'<line x1="{bx:.1f}" y1="{top-12}" x2="{bx:.1f}" y2="{height-bottom+8}" class="baseline"/>\n')
        parts.append(f'<text x="{bx+4:.1f}" y="{top-18}" class="label">baseline {baseline:.3f}</text>\n')

    for i, (name, value) in enumerate(data):
        y = top + i * row_h
        x0 = xscale(0.0)
        x1 = xscale(value)
        bar_x = min(x0, x1)
        bar_w = abs(x1 - x0)
        parts.append(f'<text x="{left-8}" y="{y+15}" text-anchor="end" class="label">{esc(name)}</text>\n')
        parts.append(f'<rect x="{bar_x:.1f}" y="{y+3}" width="{bar_w:.1f}" height="16" class="bar"/>\n')
        parts.append(f'<text x="{x1+4:.1f}" y="{y+15}" class="label">{value:.3f}</text>\n')

    parts.append(f'<text x="{left + plot_w/2:.1f}" y="{height-18}" text-anchor="middle">{esc(xlabel)}</text>\n')
    parts.append(svg_footer())
    write(out_path, "".join(parts))


def parse_phase_ranking(text: str) -> Dict[int, float]:
    out: Dict[int, float] = {}
    for item in str(text).split(","):
        if ":" not in item:
            continue
        p, s = item.split(":", 1)
        try:
            out[int(p.strip())] = float(s.strip())
        except Exception:
            pass
    return out


def heat_color(value: float, vmin: float, vmax: float) -> str:
    if math.isnan(value):
        return "#f0f0f0"
    t = 0.0 if vmax == vmin else (value - vmin) / (vmax - vmin)
    t = max(0.0, min(1.0, t))
    r = int(240 - 150 * t)
    g = int(248 - 120 * t)
    b = int(255 - 55 * t)
    return f"#{r:02x}{g:02x}{b:02x}"


def phase_heatmap(rows: List[Dict[str, str]], out_path: Path) -> None:
    ordered = sorted(rows, key=lambda r: to_float(r["full_best_score"]), reverse=True)
    vals: List[List[float]] = []
    for r in ordered:
        parsed = parse_phase_ranking(r["full_phase_ranking"])
        vals.append([parsed.get(i, math.nan) for i in range(7)])
    flat = [v for row in vals for v in row if not math.isnan(v)]
    vmin, vmax = min(flat), max(flat)

    left, top = 220, 60
    cell_w, cell_h = 64, 22
    width = left + cell_w * 7 + 70
    height = top + cell_h * len(ordered) + 70
    parts = [svg_header(width, height)]
    parts.append(f'<text x="{left}" y="28" class="title">HARP phase-score landscape</text>\n')
    for j in range(7):
        parts.append(f'<text x="{left+j*cell_w+cell_w/2}" y="{top-12}" text-anchor="middle">phase {j}</text>\n')
    for i, r in enumerate(ordered):
        y = top + i * cell_h
        parts.append(f'<text x="{left-8}" y="{y+15}" text-anchor="end" class="label">{esc(r["taxon"])}</text>\n')
        for j, value in enumerate(vals[i]):
            x = left + j * cell_w
            fill = heat_color(value, vmin, vmax)
            label = "" if math.isnan(value) else f"{value:.2f}"
            parts.append(f'<rect x="{x}" y="{y}" width="{cell_w}" height="{cell_h}" fill="{fill}" stroke="#fff"/>\n')
            parts.append(f'<text x="{x+cell_w/2}" y="{y+15}" text-anchor="middle" class="label">{label}</text>\n')
    parts.append(svg_footer())
    write(out_path, "".join(parts))


def scatter_score_margin(
    rows: List[Dict[str, str]],
    out_path: Path,
    margin_threshold: float = 0.05,
    show_guides: bool = True,
) -> None:
    points = []
    for r in rows:
        x = to_float(r["full_best_score"])
        y = to_float(r["full_margin"])
        if not math.isnan(x) and not math.isnan(y):
            points.append((r["taxon"], x, y))

    width, height = 980, 640
    left, right, top, bottom = 95, 55, 70, 85
    plot_w, plot_h = width - left - right, height - top - bottom
    xs = [p[1] for p in points]
    ys = [p[2] for p in points]
    xmin, xmax = min(xs), max(xs)
    xpad = (xmax - xmin) * 0.04 if xmax > xmin else 0.05
    xmin -= xpad
    xmax += xpad
    ymin = 0.0
    ymax = max(ys) * 1.20 if max(ys) > 0 else 1.0

    def sx(x: float) -> float:
        return left + (x - xmin) / (xmax - xmin) * plot_w if xmax != xmin else left + plot_w / 2

    def sy(y: float) -> float:
        return top + plot_h - (y - ymin) / (ymax - ymin) * plot_h if ymax != ymin else top + plot_h / 2

    parts = [svg_header(width, height)]
    parts.append(f'<text x="{left}" y="28" class="title">HARP best score vs phase margin</text>\n')
    parts.append(
        f'<text x="{left}" y="47" class="subtitle">Score reflects phase-aligned contact enrichment; margin reflects register decisiveness.</text>\n'
    )
    parts.append(f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top+plot_h}" class="axis"/>\n')
    parts.append(f'<line x1="{left}" y1="{top+plot_h}" x2="{left+plot_w}" y2="{top+plot_h}" class="axis"/>\n')

    # Flat a/d expectation for context. This is a visual reference only.
    flat_x = 2 / 7
    if xmin <= flat_x <= xmax:
        bx = sx(flat_x)
        parts.append(f'<line x1="{bx:.1f}" y1="{top}" x2="{bx:.1f}" y2="{top+plot_h}" class="baseline"/>\n')
        parts.append(f'<text x="{bx+5:.1f}" y="{top+14}" class="label">flat a/d expectation</text>\n')

    if show_guides:
        gy = sy(margin_threshold)
        if top <= gy <= top + plot_h:
            parts.append(f'<line x1="{left}" y1="{gy:.1f}" x2="{left+plot_w}" y2="{gy:.1f}" class="guide"/>\n')
            parts.append(
                f'<text x="{left+8}" y="{gy-6:.1f}" class="annotation">visual guide: margin {margin_threshold:.2f}</text>\n'
            )
        parts.append(f'<text x="{left+plot_w-170}" y="{top+28}" class="annotation">strong + specific</text>\n')
        parts.append(f'<text x="{left+plot_w-210}" y="{top+plot_h-18}" class="annotation">strong + multi-register</text>\n')
        parts.append(f'<text x="{left+10}" y="{top+plot_h-18}" class="annotation">weak / ambiguous</text>\n')

    for name, x, y in points:
        parts.append(f'<circle cx="{sx(x):.1f}" cy="{sy(y):.1f}" r="4" class="dot"/>\n')
        # Keep labels small and close; downstream figure editing can adjust if needed.
        display = "Myxine" if name == "Myxine_glutinosa" else name
        parts.append(f'<text x="{sx(x)+6:.1f}" y="{sy(y)-5:.1f}" class="label">{esc(display)}</text>\n')

    parts.append(
        f'<text x="{left+plot_w/2}" y="{height-24}" text-anchor="middle">'
        'Best HARP score (fraction of phase-aligned contacts)</text>\n'
    )
    parts.append(
        f'<text x="24" y="{top+plot_h/2}" transform="rotate(-90 24 {top+plot_h/2})" text-anchor="middle">'
        'Phase margin (Δ score: best − second-best phase)</text>\n'
    )
    parts.append(svg_footer())
    write(out_path, "".join(parts))


def observed_vs_null(rows: List[Dict[str, str]], out_path: Path) -> None:
    data = []
    for r in rows:
        obs = to_float(r["full_best_score"])
        null = to_float(r["full_null_p95_best"])
        if not math.isnan(obs) and not math.isnan(null):
            data.append((r["taxon"], obs, null))
    data = sorted(data, key=lambda x: x[1])

    left, right, top, bottom = 220, 60, 55, 70
    row_h = 24
    width = 960
    height = top + bottom + row_h * len(data)
    plot_w = width - left - right
    xmin = 0.0
    xmax = max(max(obs, null) for _, obs, null in data) * 1.05

    def sx(v: float) -> float:
        return left + (v - xmin) / (xmax - xmin) * plot_w

    parts = [svg_header(width, height)]
    parts.append(f'<text x="{left}" y="28" class="title">Observed HARP score vs null 95th percentile</text>\n')
    parts.append(f'<text x="{width-225}" y="26">● observed   × null p95</text>\n')
    for i, (name, obs, null) in enumerate(data):
        y = top + i * row_h
        parts.append(f'<text x="{left-8}" y="{y+15}" text-anchor="end" class="label">{esc(name)}</text>\n')
        parts.append(f'<circle cx="{sx(obs):.1f}" cy="{y+11}" r="4" class="dot"/>\n')
        parts.append(f'<text x="{sx(null):.1f}" y="{y+15}" text-anchor="middle" class="null">×</text>\n')
    parts.append(
        f'<text x="{left+plot_w/2}" y="{height-18}" text-anchor="middle">'
        'Best HARP score (observed vs shuffle_occupied_span null p95)</text>\n'
    )
    parts.append(svg_footer())
    write(out_path, "".join(parts))


def phase_counts(rows: List[Dict[str, str]], out_path: Path) -> None:
    counts = {i: 0 for i in range(7)}
    for r in rows:
        p = to_int(r["full_best_phase"])
        if p in counts:
            counts[p] += 1

    width, height = 660, 430
    left, right, top, bottom = 70, 35, 55, 65
    plot_w, plot_h = width - left - right, height - top - bottom
    max_count = max(counts.values()) or 1
    bar_w = plot_w / 7 * 0.75
    parts = [svg_header(width, height)]
    parts.append(f'<text x="{left}" y="28" class="title">Distribution of best HARP phases</text>\n')
    parts.append(f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top+plot_h}" class="axis"/>\n')
    parts.append(f'<line x1="{left}" y1="{top+plot_h}" x2="{left+plot_w}" y2="{top+plot_h}" class="axis"/>\n')
    for i in range(7):
        x = left + (i + 0.125) * plot_w / 7
        h = counts[i] / max_count * plot_h
        y = top + plot_h - h
        parts.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_w:.1f}" height="{h:.1f}" class="bar"/>\n')
        parts.append(f'<text x="{x+bar_w/2:.1f}" y="{y-5:.1f}" text-anchor="middle">{counts[i]}</text>\n')
        parts.append(f'<text x="{x+bar_w/2:.1f}" y="{top+plot_h+20}" text-anchor="middle">{i}</text>\n')
    parts.append(f'<text x="{left+plot_w/2}" y="{height-18}" text-anchor="middle">Best HARP phase</text>\n')
    parts.append(svg_footer())
    write(out_path, "".join(parts))


def main() -> None:
    args = parse_args()
    summary = Path(args.summary)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    formats = parse_formats(args.formats)
    rows = load_rows(summary)
    top_n = min(args.top_n, len(rows))

    outputs = [
        ("01_observed_vs_null_p95.svg", lambda p: observed_vs_null(rows, p)),
        (
            "02_score_vs_margin.svg",
            lambda p: scatter_score_margin(
                rows,
                p,
                margin_threshold=args.margin_threshold,
                show_guides=not args.no_scatter_guides,
            ),
        ),
        ("03_phase_score_heatmap.svg", lambda p: phase_heatmap(rows, p)),
    ]

    svg_made: List[Path] = []
    all_made: List[Path] = []
    for filename, func in outputs:
        svg_path = out_dir / filename
        func(svg_path)
        svg_made.append(svg_path)
        all_made.extend(convert_svg(svg_path, formats, args.png_dpi))

    # De-duplicate while preserving order.
    seen = set()
    unique_made: List[Path] = []
    for path in all_made:
        if path not in seen:
            unique_made.append(path)
            seen.add(path)

    manifest = out_dir / "harp_plots_manifest.txt"
    manifest.write_text(
        "HARP plot outputs\n"
        "=================\n\n"
        f"Summary input: {summary}\n"
        f"Requested formats: {','.join(formats)}\n"
        f"CairoSVG available: {HAS_CAIROSVG}\n\n"
        + "\n".join(p.name for p in unique_made)
        + "\n",
        encoding="utf-8",
    )

    print(f"Loaded {len(rows)} HARP rows from: {summary}")
    print(f"Wrote {len(svg_made)} canonical SVG plot(s) to: {out_dir}")
    if any(fmt in formats for fmt in ("png", "pdf")) and not HAS_CAIROSVG:
        print("PNG/PDF conversion requested but skipped because CairoSVG is not installed.")
    else:
        non_svg = [p for p in unique_made if p.suffix.lower() != ".svg"]
        if non_svg:
            print(f"Wrote {len(non_svg)} converted PNG/PDF file(s).")
    for p in unique_made:
        print(f"  - {p}")
    print(f"Manifest: {manifest}")



if __name__ == "__main__":
    main()
