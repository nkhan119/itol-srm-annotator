#!/usr/bin/env python3
"""
itol_annotate.py
────────────────
Generate iTOL annotation files for visualising SRM (Short Repeat Motif)
counts on a Legionella phylogenetic tree.

Outputs
-------
  itol_label_replace.txt   Node label replacement
                           Format: GCA_XXXXXX|WP_XXXXXX_Species (SRM=N)
  itol_srm_gradient.txt    Heatmap colour-strip  (green → yellow → red)
  itol_srm_bars.txt        Simple bar chart with per-bar colour coding

Usage
-----
  python itol_annotate.py --tree legionella.treefile --csv summary.csv
  python itol_annotate.py --tree legionella.treefile --csv summary.csv --outdir results/
  python itol_annotate.py --help

Author : Nadeem
License: MIT
"""

import argparse
import csv
import os
import re
import sys
from pathlib import Path


# ─────────────────────────────────────────────────────────────────────────────
# Colour helpers
# ─────────────────────────────────────────────────────────────────────────────

def _lerp_colour(c1: tuple, c2: tuple, t: float) -> str:
    """Linearly interpolate between two RGB tuples and return hex colour."""
    r = int(c1[0] + t * (c2[0] - c1[0]))
    g = int(c1[1] + t * (c2[1] - c1[1]))
    b = int(c1[2] + t * (c2[2] - c1[2]))
    return f"#{r:02x}{g:02x}{b:02x}"


# Palette: green → yellow → red  (ColourBrewer RdYlGn)
_GREEN  = (0x1a, 0x98, 0x50)   # #1a9850
_YELLOW = (0xff, 0xff, 0xbf)   # #ffffbf
_RED    = (0xd7, 0x30, 0x27)   # #d73027


def srm_colour(val: int, lo: int, mid: int, hi: int) -> str:
    """Map an SRM value onto the green→yellow→red gradient."""
    if hi == lo:
        return f"#{_YELLOW[0]:02x}{_YELLOW[1]:02x}{_YELLOW[2]:02x}"
    if val <= mid:
        t = (val - lo) / max(mid - lo, 1)
        return _lerp_colour(_GREEN, _YELLOW, t)
    else:
        t = (val - mid) / max(hi - mid, 1)
        return _lerp_colour(_YELLOW, _RED, t)


# ─────────────────────────────────────────────────────────────────────────────
# Parsers
# ─────────────────────────────────────────────────────────────────────────────

def parse_tree(tree_file: str) -> dict:
    """
    Extract all leaf labels matching GCA/GCF accession format from a
    Newick tree file.

    Returns
    -------
    dict  {accession: full_original_label}
          e.g. {"GCA_000211115.2": "GCA_000211115.2|WP_038839603.1"}
    """
    with open(tree_file) as fh:
        tree = fh.read()

    raw = re.findall(r'(GC[AF]_\d+\.\d+\|[^,:();\s]+)', tree)
    tree_ids = {}
    for label in raw:
        acc = label.split("|")[0]
        tree_ids[acc] = label

    if not tree_ids:
        sys.exit(
            f"[ERROR] No labels matching GCA/GCF pattern found in '{tree_file}'.\n"
            "        Expected format: GCA_XXXXXX.X|WP_XXXXXX"
        )
    return tree_ids


def parse_csv(csv_file: str) -> dict:
    """
    Read genome metadata from CSV.

    Required columns: GENOME, Species, SRM_COUNT

    Returns
    -------
    dict  {accession: {"species": str, "srm": int}}
    """
    required = {"GENOME", "Species", "SRM_COUNT"}

    csv_data = {}
    with open(csv_file, newline="") as fh:
        reader = csv.DictReader(fh)
        missing = required - set(reader.fieldnames or [])
        if missing:
            sys.exit(
                f"[ERROR] CSV is missing required columns: {missing}\n"
                f"        Found columns: {list(reader.fieldnames)}"
            )
        for row in reader:
            acc = row.get("GENOME", "").strip()
            if not acc:
                continue
            species = row.get("Species", "").strip()
            species = re.sub(r"[\"'\s]+", "_", species).strip("_")
            try:
                srm = int(float(row.get("SRM_COUNT", "0")))
            except ValueError:
                srm = 0
            csv_data[acc] = {"species": species or "Unknown", "srm": srm}

    return csv_data


# ─────────────────────────────────────────────────────────────────────────────
# Writers
# ─────────────────────────────────────────────────────────────────────────────

def write_labels(tree_ids: dict, csv_data: dict, out_path: str) -> int:
    """
    Write LABELS annotation file.

    New label format:
        GCA_000211115.2|WP_038839603.1_Legionella_pneumophila (SRM=N)
    """
    matched = 0
    unmatched = []

    with open(out_path, "w") as fh:
        fh.write("LABELS\n")
        fh.write("SEPARATOR COMMA\n")
        fh.write("DATA\n")
        for acc, orig in tree_ids.items():
            if acc in csv_data:
                sp  = csv_data[acc]["species"]
                srm = csv_data[acc]["srm"]
                fh.write(f"{orig},{orig}_{sp} (SRM={srm})\n")
                matched += 1
            else:
                fh.write(f"{orig},{orig}\n")
                unmatched.append(acc)

    if unmatched:
        print(
            f"  [WARN] {len(unmatched)} tree accession(s) had no CSV match "
            f"(kept original label): {unmatched[:3]}"
            f"{'...' if len(unmatched) > 3 else ''}"
        )
    return matched


def write_gradient(tree_ids: dict, csv_data: dict, srm_stats: dict, out_path: str) -> None:
    """Write DATASET_GRADIENT (heatmap colour strip) annotation file."""
    lo, mid, hi = srm_stats["min"], srm_stats["mid"], srm_stats["max"]

    with open(out_path, "w") as fh:
        fh.write("DATASET_GRADIENT\n")
        fh.write("SEPARATOR COMMA\n")
        fh.write("DATASET_LABEL,SRM Count\n")
        fh.write("COLOR,#d73027\n")
        fh.write("COLOR_MIN,#1a9850\n")
        fh.write("COLOR_MID,#ffffbf\n")
        fh.write("COLOR_MAX,#d73027\n")
        fh.write(f"USER_MIN_VALUE,{lo}\n")
        fh.write(f"USER_MID_VALUE,{mid}\n")
        fh.write(f"USER_MAX_VALUE,{hi}\n")
        fh.write("MARGIN,15\n")
        fh.write("STRIP_WIDTH,50\n")
        fh.write("BORDER_WIDTH,0\n")
        fh.write("SHOW_INTERNAL,0\n")
        fh.write("LEGEND_TITLE,SRM Count\n")
        fh.write("LEGEND_SHAPES,1,1,1\n")
        fh.write("LEGEND_COLORS,#1a9850,#ffffbf,#d73027\n")
        fh.write(f"LEGEND_LABELS,Low ({lo}),Mid ({mid}),High ({hi})\n")
        fh.write("DATA\n")
        for acc, orig in tree_ids.items():
            if acc in csv_data:
                fh.write(f"{orig},{csv_data[acc]['srm']}\n")


def write_bars(tree_ids: dict, csv_data: dict, srm_stats: dict, out_path: str) -> None:
    """Write DATASET_SIMPLEBAR annotation file with per-bar colour coding."""
    lo, mid, hi = srm_stats["min"], srm_stats["mid"], srm_stats["max"]

    with open(out_path, "w") as fh:
        fh.write("DATASET_SIMPLEBAR\n")
        fh.write("SEPARATOR COMMA\n")
        fh.write("DATASET_LABEL,SRM Count\n")
        fh.write("COLOR,#1a9850\n")
        fh.write("MARGIN,15\n")
        fh.write("WIDTH,300\n")
        fh.write("BORDER_WIDTH,0\n")
        fh.write("SHOW_VALUE,1\n")
        fh.write("SHOW_INTERNAL,0\n")
        fh.write("LEGEND_TITLE,SRM Count\n")
        fh.write("LEGEND_SHAPES,1,1,1\n")
        fh.write("LEGEND_COLORS,#1a9850,#ffffbf,#d73027\n")
        fh.write(f"LEGEND_LABELS,Low ({lo}),Mid ({mid}),High ({hi})\n")
        fh.write("DATA\n")
        for acc, orig in tree_ids.items():
            if acc in csv_data:
                val = csv_data[acc]["srm"]
                col = srm_colour(val, lo, mid, hi)
                fh.write(f"{orig},{val},{col}\n")


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="itol_annotate",
        description=(
            "Generate iTOL annotation files for SRM count visualisation\n"
            "on a Legionella phylogenetic tree."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python itol_annotate.py --tree legionella.treefile --csv summary.csv\n"
            "  python itol_annotate.py --tree legionella.treefile --csv summary.csv --outdir results/\n"
        ),
    )
    p.add_argument(
        "--tree", "-t",
        required=True,
        metavar="TREEFILE",
        help="Newick tree file (e.g. legionella.treefile)"
    )
    p.add_argument(
        "--csv", "-c",
        required=True,
        metavar="CSV",
        help="CSV with columns: GENOME, Species, SRM_COUNT"
    )
    p.add_argument(
        "--outdir", "-o",
        default=".",
        metavar="DIR",
        help="Output directory (default: current directory)"
    )
    p.add_argument(
        "--prefix", "-p",
        default="itol",
        metavar="PREFIX",
        help="Prefix for output filenames (default: itol)"
    )
    return p


def main() -> None:
    args = build_parser().parse_args()

    # Validate inputs
    for path, name in [(args.tree, "tree"), (args.csv, "CSV")]:
        if not os.path.isfile(path):
            sys.exit(f"[ERROR] {name} file not found: '{path}'")

    # Create output directory
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    label_file = str(outdir / f"{args.prefix}_label_replace.txt")
    grad_file  = str(outdir / f"{args.prefix}_srm_gradient.txt")
    bar_file   = str(outdir / f"{args.prefix}_srm_bars.txt")

    print("─" * 55)
    print("  iTOL SRM Annotation Generator")
    print("─" * 55)
    print(f"  Tree : {args.tree}")
    print(f"  CSV  : {args.csv}")
    print(f"  Out  : {args.outdir}/")
    print("─" * 55)

    # Parse inputs
    print("\n[1/5] Parsing tree labels …")
    tree_ids = parse_tree(args.tree)
    print(f"       Found {len(tree_ids)} leaf labels")

    print("[2/5] Reading CSV …")
    csv_data = parse_csv(args.csv)
    print(f"       Loaded {len(csv_data)} genome records")

    # Compute SRM statistics
    print("[3/5] Computing SRM statistics …")
    srm_values = [csv_data[acc]["srm"] for acc in tree_ids if acc in csv_data]
    if not srm_values:
        sys.exit("[ERROR] No tree accessions matched any CSV GENOME entry.")

    srm_stats = {
        "min": min(srm_values),
        "max": max(srm_values),
        "mid": round((min(srm_values) + max(srm_values)) / 2),
    }
    print(f"       SRM range: {srm_stats['min']} – {srm_stats['max']}  "
          f"(mid = {srm_stats['mid']})")

    # Write outputs
    print("[4/5] Writing annotation files …")
    matched = write_labels(tree_ids, csv_data, label_file)
    print(f"       Labels replaced : {matched}/{len(tree_ids)}")

    write_gradient(tree_ids, csv_data, srm_stats, grad_file)
    write_bars(tree_ids, csv_data, srm_stats, bar_file)

    print("[5/5] Done!\n")
    print("  Upload these files to iTOL → Datasets → Upload annotation file:")
    print(f"    📄  {label_file}")
    print(f"    🌡️   {grad_file}")
    print(f"    📊  {bar_file}")
    print()


if __name__ == "__main__":
    main()
