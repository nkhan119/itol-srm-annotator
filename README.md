# iTOL SRM Annotator

> Generate publication-ready [iTOL](https://itol.embl.de/) annotation files for visualising **Short Repeat Motif (SRM)** counts on phylogenetic trees — purpose-built for *Legionella* genomics.

## Features

- **Replaces raw node labels** with clean, informative labels:
  `GCA_000211115.2|WP_038839603.1_Legionella_pneumophila (SRM=5)`
- **Heatmap gradient strip** — green → yellow → red, scaled dynamically to your data range
- **Colour-coded bar chart** — exact SRM counts displayed as bars with per-bar colour matching the heatmap
- **CLI interface** with `--tree`, `--csv`, `--outdir`, and `--prefix` options
- Warns you about any unmatched accessions rather than silently dropping them

---

## Repository Structure
```
itol-srm-annotator/
├── itol_annotate.py        # Main script
├── README.md
```

> **Data Availability** — The *Legionella* genome tree and SRM summary CSV used in this study are available upon reasonable request. Please contact [nadeem.khan@inrs.ca](mailto:nkhan119@uottawa.ca).

---

## Requirements

- Python **3.7+**
- Standard library only — no external dependencies

---

## Quick Start

### 1. Clone the repository
```bash
git clone https://github.com/nkhan119/itol-srm-annotator.git
cd itol-srm-annotator
```

### 2. Prepare your inputs

**Tree file** — a standard Newick file with labels in the format:
```
GCA_000211115.2|WP_038839603.1
```

**CSV file** — must contain these three columns:

| Column | Description |
|--------|-------------|
| `GENOME` | NCBI accession (e.g. `GCA_000211115.2`) |
| `Species` | Full species name (e.g. `Legionella pneumophila`) |
| `SRM_COUNT` | Number of Short Repeat Motifs found |

### 3. Run the script
```bash
# Basic usage
python itol_annotate.py --tree legionella.treefile --csv summary.csv

# Specify output directory and filename prefix
python itol_annotate.py \
  --tree legionella.treefile \
  --csv  summary.csv \
  --outdir results/ \
  --prefix legionella
```

### 4. Upload to iTOL

Upload the three generated `.txt` files to your iTOL tree:

1. Open your tree on [iTOL](https://itol.embl.de/)
2. Go to **Datasets** → **Upload annotation file**
3. Upload in this order:
   - `*_label_replace.txt` — replaces node labels
   - `*_srm_gradient.txt` — adds heatmap colour strip
   - `*_srm_bars.txt` — adds bar chart

---

## Output Files

### `*_label_replace.txt` — Node Label Replacement
Replaces raw accession labels with:
```
GCA_000211115.2|WP_038839603.1_Legionella_pneumophila (SRM=5)
```

### `*_srm_gradient.txt` — Heatmap Strip

| Colour | Meaning |
|--------|---------|
| 🟢 Green `#1a9850` | Low SRM count |
| 🟡 Yellow `#ffffbf` | Mid SRM count |
| 🔴 Red `#d73027` | High SRM count |

### `*_srm_bars.txt` — Bar Chart
Per-bar colour matches the gradient, with numeric values printed on each bar.

---

## Full CLI Options
```
usage: itol_annotate [-h] --tree TREEFILE --csv CSV [--outdir DIR] [--prefix PREFIX]

  --tree  / -t    Newick tree file  (required)
  --csv   / -c    CSV metadata file (required)
  --outdir / -o   Output directory  (default: current directory)
  --prefix / -p   Output filename prefix (default: itol)
```

---

## Background

**Short Repeat Motifs (SRMs)** are repetitive sequence elements found in bacterial genomes. Mapping their counts onto a phylogenetic tree allows rapid visual identification of lineages with elevated repeat content, which may be associated with genomic instability or host adaptation in *Legionella*.

---

## Citation / Reuse

If you use or adapt this pipeline, please credit the author.

---

## Author
---

**Nadeem Khan, PhD**
Bioinformatician — INRS–Centre Armand-Frappier Santé-Biotechnologie, Laval, QC, Canada
nkhan119@uottawa.ca
[@nkhan119](https://github.com/nkhan119)
