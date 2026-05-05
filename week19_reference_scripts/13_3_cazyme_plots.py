#!/usr/bin/env python3
"""
Visualise dbCAN output across multiple samples.

Expects the following layout:

    <dbcan_dir>/
        sample_A/overview.tsv
        sample_B/overview.tsv
        ...

Each overview.tsv is the standard run_dbcan output with columns:
    Gene ID  EC#  dbCAN_hmm  dbCAN_sub  DIAMOND  #ofTools  Recommend Results  Substrate

Produces, in the chosen output directory:

    01_cazyme_class_counts.png      Stacked-bar of CAZyme class counts per sample
    02_substrate_categories.png     Multi-panel figure (one panel per category):
                                      - HMO GH families
                                      - Plant substrates (pectin / xylan / xyloglucan)
                                      - Glycogen / alpha-glucan families

Also writes the underlying counts as CSVs alongside the plots.

Step 1: Load environment (HPC only)
    module load anaconda3/5.3.1

Step 2: Run the script
    python scripts/13_3_cazyme_plots.py /path/to/dbcan_dir -o cazyme_plot
"""

from __future__ import annotations

import argparse
import math
import re
import sys
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# --------------------------------------------------------------------------- #
# Configuration
# --------------------------------------------------------------------------- #

CAZYME_CLASSES = ["GH", "PL"]

# --- Substrate / family category definitions ------------------------------- #
# Each category renders as a panel in plot 2.
#   kind="family"     -> count occurrences of each listed family
#                        in the "Recommend Results" column
#   kind="substrate"  -> count occurrences of each listed substrate string
#                        in the "Substrate" column

# Human-milk-oligosaccharide-relevant glycoside hydrolase families.
HMO_FAMILIES = ["GH112", "GH136"]

# Plant cell-wall substrates (read from the Substrate column).
PLANT_SUBSTRATES = ["pectin", "xylan", "xyloglucan"]

# Lactobacillus / α-glucan storage-degradation families.
ALPHA_GLUCAN_FAMILIES = [
    "GH13_20", "GH13_30", "GH13_31",
    "GH31_1", "GH65"
]

SUBSTRATE_CATEGORIES = [
    {
        "name": "HMO GH families",
        "kind": "family",
        "items": HMO_FAMILIES,
        "ylabel": "# domains",
        "csv": "hmo_family_counts.csv",
    },
    {
        "name": "Plant substrates (pectin / xylan / xyloglucan)",
        "kind": "substrate",
        "items": PLANT_SUBSTRATES,
        "ylabel": "# substrate annotations",
        "csv": "plant_substrate_counts.csv",
    },
    {
        "name": "Glycogen / α-glucan families",
        "kind": "family",
        "items": ALPHA_GLUCAN_FAMILIES,
        "ylabel": "# domains",
        "csv": "alpha_glucan_family_counts.csv",
    }
]

CLASS_COLORS = {
    "GH": "#1f77b4",
    "PL": "#2ca02c"
}

# Match a CAZyme domain at the start of a token, e.g. "GH112_e85" -> ("GH", "GH112")
_CLASS_RE = re.compile(r"^(GH|PL)\d+")
# Capture the family with optional subfamily suffix:
#   GH8_e4      -> GH8
#   GH13_18     -> GH13_18   (subfamily preserved)
#   GH13_e122   -> GH13      (e-cluster suffix dropped)
#   CBM32_e14   -> CBM32
_FAMILY_RE = re.compile(
    r"^((?:GH|PL)\d+(?:_\d+)?)"
)


# --------------------------------------------------------------------------- #
# Parsing helpers
# --------------------------------------------------------------------------- #

def parse_recommend_field(field: str):
    """Yield individual domain entries from a 'Recommend Results' cell."""
    if not isinstance(field, str):
        return
    field = field.strip()
    if not field or field == "-":
        return
    # dbCAN can join multiple domains for one gene with '|' or '+'.
    for part in re.split(r"[|+]", field):
        part = part.strip()
        if part and part != "-":
            yield part


def domain_to_class(domain: str):
    m = _CLASS_RE.match(domain)
    return m.group(1) if m else None


def domain_to_family(domain: str):
    """e.g. GH8_e4 -> GH8, GH13_18 -> GH13_18, GH13_e122 -> GH13, CBM32_e14 -> CBM32"""
    m = _FAMILY_RE.match(domain)
    if not m:
        return None
    fam = m.group(1)
    # Drop trailing _<number> if it actually represents a dbCAN-sub e-cluster
    # (those look like "_e\d+" which the regex above excluded already, so this
    # branch keeps real subfamily numbers).
    return fam


def parse_substrate_field(field: str):
    """Yield individual substrates (lower-cased) from a Substrate cell.

    dbCAN-sub uses ';' to delimit multi-substrate annotations
    (e.g. 'alpha-glucan;sucrose') and ',' for some older outputs. We allow
    ';', ',', '|', and '/' to be safe.
    """
    if not isinstance(field, str):
        return
    field = field.strip()
    if not field or field == "-":
        return
    for part in re.split(r"[;,|/]", field):
        part = part.strip().lower()
        if part and part != "-":
            yield part


# --------------------------------------------------------------------------- #
# Sample discovery and per-sample counting
# --------------------------------------------------------------------------- #

def find_overview_files(root: Path) -> dict:
    """Return {sample_name: path} for the subfolder-per-sample layout."""
    samples = {}
    for sub in sorted(Path(root).iterdir()):
        if not sub.is_dir():
            continue
        for name in ("overview.tsv", "overview.txt"):
            cand = sub / name
            if cand.is_file():
                samples[sub.name] = cand
                break
    return samples


def load_sample(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("-")
    df.columns = [c.strip() for c in df.columns]
    return df


def count_per_sample(df: pd.DataFrame):
    """Return (class_counter, family_counter, substrate_counter) for one sample."""
    class_c, family_c, sub_c = Counter(), Counter(), Counter()

    rec_col = df.get("Recommend Results")
    if rec_col is not None:
        for cell in rec_col:
            for dom in parse_recommend_field(cell):
                cls = domain_to_class(dom)
                fam = domain_to_family(dom)
                if cls:
                    class_c[cls] += 1
                if fam:
                    family_c[fam] += 1

    sub_col = df.get("Substrate")
    if sub_col is not None:
        for cell in sub_col:
            for s in parse_substrate_field(cell):
                sub_c[s] += 1

    return class_c, family_c, sub_c


def category_matrix(category: dict, family_counts: dict, substrate_counts: dict,
                    sample_order: list) -> pd.DataFrame:
    """Build a (rows=item, cols=sample) DataFrame for one substrate category."""
    items = category["items"]
    if category["kind"] == "family":
        data = {
            s: {it: family_counts[s].get(it, 0) for it in items}
            for s in sample_order
        }
    elif category["kind"] == "substrate":
        data = {
            s: {it: substrate_counts[s].get(it, 0) for it in items}
            for s in sample_order
        }
    else:
        raise ValueError(f"Unknown category kind: {category['kind']!r}")
    df = pd.DataFrame(data).reindex(items).fillna(0).astype(int)
    return df


# --------------------------------------------------------------------------- #
# Plotting
# --------------------------------------------------------------------------- #

def plot_class_counts(class_df: pd.DataFrame, outpath: Path):
    samples = class_df.index.tolist()
    fig, ax = plt.subplots(figsize=(max(8.0, 0.7 * len(samples) + 4), 6))

    bottom = np.zeros(len(samples))
    for cls in CAZYME_CLASSES:
        vals = class_df.get(cls, pd.Series(0, index=samples)).values
        ax.bar(samples, vals, bottom=bottom, label=cls, color=CLASS_COLORS[cls])
        bottom += vals

    ax.set_ylabel("Number of CAZyme domains")
    ax.set_title("CAZyme class counts per sample")
    ax.legend(title="Class", bbox_to_anchor=(1.02, 1), loc="upper left")
    plt.xticks(rotation=45, ha="right")
    fig.tight_layout()
    fig.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close(fig)


def _grouped_bar(ax, df: pd.DataFrame, title: str, ylabel: str, show_legend: bool = True):
    """df: rows = category, cols = sample. Plots grouped bars (one group per row)."""
    if df.empty or df.values.sum() == 0:
        ax.text(0.5, 0.5, "No matches", ha="center", va="center",
                transform=ax.transAxes)
        ax.set_title(title)
        ax.set_xticks([])
        ax.set_yticks([])
        return
    df.plot(kind="bar", ax=ax, width=0.8)
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    ax.set_xlabel("")
    rotation = 0 if df.shape[0] <= 4 else 45
    ax.tick_params(axis="x", rotation=rotation)
    if rotation:
        for lbl in ax.get_xticklabels():
            lbl.set_ha("right")
    if show_legend:
        ax.legend(title="Sample", bbox_to_anchor=(1.02, 1), loc="upper left",
                  fontsize=8)
    else:
        leg = ax.get_legend()
        if leg is not None:
            leg.remove()


def plot_substrate_panels(category_dfs: dict, outpath: Path, ncols: int = 3):
    """category_dfs: ordered dict of {category_name: DataFrame(rows=item, cols=sample)}.

    Lays the panels out on a grid (ncols per row), shows the legend only on the
    first panel of each row.
    """
    if not category_dfs:
        return
    n = len(category_dfs)
    nrows = math.ceil(n / ncols)
    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(7 * ncols, 5.5 * nrows),
        squeeze=False,
    )
    flat_axes = [ax for row in axes for ax in row]

    for idx, (name, df) in enumerate(category_dfs.items()):
        ax = flat_axes[idx]
        # legend only on first panel of each row to reduce clutter
        show_legend = (idx % ncols == 0)
        # find ylabel from category config (lookup by name)
        ylabel = next(
            (c["ylabel"] for c in SUBSTRATE_CATEGORIES if c["name"] == name),
            "# hits",
        )
        _grouped_bar(ax, df, name, ylabel, show_legend=show_legend)

    # turn off any unused axes
    for ax in flat_axes[n:]:
        ax.set_visible(False)

    fig.suptitle("CAZymes by substrate category", fontsize=14, y=1.0)
    fig.tight_layout()
    fig.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #

def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Plot CAZyme summaries from a dbCAN output directory.",
    )
    parser.add_argument(
        "dbcan_dir",
        help="Root dbCAN output directory; one subfolder per sample, "
             "each containing overview.tsv.",
    )
    parser.add_argument(
        "-o", "--outdir", default="dbcan_plots",
        help="Output directory for PNG plots and CSV summaries.",
    )
    parser.add_argument(
        "--top-core", type=int, default=20,
        help="How many top core families to show in plot 3 (default 20).",
    )
    parser.add_argument(
        "--ncols", type=int, default=3,
        help="Columns in the substrate-category panel grid (default 3).",
    )
    args = parser.parse_args(argv)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    samples = find_overview_files(Path(args.dbcan_dir))
    if not samples:
        sys.exit(
            f"No overview.tsv files found under {args.dbcan_dir}/<sample>/overview.tsv"
        )
    print(f"Found {len(samples)} samples: {', '.join(samples.keys())}")

    class_counts, family_counts, substrate_counts = {}, {}, {}
    for sample, path in samples.items():
        df = load_sample(path)
        cls_c, fam_c, sub_c = count_per_sample(df)
        class_counts[sample] = cls_c
        family_counts[sample] = fam_c
        substrate_counts[sample] = sub_c

    sample_order = list(samples.keys())

    # --- Wide matrices ---------------------------------------------------------
    class_df = (
        pd.DataFrame(class_counts).T
        .reindex(index=sample_order)
        .fillna(0).astype(int)
        .reindex(columns=CAZYME_CLASSES, fill_value=0)
    )
    family_df = (
        pd.DataFrame(family_counts).T
        .reindex(index=sample_order)
        .fillna(0).astype(int)
    )

    # Per-category matrices (rows = item, cols = sample)
    category_dfs = {}
    for cat in SUBSTRATE_CATEGORIES:
        cdf = category_matrix(cat, family_counts, substrate_counts, sample_order)
        category_dfs[cat["name"]] = cdf
        cdf.to_csv(outdir / cat["csv"])

    # --- Top-level CSV outputs ------------------------------------------------
    class_df.to_csv(outdir / "cazyme_class_counts.csv")
    family_df.to_csv(outdir / "cazyme_family_counts.csv")
    pd.DataFrame(substrate_counts).fillna(0).astype(int).to_csv(
        outdir / "substrate_annotation_counts.csv"
    )

    # --- Plots -----------------------------------------------------------------
    plot_class_counts(class_df, outdir / "01_cazyme_class_counts.png")
    plot_substrate_panels(
        category_dfs, outdir / "02_substrate_categories.png", ncols=args.ncols,
    )

if __name__ == "__main__":
    main()
