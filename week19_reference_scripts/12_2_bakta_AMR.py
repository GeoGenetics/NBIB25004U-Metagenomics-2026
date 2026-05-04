#!/usr/bin/env python3
"""
bakta_amr_pipeline_advanced.py
=============================

AMR pipeline with:
✔ GFF3 extraction
✔ AMR class + mechanism classification
✔ Seaborn plots (heatmap + barplots)

Usage:
    python bakta_amr_pipeline_advanced.py /maps/projects/course_1/scratch/<group#>/<group-project-group-#>/09_annotation_bakta_ref -o /maps/projects/course_1/scratch/<group#>/<group-project-group-#>/amr_plots
"""

from __future__ import annotations
import argparse
import re
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# -----------------------------------------------------------------------------
# Classification
# -----------------------------------------------------------------------------

def classify_amr(product: str, gene: str):
    p = product.lower()
    g = gene.lower()

    # ❌ remove false positives
    if "inhibitor" in p or "pepsy" in p:
        return None, None

    # ❌ drop low-confidence domain-only hits
    if "beta-lactamase" in p and "domain" in p:
        return None, None

    # ✅ HIGH confidence
    if g == "bla":
        return "beta-lactam", "enzymatic"

    # ✅ MEDIUM-HIGH
    if p == "beta-lactamase":
        return "beta-lactam", "enzymatic"

    if "class a beta-lactamase" in p:
        return "beta-lactam", "enzymatic"

    if "beta-lactamase" in p:
        return "beta-lactam", "enzymatic"

    # tetracycline
    if "tetracycline" in p or g.startswith("tet"):
        return "tetracycline", "efflux" if "efflux" in p else "target_protection"

    # macrolide
    if "macrolide" in p or g.startswith("erm"):
        return "macrolide", "target_modification"

    # aminoglycoside
    if re.search(r"aac|aph|ant", g):
        return "aminoglycoside", "enzymatic"

    return None, None

# -----------------------------------------------------------------------------
# Extraction
# -----------------------------------------------------------------------------

def extract_amr_from_gff(genome: str, gff_path: Path) -> pd.DataFrame:
    rows = []

    with gff_path.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue

            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue

            attr = parts[8]

            product_match = re.search(r"product=([^;]+)", attr)
            if not product_match:
                continue

            product = product_match.group(1).replace("%2C", ",")

            gene_match = re.search(r"gene=([^;]+)", attr)
            gene = gene_match.group(1) if gene_match else ""

            locus_match = re.search(r"locus_tag=([^;]+)", attr)
            locus = locus_match.group(1) if locus_match else ""

            amr_class, mechanism = classify_amr(product, gene)

            if amr_class is None:
                continue

            rows.append({
                "genome": genome,
                "locus_tag": locus,
                "gene": gene,
                "product": product,
                "amr_class": amr_class,
                "mechanism": mechanism
            })

    return pd.DataFrame(rows)

# -----------------------------------------------------------------------------
# Plotting (Seaborn)
# -----------------------------------------------------------------------------

def plot_class_bar(df: pd.DataFrame, outpath: Path):
    counts = df.groupby(["genome", "amr_class"]).size().reset_index(name="count")
    plt.figure(figsize=(8, 5))
    sns.barplot(data=counts, x="genome", y="count", hue="amr_class")
    plt.title("AMR class distribution")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(outpath)
    plt.close()

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("bakta_dir", type=Path)
    ap.add_argument("-o", "--outdir", type=Path, default=Path("amr_out"))
    args = ap.parse_args()

    args.outdir.mkdir(exist_ok=True)
    (args.outdir / "per_genome").mkdir(exist_ok=True)

    all_frames = []

    for folder in sorted(args.bakta_dir.iterdir()):
        if not folder.is_dir():
            continue

        genome = folder.name
        gff = folder / f"{genome}.gff3"

        if not gff.exists():
            continue

        df = extract_amr_from_gff(genome, gff)
        df.to_csv(args.outdir / "per_genome" / f"{genome}_amr.tsv", sep="\t", index=False)

        print(genome, len(df), "AMR genes")

        if not df.empty:
            all_frames.append(df)

    if not all_frames:
        print("No AMR detected")
        return

    long_df = pd.concat(all_frames).drop_duplicates()
    long_df.to_csv(args.outdir / "all_amr.tsv", sep="\t", index=False)

    # presence absence
    pa = (
        long_df.assign(val=1)
        .pivot_table(index="locus_tag", columns="genome", values="val", fill_value=0)
    )
    pa.to_csv(args.outdir / "presence_absence.tsv", sep="\t")

    # plots
    plot_class_bar(long_df, args.outdir / "amr_class_barplot.png")

    print("Done")


if __name__ == "__main__":
    main()
