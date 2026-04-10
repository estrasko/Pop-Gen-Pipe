#!/usr/bin/env python3

from __future__ import annotations

import argparse
import logging
import subprocess
import sys
from pathlib import Path

import pandas as pd


def setup_logging(log_file: Path) -> None:
    """Configure logging to both console and file."""
    log_file.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout),
        ],
    )


def validate_file(path: Path, label: str) -> None:
    """Raise a clear error if a required file is missing."""
    if not path.exists():
        raise FileNotFoundError(f"{label} not found: {path}")


def make_output_dirs(outdir: Path) -> dict[str, Path]:
    """Create standard output directories and return them."""
    dirs = {
        "base": outdir,
        "logs": outdir / "logs",
        "amova": outdir / "amova",
        "dapc": outdir / "dapc",
        "ibd": outdir / "ibd",
        "divmigrate": outdir / "divmigrate",
    }
    for directory in dirs.values():
        directory.mkdir(parents=True, exist_ok=True)
    return dirs


def extract_genepop_sample_names(genepop_path: Path) -> list[str]:
    """
    Extract sample names from a GENEPOP file.

    Assumes sample lines contain a comma, where the sample name is the text
    before the first comma.
    """
    sample_names: list[str] = []

    with genepop_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()

            if not line:
                continue

            # Skip POP lines
            if line.upper() == "POP":
                continue

            # In genepop, individual/sample lines contain a comma
            if "," in line:
                sample_name = line.split(",", 1)[0].strip()
                sample_names.append(sample_name)

    if not sample_names:
        raise ValueError(
            f"No sample names were extracted from the genepop file: {genepop_path}"
        )

    return sample_names


def validate_and_fix_popmap(
    genepop_path: Path,
    popmap_path: Path,
    outdir: Path,
    label: str,
) -> Path:
    """
    Validate that the popmap Sample column matches the genepop sample names,
    reorder the popmap to match the genepop order, and write a corrected copy.

    Returns the path to the corrected popmap CSV.
    """
    logging.info("Validating popmap against %s genepop file: %s", label, genepop_path)

    genepop_samples = extract_genepop_sample_names(genepop_path)
    popmap_df = pd.read_csv(popmap_path)

    if "Sample" not in popmap_df.columns:
        raise ValueError(
            "Popmap must contain a 'Sample' column for sample-name validation."
        )

    if "Population" not in popmap_df.columns:
        raise ValueError(
            "Popmap must contain a 'Population' column."
        )

    popmap_samples = popmap_df["Sample"].astype(str).tolist()

    genepop_set = set(genepop_samples)
    popmap_set = set(popmap_samples)

    if genepop_set != popmap_set:
        missing_in_popmap = sorted(genepop_set - popmap_set)
        missing_in_genepop = sorted(popmap_set - genepop_set)

        raise ValueError(
            "Sample mismatch detected between genepop and popmap.\n"
            f"Missing in popmap: {missing_in_popmap}\n"
            f"Missing in genepop: {missing_in_genepop}"
        )

    # Reorder popmap to match genepop sample order
    corrected_df = (
        popmap_df.assign(Sample=popmap_df["Sample"].astype(str))
        .set_index("Sample")
        .loc[genepop_samples]
        .reset_index()
    )

    corrected_path = outdir / "popmap_corrected.csv"
    corrected_df.to_csv(corrected_path, index=False)

    logging.info("Wrote corrected popmap to: %s", corrected_path)

    return corrected_path


def run_r_script(script_path: Path, args: list[str]) -> None:
    """Run an R script with subprocess and fail loudly on error."""
    cmd = ["Rscript", str(script_path)] + args
    logging.info("Running command: %s", " ".join(cmd))
    subprocess.run(cmd, check=True)


def run_amova(haps_genepop: Path, popmap: Path, outdir: Path, scripts_dir: Path) -> None:
    """Run AMOVA using the haplotype genepop file."""
    run_r_script(
        scripts_dir / "run_amova.R",
        [str(haps_genepop), str(popmap), str(outdir)],
    )


def run_dapc(multi_snp_genepop: Path, popmap: Path, outdir: Path, scripts_dir: Path) -> None:
    """Run DAPC using the multi-SNP-per-locus genepop file."""
    run_r_script(
        scripts_dir / "run_dapc.R",
        [str(multi_snp_genepop), str(popmap), str(outdir)],
    )


def run_ibd(
    fst_csv: Path,
    geographic_distance_csv: Path,
    outdir: Path,
    scripts_dir: Path,
    summary_stats_csv: Path | None = None,
) -> None:
    """Run IBD analyses from Fst and geographic distance matrices."""
    args = [str(fst_csv), str(geographic_distance_csv), str(outdir)]
    if summary_stats_csv is not None:
        args.append(str(summary_stats_csv))
    run_r_script(scripts_dir / "run_ibd.R", args)


def run_divmigrate(
    multi_snp_genepop: Path,
    outdir: Path,
    scripts_dir: Path,
    stat: str,
    boots: int,
    node_names: str | None = None,
) -> None:
    """Run divMigrate on the multi-SNP-per-locus genepop file."""
    args = [str(multi_snp_genepop), str(outdir), stat, str(boots)]
    if node_names is not None:
        args.append(node_names)
    run_r_script(scripts_dir / "run_divmigrate.R", args)


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the pipeline."""
    parser = argparse.ArgumentParser(
        description="Population genetics pipeline connecting AMOVA, DAPC, IBD, and divMigrate."
    )

    parser.add_argument(
        "--haps-genepop",
        required=True,
        help="Path to haplotype genepop file for AMOVA.",
    )
    parser.add_argument(
        "--multi-snp-genepop",
        required=True,
        help="Path to multiple-SNP-per-locus genepop file for DAPC and divMigrate.",
    )
    parser.add_argument(
        "--popmap",
        required=True,
        help="Path to popmap CSV with required Population column and strongly recommended Sample column for name validation.",
    )
    parser.add_argument(
        "--fst-csv",
        required=True,
        help="Path to pairwise Fst matrix CSV for IBD.",
    )
    parser.add_argument(
        "--geo-csv",
        required=True,
        help="Path to pairwise geographic distance matrix CSV for IBD.",
    )
    parser.add_argument(
        "--summary-stats-csv",
        help="Optional summary statistics CSV for downstream regressions in IBD.",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory for all results.",
    )
    parser.add_argument(
        "--scripts-dir",
        default="scripts",
        help="Directory containing the R scripts. Default: scripts",
    )

    parser.add_argument("--run-amova", action="store_true", help="Run AMOVA.")
    parser.add_argument("--run-dapc", action="store_true", help="Run DAPC.")
    parser.add_argument("--run-ibd", action="store_true", help="Run IBD.")
    parser.add_argument("--run-divmigrate", action="store_true", help="Run divMigrate.")

    parser.add_argument(
        "--divmigrate-stat",
        default="gst",
        choices=["gst", "D", "Nm"],
        help="Statistic to use for divMigrate. Default: gst",
    )
    parser.add_argument(
        "--divmigrate-boots",
        type=int,
        default=1000,
        help="Bootstrap replicates for divMigrate. Default: 1000",
    )
    parser.add_argument(
        "--divmigrate-node-names",
        help="Optional comma-separated node names for divMigrate plots.",
    )

    return parser.parse_args()


def main() -> int:
    """Run the requested population genetics analyses."""
    args = parse_args()

    haps_genepop = Path(args.haps_genepop)
    multi_snp_genepop = Path(args.multi_snp_genepop)
    original_popmap = Path(args.popmap)
    fst_csv = Path(args.fst_csv)
    geo_csv = Path(args.geo_csv)
    summary_stats_csv = Path(args.summary_stats_csv) if args.summary_stats_csv else None
    outdir = Path(args.outdir)
    scripts_dir = Path(args.scripts_dir)

    dirs = make_output_dirs(outdir)
    setup_logging(dirs["logs"] / "pipeline.log")

    logging.info("Starting population genetics pipeline.")
    logging.info("Output directory: %s", outdir)

    validate_file(haps_genepop, "Haplotype genepop file")
    validate_file(multi_snp_genepop, "Multiple SNP genepop file")
    validate_file(original_popmap, "Popmap file")
    validate_file(fst_csv, "Fst matrix CSV")
    validate_file(geo_csv, "Geographic distance matrix CSV")
    validate_file(scripts_dir / "run_amova.R", "AMOVA R script")
    validate_file(scripts_dir / "run_dapc.R", "DAPC R script")
    validate_file(scripts_dir / "run_ibd.R", "IBD R script")
    validate_file(scripts_dir / "run_divmigrate.R", "divMigrate R script")

    if summary_stats_csv is not None:
        validate_file(summary_stats_csv, "Summary statistics CSV")

    if not any([args.run_amova, args.run_dapc, args.run_ibd, args.run_divmigrate]):
        raise ValueError(
            "No analyses selected. Add at least one of: "
            "--run-amova, --run-dapc, --run-ibd, --run-divmigrate"
        )

    # Validate popmap against the multi-SNP genepop first and create corrected popmap
    corrected_popmap = validate_and_fix_popmap(
        genepop_path=multi_snp_genepop,
        popmap_path=original_popmap,
        outdir=outdir,
        label="multi-SNP",
    )

    # Validate the corrected popmap against the haplotype genepop too
    corrected_popmap = validate_and_fix_popmap(
        genepop_path=haps_genepop,
        popmap_path=corrected_popmap,
        outdir=outdir,
        label="haplotype",
    )

    try:
        if args.run_amova:
            logging.info("Running AMOVA...")
            run_amova(haps_genepop, corrected_popmap, dirs["amova"], scripts_dir)

        if args.run_dapc:
            logging.info("Running DAPC...")
            run_dapc(multi_snp_genepop, corrected_popmap, dirs["dapc"], scripts_dir)

        if args.run_ibd:
            logging.info("Running IBD...")
            run_ibd(fst_csv, geo_csv, dirs["ibd"], scripts_dir, summary_stats_csv)

        if args.run_divmigrate:
            logging.info("Running divMigrate...")
            run_divmigrate(
                multi_snp_genepop,
                dirs["divmigrate"],
                scripts_dir,
                args.divmigrate_stat,
                args.divmigrate_boots,
                args.divmigrate_node_names,
            )

    except subprocess.CalledProcessError as error:
        logging.exception("Analysis failed while running an R script.")
        return error.returncode
    except Exception:
        logging.exception("Pipeline failed.")
        return 1

    logging.info("Pipeline completed successfully.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
