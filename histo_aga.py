#!/usr/bin/env python3
"""Generate redshift histograms for all Excel files in AGAtables."""

from __future__ import annotations

import argparse
import pathlib

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def save_histograms(input_dir: pathlib.Path, output_dir: pathlib.Path) -> int:
    input_dir = input_dir.expanduser().resolve()
    output_dir = output_dir.expanduser().resolve()
    if not input_dir.is_dir():
        raise FileNotFoundError(f"Input folder not found: {input_dir}")

    output_dir.mkdir(parents=True, exist_ok=True)
    files = sorted(input_dir.glob("*.xlsx"))
    if not files:
        print(f"No .xlsx files found in {input_dir}")
        return 0

    count = 0
    for file_path in files:
        try:
            df = pd.read_excel(file_path)
        except Exception as exc:
            print(f"Failed to read {file_path.name}: {exc}")
            continue

        if "Redshift" not in df.columns:
            print(f"Skipping {file_path.name}: no Redshift column")
            continue

        plt.figure(figsize=(10, 6))
        sns.histplot(df, x="Redshift", binwidth=0.2, color="#0077CC", alpha=0.6)
        plt.title(f"Redshift distribution for {file_path.stem}")
        plt.xlabel("Redshift")
        plt.ylabel("Count")

        output_file = output_dir / f"{file_path.stem}_redshift_histogram.png"
        plt.tight_layout()
        plt.savefig(output_file)
        plt.close()
        print(f"Saved histogram to {output_file}")
        count += 1

    return count


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate redshift histogram PNG files for all Excel tables in AGAtables.",
    )
    parser.add_argument(
        "input_dir",
        nargs="?",
        default="AGAtables",
        help="Path to the AGAtables folder. Defaults to AGAtables.",
    )
    parser.add_argument(
        "--output-dir",
        default="Histo/AGA",
        help="Path to save histogram PNG files. Defaults to Histo/AGA.",
    )
    args = parser.parse_args()

    total = save_histograms(pathlib.Path(args.input_dir), pathlib.Path(args.output_dir))
    if total:
        print(f"Done: saved {total} histograms to {args.output_dir}")
    else:
        print("Done: no histograms were saved.")


if __name__ == "__main__":
    main()
