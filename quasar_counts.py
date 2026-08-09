#!/usr/bin/env python3
"""Create an Excel summary of QSO objects and quasar counts."""

from __future__ import annotations

import argparse
import os
import pathlib

import pandas as pd

OBJECT_TABLES = ["Objects Tables/Quasars.xlsx", "Objects Tables/Radio_Galaxies.xlsx"]
EXTRAGALACTIC_TYPES = ["QSO", "G", "GPair", "GTrpl", "GGroup", "GClstr"]


def clean_object_name(name: str) -> str:
    if name.startswith("QSO "):
        return name[len("QSO "):]
    if name.startswith("AGA "):
        return name[len("AGA "):]
    return name


def normalize_name(name: str) -> str:
    return "".join(ch for ch in str(name).upper() if ch.isalnum() or ch in "+.")


def get_object_info(object_name: str) -> dict[str, object]:
    info = {"z": None, "b": None, "RA": None, "DEC": None, "Object Type": None}
    normalized_target = normalize_name(object_name)

    for table_path in OBJECT_TABLES:
        if not os.path.exists(table_path):
            continue
        df = pd.read_excel(table_path)
        if {"Name", "z", "b", "Ra", "Dec"}.issubset(df.columns):
            df = df.copy()
            df["_normalized"] = df["Name"].astype(str).apply(normalize_name)

            exact_matches = df[df["_normalized"] == normalized_target]
            if not exact_matches.empty:
                row = exact_matches.iloc[0]
            else:
                partial_matches = df[df["_normalized"].str.contains(normalized_target, na=False)]
                if not partial_matches.empty:
                    row = partial_matches.iloc[0]
                else:
                    continue

            info["z"] = row.get("z")
            info["b"] = row.get("b")
            info["RA"] = row.get("Ra")
            info["DEC"] = row.get("Dec")
            info["Object Type"] = "Quasar" if "Quasars.xlsx" in table_path else "Radio Galaxy"
            return info

    return info


def count_quasars_in_file(file_path: pathlib.Path) -> int:
    df = pd.read_excel(file_path)
    if "Type" in df.columns:
        return int(df["Type"].eq("QSO").sum())
    return len(df)


def count_extragalactic_objects(aga_folder: pathlib.Path, object_name: str) -> int | None:
    aga_name = f"AGA {object_name}.xlsx"
    aga_path = aga_folder / aga_name
    if not aga_path.exists():
        return None

    df = pd.read_excel(aga_path)
    if "Type" in df.columns:
        return int(df["Type"].isin(EXTRAGALACTIC_TYPES).sum())
    return len(df)


def build_summary(
    folder_path: pathlib.Path,
    output_path: pathlib.Path,
    aga_folder: pathlib.Path,
) -> int:
    folder_path = folder_path.expanduser().resolve()
    output_path = output_path.expanduser().resolve()
    aga_folder = aga_folder.expanduser().resolve()

    if not folder_path.is_dir():
        raise FileNotFoundError(f"Folder not found: {folder_path}")

    xlsx_files = sorted(folder_path.glob("*.xlsx"))
    if not xlsx_files:
        print(f"No .xlsx files found in {folder_path}")
        return 0

    rows = []
    for xlsx_file in xlsx_files:
        stem = clean_object_name(xlsx_file.stem)
        try:
            quasar_count = count_quasars_in_file(xlsx_file)
        except Exception as exc:
            print(f"Warning: failed to read {xlsx_file.name}: {exc}")
            quasar_count = 0

        object_info = get_object_info(stem)
        extragalactic_count = count_extragalactic_objects(aga_folder, stem)

        rows.append(
            {
                "Object Name": stem,
                "RA": object_info["RA"],
                "DEC": object_info["DEC"],
                "Object Type": object_info["Object Type"],
                "z": object_info["z"],
                "b": object_info["b"],
                "Quasar Count": quasar_count,
                "Extragalactic Count": extragalactic_count,
            }
        )
        print(f"Processed {xlsx_file.name}: {quasar_count} quasar(s), {extragalactic_count} extragalactic object(s)")

    summary_df = pd.DataFrame(rows)
    summary_df = summary_df.sort_values(["RA"])
    summary_df.to_excel(output_path, index=False, sheet_name="Summary")

    return len(rows)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate an Excel summary of object names and quasar counts from QSO tables."
    )
    parser.add_argument(
        "folder",
        nargs="?",
        default="QSOtables",
        help="Path to the folder containing .xlsx QSO table files.",
    )
    parser.add_argument(
        "--aga-folder",
        default="AGAtables",
        help="Path to the folder containing .xlsx AGA table files.",
    )
    parser.add_argument(
        "--output",
        default="final results.xlsx",
        help="Path to the output Excel summary file.",
    )
    args = parser.parse_args()

    total = build_summary(
        pathlib.Path(args.folder),
        pathlib.Path(args.output),
        pathlib.Path(args.aga_folder),
    )
    if total:
        print(f"Done: wrote summary for {total} object(s) to {args.output}")
    else:
        print("Done: no summary rows written.")


if __name__ == "__main__":
    main()
