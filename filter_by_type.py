import os
import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u

RAW_DATA_DIR = "Raw data"
OBJECT_TABLES = ["Objects Tables/Quasars.xlsx", "Objects Tables/Radio_Galaxies.xlsx"]
QSO_DIR = "QSOtables"
AGA_DIR = "AGAtables"

BAD_FLAGS = {"SUN", "SSN"}
TYPE_QSO = ["QSO"]
TYPE_AGA = ["QSO", "G", "GPair", "GTrpl", "GGroup", "GClstr"]

os.makedirs(RAW_DATA_DIR, exist_ok=True)
os.makedirs(QSO_DIR, exist_ok=True)
os.makedirs(AGA_DIR, exist_ok=True)


def clean(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return df
    if "Redshift Flag" not in df.columns:
        return df

    return df[
        df["Redshift Flag"].astype(str).str.startswith("S") &
        ~df["Redshift Flag"].astype(str).str.contains("?", regex=False) &
        ~df["Redshift Flag"].isin(BAD_FLAGS)
    ]


def load_object_table() -> pd.DataFrame:
    frames = []
    for table_path in OBJECT_TABLES:
        if not os.path.exists(table_path):
            print(f"Warning: object table not found: {table_path}")
            continue
        df = pd.read_excel(table_path)
        if {"Name", "Ra", "Dec"}.issubset(df.columns):
            df = df[["Name", "Ra", "Dec"]].copy()
            frames.append(df)
        else:
            print(f"Warning: missing columns in {table_path}: {df.columns.tolist()}")

    if not frames:
        return pd.DataFrame(columns=["Name", "Ra", "Dec"])

    all_objects = pd.concat(frames, ignore_index=True)
    all_objects["Ra_str"] = all_objects["Ra"].astype(str)
    all_objects["Dec_str"] = all_objects["Dec"].astype(str)
    return all_objects


def parse_filename_to_coord(filename: str):
    base = os.path.splitext(filename)[0]
    if "_" not in base:
        raise ValueError(f"Filename does not contain coordinate separator '_': {filename}")

    ra_str, dec_str = base.split("_", 1)
    try:
        ra = float(ra_str)
        dec = float(dec_str)
    except ValueError as exc:
        raise ValueError(f"Invalid coordinate values in filename {filename}: {exc}") from exc

    return ra, dec


def find_object_name(ra_deg: float, dec_deg: float, objects_df: pd.DataFrame, tolerance_deg: float = 0.25):
    if objects_df.empty:
        return None

    target = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg, frame="icrs")
    coords = SkyCoord(
        objects_df["Ra_str"].values,
        objects_df["Dec_str"].values,
        unit=(u.hourangle, u.deg), frame="icrs"
    )
    separations = coords.separation(target).deg
    best_idx = separations.argmin()
    best_sep = separations[best_idx]

    if best_sep <= tolerance_deg:
        return objects_df.iloc[best_idx]["Name"]
    return None


def process_file(input_path: str, objects_df: pd.DataFrame):
    filename = os.path.basename(input_path)
    try:
        ra_deg, dec_deg = parse_filename_to_coord(filename)
    except ValueError as exc:
        print(f"Skipping file with invalid name: {filename} ({exc})")
        return

    object_name = find_object_name(ra_deg, dec_deg, objects_df)
    if object_name is None:
        object_name = os.path.splitext(filename)[0]
        print(f"Object name not found in tables for {filename}; using fallback name '{object_name}'")
    else:
        print(f"Processing {filename} as object '{object_name}'")

    try:
        df = pd.read_excel(input_path)
    except Exception as exc:
        print(f"Failed to read {filename}: {exc}")
        return

    if df.empty:
        print(f"Skipping empty file: {filename}")
        return

    qso_df = df[df["Type"].isin(TYPE_QSO)]
    qso_df = clean(qso_df)

    aga_df = df[df["Type"].isin(TYPE_AGA)]
    aga_df = clean(aga_df)

    qso_output = os.path.join(QSO_DIR, f"QSO {object_name}.xlsx")
    aga_output = os.path.join(AGA_DIR, f"AGA {object_name}.xlsx")

    qso_df.to_excel(qso_output, index=False)
    aga_df.to_excel(aga_output, index=False)

    print(f"Saved: {qso_output} ({len(qso_df)} rows)")
    print(f"Saved: {aga_output} ({len(aga_df)} rows)")


def main():
    objects_df = load_object_table()
    raw_files = sorted(
        f for f in os.listdir(RAW_DATA_DIR)
        if f.lower().endswith(".xlsx")
    )

    if not raw_files:
        print(f"No Excel files found in {RAW_DATA_DIR}")
        return

    for raw_file in raw_files:
        process_file(os.path.join(RAW_DATA_DIR, raw_file), objects_df)


if __name__ == "__main__":
    main()
