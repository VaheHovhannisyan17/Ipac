import pandas as pd
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
import os

# ── Configuration ─────────────────────────────────────────────────────────────

# Path to your existing 6-degree Excel file
INPUT_EXCEL = "QSOtables/QSO 1555+00.xlsx"   # ← change to your file name/path

# Central coordinate (provide as strings or decimal degrees)
CENTRAL_RA_HMS  = "15h57m51.4339s"   # ← change if needed
CENTRAL_DEC_DMS = "-00d01m50.413s"   # ← change if needed

# Column names in your Excel file
RA_COL  = "RA"
DEC_COL = "DEC"

# Output directory (leave as "" to save in the same folder as this script)
OUTPUT_DIR = ""

# ── Helpers ───────────────────────────────────────────────────────────────────

def angular_distance_deg(ra, dec, ra0, dec0):
    """Vectorised great-circle distance in degrees."""
    ra_rad   = np.radians(ra)
    dec_rad  = np.radians(dec)
    ra0_rad  = np.radians(ra0)
    dec0_rad = np.radians(dec0)

    cos_d = (
        np.sin(dec_rad) * np.sin(dec0_rad) +
        np.cos(dec_rad) * np.cos(dec0_rad) * np.cos(ra_rad - ra0_rad)
    )
    cos_d = np.clip(cos_d, -1, 1)
    return np.degrees(np.arccos(cos_d))


def save_excel(dataframe, path):
    """Save a DataFrame to Excel, flushing to disk."""
    with open(path, "wb") as f:
        with pd.ExcelWriter(f, engine="openpyxl") as writer:
            dataframe.to_excel(writer, index=False)
            writer.book.save(f)
            f.flush()
            os.fsync(f.fileno())

# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    # Parse central coordinate
    central = SkyCoord(CENTRAL_RA_HMS, CENTRAL_DEC_DMS)
    ra0  = central.ra.deg
    dec0 = central.dec.deg
    print(f"Central coordinate: RA={ra0:.4f}°  DEC={dec0:.4f}°")

    # Load the 6-degree catalogue
    print(f"Loading {INPUT_EXCEL} …")
    try:
        df = pd.read_csv(INPUT_EXCEL)
    except Exception:
        df = pd.read_excel(INPUT_EXCEL, engine="openpyxl")
    print(f"  {len(df)} objects loaded.")

    # Compute angular distance for every object
    df["_dist_deg"] = angular_distance_deg(
        df[RA_COL].values, df[DEC_COL].values, ra0, dec0
    )

    base = os.path.splitext(os.path.basename(INPUT_EXCEL))[0]
    out_dir = OUTPUT_DIR or os.path.dirname(os.path.abspath(INPUT_EXCEL))

    # Save one file per radius bin
    for radius in range(1, 6):          # 1, 2, 3, 4, 5
        subset = df[df["_dist_deg"] <= radius].drop(columns=["_dist_deg"])
        out_path = os.path.join(out_dir, f"{base}_{radius}deg.xlsx")
        save_excel(subset, out_path)
        print(f"  ≤{radius}°  →  {len(subset):>5} objects  →  {out_path}")

    print("Done.")


if __name__ == "__main__":
    main()