import os
import sys
import time

import numpy as np
import pandas as pd
from astroquery.ipac.ned import Ned
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.table import vstack

Ned.TIMEOUT = 1200
Ned.ROW_LIMIT = -1

RAW_DATA_DIR = "Raw data"
CHECKPOINT_FILE = "checkpoint.txt"
OBJECT_TABLES = ["Objects Tables/Quasars.xlsx", "Objects Tables/Radio_Galaxies.xlsx"]
OBJECT_NAME = "0422+00"
os.makedirs(RAW_DATA_DIR, exist_ok=True)


def load_object_tables():
    frames = []
    for table_path in OBJECT_TABLES:
        if not os.path.exists(table_path):
            print(f"Warning: object table not found: {table_path}")
            continue
        table = pd.read_excel(table_path)
        if {"Name", "Ra", "Dec"}.issubset(table.columns):
            frames.append(table[["Name", "Ra", "Dec"]].copy())
        else:
            print(f"Warning: missing Name/Ra/Dec columns in {table_path}: {table.columns.tolist()}")

    if not frames:
        return pd.DataFrame(columns=["Name", "Ra", "Dec"])
    return pd.concat(frames, ignore_index=True)


def find_object_coord(object_name: str, objects_df: pd.DataFrame):
    if objects_df.empty:
        raise ValueError("No object table data available.")

    normalized = objects_df["Name"].astype(str).str.strip().str.casefold()
    target_key = object_name.strip().casefold()
    matches = objects_df[normalized == target_key]
    if not matches.empty:
        row = matches.iloc[0]
        return SkyCoord(row["Ra"], row["Dec"], unit=(u.hourangle, u.deg), frame="icrs")

    raise ValueError(f"Object '{object_name}' not found in object tables.")


# ── Resume logic ──────────────────────────────────────────────────────────────
try:
    with open(CHECKPOINT_FILE, "r") as f:
        saved = int(f.read().strip())
except Exception:
    saved = 0

notSaved = 0

object_name = OBJECT_NAME
if len(sys.argv) > 1:
    object_name = sys.argv[1]

objects_df = load_object_tables()
centralCoord = find_object_coord(object_name, objects_df)
centralRA = centralCoord.ra.deg
centralDEC = centralCoord.dec.deg
DBname = os.path.join(RAW_DATA_DIR, f"{centralRA:.2f}_{centralDEC:.2f}.xlsx")

try:
    df = pd.read_excel(DBname, engine="openpyxl")
except Exception:
    df = pd.DataFrame()

print(f"Using object: {object_name} -> RA={centralRA:.6f}, DEC={centralDEC:.6f}")

# ── Helpers ───────────────────────────────────────────────────────────────────
def filterByAngularDistance(df, ra0=centralRA, dec0=centralDEC,
                             max_distance_deg=6, raCol='RA', decCol='DEC'):
    ra  = np.radians(df[raCol].values)
    dec = np.radians(df[decCol].values)
    ra0_rad  = np.radians(ra0)
    dec0_rad = np.radians(dec0)

    cos_d = (
        np.sin(dec) * np.sin(dec0_rad) +
        np.cos(dec) * np.cos(dec0_rad) * np.cos(ra - ra0_rad)
    )
    cos_d = np.clip(cos_d, -1, 1)
    distances = np.degrees(np.arccos(cos_d))

    return df[distances <= max_distance_deg]


def save_checkpoint(value: int):
    with open(CHECKPOINT_FILE, "w") as f:
        f.write(str(value))


def save_df():
    global df, DBname
    os.makedirs(os.path.dirname(DBname), exist_ok=True)
    df.to_excel(DBname, index=False, engine="openpyxl")


def retry_call(func, *args, attempts=10, delay=10, **kwargs):
    last_exc = None
    for attempt in range(1, attempts + 1):
        try:
            return func(*args, **kwargs)
        except Exception as exc:
            last_exc = exc
            print(f"Attempt {attempt}/{attempts} failed: {type(exc).__name__}: {exc}")
            if attempt < attempts:
                time.sleep(delay)
    raise last_exc


def get60Arcmin(coord):
    ra  = coord.ra.deg
    dec = coord.dec.deg
    radius = 60 * u.arcmin
    results = []

    for raOffset in [-radius/1.5, 0, radius/1.5]:
        for decOffset in [-radius/1.5, 0, radius/1.5]:
            patch_center = SkyCoord(
                coord.ra  + raOffset,
                coord.dec + decOffset,
                frame=coord.frame
            )
            patch_radius = radius / 3 * 1.415

            for attempt in range(3):
                try:
                    result = Ned.query_region(patch_center, radius=patch_radius)
                    results.append(result)
                    print(f"✓ Success: Patch at {patch_center.ra.deg:.2f}, {patch_center.dec.deg:.2f}")
                    break
                except Exception as e:
                    print(f"  Patch attempt {attempt + 1}/3 failed at {patch_center.ra.deg:.2f}, {patch_center.dec.deg:.2f}: {type(e).__name__}: {e}")
                    if attempt < 2:
                        time.sleep(5)
                    else:
                        print(f"  Skipping patch after 3 failed attempts: {patch_center.ra.deg:.2f}, {patch_center.dec.deg:.2f}")

    if not results:
        return pd.DataFrame()

    for t in results:
        t.meta.clear()
    try:
        combined = vstack(results)
    except Exception as exc:
        print(f"Failed to combine patch results: {exc}")
        return pd.DataFrame()

    try:
        result_df = combined.to_pandas()
    except Exception as exc:
        print(f"Failed to convert patch results to pandas DataFrame: {exc}")
        return pd.DataFrame()

    if not {"RA", "DEC", "Redshift Flag", "Object Name"}.issubset(result_df.columns):
        print(f"Missing required columns in query result: {result_df.columns.tolist()}")
        return pd.DataFrame()

    result_df = result_df[result_df["Redshift Flag"] != ""]
    true_df = filterByAngularDistance(result_df, ra0=ra, dec0=dec, max_distance_deg=1)
    return true_df.drop_duplicates(subset=["Object Name"])


def get6degree(coord):
    global saved, df, DBname, notSaved

    ra  = coord.ra.deg
    dec = coord.dec.deg

    # Build the full list of patch offsets (arcmin) in one place
    offsets = []

    for raOff in range(-240, 241, 120):
        for decOff in range(-240, 241, 120):
            offsets.append((raOff, decOff))

    for raOff in range(-180, 181, 120):
        for decOff in range(-180, 181, 120):
            offsets.append((raOff, decOff))

    for raOff in [-300, 300]:
        for decOff in range(-180, 181, 120):
            offsets.append((raOff, decOff))

    for decOff in [-300, 300]:
        for raOff in range(-180, 181, 120):
            offsets.append((raOff, decOff))

    for raOff in [-360, 360]:
        for decOff in [-120, 0, 120]:
            offsets.append((raOff, decOff))

    for decOff in [-360, 360]:
        for raOff in [-120, 0, 120]:
            offsets.append((raOff, decOff))

    for raOff, decOff in offsets:
        notSaved += 1
        if notSaved <= saved:
            continue

        patchCenter = SkyCoord(
            coord.ra  + raOff  * u.arcmin,
            coord.dec + decOff * u.arcmin,
            frame=coord.frame
        )

        try:
            result = retry_call(get60Arcmin, patchCenter, attempts=10, delay=10)
        except Exception as exc:
            print(f"Failed after retries for patch {notSaved} at {patchCenter.ra.deg:.2f}, {patchCenter.dec.deg:.2f}: {type(exc).__name__}: {exc}")
            result = pd.DataFrame()

        try:
            true_df = filterByAngularDistance(result, ra0=ra, dec0=dec, max_distance_deg=6)
        except Exception as exc:
            print(f"Distance filtering failed for patch {notSaved}: {type(exc).__name__}: {exc}")
            true_df = pd.DataFrame()

        df = pd.concat([df, true_df], axis=0, ignore_index=True)
        if "Object Name" in df.columns:
            df = df.drop_duplicates(subset=["Object Name"])

        save_df()
        print(f"{notSaved} done — {len(df)} objects so far")

    print("All patches complete.")
    save_df()
    save_checkpoint(0)
    print("Checkpoint reset to 0.")
    print(f"Data saved: {DBname}")

try:
    get6degree(centralCoord)
except KeyboardInterrupt:
    print(f"\nInterrupted at patch {notSaved}")
    save_checkpoint(max(0, notSaved - 1))
    save_df()
    print(f"Checkpoint saved: {max(0, notSaved - 1)}")
    print(f"Data saved: {DBname}")
except Exception as e:
    print(f"\nCrashed at patch {notSaved}: {type(e).__name__}: {e}")
    save_checkpoint(max(0, notSaved - 1))
    save_df()
    print(f"Checkpoint saved: {max(0, notSaved - 1)}")
    print(f"Data saved: {DBname}")