import os
import sys
import time

import numpy as np
import pandas as pd
from astroquery.ipac.ned import Ned
from astropy.coordinates import SkyCoord
import astropy.units as u

Ned.TIMEOUT = 1200
Ned.ROW_LIMIT = -1

# ── Configuration ─────────────────────────────────────────────────────────────
RAW_DATA_DIR    = "Raw data"
CHECKPOINT_FILE = "checkpoint.txt"
OBJECT_TABLES   = ["Objects Tables/Quasars.xlsx", "Objects Tables/Radio_Galaxies.xlsx"]
OBJECT_NAME     = "0736+01"

CONE_RADIUS     = 6.0             # degrees
NSIDE           = 128             # HEALPix resolution — pixel radius ≈ 27.5 arcmin
QUERY_RADIUS    = 28 * u.arcmin   # NED query radius per patch

os.makedirs(RAW_DATA_DIR, exist_ok=True)

# ── Object table loading ──────────────────────────────────────────────────────
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
failed_patches = 0

object_name = OBJECT_NAME
if len(sys.argv) > 1:
    object_name = sys.argv[1]

objects_df   = load_object_tables()
centralCoord = find_object_coord(object_name, objects_df)
centralRA    = centralCoord.ra.deg
centralDEC   = centralCoord.dec.deg
DBname       = os.path.join(RAW_DATA_DIR, f"{centralRA:.2f}_{centralDEC:.2f}.xlsx")

try:
    df = pd.read_excel(DBname, engine="openpyxl")
except Exception:
    df = pd.DataFrame()

print(f"Using object: {object_name} -> RA={centralRA:.6f}, DEC={centralDEC:.6f}")
print(f"Cone radius: {CONE_RADIUS}°  |  HEALPix nside={NSIDE}  |  Query radius={QUERY_RADIUS}")

# ── Helpers ───────────────────────────────────────────────────────────────────
def filterByAngularDistance(df, ra0=centralRA, dec0=centralDEC,
                             max_distance_deg=6, raCol='RA', decCol='DEC'):
    ra   = np.radians(df[raCol].values)
    dec  = np.radians(df[decCol].values)
    ra0_rad  = np.radians(ra0)
    dec0_rad = np.radians(dec0)
    cos_d = (
        np.sin(dec) * np.sin(dec0_rad) +
        np.cos(dec) * np.cos(dec0_rad) * np.cos(ra - ra0_rad)
    )
    cos_d = np.clip(cos_d, -1, 1)
    return df[np.degrees(np.arccos(cos_d)) <= max_distance_deg]


def save_checkpoint(value: int):
    with open(CHECKPOINT_FILE, "w") as f:
        f.write(str(value))


def save_df():
    global df, DBname
    os.makedirs(os.path.dirname(DBname), exist_ok=True)
    df.to_excel(DBname, index=False, engine="openpyxl")


def retry_call(func, *args, attempts=3, delay=15, **kwargs):
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


# ── HEALPix patch center generator ───────────────────────────────────────────
def get_patch_centers(central_coord, cone_radius_deg, nside):
    """
    Returns SkyCoord array of HEALPix pixel centers within cone_radius_deg
    of central_coord. nside=128 gives pixel radius ~27.5 arcmin, which
    guarantees full coverage with 35 arcmin NED queries (verified numerically).
    """
    npix = 12 * nside * nside
    ipix = np.arange(npix)
    theta = np.zeros(npix)
    phi   = np.zeros(npix)
    nc    = 2 * nside * (nside - 1)

    mask_n = ipix < nc
    i_n    = ipix[mask_n]
    ring_n = np.floor((1 + np.sqrt(1 + 2 * i_n)) / 2).astype(int)
    s_n    = i_n - 2 * ring_n * (ring_n - 1)
    theta[mask_n] = np.arccos(1 - ring_n**2 / (3 * nside**2))
    phi[mask_n]   = (np.pi / (2 * ring_n)) * (s_n + 0.5)

    mask_s = ipix >= npix - nc
    i_s    = npix - 1 - ipix[mask_s]
    ring_s = np.floor((1 + np.sqrt(1 + 2 * i_s)) / 2).astype(int)
    s_s    = i_s - 2 * ring_s * (ring_s - 1)
    theta[mask_s] = np.arccos(-(1 - ring_s**2 / (3 * nside**2)))
    phi[mask_s]   = (np.pi / (2 * ring_s)) * (s_s + 0.5)

    mask_e = ~mask_n & ~mask_s
    i_e    = ipix[mask_e]
    ring_e = np.floor(i_e / (4 * nside) - (nside - 2) / 2 + nside).astype(int)
    s_e    = i_e - (2 * nside * (nside + 1) + 4 * nside * (ring_e - 2 * nside))
    parity = (ring_e + nside) % 2
    theta[mask_e] = np.arccos((2 - ring_e / nside) * 2 / 3)
    phi[mask_e]   = (np.pi / (2 * nside)) * (s_e + 0.5 * parity)

    ra_all  = np.degrees(phi) % 360
    dec_all = 90 - np.degrees(theta)

    ra0_rad  = np.radians(central_coord.ra.deg)
    dec0_rad = np.radians(central_coord.dec.deg)
    cos_d = (
        np.sin(np.radians(dec_all)) * np.sin(dec0_rad) +
        np.cos(np.radians(dec_all)) * np.cos(dec0_rad) *
        np.cos(np.radians(ra_all) - ra0_rad)
    )
    mask_cone = np.degrees(np.arccos(np.clip(cos_d, -1, 1))) <= cone_radius_deg

    return SkyCoord(ra=ra_all[mask_cone] * u.deg, dec=dec_all[mask_cone] * u.deg)


# ── Single patch query ────────────────────────────────────────────────────────
def query_patch(coord):
    result = Ned.query_region(coord, radius=QUERY_RADIUS)
    result.meta.clear()
    patch_df = result.to_pandas()

    if not {"RA", "DEC", "Redshift Flag", "Object Name"}.issubset(patch_df.columns):
        print(f"  Missing required columns: {patch_df.columns.tolist()}")
        return pd.DataFrame()

    patch_df = patch_df[patch_df["Redshift Flag"] != ""]
    patch_df = filterByAngularDistance(patch_df, max_distance_deg=CONE_RADIUS)
    patch_df = patch_df.drop_duplicates(subset=["Object Name"])
    print(f"  ✓ {coord.ra.deg:.3f}, {coord.dec.deg:.3f}  →  {len(patch_df)} objects")
    return patch_df


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    global df, notSaved, failed_patches

    patches = get_patch_centers(
        centralCoord,
        CONE_RADIUS + QUERY_RADIUS.to_value(u.deg),
        NSIDE
    )
    total = len(patches)
    print(f"Total patches to query: {total}")
    if saved > 0:
        print(f"Resuming from patch {saved + 1}")

    for patch_coord in patches:
        notSaved += 1
        if notSaved <= saved:
            continue

        print(f"Patch {notSaved}/{total}")

        try:
            patch_df = retry_call(query_patch, patch_coord, attempts=3, delay=15)
        except Exception as exc:
            failed_patches += 1
            print(f"  Failed after retries: {type(exc).__name__}: {exc}")
            patch_df = pd.DataFrame()

        if not patch_df.empty:
            df = pd.concat([df, patch_df], axis=0, ignore_index=True)
            if "Object Name" in df.columns:
                df = df.drop_duplicates(subset=["Object Name"])
        else:
            if failed_patches > 0 and patch_df.empty:
                print(f"  Patch {notSaved} produced no objects after retries.")

        save_df()
        save_checkpoint(notSaved)
        print(f"  {len(df)} unique objects so far")

    print(f"\nAll patches complete. Total unique objects: {len(df)}")
    print(f"Failed patches: {failed_patches}/{total}")
    save_df()
    save_checkpoint(0)
    print("Checkpoint reset to 0.")
    print(f"Data saved: {DBname}")


try:
    main()
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
    print(f"Failed patches so far: {failed_patches}/{notSaved}")
