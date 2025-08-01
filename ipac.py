from astroquery.ipac.ned import Ned
from astropy.coordinates import SkyCoord
import astropy.units as u
import pandas as pd
from astropy.table import vstack
import time
import os
import numpy as np

Ned.TIMEOUT = 600
Ned.ROW_LIMIT = -1

saved = 0
notSaved = 0

centralCoord = SkyCoord("07h46m07.3937s", "17d42m18.998s")
centralRA = centralCoord.ra.deg
centralDEC = centralCoord.dec.deg
DBname = f"{centralRA:.2f}_{centralDEC:.2f}.xlsx"

try:
    df = pd.read_excel(DBname)
except:
    df = pd.DataFrame()
    df.to_excel(DBname)

def filterByAngularDistance(df, theta0 = centralRA, phi0 = centralDEC, max_distance_deg = 6, thetaCol='RA', phiCol='DEC'):
    theta = np.radians(df[thetaCol].values)
    phi = np.radians(df[phiCol].values)
    theta0Rad = np.radians(theta0)
    phi0Rad = np.radians(phi0)

    cos_d = (
        np.sin(theta) * np.sin(theta0Rad) +
        np.cos(theta) * np.cos(theta0Rad) * np.cos(phi - phi0Rad)
    )
    cos_d = np.clip(cos_d, -1, 1)
    distances = np.degrees(np.arccos(cos_d))

    mask = distances <= max_distance_deg
    return df[mask]

def get60Arcmin(coord):
    ra = coord.ra.deg
    dec = coord.dec.deg
    radius = 60 * u.arcmin
    results = []
    
    for raOffset in [-radius/1.5, 0, radius/1.5]:
        for decOffset in [-radius/1.5, 0, radius/1.5]:
            patch_center = SkyCoord(
                coord.ra + raOffset,
                coord.dec + decOffset,
                frame=coord.frame
            )
            
            patch_radius = radius/3 * 1.415  
            
            for attempt in range(3):
                try:
                    result = Ned.query_region(patch_center, radius=patch_radius)
                    results.append(result)
                    print(f"✓ Success: Patch at {patch_center.ra.deg:.2f}, {patch_center.dec.deg:.2f}")
                    break
                except Exception as e:
                    if attempt == 2:
                        print(f"✗ Failed: Patch at {patch_center.ra.deg:.2f}, {patch_center.dec.deg:.2f} - {str(e)}")
                    time.sleep(5)
    
    for t in results:
        t.meta.clear()
    combined = vstack(results) if results else None
    if combined:
        df = combined.to_pandas()
        df = df[df["Redshift Flag"] != ""]
        true_df = filterByAngularDistance(df, theta0=ra, phi0=dec, max_distance_deg=1)
        return true_df.drop_duplicates(subset=["Object Name"])
    return pd.DataFrame()

def get6degree(coord):
    global saved
    global df
    global DBname
    global notSaved
    ra = coord.ra.deg
    dec = coord.dec.deg
    for raOffset in range(-240, 241, 120):
        for decOffset in range(-240, 241, 120):
            patchCenter = SkyCoord(
                coord.ra + raOffset * u.arcmin,
                coord.dec + decOffset * u.arcmin,
                frame=coord.frame
            )

            notSaved += 1
            if notSaved <= saved:
                continue

            result = get60Arcmin(patchCenter)
            true_df = filterByAngularDistance(result, theta0=ra, phi0=dec, max_distance_deg=6)
            df = pd.concat([df, true_df], axis=0, ignore_index=True)
            df = df.drop_duplicates(subset=["Object Name"])

            print(notSaved, "done")
            print("converted to excel")

    for raOffset in range(-180, 181, 120):
        for decOffset in range(-180, 181, 120):
            patchCenter = SkyCoord(
                coord.ra + raOffset * u.arcmin,
                coord.dec + decOffset * u.arcmin,
                frame=coord.frame
            )

            notSaved += 1
            if notSaved <= saved:
                continue

            result = get60Arcmin(patchCenter)
            true_df = filterByAngularDistance(result, theta0=ra, phi0=dec, max_distance_deg=6)
            df = pd.concat([df, true_df], axis=0, ignore_index=True)
            df = df.drop_duplicates(subset=["Object Name"])

            print(notSaved, "done")
            print("converted to excel")

    for raOffset in [-300, 300]:
        for decOffset in range(-180, 181, 120):
            patchCenter = SkyCoord(
                coord.ra + raOffset * u.arcmin,
                coord.dec + decOffset * u.arcmin,
                frame=coord.frame
            )

            notSaved += 1
            if notSaved <= saved:
                continue

            result = get60Arcmin(patchCenter)
            true_df = filterByAngularDistance(result, theta0=ra, phi0=dec, max_distance_deg=6)
            df = pd.concat([df, true_df], axis=0, ignore_index=True)
            df = df.drop_duplicates(subset=["Object Name"])

            print(notSaved, "done")
            print("converted to excel")

    for decOffset in [-300, 300]:
        for raOffset in range(-180, 181, 120):
            patchCenter = SkyCoord(
                coord.ra + raOffset * u.arcmin,
                coord.dec + decOffset * u.arcmin,
                frame=coord.frame
            )

            notSaved += 1
            if notSaved <= saved:
                continue

            result = get60Arcmin(patchCenter)
            true_df = filterByAngularDistance(result, theta0=ra, phi0=dec, max_distance_deg=6)
            df = pd.concat([df, true_df], axis=0, ignore_index=True)
            df = df.drop_duplicates(subset=["Object Name"])

            print(notSaved, "done")
            print("converted to excel")

    for raOffset in [-360, 360]:
        for decOffset in [-120, 0, 120]:
            patchCenter = SkyCoord(
                coord.ra + raOffset * u.arcmin,
                coord.dec + decOffset * u.arcmin,
                frame=coord.frame
            )

            notSaved += 1
            if notSaved <= saved:
                continue

            result = get60Arcmin(patchCenter)
            true_df = filterByAngularDistance(result, theta0=ra, phi0=dec, max_distance_deg=6)
            df = pd.concat([df, true_df], axis=0, ignore_index=True)
            df = df.drop_duplicates(subset=["Object Name"])

            print(notSaved, "done")
            print("converted to excel")

    for decOffset in [-360, 360]:
        for raOffset in [-120, 0, 120]:
            patchCenter = SkyCoord(
                coord.ra + raOffset * u.arcmin,
                coord.dec + decOffset * u.arcmin,
                frame=coord.frame
            )

            notSaved += 1
            if notSaved <= saved:
                continue

            result = get60Arcmin(patchCenter)
            true_df = filterByAngularDistance(result, theta0=ra, phi0=dec, max_distance_deg=6)
            df = pd.concat([df, true_df], axis=0, ignore_index=True)
            df = df.drop_duplicates(subset=["Object Name"])

            print(notSaved, "done")
            print("converted to excel")

    with open(DBname, "wb") as f:
        with pd.ExcelWriter(f, engine="openpyxl") as writer:
            df.to_excel(writer, index=False)
            writer.book.save(f)
            f.flush()
            os.fsync(f.fileno())

try:
    get6degree(centralCoord)
except KeyboardInterrupt:
    codename = "ipac5.py"
    with open(codename, "r") as f:
        lines = f.readlines()

    lines[12] = f"saced = {notSaved - 1}\n"

    with open(codename, "w") as f:
        f.writelines(lines)

    with open(DBname, "wb") as f:
        with pd.ExcelWriter(f, engine="openpyxl") as writer:
            df.to_excel(writer, index=False)
            writer.book.save(f)
            f.flush()
            os.fsync(f.fileno())