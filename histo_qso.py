import pathlib

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

INPUT_DIR = pathlib.Path("QSOtables")
OUTPUT_DIR = pathlib.Path("Histo")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

for excel_file in sorted(INPUT_DIR.glob("*.xlsx")):
    try:
        df = pd.read_excel(excel_file)
    except Exception as exc:
        print(f"Failed to read {excel_file.name}: {exc}")
        continue

    if "Redshift" not in df.columns:
        print(f"Skipping {excel_file.name}: no Redshift column")
        continue

    plt.figure(figsize=(10, 6))
    sns.histplot(df, x="Redshift", binwidth=0.2, color="#00FF00", alpha=0.6)
    plt.title(f"Redshift distribution for {excel_file.stem}")
    plt.xlabel("Redshift")
    plt.ylabel("Count")

    output_file = OUTPUT_DIR / f"{excel_file.stem}_redshift_histogram.png"
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()
    print(f"Saved histogram to {output_file}")
