# HARV weighing-gauge precipitation.
# Select only daily data; each CSV contains one bulk precipitation value per day.

from pathlib import Path
import pandas as pd

PRJ_DIR    = Path(__file__).parent
DATA_DIR   = PRJ_DIR / "data" / "NEON_precip-weighing"
OUTPUT_DIR = PRJ_DIR / "output"

# False: Include nonmissing readings even if a potential QC issue is present.
# True:  Use only readings with a passing quality flag.
EXCLUDE_QC_ISSUES = True


# Helper function to read and combine matching files in a single directory.
#
def readNeonFiles(path: Path, pattern: str):
    files = sorted(path.glob(pattern))
    if not files: raise FileNotFoundError(f"No files found in {path} matching: {pattern}")

    data = pd.concat( [pd.read_csv(file) for file in files],
                      ignore_index=True )
    return data


def getDailySummary(path: Path):
    precipData = readNeonFiles(
        path, "NEON.D01.HARV.DP1.00044.001.900.000.01D.WEIPRE_daily.*.csv")

    if precipData.empty: raise ValueError("The matching files contain no observations.")

    # Parse the daily date.
    precipData["date"] = pd.to_datetime(precipData["date"])

    if precipData["date"].isna().any():
        raise ValueError("Missing or invalid dates found.")

    # Remove identical records. Stop if different records describe the same day,
    # for example when multiple versions of a monthly file are in this directory.
    precipData = precipData.drop_duplicates()
    if precipData["date"].duplicated().any():
        raise ValueError("Conflicting duplicate days found. Keep one version of each input file.")

    # Exclude suspect daily precipitation values.
    if EXCLUDE_QC_ISSUES:
        qcPass = precipData["finalQF"].eq(0)

        # Expanded files may also contain a science-review flag.
        # Blank or 0: no unresolved science-review failure. 1 or 2: fail.
        if "finalQFSciRvw" in precipData.columns:
            qcPass = qcPass & ( precipData["finalQFSciRvw"].isna()
                              | precipData["finalQFSciRvw"].eq(0) )

        precipData.loc[~qcPass, "precipBulk"] = float("nan")

    # Each input row already represents total precipitation over one day.
    summaryData = precipData.loc[
        :, ["date", "precipBulk"]
    ].copy()

    # Retain days with no records between the first and last dates in the data.
    allDates = pd.date_range(precipData["date"].min(), precipData["date"].max(),
                             freq="D", name="date")
    summaryData = summaryData.set_index("date").reindex(allDates).reset_index()

    # Format the output table: precipitation in mm.
    summaryData["date"] = summaryData["date"].dt.date
    summaryData["precipBulk"] = summaryData["precipBulk"].round(2)

    return summaryData


if __name__ == "__main__":
    targetDir = DATA_DIR / "NEON.D01.HARV.DP1.00044.001.2024-07.expanded.20260123T000749Z.RELEASE-2026"
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    summaryDf = getDailySummary(targetDir)

    print( summaryDf.to_string(index=False) )
    summaryDf.to_csv(OUTPUT_DIR / "neon-precip-by-day.csv", index=False)