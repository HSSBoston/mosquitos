# HARV soil-array sensor: HOR.VER = 003.000, recorded height = 0.65 m.
# Select only 30-minute data; each CSV contains both RH and air temperature.

from pathlib import Path
import pandas as pd

PRJ_DIR    = Path(__file__).parent
DATA_DIR   = PRJ_DIR / "data" / "NEON_rel-humidity"
OUTPUT_DIR = PRJ_DIR / "output"

# False: Include nonmissing readings even if a potential QC issue is present.
# True:  Use only readings with a passing quality flag, separately for RH and temperature.
EXCLUDE_QC_ISSUES = True

# A complete UTC day contains 48 half-hour intervals.
EXPECTED_INTERVALS = 48


# Helper function to read and combine matching files in a single directory.
#
def readNeonFiles(path: Path, pattern: str):
    files = sorted(path.glob(pattern))
    if not files: raise FileNotFoundError(f"No files found in {path} matching: {pattern}")

    data = pd.concat( [pd.read_csv(file) for file in files],
                      ignore_index=True )
    return data


def getDailySummary(path:Path):
    weatherData = readNeonFiles(
        path, "NEON.D01.HARV.DP1.00098.001.003.000.030.RH_30min.*.csv")

    if weatherData.empty: raise ValueError("The matching files contain no observations.")

    # Parse timestamps as UTC. Assign each interval to its start date.
    weatherData["startDateTime"] = pd.to_datetime(weatherData["startDateTime"], utc=True)
    weatherData["endDateTime"]   = pd.to_datetime(weatherData["endDateTime"],   utc=True)

    invalidTime = ( weatherData["startDateTime"].isna()
                  | weatherData["endDateTime"].isna()
                  | (weatherData["endDateTime"] - weatherData["startDateTime"])
                    .ne(pd.Timedelta(minutes=30)) )
    if invalidTime.any(): raise ValueError("Missing timestamps or intervals other than 30 minutes.")

    # Remove identical records. Stop if different records describe the same interval,
    # for example when multiple versions of a monthly file are in this directory.
    weatherData = weatherData.drop_duplicates()
    if weatherData["startDateTime"].duplicated().any():
        raise ValueError("Conflicting duplicate intervals found. Keep one version of each input file.")

    # Exclude suspect readings separately so a failed RH reading does not remove
    # a valid temperature reading from the same interval (or vice versa).
    if EXCLUDE_QC_ISSUES:
        for valueColumn, flagColumn in [("RHMean", "RHFinalQF"),
                                        ("tempRHMean", "tempRHFinalQF")]:
            qcPass = weatherData[flagColumn].eq(0)

            # Expanded files may also contain a science-review flag.
            # Blank or 0: no unresolved science-review failure. 1 or 2: fail.
            reviewColumn = flagColumn + "SciRvw"
            if reviewColumn in weatherData.columns:
                qcPass = qcPass & ( weatherData[reviewColumn].isna()
                                 | weatherData[reviewColumn].eq(0) )

            weatherData.loc[~qcPass, valueColumn] = float("nan")

    weatherData["date"] = weatherData["startDateTime"].dt.floor("D")

    # Average the available, accepted half-hour means with equal weight.
    # mean() ignores missing values; count() counts only contributing readings.
    summaryData = weatherData.groupby("date").agg(
        tempMean      =("tempRHMean", "mean"),
        rhMean        =("RHMean", "mean"),
        tempIntervals =("tempRHMean", "count"),
        rhIntervals   =("RHMean", "count") )

    # Retain days with no records between the first and last dates in the data.
    allDates = pd.date_range(weatherData["date"].min(), weatherData["date"].max(),
                            freq="D", name="date")
    summaryData = summaryData.reindex(allDates)
    for columnName in ["tempIntervals", "rhIntervals"]:
        summaryData[columnName] = summaryData[columnName].fillna(0).astype(int)

    # Coverage is the percentage of expected half-hour intervals contributing
    # to each mean. No minimum daily coverage threshold is imposed here.
    summaryData["tempCoverage"] = summaryData["tempIntervals"] / EXPECTED_INTERVALS * 100
    summaryData["rhCoverage"]   = summaryData["rhIntervals"] / EXPECTED_INTERVALS * 100

    # Format the output table: temperature in degrees C; RH and coverage in percent.
    summaryData = summaryData.reset_index()
    summaryData["date"] = summaryData["date"].dt.date
    summaryData = summaryData.round({"tempMean": 3, "rhMean": 2,
                                    "tempCoverage": 1, "rhCoverage": 1})
    return summaryData


if __name__ == "__main__":
    targetDir = DATA_DIR / "NEON.D01.HARV.DP1.00098.001.2024-07.expanded.20260123T000749Z.RELEASE-2026"

    summaryDf = getDailySummary(targetDir)

    print( summaryDf.to_string(index=False) )
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    summaryDf.to_csv(OUTPUT_DIR / "neon-rel-humidity-by-day.csv", index=False)
