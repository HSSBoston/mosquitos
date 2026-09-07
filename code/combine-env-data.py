# Combine daily NEON temperature, relative humidity, and precipitation data.

from pathlib import Path
import pandas as pd

START_DATE = "2017-01-01"
END_DATE   = "2024-12-31"

PRJ_DIR    = Path(__file__).parent
OUTPUT_DIR = PRJ_DIR / "output"

TEMP_RH_FILE = OUTPUT_DIR / "neon-rel-humidity-by-day-multiple-mo.csv"
PRECIP_FILE  = OUTPUT_DIR / "neon-precip-by-day-multiple-mo.csv"


def getDailyEnvironmentData():
    tempRhData = pd.read_csv(TEMP_RH_FILE)
    precipData = pd.read_csv(PRECIP_FILE)

    # Parse dates.
    tempRhData["date"] = pd.to_datetime(tempRhData["date"])
    precipData["date"] = pd.to_datetime(precipData["date"])

    # Stop if either input file contains more than one row for the same day.
    if tempRhData["date"].duplicated().any():
        raise ValueError("Duplicate dates found in temperature/humidity data.")

    if precipData["date"].duplicated().any():
        raise ValueError("Duplicate dates found in precipitation data.")

    # Keep only variables needed for the daily environmental table.
    tempRhData = tempRhData.loc[ :, ["date", "tempMean", "rhMean"] ].copy()
    precipData = precipData.loc[ :, ["date", "precipBulk"] ].copy()

    # Preserve all dates available from either dataset.
    dailyData = tempRhData.merge(
        precipData,
        on="date",
        how="outer") # Keep every "date" value in either DataFrame

    dailyData = dailyData.sort_values("date").reset_index(drop=True)

    # Restrict to the period with START_DATE and END_DATE
    dailyData = dailyData.loc[
        (dailyData["date"] >= START_DATE)
        & (dailyData["date"] <= END_DATE)
    ].copy()

    dailyData["date"] = dailyData["date"].dt.date

    return dailyData


if __name__ == "__main__":
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    dailyDf = getDailyEnvironmentData()
    
    print(dailyDf.to_string(index=False))
    dailyDf.to_csv( OUTPUT_DIR / "neon-environment-by-day.csv", index=False)