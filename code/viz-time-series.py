from pathlib import Path
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import pandas as pd

START_DATE = "2017-01-01"
END_DATE   = "2024-12-31"

TEMP = True
PRECIP = True
HUMIDITY = True

SHOW_TEMP_ROLLING_MEAN     = True
SHOW_HUMIDITY_ROLLING_MEAN = True
ROLLING_WINDOW             = 7

USE_LOG_ABUNDANCE     = False
SHADE_ALTERNATE_YEARS = False

PRJ_DIR    = Path(__file__).parent
OUTPUT_DIR = PRJ_DIR / "output"

ENV_FILE = OUTPUT_DIR / "neon-environment-by-day.csv"
MOS_FILE = OUTPUT_DIR / "neon-mos-abundance-by-event-multiple-mo.csv"


def requireColumns(data: pd.DataFrame, requiredColumns: list[str], dataName: str):
    missingColumns = [
        column for column in requiredColumns
        if column not in data.columns
    ]

    if missingColumns:
        raise ValueError(
            f"Missing required columns in {dataName}: {missingColumns}"
        )

def readInputData():
    envData = pd.read_csv(ENV_FILE)
    mosData = pd.read_csv(MOS_FILE)

    requiredEnvColumns = ["date"]

    if TEMP:
        requiredEnvColumns.append("tempMean")

    if HUMIDITY:
        requiredEnvColumns.append("rhMean")

    if PRECIP:
        requiredEnvColumns.append("precipBulk")

    requireColumns(
        envData,
        requiredEnvColumns,
        "environment data"
    )

    requireColumns(
        mosData,
        ["eventID", "eventStart", "completePlots", "abundance24h"],
        "mosquito abundance data"
    )

    envData["date"] = pd.to_datetime(envData["date"])
    mosData["eventStart"] = pd.to_datetime(mosData["eventStart"])

    envData = envData.loc[
        (envData["date"] >= START_DATE)
        & (envData["date"] <= END_DATE)
    ].copy()

    mosData = mosData.loc[
        (mosData["eventStart"] >= START_DATE)
        & (mosData["eventStart"] <= END_DATE)
    ].copy()

    envData = envData.sort_values("date").reset_index(drop=True)
    mosData = mosData.sort_values("eventStart").reset_index(drop=True)

    mosData["year"] = mosData["eventStart"].dt.year

    if TEMP and SHOW_TEMP_ROLLING_MEAN:
        envData["tempMeanRolling"] = (
            envData["tempMean"]
            .rolling(ROLLING_WINDOW, min_periods=1)
            .mean()
        )

    if HUMIDITY and SHOW_HUMIDITY_ROLLING_MEAN:
        envData["rhMeanRolling"] = (
            envData["rhMean"]
            .rolling(ROLLING_WINDOW, min_periods=1)
            .mean()
        )

    return envData, mosData


def addYearGuides(ax, startYear: int, endYear: int):
    if SHADE_ALTERNATE_YEARS:
        for year in range(startYear, endYear + 1):
            if (year - startYear) % 2 == 1:
                ax.axvspan(
                    pd.Timestamp(f"{year}-01-01"),
                    pd.Timestamp(f"{year + 1}-01-01"),
                    alpha=0.06
                )

    for year in range(startYear + 1, endYear + 1):
        ax.axvline(
            pd.Timestamp(f"{year}-01-01"),
            linewidth=0.8,
            alpha=0.45
        )


def makeFigure(envData: pd.DataFrame, mosData: pd.DataFrame):
    startYear = pd.Timestamp(START_DATE).year
    endYear   = pd.Timestamp(END_DATE).year

    panelCount = 1 + TEMP + HUMIDITY + PRECIP

    fig, axes = plt.subplots(
        panelCount,
        1,
        figsize=(14, 2.3 * panelCount + 1.5),
        sharex=True,
        gridspec_kw={"hspace": 0.08}
    )

    if panelCount == 1:
        axes = [axes]
    else:
        axes = list(axes)

    panelIndex = 0

    # ------------------------------------------------------------
    # Panel A: Mosquito abundance
    # ------------------------------------------------------------

    axMos = axes[panelIndex]
    panelIndex += 1

    for _, yearData in mosData.groupby("year"):
        yearData = yearData.sort_values("eventStart")

        axMos.plot(
            yearData["eventStart"],
            yearData["abundance24h"],
            linewidth=1.0,
            alpha=0.7,
            color="0.55"
        )

    axMos.scatter(
        mosData["eventStart"],
        mosData["abundance24h"],
        s=26,
        zorder=3,
        color="0.10"
    )

    if USE_LOG_ABUNDANCE:
        axMos.set_yscale("symlog", linthresh=10)

    axMos.set_ylabel(
        "Mosquito\nabundance24h"
    )

    axMos.set_title(
        "Seasonal Mosquito Abundance and Meteorological Conditions "
        "at HARV, 2017–2024"
    )

    # ------------------------------------------------------------
    # Temperature
    # ------------------------------------------------------------

    if TEMP:
        axTemp = axes[panelIndex]
        panelIndex += 1

        validTemp = envData["tempMean"].notna()

        axTemp.plot(
            envData.loc[validTemp, "date"],
            envData.loc[validTemp, "tempMean"],
            linewidth=0.8,
            alpha=0.45,
            color="0.55",
            label="Daily mean"
        )

        if SHOW_TEMP_ROLLING_MEAN:
            validRolling = envData["tempMeanRolling"].notna()

            axTemp.plot(
                envData.loc[validRolling, "date"],
                envData.loc[validRolling, "tempMeanRolling"],
                linewidth=1.5,
                color="0.10",
                label=f"{ROLLING_WINDOW}-day mean"
            )

            axTemp.legend(
                loc="upper left",
                frameon=False
            )

        axTemp.set_ylabel(
            "Temp\n(°C)"
        )

    # ------------------------------------------------------------
    # Relative humidity
    # ------------------------------------------------------------

    if HUMIDITY:
        axHumidity = axes[panelIndex]
        panelIndex += 1

        validHumidity = envData["rhMean"].notna()

        axHumidity.plot(
            envData.loc[validHumidity, "date"],
            envData.loc[validHumidity, "rhMean"],
            linewidth=0.8,
            alpha=0.45,
            color="0.55",
            label="Daily mean"
        )

        if SHOW_HUMIDITY_ROLLING_MEAN:
            validRolling = envData["rhMeanRolling"].notna()

            axHumidity.plot(
                envData.loc[validRolling, "date"],
                envData.loc[validRolling, "rhMeanRolling"],
                linewidth=1.5,
                color="0.10",
                label=f"{ROLLING_WINDOW}-day mean"
            )

            axHumidity.legend(
                loc="upper left",
                frameon=False
            )

        axHumidity.set_ylabel(
            "Relative\nhumidity (%)"
        )

    # ------------------------------------------------------------
    # Precipitation
    # ------------------------------------------------------------

    if PRECIP:
        axPrecip = axes[panelIndex]
        panelIndex += 1

        validPrecip = envData["precipBulk"].notna()

        axPrecip.bar(
            envData.loc[validPrecip, "date"],
            envData.loc[validPrecip, "precipBulk"],
            width=1.0,
            align="center",
            color="0.55",
            alpha=0.9
        )

        axPrecip.set_ylabel(
            "Precip\n(mm)"
        )

    # ------------------------------------------------------------
    # Shared formatting
    # ------------------------------------------------------------

    for ax in axes:
        addYearGuides(
            ax,
            startYear,
            endYear
        )

        ax.set_xlim(
            pd.Timestamp(START_DATE),
            pd.Timestamp(END_DATE)
        )

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    axes[-1].set_xlabel("Date")

    axes[-1].xaxis.set_major_locator(
        mdates.YearLocator()
    )

    axes[-1].xaxis.set_major_formatter(
        mdates.DateFormatter("%Y")
    )

    fig.align_ylabels(axes)
    fig.tight_layout()

    return fig


if __name__ == "__main__":
    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True
    )

    envDf, mosDf = readInputData()

    fig = makeFigure(
        envDf,
        mosDf
    )

    fig.savefig(
        OUTPUT_DIR / "neon-seasonal-mosquito-weather-timeseries.png",
        dpi=300,
        bbox_inches="tight"
    )

    fig.savefig(
        OUTPUT_DIR / "neon-seasonal-mosquito-weather-timeseries.pdf",
        bbox_inches="tight"
    )

    plt.show()
