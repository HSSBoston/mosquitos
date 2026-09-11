from pathlib import Path
import matplotlib.pyplot as plt, numpy as np, pandas as pd
from scipy.stats import spearmanr

START_DATE = "2017-01-01"
END_DATE   = "2024-12-31"

PRJ_DIR    = Path(__file__).parent
OUTPUT_DIR = PRJ_DIR / "output"

ENV_FILE = OUTPUT_DIR / "neon-environment-by-day.csv"
MOS_FILE = OUTPUT_DIR / "neon-mos-abundance-by-event-multiple-mo.csv"


def requireColumns(data: pd.DataFrame, requiredColumns: list[str], dataName: str):
    missingColumns = [ column for column in requiredColumns
                       if column not in data.columns ]
    if missingColumns:
        raise ValueError(f"Missing required columns in {dataName}: {missingColumns}")


def readInputData():
    envData = pd.read_csv(ENV_FILE)
    mosData = pd.read_csv(MOS_FILE)

    requireColumns(envData,
                   ["date", "tempMean", "rhMean", "precipBulk"],
                   "environment data")

    requireColumns(mosData,
                   ["eventID", "eventStart", "completePlots", "abundance24h"],
                   "mosquito abundance data")

    envData["date"]       = pd.to_datetime(envData["date"])
    mosData["eventStart"] = pd.to_datetime(mosData["eventStart"])

    if envData["date"].duplicated().any():
        raise ValueError("Duplicate dates found in environment data.")

    if mosData["eventID"].duplicated().any():
        raise ValueError("Duplicate event IDs found in mosquito data.")

    envData = envData.loc[
        (envData["date"] >= START_DATE)
        & (envData["date"] <= END_DATE) ].copy()

    mosData = mosData.loc[
        (mosData["eventStart"] >= START_DATE)
        & (mosData["eventStart"] <= END_DATE) ].copy()

    envData = envData.sort_values("date").reset_index(drop=True)
    mosData = mosData.sort_values("eventStart").reset_index(drop=True)

    return envData, mosData


def getLagValue(
    envData: pd.DataFrame,
    eventStart: pd.Timestamp,
    columnName: str,
    days: int,
    summaryType: str):
    # Example for a 7-day window and eventStart = July 9:
    # use July 2 through July 8.
    startDate = eventStart - pd.Timedelta(days=days)
    endDate   = eventStart - pd.Timedelta(days=1)

    windowDates = pd.date_range(startDate, endDate, freq="D")
    windowData = envData.set_index("date").reindex(windowDates)[columnName]

    # Require a valid value for every day in the lag window.
    if windowData.notna().sum() != days: return np.nan
    if summaryType == "mean":            return windowData.mean()
    if summaryType == "sum":             return windowData.sum()

    raise ValueError(f"Unknown summary type: {summaryType}")


def addLagVariables(envData: pd.DataFrame, mosData: pd.DataFrame):
    mosData = mosData.copy()

    mosData["tempMean7d"] = mosData["eventStart"].apply(
        lambda eventStart:
            getLagValue(
                envData,
                eventStart,
                "tempMean",
                7,
                "mean"
            )
    )
    mosData["tempMean14d"] = mosData["eventStart"].apply(
        lambda eventStart:
            getLagValue(
                envData,
                eventStart,
                "tempMean",
                14,
                "mean"
            )
    )
    
    mosData["precipSum7d"] = mosData["eventStart"].apply(
        lambda eventStart:
            getLagValue(
                envData,
                eventStart,
                "precipBulk",
                7,
                "sum"
            )
    )
    mosData["precipSum14d"] = mosData["eventStart"].apply(
        lambda eventStart:
            getLagValue(
                envData,
                eventStart,
                "precipBulk",
                14,
                "sum"
            )
    )

    mosData["rhMean7d"] = mosData["eventStart"].apply(
        lambda eventStart:
            getLagValue(
                envData,
                eventStart,
                "rhMean",
                7,
                "mean"
            )
    )
    mosData["rhMean14d"] = mosData["eventStart"].apply(
        lambda eventStart:
            getLagValue(
                envData,
                eventStart,
                "rhMean",
                14,
                "mean"
            )
    )
    return mosData


def getAnalysisData(mosData: pd.DataFrame):
    # Use the same mosquito events in all three panels.
    analysisData = mosData.dropna(
        subset=["abundance24h",
                "tempMean7d", "tempMean14d",
                "precipSum7d", "precipSum14d",
                "rhMean7d", "rhMean14d"]).copy()

    analysisData["logAbundance"] = np.log10( analysisData["abundance24h"] + 1 )

    return analysisData


def addScatterPanel(
    ax,
    data: pd.DataFrame,
    xColumn: str,
    xLabel: str,
    panelLabel: str
):
    ax.scatter(
        data[xColumn],
        data["logAbundance"],
        s=32,
        alpha=0.7,
        color="0.15")

    # Spearman correlation is calculated using the original abundance
    # values. The log transformation is only for visualization.
    rho, _ = spearmanr(
        data[xColumn],
        data["abundance24h"])

    ax.text(
        0.05,
        0.95,
        f"{panelLabel}\n"
        f"Spearman ρ = {rho:.2f}\n"
        f"n = {len(data)}",
        transform=ax.transAxes,
        ha="left",
        va="top")

    ax.set_xlabel(xLabel)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def makeFigure(analysisData: pd.DataFrame):
    fig, axes = plt.subplots(
        1,
        3,
        figsize=(13, 4.3),
        sharey=True)

    addScatterPanel(
        axes[0],
        analysisData,
        "tempMean7d",
        "7-day mean temperature (°C)",
        "A")

    addScatterPanel(
        axes[1],
        analysisData,
        "precipSum14d",
        "14-day cumulative precipitation (mm)",
        "B")

    addScatterPanel(
        axes[2],
        analysisData,
        "rhMean7d",
        "7-day mean relative humidity (%)",
        "C")

    axes[0].set_ylabel("Mosquito abundance\nlog10(abundance24h + 1)")

    fig.tight_layout( rect=[0, 0, 1, 0.93] )
    return fig


if __name__ == "__main__":
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    envDf, mosDf = readInputData()
    mosDf        = addLagVariables(envDf, mosDf)
    analysisDf   = getAnalysisData(mosDf)

    print(f"Mosquito events from 2017–2024: {len(mosDf)}")
    print(f"Events with complete data for all three panels: {len(analysisDf)}")

    fig = makeFigure(analysisDf)
    fig.savefig(
        OUTPUT_DIR / "neon-mosquito-environment-scatterplots.png",
        dpi=300,
        bbox_inches="tight")

    plt.show()
