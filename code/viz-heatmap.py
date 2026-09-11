from pathlib import Path
import matplotlib.pyplot as plt, numpy as np, pandas as pd
from scipy.stats import spearmanr


PRJ_DIR    = Path(__file__).parent
OUTPUT_DIR = PRJ_DIR / "output"

ENV_FILE = OUTPUT_DIR / "neon-environment-by-day.csv"
MOS_FILE = OUTPUT_DIR / "neon-mos-abundance-by-event-multiple-mo.csv"

START_DATE = "2018-01-01"
END_DATE   = "2024-12-31"


def requireColumns(data: pd.DataFrame, requiredColumns: list[str], dataName: str):
    missingColumns = [ column for column in requiredColumns
                       if column not in data.columns ]
    if missingColumns:
        raise ValueError(f"Missing required columns in {dataName}: {missingColumns}")


def readInputData():
    envData = pd.read_csv(ENV_FILE)
    mosData = pd.read_csv(MOS_FILE)

    requireColumns(envData, ["date", "tempMean", "rhMean", "precipBulk"],
                   "environment data")
    requireColumns(mosData, ["eventID", "eventStart", "completePlots", "abundance24h"],
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
    firstDay: int,
    lastDay: int,
    summaryType: str
):
    # Example:
    # firstDay = 1, lastDay = 7
    # eventStart = July 22
    # window = July 15 through July 21.

    startDate = eventStart - pd.Timedelta(days=lastDay)
    endDate   = eventStart - pd.Timedelta(days=firstDay)

    windowDates = pd.date_range(startDate, endDate, freq="D")
    windowData = (envData.set_index("date").reindex(windowDates)[columnName] )
    expectedDays = lastDay - firstDay + 1

    # Require a valid environmental value for every day in the window.
    if windowData.notna().sum() != expectedDays:
        return np.nan

    if summaryType == "mean": return windowData.mean()
    if summaryType == "sum":  return windowData.sum()

    raise ValueError(f"Unknown summary type: {summaryType}")


def addLagVariables(envData: pd.DataFrame, mosData: pd.DataFrame):
    mosData = mosData.copy()

    lagWindows = [
        (1, 7, "Lag1to7d"),
        (8, 14, "Lag8to14d"),
        (15, 21, "Lag15to21d")
    ]

    for firstDay, lastDay, lagName in lagWindows:

        mosData[f"tempMean{lagName}"] = mosData["eventStart"].apply(
            lambda eventStart:
                getLagValue(
                    envData,
                    eventStart,
                    "tempMean",
                    firstDay,
                    lastDay,
                    "mean"
                )
        )

        mosData[f"precipSum{lagName}"] = mosData["eventStart"].apply(
            lambda eventStart:
                getLagValue(
                    envData,
                    eventStart,
                    "precipBulk",
                    firstDay,
                    lastDay,
                    "sum"
                )
        )

        mosData[f"rhMean{lagName}"] = mosData["eventStart"].apply(
            lambda eventStart:
                getLagValue(
                    envData,
                    eventStart,
                    "rhMean",
                    firstDay,
                    lastDay,
                    "mean"
                )
        )

    return mosData


def getAnalysisData(mosData: pd.DataFrame):
    lagColumns = [
        "tempMeanLag1to7d",
        "tempMeanLag8to14d",
        "tempMeanLag15to21d",

        "precipSumLag1to7d",
        "precipSumLag8to14d",
        "precipSumLag15to21d",

        "rhMeanLag1to7d",
        "rhMeanLag8to14d",
        "rhMeanLag15to21d"]

    # Use exactly the same mosquito events for all nine correlations.
    analysisData = mosData.dropna( subset=["abundance24h"] + lagColumns ).copy()

    return analysisData


def getCorrelationData(analysisData: pd.DataFrame):
    predictorColumns = {
        "Temperature":       ["tempMeanLag1to7d", "tempMeanLag8to14d", "tempMeanLag15to21d"],
        "Precipitation":     ["precipSumLag1to7d", "precipSumLag8to14d", "precipSumLag15to21d"],
        "Relative humidity": ["rhMeanLag1to7d", "rhMeanLag8to14d", "rhMeanLag15to21d"] }

    correlationData = pd.DataFrame(
        index=predictorColumns.keys(),
        columns=["1–7 days", "8–14 days", "15–21 days"], 
        dtype=float)

    for predictorName, columns in predictorColumns.items():
        for lagLabel, columnName in zip(correlationData.columns, columns):
            rho, _ = spearmanr(
                analysisData[columnName],
                analysisData["abundance24h"]
            )

            correlationData.loc[
                predictorName,
                lagLabel
            ] = rho

    return correlationData


def makeFigure(correlationData: pd.DataFrame, sampleSize: int):
    fig, ax = plt.subplots( figsize=(8.5, 4.8) )

    image = ax.imshow(
        correlationData.to_numpy(),
        cmap="coolwarm",
        vmin=-1,
        vmax=1,
        aspect="auto")

    ax.set_xticks(np.arange(len(correlationData.columns)))
    ax.set_xticklabels(correlationData.columns)

    ax.set_yticks(np.arange(len(correlationData.index)))
    ax.set_yticklabels(correlationData.index)
    
    ax.set_xlabel("Days before mosquito sampling")

    # Display Spearman rho in each cell.
    for rowIndex in range(len(correlationData.index)):
        for columnIndex in range(len(correlationData.columns)):

            rho = correlationData.iloc[ rowIndex, columnIndex ]
            textColor = ( "white" if abs(rho) >= 0.5 else "black" )
            ax.text(columnIndex,
                    rowIndex,
                    f"{rho:.2f}",
                    ha="center",
                    va="center",
                    fontsize=12,
                    color=textColor)

    colorBar = fig.colorbar(image, ax=ax, pad=0.03)
    colorBar.set_label("Spearman ρ")

#     ax.set_title(
#         "Lagged Associations Between Mosquito Abundance and "
#         "Environmental Conditions at HARV, 2017–2024\n"
#         f"Spearman correlations; n = {sampleSize}"
#     )

    fig.tight_layout()

    return fig


if __name__ == "__main__":
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    envDf, mosDf = readInputData()
    mosDf        = addLagVariables(envDf, mosDf)

    analysisDf    = getAnalysisData(mosDf)
    correlationDf = getCorrelationData(analysisDf)

    print(f"Mosquito events from 2017–2024: {len(mosDf)}")
    print(f"Events with complete data for all nine lag variables: {len(analysisDf)}")
    print()
    print(correlationDf)

    correlationDf.to_csv(OUTPUT_DIR / "neon-mosquito-environment-lag-correlations.csv")

    fig = makeFigure(correlationDf, len(analysisDf))
    fig.savefig(
        OUTPUT_DIR / "neon-mosquito-environment-lag-heatmap.png",
        dpi=300,
        bbox_inches="tight")

    plt.show()
