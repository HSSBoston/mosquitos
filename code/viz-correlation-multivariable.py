from pathlib import Path
import matplotlib.pyplot as plt, numpy as np, pandas as pd
from scipy.stats import spearmanr
import statsmodels.api as sm


START_DATE = "2017-01-01"
END_DATE   = "2024-12-31"

PRJ_DIR    = Path(__file__).parent
OUTPUT_DIR = PRJ_DIR / "output"

ENV_FILE = OUTPUT_DIR / "neon-environment-by-day.csv"
MOS_FILE = OUTPUT_DIR / "neon-mos-abundance-by-event-multiple-mo.csv"


# Model predictors are listed in coefficient order after the intercept.
#
MODEL_SPECS = [
    ("T",
     ["zTemp"]),

    ("H",
     ["zHumidity"]),

    ("P",
     ["zPrecip"]),

    ("T + H",
     ["zTemp", "zHumidity", "tempHumidity"]),

    ("T + P",
     ["zTemp", "zPrecip", "tempPrecip"]),

    ("T + H + P",
     ["zTemp", "zHumidity", "zPrecip",
      "tempHumidity", "tempPrecip"])
]


TERM_LABELS = {
    "zTemp":        "T",
    "zHumidity":    "H",
    "zPrecip":      "P",
    "tempHumidity": "T×H",
    "tempPrecip":   "T×P"
}


def requireColumns(data: pd.DataFrame, requiredColumns: list[str], dataName: str):
    missingColumns = [ column for column in requiredColumns
                       if column not in data.columns ]

    if missingColumns:
        raise ValueError(
            f"Missing required columns in {dataName}: {missingColumns}")


def readInputData():
    envData = pd.read_csv(ENV_FILE)
    mosData = pd.read_csv(MOS_FILE)

    requireColumns(
        envData,
        ["date", "tempMean", "rhMean", "precipBulk"],
        "environment data")

    requireColumns(
        mosData,
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
    if windowData.notna().sum() != days:
        return np.nan

    if summaryType == "mean": return windowData.mean()
    if summaryType == "sum":  return windowData.sum()

    raise ValueError(f"Unknown summary type: {summaryType}")


def addLagVariables(envData: pd.DataFrame, mosData: pd.DataFrame):
    mosData = mosData.copy()

    # T: preceding 7-day mean temperature.
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

    # H: preceding 7-day mean relative humidity.
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

    # P: preceding 14-day cumulative precipitation.
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

    return mosData


def standardizeColumn(data: pd.Series):
    standardDeviation = data.std()

    if pd.isna(standardDeviation) or standardDeviation == 0:
        raise ValueError(
            f"Cannot standardize {data.name}: "
            "standard deviation is zero or missing.")

    return (data - data.mean()) / standardDeviation


def getAnalysisData(mosData: pd.DataFrame):
    # Use exactly the same mosquito events for all six models.
    analysisData = mosData.dropna(
        subset=[
            "abundance24h",
            "tempMean7d",
            "rhMean7d",
            "precipSum14d"
        ]
    ).copy()

    # Response transformation.
    analysisData["logAbundance"] = np.log10(
        analysisData["abundance24h"] + 1
    )

    # Precipitation transformation before standardization.
    if (analysisData["precipSum14d"] < 0).any():
        raise ValueError("Negative precipitation values found.")

    analysisData["logPrecip"] = np.log1p(
        analysisData["precipSum14d"]
    )

    # Standardize the three environmental predictors.
    analysisData["zTemp"] = standardizeColumn(
        analysisData["tempMean7d"]
    )

    analysisData["zHumidity"] = standardizeColumn(
        analysisData["rhMean7d"]
    )

    analysisData["zPrecip"] = standardizeColumn(
        analysisData["logPrecip"]
    )

    # Interaction terms are formed after standardization.
    analysisData["tempHumidity"] = (
        analysisData["zTemp"]
        * analysisData["zHumidity"]
    )

    analysisData["tempPrecip"] = (
        analysisData["zTemp"]
        * analysisData["zPrecip"]
    )

    return analysisData


def fitModels(analysisData: pd.DataFrame):
    summaryRows     = []
    coefficientRows = []
    fittedData      = {}

    observedData = analysisData["abundance24h"]
    responseData = analysisData["logAbundance"]

    for modelName, predictorColumns in MODEL_SPECS:

        xData = sm.add_constant(
            analysisData[predictorColumns],
            has_constant="add"
        )

        model = sm.OLS(
            responseData,
            xData
        ).fit()

        fittedValues = model.fittedvalues

        # Spearman correlation between observed mosquito abundance
        # and the model's fitted score.
        rho, _ = spearmanr(
            observedData,
            fittedValues
        )

        rmse = np.sqrt(
            np.mean(
                (responseData - fittedValues) ** 2
            )
        )

        confidenceIntervals = model.conf_int()

        summaryRow = {
            "model":            modelName,
            "terms":            " + ".join(
                                    TERM_LABELS[column]
                                    for column in predictorColumns),
            "n":                len(analysisData),
            "spearmanRho":      rho,
            "rSquared":         model.rsquared,
            "adjustedRSquared": model.rsquared_adj,
            "rmse":             rmse,
            "aic":              model.aic,
            "bic":              model.bic
        }

        # beta0 is the intercept. beta1, beta2, ... follow the order
        # specified in MODEL_SPECS.
        coefficientNames = ["const"] + predictorColumns

        for betaIndex, coefficientName in enumerate(coefficientNames):
            betaName = f"beta{betaIndex}"

            if coefficientName == "const":
                termName = "Intercept"
            else:
                termName = TERM_LABELS[coefficientName]

            summaryRow[f"{betaName}Term"] = termName
            summaryRow[betaName] = model.params[coefficientName]

            summaryRow[f"{betaName}CiLow"] = (
                confidenceIntervals.loc[coefficientName, 0]
            )

            summaryRow[f"{betaName}CiHigh"] = (
                confidenceIntervals.loc[coefficientName, 1]
            )

            summaryRow[f"{betaName}PValue"] = (
                model.pvalues[coefficientName]
            )

            coefficientRows.append({
                "model":       modelName,
                "beta":        betaName,
                "term":        termName,
                "coefficient": model.params[coefficientName],
                "ciLow":       confidenceIntervals.loc[coefficientName, 0],
                "ciHigh":      confidenceIntervals.loc[coefficientName, 1],
                "pValue":      model.pvalues[coefficientName]
            })

        # The largest model has beta0 through beta5.
        # Leave unused coefficient columns blank for smaller models.
        for betaIndex in range(len(coefficientNames), 6):
            betaName = f"beta{betaIndex}"

            summaryRow[f"{betaName}Term"]   = np.nan
            summaryRow[betaName]            = np.nan
            summaryRow[f"{betaName}CiLow"]  = np.nan
            summaryRow[f"{betaName}CiHigh"] = np.nan
            summaryRow[f"{betaName}PValue"] = np.nan

        summaryRows.append(summaryRow)
        fittedData[modelName] = fittedValues

    summaryData = pd.DataFrame(summaryRows)
    coefficientData = pd.DataFrame(coefficientRows)

    return summaryData, coefficientData, fittedData


def makeFigure(
    analysisData: pd.DataFrame,
    summaryData: pd.DataFrame,
    fittedData: dict):

    fig, axes = plt.subplots(
        2,
        3,
        figsize=(13, 8),
        sharex=True,
        sharey=True
    )

    axes = axes.ravel()

    observedLog = analysisData["logAbundance"]

    # Use the same axes in every panel.
    allValues = [ observedLog.to_numpy() ]

    for fittedValues in fittedData.values():
        allValues.append(
            fittedValues.to_numpy()
        )

    axisMin = min(
        values.min()
        for values in allValues
    )

    axisMax = max(
        values.max()
        for values in allValues
    )

    axisPadding = (axisMax - axisMin) * 0.05
    axisMin -= axisPadding
    axisMax += axisPadding

    panelLabels = ["A", "B", "C", "D", "E", "F"]

    for ax, panelLabel, (modelName, _) in zip(
        axes,
        panelLabels,
        MODEL_SPECS):

        fittedValues = fittedData[modelName]

        resultRow = summaryData.loc[
            summaryData["model"] == modelName
        ].iloc[0]

        ax.scatter(
            fittedValues,
            observedLog,
            s=32,
            alpha=0.7,
            color="0.15"
        )

        # Fitted regression line.
#         slope, intercept = np.polyfit(
#             fittedValues,
#             observedLog,
#             1
#         )
# 
#         xLine = np.linspace(
#             fittedValues.min(),
#             fittedValues.max(),
#             100
#         )
# 
#         yLine = intercept + slope * xLine
# 
#         ax.plot(
#             xLine,
#             yLine,
#             linewidth=1.2,
#             color="0.15",
#             linestyle="--"
#         )
        
        # 1:1 reference line, not a fitted regression line.
        ax.plot([axisMin, axisMax], [axisMin, axisMax],
                linewidth=1.0, linestyle="--", color="0.65")

        ax.text(0.05, 0.95,
#                 f"{panelLabel}  {modelName}\n"
                f"{modelName}\n"
                f"Spearman ρ = {resultRow['spearmanRho']:.3f}\n"
#                 f"R² = {resultRow['rSquared']:.2f}\n"
#                 f"Adjusted R² = {resultRow['adjustedRSquared']:.2f}\n"
                f"n = {int(resultRow['n'])}",
                transform=ax.transAxes,
                ha="left", va="top")

        ax.set_xlim(axisMin, axisMax)
        ax.set_ylim(axisMin, axisMax)

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    fig.supxlabel("Fitted mosquito abundance\n"
                  "log10(abundance24h + 1)")

    fig.supylabel("Observed mosquito abundance\n"
                  "log10(abundance24h + 1)")

    fig.tight_layout()

    return fig


if __name__ == "__main__":
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    envDf, mosDf = readInputData()
    mosDf = addLagVariables(envDf, mosDf)
    analysisDf = getAnalysisData(mosDf)
    summaryDf, coefficientDf, fittedData = fitModels(analysisDf)

    print(f"Mosquito events from 2017–2024: {len(mosDf)}")
    print(f"Events used in all six models: {len(analysisDf)}")
    print()

    displayColumns = [
        "model",
        "terms",
        "n",
        "spearmanRho",
        "rSquared", "adjustedRSquared",
        "beta0", "beta1", "beta2", "beta3", "beta4", "beta5" ]

    print(summaryDf.loc[:, displayColumns].round(4).to_string(index=False))

    summaryDf.to_csv(     OUTPUT_DIR / "neon-mosquito-model-comparison.csv", index=False )
    coefficientDf.to_csv( OUTPUT_DIR / "neon-mosquito-model-coefficients.csv", index=False )

    fig = makeFigure(analysisDf, summaryDf, fittedData)

    plt.show()
