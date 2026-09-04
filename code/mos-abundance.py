from pathlib import Path
import numpy as np, pandas as pd

PRJ_DIR    = Path(__file__).parent
DATA_DIR   = PRJ_DIR / "data" / "NEON_count-mosquitoes"
OUTPUT_DIR = PRJ_DIR / "output"

# False: Keep sampled traps even if a potential QC issue is present.
#
# True:
#   Exclude traps with obvious equipment/sample problems.
EXCLUDE_QC_ISSUES = False


# Helper function to read and combine all monthly files that match a given pattern.
#
def readNeonFile(path: Path, pattern: str, *, required: bool):
    if not path.is_dir(): raise NotADirectoryError(f"Data directory not found: {path}")

    files = [file for file in path.glob(pattern) if file.is_file()]
    if len(files) == 0:
        if required:
            raise FileNotFoundError(f"No file found matching: {path / pattern}")
        return None
    if len(files) > 1:
        raise RuntimeError(f"Expected one file matching {path / pattern}, but found {len(files)}")

    return pd.read_csv(files[0])

def getMonthlySummary(path:Path):
    trappingData  = readNeonFile(path, "*mos_trapping*.csv", required=True)
    sortingData   = readNeonFile(path, "*mos_sorting*.csv",  required=False)
    idData        = readNeonFile(path, "*mos_expertTaxonomistIDProcessed*.csv", required=False)

    # Remove duplicate records, if any
    if "uid" in trappingData.columns: trappingData = trappingData.drop_duplicates(subset="uid")
    
    if sortingData is not None and "uid" in sortingData.columns:
        sortingData = sortingData.drop_duplicates(subset="uid")
    if idData is not None and "uid" in idData.columns:
        idData = idData.drop_duplicates(subset="uid")

    trappingData["trapHours"] = pd.to_numeric(trappingData["trapHours"], errors="coerce")
    positiveCatch = trappingData["targetTaxaPresent"].eq("Y")

    if sortingData is None:
        # idData cannot be interpreted without its parent sorting table (sortingData).
        if idData is not None: raise ValueError("ID data exist, but sorting data are missing.")
        # sortingData is missing when positive catch exists.
        if positiveCatch.any():
            affectedRows = trappingData.loc[
                positiveCatch,
                ["eventID", "plotID", "sampleID", "targetTaxaPresent"] ]
            raise ValueError("Some traps contained mosquitoes, but no sorting file was found:\n"
                             + affectedRows.to_string(index=False) )

        trappingData["estimatedCount"] = np.nan

    else:
        if idData is None: raise FileNotFoundError("Sorting data exist, but the ID is missing.")

        # Sum up identified mosquitoes within each subsample.
        # One subsample can have multiple rows because mosquitoes are separated by
        # scientificName, sex, etc. in idData. Generates a DataFrame that shows the
        # total nubmer of identified mosquitos per subsampleID; for example:
        #    subsampleID  identifiedCount
        #  0          S1               55
        #  1          S2               15
        identifiedCounts = idData.groupby("subsampleID").agg(
                identifiedCount=("individualCount", "sum")
            ).reset_index()

        # Merge two DataFrames: sortingData and identifiedCounts.
        sampleCounts = sortingData[ ["sampleID", "subsampleID", "proportionIdentified"] ].merge(
            identifiedCounts,
            on="subsampleID",
            how="left")
            # how="left": Keep all rows in the left DF (sortingData) and attach matching info
            # from the right DF (identifiedCounts)

        # Check for invalid proportions.
        invalidProportion = ( (sampleCounts["proportionIdentified"].isna())
                            | (sampleCounts["proportionIdentified"] <= 0)
                            | (sampleCounts["proportionIdentified"] > 1) )

        if invalidProportion.any():
            invalidRows = sampleCounts.loc[ invalidProportion,
                                            ["sampleID","subsampleID","proportionIdentified"]]
            raise ValueError(
                "Invalid proportionIdentified values found:\n"
                + invalidRows.to_string(index=False))

        # Estimate the total sample count based on identifiedCount and proportionIdentified
        # columns in sortingData: 
        #   estimatedCount = identifiedCount / proportionIdentified
        sampleCounts["estimatedCount"] = (
            sampleCounts["identifiedCount"] / sampleCounts["proportionIdentified"] )

        # Add estimated mosquito counts to trappingData
        trappingData = trappingData.merge(
            sampleCounts[["sampleID","estimatedCount"]],
            on="sampleID",
            how="left")

    # targetTaxaPresent == N is the zero-catch indicator.
    trappingData.loc[
        trappingData["targetTaxaPresent"].eq("N"),
        "estimatedCount"] = 0

    # Positive trapHours indicates that the trap actually operated.
    # Do not require sampleID: valid zero catches may not have one.
    sampledData = trappingData.loc[
        trappingData["trapHours"] > 0
    ].copy()

    # Basic quality control flag
    sampledData["qcPass"] = True

    def requireValue(columnName, acceptableValue):
        if columnName not in sampledData.columns: return
        sampledData.loc[
            sampledData[columnName].notna() & ~sampledData[columnName].eq(acceptableValue),
            "qcPass"] = False

    requireValue("samplingImpractical", "OK")
    requireValue("fanStatus",           "On")
    requireValue("catchCupStatus",      "OK")
    requireValue("sampleCondition",     "No known compromise")
    requireValue("CO2Status",           "Present")

    if "dataQF" in sampledData.columns:
        sampledData.loc[
            sampledData["dataQF"].notna(),
            "qcPass"] = False

    if EXCLUDE_QC_ISSUES:
        sampledData = sampledData.loc[
            sampledData["qcPass"]==True].copy()

    # Verify that every sampled trap has an abundance estimate
    missingAbundance = sampledData["estimatedCount"].isna()

    if missingAbundance.any():
        missingRows = sampledData.loc[
            missingAbundance,
            ["eventID", "plotID", "sampleID", "targetTaxaPresent"] ]
        raise ValueError(
            "Some sampled traps have no estimated mosquito count.\n"
            "Check the sorting/identification data for these samples:\n"
            + missingRows.to_string(index=False) )

    # Parse dates
    sampledData["setDate"]     = pd.to_datetime(sampledData["setDate"],     utc=True)
    sampledData["collectDate"] = pd.to_datetime(sampledData["collectDate"], utc=True)

    sampledData["isDay"]   = sampledData["nightOrDay"].eq("day")
    sampledData["isNight"] = sampledData["nightOrDay"].eq("night")

    # Combine day + night samples for each plot within each event.
    # Example output:
    #     eventID  plotID  eventStart estimatedMosquitoes totalTrapHours dayCount nightCount 
    #  0  HARV...  HARV...      
    #  1  HARV...  HARV...      
    plotEventData = sampledData.groupby(["eventID","plotID"], as_index=False).agg(
        eventStart          =("setDate", "min"),
        estimatedMosquitoes =("estimatedCount", "sum"),
        totalTrapHours      =("trapHours", "sum"),
        dayCount=("isDay", "sum"),
        nightCount=("isNight", "sum") )

    # Keep only complete plots. In 2018 and after, a complete plot has both:
    #   - one daytime trapping interval
    #   - one nighttime trapping interval
    # In 2017 and before, a complete plot has both:
    #   - one daytime trapping interval
    #   - two nighttime trapping interval
    pre2018 = plotEventData["eventStart"].dt.year < 2018
    completePlot = (
        # Before 2018: one day and two nights
        (pre2018
         & plotEventData["dayCount"].eq(1)
         & plotEventData["nightCount"].eq(2)
        )
        |
        # Beginning in 2018: one day and one night
        (~pre2018
         & plotEventData["dayCount"].eq(1)
         & plotEventData["nightCount"].eq(1)
        )
    )
    completePlotData = plotEventData.loc[completePlot].copy()

    # Normalize each plot to 24 trap-hours
    completePlotData["abundance24hPlot"] = (
        completePlotData["estimatedMosquitoes"] / completePlotData["totalTrapHours"] * 24 )

    # Calculate site-level abundance for each sampling event
    summaryData = completePlotData.groupby("eventID", as_index=False).agg(
        eventStart    =("eventStart", "min"),
        completePlots =("plotID", "nunique"),
        abundance24h  =("abundance24hPlot", "mean") )

    # Format the output table
    summaryData["eventStart"]   = summaryData["eventStart"].dt.date
    summaryData["abundance24h"] = summaryData["abundance24h"].round(1)
    summaryData = summaryData.sort_values("eventStart").reset_index(drop=True)
    
    return summaryData


def getSummary(path:Path):
    if not path.is_dir(): raise NotADirectoryError(f"Data directory not found: {path}")
    
    directories = sorted(
        directory for directory in path.glob("NEON.D01.HARV.DP1.10043.001.*")
        if directory.is_dir())
    
    if not directories: raise FileNotFoundError(f"No NEON monthly directories found in: {path}")
    
    summaries = [getMonthlySummary(directory) for directory in directories]
    print( f"{len(summaries)}-month data were inspected.")
    
    return pd.concat(summaries, ignore_index=True).sort_values("eventStart").reset_index(drop=True)


if __name__ == "__main__":
#     targetDir = DATA_DIR / "NEON.D01.HARV.DP1.10043.001.2024-07.expanded.20260123T000749Z.RELEASE-2026"
#     targetDir = DATA_DIR / "NEON.D01.HARV.DP1.10043.001.2018-04.expanded.20260123T000749Z.RELEASE-2026"
    targetDir = DATA_DIR / "NEON.D01.HARV.DP1.10043.001.2018-05.expanded.20260123T000749Z.RELEASE-2026"
    
    summaryDf = getMonthlySummary(targetDir)
    print( summaryDf.to_string(index=False) )
    summaryDf.to_csv(OUTPUT_DIR / "mos-abundance-by-event-single-mo.csv", index=False)
    
    combinedDf = getSummary(DATA_DIR)
    print( combinedDf.to_string(index=False) )
    combinedDf.to_csv(OUTPUT_DIR / "mos-abundance-by-event-multiple-mo.csv", index=False)
    
    
