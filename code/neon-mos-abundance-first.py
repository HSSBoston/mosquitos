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
def readNeonFiles(path: Path, pattern: str):
    files = list(path.glob(pattern))
    if not files: raise FileNotFoundError(f"No files found matching: {pattern}")

    data = pd.concat( [pd.read_csv(file) for file in files],
                      ignore_index=True )
    return data

def getMonthlySummary(path:Path):
    trappingData  = readNeonFiles(path, "*mos_trapping*.csv")
    sortingData   = readNeonFiles(path, "*mos_sorting*.csv")
    idData        = readNeonFiles(path, "*mos_expertTaxonomistIDProcessed*.csv")

    # Remove duplicate records, if any
    if "uid" in trappingData.columns: trappingData = trappingData.drop_duplicates(subset="uid")
    if "uid" in sortingData.columns:  sortingData = sortingData.drop_duplicates(subset="uid")
    if "uid" in idData.columns:       idData = idData.drop_duplicates(subset="uid")


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
    sampleCounts["estimatedCount"] = sampleCounts["identifiedCount"] / sampleCounts["proportionIdentified"]


    # Add estimated mosquito counts to trappingData
    trappingData = trappingData.merge(
        sampleCounts[["sampleID","estimatedCount"]],
        on="sampleID",
        how="left")

    # Distinguish true zero catches from unsampled events. 
    # true zero catch:
    #   → trapping record exists (sampleID existis in trappingData)
    #   → targetTaxaPresent = "N" in trappingData
    #   → no row exist in sortingData
    #   → no proportionIdentified value
    #   → estimatedCount = NaN
    # Unsampled records are not treated as zeros.
    zeroCatch = (trappingData["sampleID"].notna() &
                 trappingData["targetTaxaPresent"].eq("N") )

    trappingData.loc[zeroCatch & trappingData["estimatedCount"].isna(),
                     "estimatedCount"] = 0

    # Keep records where trapping actually occurred
    trappingData["trapHours"] = pd.to_numeric(trappingData["trapHours"], errors="coerce")

    sampledData = trappingData.loc[
        trappingData["sampleID"].notna() & (trappingData["trapHours"] > 0) ].copy()


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

    # Combine day + night samples for each plot within each event.
    # Example output:
    #     eventID  plotID  eventStart estimatedMosquitoes totalTrapHours intervalCount 
    #  0  HARV...  HARV...      
    #  1  HARV...  HARV...      
    plotEventData = sampledData.groupby(["eventID","plotID"], as_index=False).agg(
        eventStart          =("setDate", "min"),
        estimatedMosquitoes =("estimatedCount", "sum"),
        totalTrapHours      =("trapHours", "sum"),
        intervalCount       =("nightOrDay", "nunique") )

    # Keep only complete plots. A complete plot has both:
    #   - one daytime trapping interval
    #   - one nighttime trapping interval
    completePlotData = plotEventData.loc[ plotEventData["intervalCount"]==2 ].copy()

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



if __name__ == "__main__":
    targetDir = DATA_DIR / "NEON.D01.HARV.DP1.10043.001.2024-07.expanded.20260123T000749Z.RELEASE-2026"
    
    summaryDf = getMonthlySummary(targetDir)

    print( summaryDf.to_string(index=False) )
    summaryDf.to_csv(OUTPUT_DIR / "neon-mos-abundance-by-event.csv", index=False)
