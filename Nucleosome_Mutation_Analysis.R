#Ensure that BiocManager is installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("lomb")
BiocManager::install("quantmod")

library(lomb)
library(data.table)

source("Parse_Nucleosome_Mutation_Data.R")
source("Assymetry_Analysis.R")

# Generate a list of prefixes for the data files so that iterating through them is easier later.
filePrefixes = c("","CtoA_omitted_","CtoG_omitted_","CtoT_omitted_",
                 "TtoA_omitted_","TtoC_omitted_","TtoG_omitted_")

# Generate some lists that will store all of our output data.
dataSetNames=c("No Omissions",filePrefixes[-1])

peakPeriodicities = numeric(length(filePrefixes))
periodictyPValues = numeric(length(filePrefixes))

peakAssymetryTValue = numeric(length(filePrefixes))
peakAssymetryPValue = numeric(length(filePrefixes))

valleyAssymetryTValue = numeric(length(filePrefixes))
valleyAssymetryPValue = numeric(length(filePrefixes))

for (i in 1:length(filePrefixes)) {
  
  ##### Parse Data #####
  
  # Parse the data into a normalized format
  normalizedData = parseWyrickNucleosomeMutationData(
    paste0("Data/",filePrefixes[i],"ESAD-UK_assymetry_analysis.txt"),
    paste0("Data/",filePrefixes[i],"ESAD-UK_assymetry_analysis_strand_background.txt"))
  
  ##### Periodicity Analysis #####
  
  # Calculate the periodicity of the data using a Lomb-Scargle periodiagram.
  lombResult = lsp(normalizedData[,.(Dyad_Position,Normalized_Both_Strands)], type = "period", plot = FALSE)
  # Store the relevant results!
  peakPeriodicities[i] = lombResult$peak.at[1]
  pValues[i] = lombResult$p.value
  
  
  ##### Assymetry Analysis #####
  
  # BEWARE THE TEXAS SHARPSHOOTER!
  
  # # Run assymetry analysis on the peaks and valleys of each strand.
  # runExtremeAnalysisSuite(normalizedData$Normalized_Plus_Strand, dyadPosCutoff = 68)
  # runExtremeAnalysisSuite(normalizedData$Normalized_Minus_Strand, dyadPosCutoff = 68)
  # runExtremeAnalysisSuite(normalizedData$Normalized_Plus_Strand, maxes = FALSE, dyadPosCutoff = 68)
  # runExtremeAnalysisSuite(normalizedData$Normalized_Minus_Strand, maxes = FALSE, dyadPosCutoff = 68)
  
  # Run assymetry analysis on the peaks and valleys of the aligned strands.
  peakAssymetryResult = 
    runExtremeAnalysisSuite(normalizedData$Normalized_Aligned_Strands, dyadPosCutoff = 68)
  peakAssymetryTValue[i] = peakAssymetryResult$statistic
  peakAssymetryPValue[i] = peakAssymetryResult$p.value
  
  
  valleyAssymetryResult = 
    runExtremeAnalysisSuite(normalizedData$Normalized_Aligned_Strands, maxes = FALSE, dyadPosCutoff = 68)
  valleyAssymetryTValue[i] = valleyAssymetryResult$statistic
  valleyAssymetryPValue[i] = valleyAssymetryResult$p.value
  
  # # Negative controls for assymetry analysis on peaks and valleys
  # runExtremeAnalysisSuite(normalizedData$Normalized_Both_Strands, dyadPosCutof = 68)
  # runExtremeAnalysisSuite(normalizedData$Normalized_Both_Strands, maxes = FALSE, dyadPosCutoff = 68)
  
}

# Create data.tables for all the results.
periodicityResults = data.table(Data_Set=dataSetNames,Peak_Periodicity=peakPeriodicities,PValue=periodicityPValues)

extremeAssymetryResults = data.table(Data_Set=dataSetNames, 
                                     Peak_Assymetry_tValue = peakAssymetryTValue, 
                                     Peak_Assymetry_PValue = peakAssymetryPValue,
                                     Valley_Assymetry_tValue = valleyAssymetryTValue, 
                                     Valley_Assymetry_PValue = valleyAssymetryPValue,)