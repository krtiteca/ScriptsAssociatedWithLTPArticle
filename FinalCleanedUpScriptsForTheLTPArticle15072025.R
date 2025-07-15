################ Document with the clean-up and calculations for figures of LTP article (updated last version after reviewer comments)
# Note: The current version of the script is a reduced version derived from larger scripts, 
# and only contains the parts focusing on the main lines of code for the analysis and visualization, 
# as well as some side graphs and side checks. Hence, it might contain some historic intermediate variable assignments 
# that seem unneccessary but that I kept in to maintain consistency.

## Legend for figure files and some other parts:
# # Done 
# (#) Done in other context, or consequence of other files
#  #X Extra input
#  #. Relevant files for this part

################ Structure of this document

######## General functions 

######## MS-data input & clean-up
#### Extended Data Table 4: With screen data combined

######## Figure 1 & supplementary material: quality tests
#### Extended Data Figure 1c: Quality control by mean retention time vs. lipid carbon chain length and unsaturation
#### Extended Data Figure 1d: General adducts distribution
#### Fig.1b-c # See further part
#### Fig.1d # See further part

######## Figure 5 
#### Fig.5a_a barchart on the top 
#### Extended Data Figure 2b_b: Donut-diagram
#### Fig.5a: circle-legend
#### Extended Data Figure 2b_a: Side-histogram
#### Fig.5a: circle heatmap visualization
#### #X Extra interesting visualization, not part of included figures anymore. (Previously Fig.2c: unsaturation distribution)
#### Figs.5a&5b extra layers: see further below

######## Figure 2
#### Fig.2a
#### Extended Data Figure 2c: the protein domain based ordering of LTPs 
#### Fig.6c & Extended Data Fig.6a
#### Fig.2b: see further, and external input for Figs.2c-d 

######## Figure 3
#### Fig.3a: Species version
#### (Fig. 3a: Subclasses version) 
#### Fig.3b

######## Figure 4
#### Panel 4A
#### Panel 4B
#### Extended Data Table 9A: subcellular-colocalization lipid species # See further
#### Panel 4C
#### (Alternatives for panels of figure 4)
#### Extended Data Table 9B: Results of the Fisher's exact tests

######## Figure 6
#### Fig.6a (see also further)
#### Fig.6c (see also further and under chapter for Figure 2)
#### (Extra information on the effect of lipid subclasses upon CERT overexpression in HEK293 cells)
#### Extended Data Figure 5a: Extra information on the effect of lipid subclasses upon CERT overexpression in HeLa cells # See further
#### Extended Data Figs.5c: Observations of very long saturated DHCer and tCer in METASPACE

####### Additional entries
#### Ext. Data Table X1/7B #!
#### Ext. Data Table X2/7C #!
#### Fig.2b
#### Ext. Data Table 7D #!
#### Ext. Data Table 10A
#### Liposome layer for Fig.6c
#### Liposome layer for Fig.6a
#### Fig.5a layer
#### Extended Data Table 9A
#### Upper part of Fig.5b_b
#### Lower part of Fig.5b_b & legend
#### Upper part of Fig.5b_a
#### Lower part of Fig.5b_a & legend
#### Fig.1b&c basis
#### Fig.1d parts
#### Parts Ext. Data Table 6
#### Extended Data Table 5A
#### Extended Data Table 5B
#### Extended Data Table 5C
#### Basis for Extended Data Table 11
#### Extended Data Table 10B
#### Extended Data Table 7A
#### Extended Data Table 8
#### Extended Data Fig.5a


# Capture previous options and set new general options
OriginalOptions <- options()

options(warn = 0)
options(stringsAsFactors = TRUE)

######## General functions
#### Function to substitute rownames by first column and then remove the first column

Col1ToRowNames <- function(x){
  y <- x[,-1]
  
  rownames(y) <- x[,1]
  return(y)}


#### Function for min-max normalization (for wide format and for vectors)

MinMaxNormMatrixFunc <- function(inputdata){
  Int <- log10(inputdata)
  
  MaxInt <- max(Int[is.finite(as.matrix(Int))])
  MinInt <- min(Int[is.finite(as.matrix(Int))])
  
  return((Int - MinInt)/(MaxInt - MinInt))
}


#### Function for min-max normalization with integrated log10 and separate handling of ion modes (for long dataframe format)

MinMaxNormFuncn <- function(inputint, ionmode, inputdata){
  MaxInt <- max(log10(inputdata[inputdata$IonMode == ionmode, "Intensity"]))
  
  MinInt <- min(log10(inputdata[inputdata$IonMode == ionmode, "Intensity"]))
  return((log10(inputint) - MinInt)/(MaxInt - MinInt))
  
}


#### Function to convert the Zeros to NAs, and function to do the inverse
ZerosToNAsConverter <- function(x){x[x == 0] <- NA; return(x)}

NAsToZerosConverter <- function(x){x[is.na(x)] <- 0; return(x)}


#### NAs from vector removal function
RemoveNAsFromVector <- function(x){x[!is.na(x)]}

#### Function to determine Fisher exact tests for negative vs. positive parts (especially for the co-regulation comparison, and not in the final figures)
FishersExactTestNegPos2 <- function(SubsetX,SupersetY){
  
  PosX <- sum(SubsetX >= 0, na.rm = TRUE)
  NegX <- sum(SubsetX < 0, na.rm = TRUE)
  
  PosY <- sum(SupersetY >= 0, na.rm = TRUE)
  NegY <- sum(SupersetY < 0, na.rm = TRUE)
  
  InputMatrixForFisherTest <- rbind("No" = c(NegY-NegX,PosY-PosX), "Yes" = c(NegX,PosX))
  colnames(InputMatrixForFisherTest) <- c("Neg", "Pos")
  
  OutputOfFisherTest <- fisher.test(InputMatrixForFisherTest)
  return(c(InputMatrixForFisherTest["No", "Neg"],
           
           InputMatrixForFisherTest["No", "Pos"],
           InputMatrixForFisherTest["Yes", "Neg"],
           
           InputMatrixForFisherTest["Yes", "Pos"],
           OutputOfFisherTest$conf.int[1],
           
           OutputOfFisherTest$conf.int[2],
           OutputOfFisherTest$estimate,
           
           OutputOfFisherTest$p.value
  ))
  
}


#### Function to determine Fisher exact tests for below or above 0.5 (especially for non-co-regulation comparisons, and not in the final figures)
FishersExactTestBelowAbove05 <- function(SubsetX,SupersetY){
  
  PosX <- sum(SubsetX >= 0.5, na.rm = TRUE)
  NegX <- sum(SubsetX < 0.5, na.rm = TRUE)
  
  PosY <- sum(SupersetY >= 0.5, na.rm = TRUE)
  NegY <- sum(SupersetY < 0.5, na.rm = TRUE)
  
  InputMatrixForFisherTest <- rbind("No" = c(NegY-NegX,PosY-PosX), "Yes" = c(NegX,PosX))
  colnames(InputMatrixForFisherTest) <- c("Neg", "Pos")
  
  OutputOfFisherTest <- fisher.test(InputMatrixForFisherTest)
  return(c(InputMatrixForFisherTest["No", "Neg"],
           
           InputMatrixForFisherTest["No", "Pos"],
           InputMatrixForFisherTest["Yes", "Neg"],
           
           InputMatrixForFisherTest["Yes", "Pos"],
           OutputOfFisherTest$conf.int[1],
           
           OutputOfFisherTest$conf.int[2],
           OutputOfFisherTest$estimate,
           
           OutputOfFisherTest$p.value
  ))
  
}


#### Function to determine Fisher exact tests based on the medians (for all comparisons & in final figures)
FishersExactTestNegPosEqualBackground2 <- function(SubsetX,SupersetY){
  
  PosXEV <- sum(SubsetX >= median(SupersetY), na.rm = TRUE)
  NegXEV <- sum(SubsetX < median(SupersetY), na.rm = TRUE)
  
  PosYEV <- ceiling(length(SubsetX)/2)
  NegYEV <- floor(length(SubsetX)/2)
  
  InputMatrixForFisherTestEV <- rbind("No" = c(NegYEV, PosYEV), "Yes" = c(NegXEV, PosXEV))
  colnames(InputMatrixForFisherTestEV) <- c("Neg", "Pos")
  
  OutputFisherTestEV <- fisher.test(InputMatrixForFisherTestEV)
  return(c(InputMatrixForFisherTestEV["No", "Neg"],
           
           InputMatrixForFisherTestEV["No", "Pos"],
           InputMatrixForFisherTestEV["Yes", "Neg"],
           
           InputMatrixForFisherTestEV["Yes", "Pos"],
           OutputFisherTestEV$conf.int[1],
           
           OutputFisherTestEV$conf.int[2],
           OutputFisherTestEV$estimate,
           
           OutputFisherTestEV$p.value
  ))
  
}


#### Function to extract data between brackets in strings 
GetStuffBetweenBrackets <- function(x){regmatches(x, gregexpr("(?<=\\().*?(?=\\))", x, perl=T))} 

#### Function to calculate the Manders' Overlap Coefficient
MOCFunction <- function(VectorA, VectorB){sum(VectorA*VectorB)/sqrt(sum(VectorA^2)*sum(VectorB^2))}

#### Function to convert matrices from a wider version to a longer version (simple version of melt)
# Not advised for dataframes (rownames are taken as numbers, and it forces everything into a factor for the value column)

meltmatrix <- function(data){do.call("rbind", lapply(colnames(data), function(x){
  data.frame(Var1 = factor(rownames(data), levels = rownames(data)), Var2 = x, value = as.vector(data[,x]))}))}


#### Function to rbind multiple dataframes even if they do not have fully fitting columns, and fill the non-fitting parts with a filler value (FillerValue) (default NA)

RbindMultipleDataframesFillNonmatches <- function(dflist, FillerValue = NA){
  do.call(rbind.data.frame,
          
          lapply(dflist,
                 function(x){do.call("cbind.data.frame", c(x, sapply(setdiff(unique(unlist(lapply(dflist, names))), names(x)), function(y){FillerValue})))}))}


#### Function to covert dataframes to a wider dataframe based on one column and with option to keep multiple columns intact & aggregate similar data

# Has default inputs and more parameters that can be optimized including data type of output and rownames.
# Takes original order without accounting for factor levels, but does not affect the levels in final dataframe.

# inputdf = input dataframe; ColumnsLong = column names (character) to keep intact (can be multiple ones: character vector)
# ColumnWide = column name (character) to widen (can only be one)

# Does not drop columns from a factor-column: uses all levels from it, even if they have no value associated
# ColumnValue = column name (character) to use for values as input to function to populate new columns; depending on input and aggregation function might create 0s where NAs expected.

# AggregatingFunction = user-defined function to aggregate data that correspond to the same new dataframe cells (without definition it defaults to the length function)
# FunctionOutputValueType = expected value type as output from the function (similar to how vapply works) (without definition it defaults to integer)

# RowNamesColumnsLong = logical determining if rownames are textual aggregations of content of ColumnLong-defined columns (TRUE) (the default)
# or if the rownames are just a number running from the first to the last row (1:nrow) (FALSE).

# SortRows = logical determining if the rows need to be sorted alphabethically (TRUE) (the default)
# or if the original order of the rows needs to be maintained as they occur in the input data (FALSE).

# DropUnusedRowsLevels = logical indicating if unused levels of the column indicated by ColumnLong need to be dropped (TRUE) (the default)
# if ColumnLong only defines one column and if that column is a factor, then DropUnusedRowsLevels = FALSE will allow to include all levels, 

# and if SortRows is TRUE the levels will be sorted, but if SortRows is FALSE the original order of the levels will be used.
# FillerValue = filler value for filling cells that did not exist before, but were created by the conversion.

# SortWide = logical indicating if the columns have to be reordered by having the ColumnsLong first followed by the sorted rest of the columns (TRUE) (the default).
# If SortWide is FALSE the columns will be ordered in the order that the function encounters them.

widen5 <- function(inputdf, ColumnsLong, ColumnWide, ColumnValue, AggregatingFunction = length, FunctionOutputValueType = integer(1), RowNamesColumnsLong = TRUE, SortRows = TRUE, DropUnusedRowsLevels = TRUE, FillerValue = NA, SortWide = TRUE){
  if((length(ColumnsLong) == 1) && is.factor(inputdf[,ColumnsLong]) && (DropUnusedRowsLevels == FALSE)){
    
    outputdf <- RbindMultipleDataframesFillNonmatches(
      lapply(split(inputdf, f = inputdf[, ColumnsLong], drop = FALSE)[if(SortRows){sort(levels(inputdf[, ColumnsLong]))}else{levels(inputdf[, ColumnsLong])}],
             
             function(y){cbind(y[1, ColumnsLong, drop = FALSE],t(vapply(split(y, f = y[, ColumnWide]),
                                                                        FUN = function(z){AggregatingFunction(z[,ColumnValue])},
                                                                        
                                                                        FUN.VALUE = FunctionOutputValueType,
                                                                        USE.NAMES = TRUE)))}), FillerValue = FillerValue)
    
    outputdf[,ColumnsLong] <- rownames(outputdf)
  }else{
    
    
    outputdf <- RbindMultipleDataframesFillNonmatches(
      
      lapply(split(inputdf, f = inputdf[, ColumnsLong], drop = TRUE)[if(length(ColumnsLong)>1){if(SortRows){sort(do.call(paste, c(unique(inputdf[, ColumnsLong]), sep=".")))}else{do.call(paste, c(unique(inputdf[, ColumnsLong]), sep="."))}
      }else{if(SortRows){sort(as.character(unique(inputdf[, ColumnsLong])))}else{as.character(unique(inputdf[, ColumnsLong]))}}], # both have option to be sorted or be in in order they original occur
      
      function(y){cbind(y[1, ColumnsLong, drop = FALSE],t(vapply(split(y, f = y[, ColumnWide]),
                                                                 FUN = function(z){AggregatingFunction(z[,ColumnValue])},
                                                                 
                                                                 FUN.VALUE = FunctionOutputValueType,
                                                                 USE.NAMES = TRUE)))}), FillerValue = FillerValue)}
  
  if(RowNamesColumnsLong == FALSE){rownames(outputdf) <- as.character(1:dim(outputdf)[1])}
  if(SortWide == TRUE){outputdf <- outputdf[,c(ColumnsLong, sort(colnames(outputdf)[!(colnames(outputdf) %in% ColumnsLong)]))]}
  
  return(outputdf)}


# Function to take a rowsums if x is a dataframe and otherwise just return x itself.
rowSumsDFS <- function(x){if(class(x) == "data.frame"){rowSums(x, na.rm = TRUE)}else{x}}

# Function to find percentage fraction of the maximum for a vector
DivMaxToPercent <- function(x){x*100/max(x, na.rm = TRUE)}

######## MS-data input 
#### MS results: in cellulo and in vitro: clean-up of input data 

# Note: although SEC14L1 data might be present at different stages, all entries associated with these were removed from the final results as a conservative precaution, 
# because of uncertainties about a potential experimental mix-up (in some cases) of CERT into the SEC14L1 entry (but not the inverse mix-up).

AntonellaData18062019 <- read.csv(file = "./InputData/invivo_postproc_input_final_18062019.tsv", header = TRUE, sep = "\t", as.is = TRUE)
EnricData18062019 <- read.csv(file = "./InputData/invitro_postproc_input_final_18062019.tsv", header = TRUE, sep = "\t", as.is = TRUE)

PureVersionAntonellaData18062019 <- AntonellaData18062019[AntonellaData18062019$cls == "I" | AntonellaData18062019$cls == "IV.7",]
PureVersionEnricData18062019 <- EnricData18062019[EnricData18062019$cls == "I" | EnricData18062019$cls == "IV.7",]

CombinedDataWithPureClasses18062019 <- rbind(cbind(PureVersionAntonellaData18062019, ScreenType = "in vivo"), cbind(PureVersionEnricData18062019 , ScreenType = "in vitro"))
write.table(CombinedDataWithPureClasses18062019, file="./Output/CombinedDataWithPureClasses18062019.csv", sep="\t", row.names = TRUE, quote = FALSE)

# Updated to avoid use of needless extra packages, and increase robustness
CombinedDataWithPureClasses180620192x <- cbind.data.frame(CompletelyUnique = as.character(do.call(paste, c(CombinedDataWithPureClasses18062019[,c("protein","ionm","uid","mz","intensity")], sep = "_"))), CombinedDataWithPureClasses18062019)

CombinedDataWithPureClasses180620192x$CompletelyUnique <- as.character(CombinedDataWithPureClasses180620192x$CompletelyUnique)
CombinedDataWithPureClasses180620192 <- CombinedDataWithPureClasses180620192x

CombinedDataWithPureClasses180620194 <- cbind(CombinedDataWithPureClasses180620192, LipidGroup = sapply(CombinedDataWithPureClasses180620192$CompletelyUnique, function(x){paste(CombinedDataWithPureClasses180620192[CombinedDataWithPureClasses180620192$CompletelyUnique == x, "uhgroup"], collapse = "; ")}))
PrePieceHeadgroupConversion <- cbind(unique(CombinedDataWithPureClasses180620194$uhgroup), c("PC(","PC(O-", "SM(t*", "SM(t*", "SM(t*", "PE(", "PS(", "PI(", "BMP(", "LPC(", "PG(", "HexCer(t*", "HexCer(t*", "HexCer(t*", "HexCer(d*", "HexCer(d*", "Hex2Cer(t*", "Hex2Cer(t*", "Hex2Cer(t*", "SM(d*", "SM(d*", "CerP(d*", "CerP(d*", NA, "VA", "Cer(d*", "Cer(d*", "DAG(", "PE(O-", "TAG(", "Cer(t*", "Cer(t*", "FA(", "LPE(", "LPG(", "FAL(", "SHexCer(d*", "SHexCer(d*", "LPE(O-", "PA(", "PGP(", "CL(", "PA(O-"))

DataRow <- CombinedDataWithPureClasses180620194[2,]
lncfA <- function(DataRow, PrePieceHeadgroupConversion){paste0(PrePieceHeadgroupConversion[match(DataRow[,"uhgroup"], PrePieceHeadgroupConversion[,1]),2], DataRow[,"carb"], ":", DataRow[,"unsat"], ")")}

lncfD1 <- function(DataRow){paste0("PC(O-", DataRow[,"carb"], ":", DataRow[,"unsat"], ")/LPC(", DataRow[,"carb"], ":", DataRow[,"unsat"], ")")}
lncfE <- function(DataRow){paste0("PG/BMP(", DataRow[,"carb"], ":", DataRow[,"unsat"], ")")}

lncfH1 <- function(DataRow){paste0("SM(DH", DataRow[,"carb"], ":", DataRow[,"unsat"], ")")}
lncfJ <- function(DataRow){paste0("Cer(d", DataRow[,"carb"], ":", DataRow[,"unsat"], ")")}

lncfH2 <- function(DataRow){paste0("Cer(DH", DataRow[,"carb"], ":", DataRow[,"unsat"], ")")}
lncfK <- function(DataRow){paste0("Cer(t", DataRow[,"carb"], ":", DataRow[,"unsat"], ")")}

lncfL <- function(DataRow){paste0("Cer(DH", DataRow[,"carb"], ":", DataRow[,"unsat"], "-OH)/Cer(t",  DataRow[,"carb"], ":", DataRow[,"unsat"], ")")}
lncfD2 <- function(DataRow){paste0("PE(O-", DataRow[,"carb"], ":", DataRow[,"unsat"], ")/LPE(", DataRow[,"carb"], ":", DataRow[,"unsat"], ")")}

ListTypesOfDifferentClasses <- unique(lapply(unique(CombinedDataWithPureClasses180620192$CompletelyUnique), function(x){CombinedDataWithPureClasses180620192[CombinedDataWithPureClasses180620192$CompletelyUnique == x, "uhgroup"]}))
VectorTypesOfDifferentClasses <- sapply(ListTypesOfDifferentClasses, function(x){paste(x, collapse = "; ")})

lncfList <- c(lncfA, lncfA, NA, lncfA, NA, lncfA, NA, lncfA, lncfA, lncfD1, lncfE, rep(list(lncfA), 15), NA, rep(list(lncfA), 9), lncfH1, NA, rep(list(lncfA), 3), NA, NA, NA, "VA", lncfJ, rep(list(lncfA), 3), NA, lncfD1, NA, NA, lncfA, lncfA, lncfH2, lncfK, lncfL, rep(list(lncfA), 5), lncfD2, rep(list(lncfA), 6), lncfE, rep(list(lncfA), 4), NA, NA)
names(lncfList) <- VectorTypesOfDifferentClasses

CombinedDataWithPureClasses180620195 <- cbind(CombinedDataWithPureClasses180620194, ConsensusLipid = sapply(1:nrow(CombinedDataWithPureClasses180620194), function(i){
  if(is.function(lncfList[[as.character(CombinedDataWithPureClasses180620194[i,"LipidGroup"])]])){
    
    if(identical(lncfList[[as.character(CombinedDataWithPureClasses180620194[i,"LipidGroup"])]], lncfA)){
      lncfList[[as.character(CombinedDataWithPureClasses180620194[i,"LipidGroup"])]](CombinedDataWithPureClasses180620194[i,], PrePieceHeadgroupConversion)
      
    }else{lncfList[[as.character(CombinedDataWithPureClasses180620194[i,"LipidGroup"])]](CombinedDataWithPureClasses180620194[i,])
    }
    
  }else{lncfList[[as.character(CombinedDataWithPureClasses180620194[i,"LipidGroup"])]]}
}))

SelectedColumnsAndNames <- cbind(c("protein","ionm", "uid", "mz", "intensity", "rtmean", "rtmin", "rtmax", "carb", "unsat", "fa1c", "fa1u", "fa2c", "fa2u", "fa3c", "fa3u", "fa4c", "fa4u", "adduct", "domain", "lit", "manual_id1", "ScreenType", "LipidGroup", "ConsensusLipid"),
                                 c("LTPProtein", "IonMode", "UniqueEntry", "mz", "Intensity", "MeanOfRetentionTime", "MinimumOfRetentionTime", "MaximumOfRetentionTime", "TotalCarbonChainLength", "TotalCarbonChainUnsaturations", "CarbonsOfFattyAcidA", 
                                   
                                   "UnsaturationsOfFattyAcidA", "CarbonsOfFattyAcidB", "UnsaturationsOfFattyAcidB", "CarbonsOfFattyAcidC", "UnsaturationsOfFattyAcidC", "CarbonsOfFattyAcidD", 
                                   "UnsaturationsOfFattyAcidD", "Adduct", "ProteinDomain", "LiteraturePresence", "FullIdentityOfLipid", "Screen", "LipidAmbiguity", "Lipid"))

SelectedColumnsAndNames2 <- SelectedColumnsAndNames[match(c("protein","domain","manual_id1","ConsensusLipid", "ScreenType", "intensity", "lit", "adduct", "carb", "unsat", "ionm", "mz", "rtmean", "rtmin", "rtmax", "fa1c", "fa1u", "fa2c", "fa2u", "fa3c", "fa3u", "fa4c", "fa4u", "LipidGroup", "uid"), SelectedColumnsAndNames[,1]),]


CombinedDataWithPureClasses180620197 <- CombinedDataWithPureClasses180620195[!as.logical(CombinedDataWithPureClasses180620195$set2),SelectedColumnsAndNames2[,1]]
colnames(CombinedDataWithPureClasses180620197) <- SelectedColumnsAndNames2[,2]

CombinedDataWithPureClasses180620198 <- CombinedDataWithPureClasses180620197[!is.na(CombinedDataWithPureClasses180620197$Lipid),]
CombinedDataWithPureClasses1806201910 <- unique(CombinedDataWithPureClasses180620198)

CombinedDataWithPureClasses1806201910$Lipid <- as.character(CombinedDataWithPureClasses1806201910$Lipid)
CombinedDataWithPureClasses1806201910$Screen <- as.character(CombinedDataWithPureClasses1806201910$Screen)

CombinedDataWithPureClasses1806201910$LipidAmbiguity <- as.character(CombinedDataWithPureClasses1806201910$LipidAmbiguity)


CombinedDataWithPureClasses1806201911 <- CombinedDataWithPureClasses1806201910[!((CombinedDataWithPureClasses1806201910$LipidAmbiguity == "PG; BMP" | CombinedDataWithPureClasses1806201910$LipidAmbiguity == "BMP; PG") & CombinedDataWithPureClasses1806201910$IonMode == "neg" & is.na(CombinedDataWithPureClasses1806201910$CarbonsOfFattyAcidA)), ]
CombinedDataWithPureClasses1806201913 <- CombinedDataWithPureClasses1806201911[!((CombinedDataWithPureClasses1806201911$LipidAmbiguity %in% c("PC-O; Lyso-PC", "Lyso-PC; PC-O", "Lyso-PE; PE-O")) & !is.na(CombinedDataWithPureClasses1806201911$CarbonsOfFattyAcidA)),]

write.table(CombinedDataWithPureClasses1806201913, file="./Output/CleanConservativeDataWithoutFilters17072019.csv", sep="\t", row.names = FALSE, quote = FALSE)
# Basis for Extended Data Table 4


library(stringr)

CombinedDataWithPureClasses1806201914 <- cbind(CombinedDataWithPureClasses1806201913, do.call("rbind", lapply(CombinedDataWithPureClasses1806201913$Lipid, function(x){
  if(str_count(x, "\\)") == 0){c(x,NA,NA)
    
  }else{if(str_count(x, "\\)") == 1){
    if(grepl(pattern = "(O-", x, fixed = TRUE)){paste0(strsplit(gsub("\\)","",x), "\\(O-|:")[[1]], c("-O","",""))}
    
    else if(grepl(pattern = "(d*", x, fixed = TRUE)){paste0(c("d*","",""), strsplit(gsub("\\)","",x), "\\(d\\*|:")[[1]])}
    else if(grepl(pattern = "(t*", x, fixed = TRUE)){paste0(c("t*","",""), strsplit(gsub("\\)","",x), "\\(t\\*|:")[[1]])}
    
    else if(grepl(pattern = "\\(d\\d", x, fixed = FALSE)){paste0(c("d","",""), strsplit(gsub("\\)","",x), "\\(d|:")[[1]])}
    else if(grepl(pattern = "\\(t\\d", x, fixed = FALSE)){paste0(c("t","",""), strsplit(gsub("\\)","",x), "\\(t|:")[[1]])}
    
    else if(grepl(pattern = "\\(DH\\d", x, fixed = FALSE)){paste0(c("DH","",""), strsplit(gsub("\\)","",x), "\\(DH|:")[[1]])}
    else if(grepl(pattern = "\\(\\d", x, fixed = FALSE)){strsplit(gsub("\\)","",x), "\\(|:")[[1]]}
    
    
  }else{
    
    if(grepl(pattern = "/C", x, fixed = TRUE)){
      c(paste0("DHOH*", strsplit(gsub("\\)","", strsplit(x,"/")[[1]][2]),"\\(t")[[1]][1]), strsplit(strsplit(gsub("\\)","", strsplit(x,"/")[[1]][2]),"\\(t")[[1]][2],":")[[1]])
      
    }else{if(as.numeric(strsplit(strsplit(gsub("\\)","", strsplit(x,"/")[[1]][2]), "\\(")[[1]][2], ":")[[1]][1]) < 26){
      strsplit(gsub("\\)","", strsplit(x,"/")[[1]][2]), "\\(|:")[[1]]
      
    }else{
      iv <- strsplit(gsub("\\)","", strsplit(x,"/")[[1]][2]), "\\(|:")[[1]]
      
      c(paste0(str_sub(iv[1],-2),"-O"), iv[-1])
    }
      
    }
  }}})))

colnames(CombinedDataWithPureClasses1806201914)[26:28] <- c("LikelySubclass","ChainLength", "TotalUnsaturation")


CombinedDataWithPureClasses1806201914$ChainLength <- as.numeric(as.character(CombinedDataWithPureClasses1806201914$ChainLength))
CombinedDataWithPureClasses1806201914$TotalUnsaturation <- as.numeric(as.character(CombinedDataWithPureClasses1806201914$TotalUnsaturation))

all(CombinedDataWithPureClasses1806201914$ChainLength == CombinedDataWithPureClasses1806201914$TotalCarbonChainLength, na.rm = TRUE) # TRUE --> only difference in NA and NaN
all(CombinedDataWithPureClasses1806201914$TotalUnsaturation == CombinedDataWithPureClasses1806201914$TotalCarbonChainUnsaturations, na.rm = TRUE) # TRUE --> only difference in NA and NaN

CombinedDataWithPureClasses1806201915 <- CombinedDataWithPureClasses1806201914[,1:26]
# Ready for use with heatmaps and other visualizations

PureAntonella32 <- CombinedDataWithPureClasses1806201915[CombinedDataWithPureClasses1806201915[,"Screen"] == "in vivo",]
PureEnric32 <- CombinedDataWithPureClasses1806201915[CombinedDataWithPureClasses1806201915[,"Screen"] == "in vitro",]

# Remove SEC14L1 from in vivo data because not trustworthy; in vitro data stay same
PureAntonella32b <- PureAntonella32[PureAntonella32$LTPProtein != "SEC14L1",]

# Previous code results in: core MS-data that will be used everywhere: in cellulo = PureAntonella32b ; in vitro = PureEnric32


######## Supplementary material: quality tests
#### Supplementary material: ECN test (Extended Data Figure 1c)


library(RColorBrewer)

BlueColorRangeUnsaturations <- setNames(rev(paste0(colorRampPalette(brewer.pal(9,"Blues")[-(1:2)])(11), "52")), as.character(0:10))
OrangeColorRangeUnsaturations <- setNames(rev(paste0(colorRampPalette(brewer.pal(9,"Oranges")[-1])(11), "52")), as.character(0:10))

BlueColorRangeUnsaturationsOpaque <- setNames(rev(colorRampPalette(brewer.pal(9,"Blues")[-(1:2)])(11)), as.character(0:10))
OrangeColorRangeUnsaturationsOpaque <- setNames(rev(colorRampPalette(brewer.pal(9,"Oranges")[-1])(11)), as.character(0:10))

#. LipidomicsQualityControlECNTestsUnsaturationInfoColored21012022.pdf #
pdf("./Output/LipidomicsQualityControlECNTestsUnsaturationInfoColored21012022.pdf")

plot(PureAntonella32b$MeanOfRetentionTime, PureAntonella32b$TotalCarbonChainLength, main = "Quality control by ECN test: in cellulo", xlab = "Mean of retention time", ylab = "Total Carbon Chain Length", pch = 16, 
     col = BlueColorRangeUnsaturations[as.character(PureAntonella32b$TotalCarbonChainUnsaturations)])

legend("topleft", legend = names(BlueColorRangeUnsaturations), col = BlueColorRangeUnsaturations, pch = 16, cex = 0.8, title = "Unsat.")


plot(PureAntonella32b$MeanOfRetentionTime, PureAntonella32b$TotalCarbonChainLength, main = "Quality control by ECN test: in cellulo", xlab = "Mean of retention time", ylab = "Total Carbon Chain Length", pch = 16, 
     col = BlueColorRangeUnsaturationsOpaque[as.character(PureAntonella32b$TotalCarbonChainUnsaturations)])

legend("topleft", legend = names(BlueColorRangeUnsaturationsOpaque), col = BlueColorRangeUnsaturationsOpaque, pch = 16, cex = 0.8, title = "Unsat.")
# Previous: selected plot for inclusion in the figures

plot(PureEnric32$MeanOfRetentionTime, PureEnric32$TotalCarbonChainLength, main = "Quality control by ECN test: in cellulo", xlab = "Mean of retention time", ylab = "Total Carbon Chain Length", pch = 16, 
     col = OrangeColorRangeUnsaturations[as.character(PureEnric32$TotalCarbonChainUnsaturations)])

legend("topleft", legend = names(OrangeColorRangeUnsaturations), col = OrangeColorRangeUnsaturations, pch = 16, cex = 0.8, title = "Unsat.")


plot(PureEnric32$MeanOfRetentionTime, PureEnric32$TotalCarbonChainLength, main = "Quality control by ECN test: in cellulo", xlab = "Mean of retention time", ylab = "Total Carbon Chain Length", pch = 16, 
     col = OrangeColorRangeUnsaturationsOpaque[as.character(PureEnric32$TotalCarbonChainUnsaturations)])

legend("topleft", legend = names(OrangeColorRangeUnsaturationsOpaque), col = OrangeColorRangeUnsaturationsOpaque, pch = 16, cex = 0.8, title = "Unsat.")
# Previous: selected plot for inclusion in the figures

dev.off()


#### Supplementary material: adducts distribution (Extended Data Figure 1d)
#. ObservedNumberOfAdductsPerEachScreen15032022.pdf #X

#. LipidomicsQualityControlLTPLipidSpeciesCombinationsDistributionOfAdductsForTheScreensIndependentlyInOneGraph110320221g2.pdf (#)
#. (LipidomicsQualityControlLTPLipidSpeciesCombinationsDistribution20012022.pdf) #X

#. SupplementaryFigure1BLipidomicsQualityControlLTPLipidSpeciesCombinationsDistributionOfAdductsForTheScreensIndependentlyInOneGraph110320221g203072022 #


all(colnames(PureAntonella32b) == colnames(PureEnric32)) # TRUE
PureCombined32 <- rbind(PureAntonella32b, PureEnric32)


library(RColorBrewer)

ColorMatrixTryOut <- rbind(even = c('in cellulo' = colorRampPalette(brewer.pal(9,"Blues")[-1])(9)[7], 'in vitro' = colorRampPalette(brewer.pal(9,"Oranges")[-1])(9)[7]),
                           odd = c('in cellulo' = colorRampPalette(brewer.pal(9,"Blues")[-1])(9)[5], 'in vitro' = colorRampPalette(brewer.pal(9,"Oranges")[-1])(9)[5]))

MergedTablesForAdducts <- t(Col1ToRowNames(merge(merge(cbind(1:10,1:10), table(table(paste(PureCombined32[PureCombined32[,"Screen"] == "in vivo","LTPProtein"], PureCombined32[PureCombined32[,"Screen"] == "in vivo","Lipid"], sep = "_"))), by = 1, all = TRUE),
                                                 merge(cbind(1:10,1:10), table(table(paste(PureCombined32[PureCombined32[,"Screen"] == "in vitro","LTPProtein"], PureCombined32[PureCombined32[,"Screen"] == "in vitro","Lipid"], sep = "_"))), by = 1, all = TRUE), by = 1, all = TRUE)[,c(1,3,5)]))

rownames(MergedTablesForAdducts) <- c("in cellulo", "in vitro")


# Adduct distribution
AdductTabelScreensCombined <- table(table(paste(PureCombined32[,"LTPProtein"], PureCombined32[,"Lipid"], sep = "_")))


pdf("./Output/LipidomicsQualityControlLTPLipidSpeciesCombinationsDistributionOfAdductsForTheScreensIndependentlyInOneGraph11032022.pdf")

barplot(MergedTablesForAdducts, beside = TRUE, las = 2, col = ColorMatrixTryOut["odd",], main = "Number of LTP - lipid species as different adducts", las = 1, xlab = "Number of adducts observed", ylab = "Number of LTPs - lipid species")
abline(h = seq(from = 50, to = 200, by = 50), col = "#FFFFFF33", lwd = 3.2)

dev.off()
# Previous: selected plot for inclusion in the figures

# What percentages covered with multiple adducts per screen
100-((MergedTablesForAdducts[,1]*100)/rowSums(MergedTablesForAdducts, na.rm = TRUE))

# in cellulo   in vitro 
# 26.07004   44.62617 

# What percentage covered with multiple adducts in both screens combined
100-(AdductTabelScreensCombined[1]*100)/sum(AdductTabelScreensCombined, na.rm = TRUE) # 40.94488%


# Independent visualization of previous figures

par()$mar # 5.1 4.1 4.1 2.1
pdf("./Output/ObservedNumberOfAdductsPerEachScreen15032022.pdf")


par(mar = c(12.8, 4.1, 4.1, 2))

barplot(sort(table(PureCombined32[PureCombined32[,"Screen"] == "in vivo","Adduct"]), decreasing = TRUE), col = ColorMatrixTryOut["odd",1], las = 2, ylab = "Observed number of LTP - lipid species pairs", main = "Distribution of adducts in cellulo")
barplot(sort(table(PureCombined32[PureCombined32[,"Screen"] == "in vitro","Adduct"]), decreasing = TRUE), col = ColorMatrixTryOut["odd",2], las = 2, ylab = "Observed number of LTP - lipid species pairs", main = "Distribution of adducts in vitro")

dev.off()
par(mar = c(5.3, 4.1, 4.1, 2))

# Combined visualization of lipid species as different adduct or in different screen
pdf("./Output/LipidomicsQualityControlLTPLipidSpeciesCombinationsDistribution20012022.pdf")

barplot(AdductTabelScreensCombined, col = "Grey", main = "Number of LTP - lipid species observations \n as different adduct or in different screen", las = 1, xlab = "Number of observations as adducts or in screens", ylab = "Number of LTPs - lipid species")
abline(h = seq(from = 50, to = 350, by = 50), col = "#FFFFFF33", lwd = 3.2)

dev.off()


######## Figure 5
#### Fig.5a barchart on the top

#. TemporaryNewBarplotsForLipidomeOverviewGraph26042020.pdf #



# Carb distribution

CarbDistribution2403220 <- Col1ToRowNames(merge(aggregate(PureAntonella32b$Intensity, by = list(PureAntonella32b$TotalCarbonChainLength), FUN = sum), 
                                                aggregate(PureEnric32$Intensity, by = list(PureEnric32$TotalCarbonChainLength), FUN = sum), by = "Group.1", all = TRUE))

CarbDistribution2403220[setdiff(as.character(16:70), rownames(CarbDistribution2403220)),] <- 0
CarbDistribution2403220[is.na(CarbDistribution2403220)] <- 0

CarbDistribution2403220 <- CarbDistribution2403220[as.character(16:70),]
colnames(CarbDistribution2403220) <- c("in vivo", "in vitro")

ColorVectora <- c("white", brewer.pal(9,"Blues"))
ColorVectore <- c("white", brewer.pal(9,"Oranges"))


pdf("./Output/TemporaryNewBarplotsForLipidomeOverviewGraph26042020b.pdf")

barplot(t(CarbDistribution2403220), main="", col=c(ColorVectora[8],ColorVectore[8]), beside= TRUE, las = 3, ylim = c(0,50000000))
abline(h = 40000000)

barplot(t(CarbDistribution2403220), main="", col=c(ColorVectora[8],ColorVectore[8]), beside=TRUE, las = 3, ylim = c(0,900000000))
abline(h = 40000000)

abline(h = 120000000, col = "red")


barplot(t(CarbDistribution2403220), main="", col=c(ColorVectora[8],ColorVectore[8]), beside=TRUE, las = 3, ylim = c(0,4000000000), border = NA)
abline(h = 200000000, col = "red")

barplot(t(CarbDistribution2403220), main="", col=c(ColorVectora[8],ColorVectore[8]), beside=TRUE, las = 3, ylim = c(0,4000000000), border = NA)
abline(h = 120000000, col = "red")

barplot(t(CarbDistribution2403220), main="", col=c(ColorVectora[8],ColorVectore[8]), beside=TRUE, las = 3, ylim = c(0,2000000000))
abline(h = 200000000, col = "red")

barplot(t(CarbDistribution2403220), main="", col=c(ColorVectora[8],ColorVectore[8]), beside=TRUE, las = 3, ylim = c(0,2000000000))
abline(h = 120000000, col = "red")

barplot(t(CarbDistribution2403220), main="", col=c(ColorVectora[8],ColorVectore[8]), beside=TRUE, las = 3, ylim = c(0,200000000), border = NA)
abline(h = 120000000, col = "red")

barplot(t(CarbDistribution2403220), main="", col=c(ColorVectora[8],ColorVectore[8]), beside=TRUE, las = 3, ylim = c(0,200000000), border = FALSE)
abline(h = 120000000, col = "red")

barplot(t(CarbDistribution2403220), main="", col=c(ColorVectora[8],ColorVectore[8]), beside=TRUE, las = 3, ylim = c(0,200000000))
abline(h = 120000000, col = "red")

dev.off()
# Previous figures: Different variants used to add upper part to figure 5a further in Adobe Illustrator

#### Extended Data Figure 2b: Donut-diagram
#. Panel2CColorVsColorSidesFlippedOver09122021CleanedRightText10122021WithBiggerText.pdf (#)

#. Panel2CColorVsColorSidesFlippedOver09122021.pdf #
# InVivoDataSetslc 1195; InVitroDataSetslc 1196; ReorderedLTPsByManualSeriation 1448; LipidClassLiteratureDataSetslc 1258


ScreenColor <- cbind(c("InCellulo", "InVitro", "Cellular"), c("#9ECAE1", "#FDAE6B", "Grey"))

# LipidClassLiteratureDataSetslc2 <- LipidClassLiteratureDataSetslc
# LipidClassLiteratureDataSetslc2["PIPs","OSBPL2"] <- FALSE
# 
# NovelAmountsSubclassLTPConnections2 <- cbind(InCellulo = c(NovelAmount = sum((InVivoDataSetslc[,ReorderedLTPsByManualSeriation] != 0) & (LipidClassLiteratureDataSetslc2[,ReorderedLTPsByManualSeriation])),
#                                                            TotalAmount = sum(InVivoDataSetslc[,ReorderedLTPsByManualSeriation] != 0)),
#                                              
#                                              InVitro = c(NovelAmount = sum((InVitroDataSetslc[,ReorderedLTPsByManualSeriation] != 0) & (LipidClassLiteratureDataSetslc2[,ReorderedLTPsByManualSeriation])),
#                                                          TotalAmount = sum(InVitroDataSetslc[,ReorderedLTPsByManualSeriation] != 0)))
# 
# CompiledCharacteristicsOfNovelAmountsSubclassLTPConnections2 <- as.data.frame(do.call("cbind", list(Screen = c("in cellulo", "in vitro"), 
#                                                                                                     t(NovelAmountsSubclassLTPConnections2), 
#                                                                                                     
#                                                                                                     Percentage = round(NovelAmountsSubclassLTPConnections2["NovelAmount",]*100/NovelAmountsSubclassLTPConnections2["TotalAmount",]),
#                                                                                                     Color = ScreenColor[1:2,2])))

SumIntensitiesEvenVsOdd <- cbind(c(sum(PureAntonella32b[PureAntonella32b$TotalCarbonChainLength %% 2 == 0,"Intensity"], na.rm = TRUE),
                                   sum(PureAntonella32b[PureAntonella32b$TotalCarbonChainLength %% 2 == 1,"Intensity"], na.rm = TRUE)),
                                 
                                 c(sum(PureEnric32[PureEnric32$TotalCarbonChainLength %% 2 == 0,"Intensity"], na.rm = TRUE),
                                   sum(PureEnric32[PureEnric32$TotalCarbonChainLength %% 2 == 1,"Intensity"], na.rm = TRUE)))

colnames(SumIntensitiesEvenVsOdd) <- c("in cellulo", "in vitro")
rownames(SumIntensitiesEvenVsOdd) <- c("even", "odd")

ProcentIntensitiesEvenVsOdd <- cbind(SumIntensitiesEvenVsOdd[,1]*100 / sum(SumIntensitiesEvenVsOdd[,1]), SumIntensitiesEvenVsOdd[,2]*100 / sum(SumIntensitiesEvenVsOdd[,2]))
colnames(ProcentIntensitiesEvenVsOdd) <- c("in cellulo", "in vitro")

ProcentIntensitiesEvenVsOdd <- rbind(ProcentIntensitiesEvenVsOdd, color = as.character(ScreenColor[1:2,2])) # Simplified to avoid jumping back-and-forth. # CompiledCharacteristicsOfNovelAmountsSubclassLTPConnections2$Color


library(RColorBrewer)
pdf("./Output/Panel2CColorVsColorSidesFlippedOver09122021.pdf")

ColorMatrixTryOut <- rbind(even = c('in cellulo' = colorRampPalette(brewer.pal(9,"Blues")[-1])(9)[7], 'in vitro' = colorRampPalette(brewer.pal(9,"Oranges")[-1])(9)[7]),
                           odd = c('in cellulo' = colorRampPalette(brewer.pal(9,"Blues")[-1])(9)[5], 'in vitro' = colorRampPalette(brewer.pal(9,"Oranges")[-1])(9)[5]))

library(circlize)
circos.clear()

circos.par("start.degree" = 90, cell.padding = c(0, 0, 0, 0))
circos.initialize("a", xlim = c(0, 100))

circos.track(ylim = c(0.5, dim(ProcentIntensitiesEvenVsOdd)[2]+0.5), track.height = 0.8, 
             bg.border = NA, panel.fun = function(x, y) {
               
               xlim = CELL_META$xlim
               
               
               circos.rect(rep(0, 2), 2:1 - 0.45, as.numeric(as.character(ProcentIntensitiesEvenVsOdd["even",])), 2:1 + 0.45,
                           col = ColorMatrixTryOut["even",], border = ColorMatrixTryOut["even",], lwd = 5)
               
               circos.rect(rep(100, 2), 2:1 - 0.45, as.numeric(as.character(ProcentIntensitiesEvenVsOdd["even",])), 2:1 + 0.45,
                           col = ColorMatrixTryOut["odd",], border = ColorMatrixTryOut["even",], lwd = 5)
               
               
               circos.text(rep(xlim[1], 2), 2:1, 
                           
                           paste0(colnames(ProcentIntensitiesEvenVsOdd), 
                                  ": ", 
                                  
                                  as.character(round(as.numeric(as.character(ProcentIntensitiesEvenVsOdd["even",])), digits = 1)), 
                                  "% even"
                                  
                           ),
                           
                           
                           facing = "downward", adj = c(1.05, 0.5), cex = 0.8) 
               breaks = seq(0, 85, by = 5)
               
             })
dev.off()

# Previous figure: basis for further clean-up in Adobe Illustrator: Especially text has big increased in size, updated in font, and white background added with borders


#### Fig.5a circle-legend
#. CirclesForTheLegendOfTheLipidomeGraph09102020.pdf #

# Make complex heatmap with the different sizes of circles to make a legend
library(ComplexHeatmap)

BlueHighlights <- structure(do.call("rbind", list(0:10,NA,NA)), dimnames = list(c("in vivo", "in vitro", "HEK293"), as.character(0:10)))
OrangeHighlights <- structure(do.call("rbind", list(NA,0:10,NA)), dimnames = list(c("in vivo", "in vitro", "HEK293"), as.character(0:10)))

GreyHighlights <- structure(do.call("rbind", list(NA,NA,0:10)), dimnames = list(c("in vivo", "in vitro", "HEK293"), as.character(0:10)))


library(circlize)
col_funb <- colorRamp2(0:10, c("white", colorRampPalette(brewer.pal(9,"Blues")[-1])(10)))

col_funo <- colorRamp2(0:10, c("white", colorRampPalette(brewer.pal(9,"Oranges")[-1])(10)))
col_fung <- colorRamp2(0:10, c("white", colorRampPalette(brewer.pal(9,"Greys")[2:7])(10)))

LegendColor <- col_fung


pdf("./Output/CirclesForTheLegendOfTheLipidomeGraph09102020.pdf",
    width = unit(23, "mm"), height = unit(8, "mm"))

Heatmap(BlueHighlights, name = "LegendCircles", col = LegendColor, rect_gp = gpar(type = "none"),  
        width = unit(100*dim(BlueHighlights)[2]/min(dim(BlueHighlights)), "mm"), height = unit(100*dim(BlueHighlights)[1]/min(dim(BlueHighlights)), "mm"), show_heatmap_legend = FALSE, column_title = "", 
        
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "LightGrey", fill = NA, lwd = 1, lwd = 10))
          
          grid.circle(x = x, y = y, r = abs(GreyHighlights[i, j])/5.5 * min(unit.c(width, height)), gp = gpar(fill = NA, col = "LightGrey", lwd = 10, lty = 1))
          grid.circle(x = x, y = y, r = abs(OrangeHighlights[i, j])/5.5 * min(unit.c(width, height)), gp = gpar(fill = NA, col = col_funo(7), lwd = 10))
          
          grid.circle(x = x, y = y, r = abs(BlueHighlights[i, j])/5.5 * min(unit.c(width, height)), gp = gpar(fill = NA, col = col_funb(7), lwd = 10))
        }, cluster_rows = FALSE, cluster_columns = FALSE)

dev.off()
# Previous figure is basis for circle legend. Adobe Illustrator used to integrate & esthetically clean up.

#### Extended Data Figure 2b: Histograms at side
#. LipidsIntensityFrationsOfEven26012022WithLysolipidsCollapsed27012022CleanedUpVersionWithTwoColorSchemesOddAndNoReducedToEmpty1403202222binsc.pdf (#)

#. LipidsIntensityFrationsOfEven26012022WithLysolipidsCollapsed27012022 #


# Function to clean up data. Updated with widen5 function.
CovertToHeatmapMatrixn <- function(InputDataStartMatrix, CompactionMatrixForTotal, HeadgroupOrder, CarbOrder, IonMode){
  
  reshapedstartofmatrixx <- widen5(inputdf = InputDataStartMatrix[InputDataStartMatrix$IonMode == IonMode,], ColumnsLong = "LikelySubclass", ColumnWide = "TotalCarbonChainLength", ColumnValue = "Intensity", AggregatingFunction = sum, FunctionOutputValueType = double(1), FillerValue = 0, RowNamesColumnsLong = FALSE)
  reshapedstartofmatrix2x <- reshapedstartofmatrixx[(reshapedstartofmatrixx$LikelySubclass != "VE" & reshapedstartofmatrixx$LikelySubclass != "VA" & reshapedstartofmatrixx$LikelySubclass != "P40"), colnames(reshapedstartofmatrixx) != "NaN"]
  
  
  GenHeadx <- CompactionMatrixForTotal[match(reshapedstartofmatrix2x$LikelySubclass, CompactionMatrixForTotal[,1]),2] # Changed CompactionMatrixForTotaln2 To CompactionMatrixForTotal #
  
  reshapedstartofmatrix4x <- rowsum(reshapedstartofmatrix2x[,-1], GenHeadx)
  reshapedstartofmatrix4x[,setdiff(CarbOrder, colnames(reshapedstartofmatrix4x))] <- 0
  
  reshapedstartofmatrix4x[setdiff(HeadgroupOrder, rownames(reshapedstartofmatrix4x)),] <- 0
  return(reshapedstartofmatrix4x[HeadgroupOrder, CarbOrder])
  
}


# Simplified and eliminated need for external package, to increase robustness.
# Also made the input cleaner by putting script snippet in function

# Local function to add a column with MinMaxRangeNorm to the dataframe (rownames have different mode (now: chr; before: numeric), but are otherwise the same as before)
AddMinMaxNormColumn <- function(y){cbind.data.frame(y, MinMaxRangeNorm = sapply(1:dim(y)[1], function(x){MinMaxNormFuncn(y$Intensity[x],y$IonMode[x],y)}))}

Normalization1OfPureAntonella32 <- AddMinMaxNormColumn(PureAntonella32b) # Changed from PureAntonella32 to PureAntonella32b 
Normalization1OfPureEnric32 <- AddMinMaxNormColumn(PureEnric32)

Normalization1OfPureAntonella32$LikelySubclass <- as.character(Normalization1OfPureAntonella32$LikelySubclass)
Normalization1OfPureEnric32$LikelySubclass <- as.character(Normalization1OfPureEnric32$LikelySubclass)

CompactionMatrixForTotaln <- cbind(unique(c(Normalization1OfPureAntonella32$LikelySubclass, Normalization1OfPureEnric32$LikelySubclass)), unique(c(Normalization1OfPureAntonella32$LikelySubclass, Normalization1OfPureEnric32$LikelySubclass)))
CompactionMatrixForTotaln[c(2,8:14,17,18,19,21:23,29,30,31),2] <- c("PC", "HexCer", "HexCer", "Hex2Cer", "SM", "CerP", "SM", "SM", "PE", "Cer", "Cer", "Cer", "Cer", "Cer", "FA", "SHexCer", "LPE") # Adapted to new order after change of y-input.

CompactionMatrixForTotaln2 <- CompactionMatrixForTotaln[CompactionMatrixForTotaln[,1] != "VA",]
HeadgroupOrdern <- c("Cer", "CerP", "HexCer", "Hex2Cer", "SHexCer", "SM", "FA", "LPC", "LPE", "LPG", "DAG", "PA", "PC", "PE", "PI", "PS", "PG", "PG/BMP", "BMP", "TAG", "PGP", "CL") 

HeatmapMatrixPosAntonellan <- CovertToHeatmapMatrixn(CompactionMatrixForTotal = CompactionMatrixForTotaln2,
                                                     HeadgroupOrder = HeadgroupOrdern,
                                                     
                                                     InputDataStartMatrix = PureAntonella32b,
                                                     CarbOrder = as.character(16:70),
                                                     
                                                     IonMode = "pos")


HeatmapMatrixNegAntonellan <- CovertToHeatmapMatrixn(CompactionMatrixForTotal = CompactionMatrixForTotaln2,
                                                     HeadgroupOrder = HeadgroupOrdern,
                                                     
                                                     InputDataStartMatrix = PureAntonella32b,
                                                     CarbOrder = as.character(16:70),
                                                     
                                                     IonMode = "neg")


HeatmapMatrixPosEnricn <- CovertToHeatmapMatrixn(CompactionMatrixForTotal = CompactionMatrixForTotaln2,
                                                 HeadgroupOrder = HeadgroupOrdern,
                                                 
                                                 InputDataStartMatrix = PureEnric32,
                                                 CarbOrder = as.character(16:70),
                                                 
                                                 IonMode = "pos")


HeatmapMatrixNegEnricn <- CovertToHeatmapMatrixn(CompactionMatrixForTotal = CompactionMatrixForTotaln2,
                                                 HeadgroupOrder = HeadgroupOrdern,
                                                 
                                                 InputDataStartMatrix = PureEnric32,
                                                 CarbOrder = as.character(16:70),
                                                 
                                                 IonMode = "neg")


HeatmapMatrixBothAntonellan <- HeatmapMatrixNegAntonellan + HeatmapMatrixPosAntonellan
HeatmapMatrixBothEnricn <- HeatmapMatrixNegEnricn + HeatmapMatrixPosEnricn

HeatmapMatrixBothAntonellancwl <- Col1ToRowNames(aggregate(x = HeatmapMatrixBothAntonellan, 
                                                           by = list(c(rownames(HeatmapMatrixBothAntonellan)[1:7], c("PC","PE","PG"), rownames(HeatmapMatrixBothAntonellan)[-(1:10)])), FUN = function(x){sum(x, na.rm = TRUE)}))[c(2,3,8,7,17,18,6,9,10,11,15,16,14,12,13,1,4,5,19),]

HeatmapMatrixBothEnricncwl <- Col1ToRowNames(aggregate(x = HeatmapMatrixBothEnricn, 
                                                       by = list(c(rownames(HeatmapMatrixBothEnricn)[1:7], c("PC","PE","PG"), rownames(HeatmapMatrixBothEnricn)[-(1:10)])), FUN = function(x){sum(x, na.rm = TRUE)}))[c(2,3,8,7,17,18,6,9,10,11,15,16,14,12,13,1,4,5,19),]


EvenPercentagesForLipidscwl <- cbind("in cellulo" = rowSums(HeatmapMatrixBothAntonellancwl[, as.numeric(colnames(HeatmapMatrixBothAntonellancwl)) %% 2 == 0], na.rm = TRUE)*100 /
                                       
                                       (rowSums(HeatmapMatrixBothAntonellancwl[, as.numeric(colnames(HeatmapMatrixBothAntonellancwl)) %% 2 == 0], na.rm = TRUE) +
                                          rowSums(HeatmapMatrixBothAntonellancwl[, as.numeric(colnames(HeatmapMatrixBothAntonellancwl)) %% 2 == 1], na.rm = TRUE)),                   
                                     
                                     
                                     "in vitro" = rowSums(HeatmapMatrixBothEnricncwl[, as.numeric(colnames(HeatmapMatrixBothEnricncwl)) %% 2 == 0], na.rm = TRUE)*100 /
                                       
                                       (rowSums(HeatmapMatrixBothEnricncwl[, as.numeric(colnames(HeatmapMatrixBothEnricncwl)) %% 2 == 0], na.rm = TRUE) +
                                          rowSums(HeatmapMatrixBothEnricncwl[, as.numeric(colnames(HeatmapMatrixBothEnricncwl)) %% 2 == 1], na.rm = TRUE))) 

EvenPercentagesForLipidsEmptycwl <- EvenPercentagesForLipidscwl
EvenPercentagesForLipidsEmptycwl[TRUE] <- 100

EvenPercentagesForLipidsAbsencescwl <- EvenPercentagesForLipidscwl
EvenPercentagesForLipidsAbsencescwl <- is.na(EvenPercentagesForLipidsAbsencescwl)*100

par(mar = c(5.1, 4.1, 4.1, 2.1))
pdf("./Output/LipidsIntensityFrationsOfEven26012022WithLysolipidsCollapsed27012022.pdf")

barplot(t(EvenPercentagesForLipidsEmptycwl[rev(rownames(EvenPercentagesForLipidsEmptycwl)),2:1]), beside = TRUE, horiz = TRUE, las = 2, col = "white", border = "grey", xlim = c(0,1000), xaxt = "n") 
barplot(t(EvenPercentagesForLipidsAbsencescwl[rev(rownames(EvenPercentagesForLipidsAbsencescwl)),2:1]), add = TRUE, beside = TRUE, horiz = TRUE, las = 2, col = "grey", border = "grey", xlim = c(0,1000), xaxt = "n", yaxt = "n") 

barplot(t(EvenPercentagesForLipidscwl[rev(rownames(EvenPercentagesForLipidscwl)),2:1]), add = TRUE, beside = TRUE, horiz = TRUE, las = 2, col = ColorMatrixTryOut["odd",2:1], xlim = c(0,1000), xaxt = "n", yaxt = "n") 
axis(side = 1, at = c(0,25,50,75,100))

dev.off()
# Previous figure: basis for further clean-up in Adobe Illustrator: especially rastering and color schemes


#### Fig.5a: circle heatmap visualization

#. CirclesToVisualizeLipidomeAndTheBackgroundEmptyCirclesForInVivoData150620205c2 #
#. CirclesToVisualizeLipidomeAndTheBackgroundEmptyCirclesForInVivoData150620205c2DifferentGreyColorsAndCircleOrderGridToBack2SameLineTicknessCombinedWorkDocument4PGPAbovePGGappedAndCollapsed4add.pdf (#)

library(RColorBrewer)
library(ComplexHeatmap)


# Import and cleaning of shotgun lipidomics data for HEK293 cells used

OverexpressionHEK <- read.csv(file = "./InputData/ExtractedDataHEK293FromKenji08082019.txt", header = TRUE, sep = "\t", as.is = TRUE)
OverexpressionHEK2 <- OverexpressionHEK[-(1:2),]

colnames(OverexpressionHEK2) <- paste(colnames(OverexpressionHEK), OverexpressionHEK[1,], OverexpressionHEK[2,], sep = "_")
rownames(OverexpressionHEK2) <- 1:nrow(OverexpressionHEK2)

OverexpressionHEK4 <- apply(OverexpressionHEK2, 2, function(x){gsub(",",".",x)})
OverexpressionHEK5 <- cbind.data.frame(OverexpressionHEK4[,1], apply(OverexpressionHEK4[,-1], 2, function(x){as.numeric(x)}))

OverexpressionHEK5[,1] <- as.character(OverexpressionHEK5[,1])
DetailsOfOverexpressionHEKData2 <- OverexpressionHEK5[1:485,]

DetailsOfOverexpressionHEK2Split <- cbind(do.call("rbind", strsplit(DetailsOfOverexpressionHEKData2[,1], split = "(?<=[a-zA-Z]|-)\\s*(?=[0-9])", perl = TRUE)), DetailsOfOverexpressionHEKData2[,-1])
colnames(DetailsOfOverexpressionHEK2Split)[1:2] <- c("LipidClass", "CarbChainLengthAndUnsat")

DetailsOfOverexpressionHEK2Split[,1] <- as.character(DetailsOfOverexpressionHEK2Split[,1])
DetailsOfOverexpressionHEK2Split[,2] <- as.character(DetailsOfOverexpressionHEK2Split[,2])

DetailsOfOverexpressionHEK2Split[DetailsOfOverexpressionHEK2Split[,2] == "Chol :", 2] <- "27:1"
DetailsOfOverexpressionHEK2Split[DetailsOfOverexpressionHEK2Split[,1] == "GM", 2] <- gsub("3 ", "", DetailsOfOverexpressionHEK2Split[DetailsOfOverexpressionHEK2Split[,1] == "GM", 2])

DetailsOfOverexpressionHEK2Split[DetailsOfOverexpressionHEK2Split[,1] == "GM", 1] <- "GM3"
DetailsOfOverexpressionHEK2Split[grep(";2", DetailsOfOverexpressionHEK2Split[,2]),2] <- gsub(";2", "", DetailsOfOverexpressionHEK2Split[grep(";2", DetailsOfOverexpressionHEK2Split[,2]),2])

DetailsOfOverexpressionHEK2Split2 <- cbind(DetailsOfOverexpressionHEK2Split[,1:2], ControlMeans = rowMeans(DetailsOfOverexpressionHEK2Split[,c(4,10,16)])) 


LipidConversionsBL <- cbind(unique(DetailsOfOverexpressionHEK2Split2$LipidClass), c("LPA", "LPC", "LPC", "LPE", "LPE", "LPG", "LPG", "LPI", "LPI", "LPS", "LPS", "GM3", "Chol", "SM", "PC", "PC", "HexCer", "CE", "DAG", "Cer", "Hex2Cer", "PE", "PE", "PG", "PI", "PI", "PS", "SHexCer", "CerP", "CL", "PA", "PA"))
DetailsOfOverexpressionHEK2Split2bl <- cbind(DetailsOfOverexpressionHEK2Split2, LipidConversionsBL[match(DetailsOfOverexpressionHEK2Split2$LipidClass, LipidConversionsBL[,1]),2])

DetailsOfOverexpressionHEK2Split2bl2 <- cbind(DetailsOfOverexpressionHEK2Split2bl, do.call("rbind", strsplit(x = DetailsOfOverexpressionHEK2Split2bl$CarbChainLengthAndUnsat, split = ":")))
colnames(DetailsOfOverexpressionHEK2Split2bl2)[4:6] <- c("Headgroup", "ChainLength", "Unsaturation")

# Updated to widen version with different rownames order, but otherwise identical
DetailsOfOverexpressionHEK2Split2bl2WideVersion <- Col1ToRowNames(widen5(inputdf = DetailsOfOverexpressionHEK2Split2bl2, ColumnsLong = "Headgroup", ColumnWide = "ChainLength", ColumnValue = "ControlMeans", AggregatingFunction = sum, FunctionOutputValueType = double(1)))

# Historical intermediate parts removed
DetailsOfOverexpressionHEK2Split2bl2WideVersion2 <- DetailsOfOverexpressionHEK2Split2bl2WideVersion

DetailsOfOverexpressionHEK2Split2bl2WideVersion2[,setdiff(as.character(14:72), colnames(DetailsOfOverexpressionHEK2Split2bl2WideVersion2))] <- 0
DetailsOfOverexpressionHEK2Split2bl2WideVersion2[setdiff(HeadgroupOrdern, rownames(DetailsOfOverexpressionHEK2Split2bl2WideVersion2)),] <- 0

DetailsOfOverexpressionHEK2Split2bl2WideVersion2 <- DetailsOfOverexpressionHEK2Split2bl2WideVersion2[HeadgroupOrdern, as.character(14:72)]
HeatmapMatrixPosBackgroundNormalizedtb <- MinMaxNormMatrixFunc(inputdata = DetailsOfOverexpressionHEK2Split2bl2WideVersion2)*9+1

HeatmapMatrixPosBackgroundNormalizedtb[is.infinite(as.matrix(HeatmapMatrixPosBackgroundNormalizedtb))] <- 0
TestMatrixb <- as.matrix(HeatmapMatrixPosBackgroundNormalizedtb)

HeatmapMatrixPosAntonellaNormalizedn <- MinMaxNormMatrixFunc(inputdata = HeatmapMatrixPosAntonellan)*9+1
HeatmapMatrixPosAntonellaNormalizedn[is.infinite(as.matrix(HeatmapMatrixPosAntonellaNormalizedn))] <- 0

HeatmapMatrixNegAntonellaNormalizedn <- MinMaxNormMatrixFunc(inputdata = HeatmapMatrixNegAntonellan)*9+1
HeatmapMatrixNegAntonellaNormalizedn[is.infinite(as.matrix(HeatmapMatrixNegAntonellaNormalizedn))] <- 0

HeatmapMatrixPosEnricNormalizedn <- MinMaxNormMatrixFunc(inputdata = HeatmapMatrixPosEnricn)*9+1
HeatmapMatrixPosEnricNormalizedn[is.infinite(as.matrix(HeatmapMatrixPosEnricNormalizedn))] <- 0

HeatmapMatrixNegEnricNormalizedn <- MinMaxNormMatrixFunc(inputdata = HeatmapMatrixNegEnricn)*9+1
HeatmapMatrixNegEnricNormalizedn[is.infinite(as.matrix(HeatmapMatrixNegEnricNormalizedn))] <- 0

HeatmapMatrixAntonellaNormalizedn <- do.call("cbind", lapply(1:dim(HeatmapMatrixNegAntonellaNormalizedn)[2], function(x){pmax(HeatmapMatrixNegAntonellaNormalizedn[,x], HeatmapMatrixPosAntonellaNormalizedn[,x])}))
HeatmapMatrixEnricNormalizedn <- do.call("cbind", lapply(1:dim(HeatmapMatrixNegEnricNormalizedn)[2], function(x){pmax(HeatmapMatrixNegEnricNormalizedn[,x], HeatmapMatrixPosEnricNormalizedn[,x])}))

rownames(HeatmapMatrixAntonellaNormalizedn) <- rownames(HeatmapMatrixPosAntonellaNormalizedn)
colnames(HeatmapMatrixAntonellaNormalizedn) <- colnames(HeatmapMatrixPosAntonellaNormalizedn)

rownames(HeatmapMatrixEnricNormalizedn) <- rownames(HeatmapMatrixPosEnricNormalizedn)
colnames(HeatmapMatrixEnricNormalizedn) <- colnames(HeatmapMatrixPosEnricNormalizedn)

HeatmapMatrixAntonellaNormalizednx2 <- as.data.frame(HeatmapMatrixAntonellaNormalizedn)
HeatmapMatrixEnricNormalizednx2 <- as.data.frame(HeatmapMatrixEnricNormalizedn)

HeatmapMatrixAntonellaNormalizednx2[,setdiff(as.character(14:72), colnames(HeatmapMatrixAntonellaNormalizednx2))] <- 0
HeatmapMatrixEnricNormalizednx2[,setdiff(as.character(14:72), colnames(HeatmapMatrixEnricNormalizednx2))] <- 0

TestMatrixa2 <- as.matrix(HeatmapMatrixAntonellaNormalizednx2[,as.character(14:72)])
TestMatrixe2 <- as.matrix(HeatmapMatrixEnricNormalizednx2[,as.character(14:72)])

library(circlize)
col_funb <- colorRamp2(0:10, c("white", colorRampPalette(brewer.pal(9,"Blues")[-1])(10)))

col_funo <- colorRamp2(0:10, c("white", colorRampPalette(brewer.pal(9,"Oranges")[-1])(10)))
col_fung <- colorRamp2(0:10, c("white", colorRampPalette(brewer.pal(9,"Greys")[2:7])(10)))

LegendName <- "Legend"
LegendColor <- col_fung

GraphName <- "Empty circles for in vivo"


pdf("./Output/CirclesToVisualizeLipidomeAndTheBackgroundEmptyCirclesForInVivoData150620205c2.pdf",
    width = unit(14, "mm"), height = unit(23, "mm"))

Heatmap(TestMatrixa2, name = LegendName, col = LegendColor, rect_gp = gpar(type = "none"),  
        width = unit(100*dim(TestMatrixa2)[2]/min(dim(TestMatrixa2)), "mm"), height = unit(100*dim(TestMatrixa2)[1]/min(dim(TestMatrixa2)), "mm"), show_heatmap_legend = FALSE, column_title = GraphName, 
        
        
        cell_fun = function(j, i, x, y, width, height, fill) {
          
          grid.circle(x = x, y = y, r = abs(TestMatrixa2[i, j])/8.5 * min(unit.c(width, height)), gp = gpar(fill = NA, col = col_funb(7), lwd = 2.3))
          grid.circle(x = x, y = y, r = abs(TestMatrixe2[i, j])/8.5 * min(unit.c(width, height)), gp = gpar(fill = NA, col = col_funo(7), lwd = 2))
          
          grid.circle(x = x, y = y, r = abs(TestMatrixb[i, j])/8.5 * min(unit.c(width, height)), gp = gpar(fill = NA, col = "LightGrey", lwd = 1.4, lty = 1))
          
          
          grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "LightGrey", fill = NA, lwd = 1))
        }, cluster_rows = FALSE, cluster_columns = FALSE)

dev.off()
# Previous figure: Basis for putting rest of figure 5a together, with shades of grey changed in Adobe Illustrator afterwards, and PGP moved in location & gaps between groups

#### #X Extra interesting visualization, not part of included figures anymore (previously: Fig.2c)
#. ComparisonOfUnsaturationPreferencesInCelluloVsInVitroVsCells140420222.pdf #

# Import of the HEK293 shotgun lipidomics data
LTPBackground1_2 <- read.csv(file = "./InputData/tableauOutput_20181023B_version2.txt", header = TRUE, sep = "\t", as.is = TRUE)

# Translate lipids manually outside R
write.table(LTPBackground1_2, file="./Output/LTPBackground1_217062021.tsv", sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Import again 
LipidSubclassesAddedToBackground17062021 <- read.csv(file = "./InputData/LipidSubclassesAddedToBackground17062021.txt", header = TRUE, sep = "\t", as.is = TRUE)

LipidSubclassesAddedToBackground170620212 <- cbind(AdaptedLipids = paste0(LipidSubclassesAddedToBackground17062021[,1], 
                                                                          "(", 
                                                                          
                                                                          gsub("O-","",sapply(strsplit(sapply(strsplit(LipidSubclassesAddedToBackground17062021[,2],";"), "[", 1)," "), "[", 2)),
                                                                          ")"),
                                                   
                                                   LipidSubclassesAddedToBackground17062021)


LipidSubclassesAddedToBackground170620212[,1] <- sapply(LipidSubclassesAddedToBackground170620212[,1], function(x){if(grepl("d\\*", x)){gsub("\\(", "\\(d\\*", gsub("d\\*", "", x))
}else{
  
  if(grepl("-O", x)){gsub("-O\\(","\\(O-",x)}else{as.character(x)}
}}) # Added as missing piece

SpeciesOverviewLipidBodyConnections <- as.data.frame(cbind(do.call("rbind", strsplit(gsub("\\)", "", LipidSubclassesAddedToBackground170620212[4:422,1]), split = "\\(")), gsub(",", ".", LipidSubclassesAddedToBackground170620212[4:422,4])))
colnames(SpeciesOverviewLipidBodyConnections) <- c("Head", "Body", "Amount")

SpeciesOverviewLipidBodyConnections$Amount <- as.numeric(as.character(SpeciesOverviewLipidBodyConnections$Amount)) # Added later as correction


IntensitiesOfObservedVsExpectedUnsatLipids <- setNames(Col1ToRowNames(merge(aggregate(PureAntonella32b$Intensity, by = list(PureAntonella32b$TotalCarbonChainUnsaturations), FUN = function(x){sum(x,na.rm=TRUE)}), 
                                                                            aggregate(SpeciesOverviewLipidBodyConnections$Amount, by = list(as.numeric(sapply(strsplit(as.character(SpeciesOverviewLipidBodyConnections$Body), ":"), "[[", 2))), FUN = function(x){sum(x,na.rm=TRUE)}), 
                                                                            
                                                                            by = 1, all = TRUE)),
                                                       nm = c("Observed", "Expected"))


InVitroSummedIntensitiesUnsaturations <- aggregate(PureEnric32$Intensity, by = list(PureEnric32$TotalCarbonChainUnsaturations), FUN = function(x){sum(x, na.rm = TRUE)})

IntensitiesUnsatLevelsScreensComparison <- setNames(Col1ToRowNames(rbind(merge(IntensitiesOfObservedVsExpectedUnsatLipids, InVitroSummedIntensitiesUnsaturations, by.x = 0, by.y = 1, all = TRUE), c(11,0,0,0)))[as.character(0:12),c(1,3,2)],
                                                    nm = c("In Cellulo", "In Vitro", "Cellular"))

IntensitiesUnsatLevelsScreensComparison[is.na(IntensitiesUnsatLevelsScreensComparison)] <- 0


IntensitiesUnsatLevelsScreensComparisonNormByMax <- apply(IntensitiesUnsatLevelsScreensComparison, 2, function(x){x*100/max(x,na.rm=TRUE)})
IntensitiesUnsatLevelsScreensComparisonNormBySum <- apply(IntensitiesUnsatLevelsScreensComparison, 2, function(x){x*100/sum(x,na.rm=TRUE)})


library(RColorBrewer)

ListLipidUnsatComp2 <- setNames(list(IntensitiesUnsatLevelsScreensComparisonNormByMax, IntensitiesUnsatLevelsScreensComparisonNormBySum), nm = c("Percentual fraction intensities of max", "Percentual fraction intensities of total"))
pdf("./Output/ComparisonOfUnsaturationPreferencesInCelluloVsInVitroVsCells140420222.pdf")

for(x in names(ListLipidUnsatComp2)){
  y <- ListLipidUnsatComp2[[x]]
  
  barplot(t(y), beside = TRUE, col = c(brewer.pal(9,"Blues")[5], brewer.pal(9,"Oranges")[5], "grey"), las = 1, xlab = "", ylab = "") #, xaxt = "n", yaxt="n"
  abline(h = seq(5, ceiling(max(y)/10)*10, 5), col = "#FFFFFF52", lwd = 3.2)
  
  title(xlab = "Total unsaturation of lipids", line = 2.5)
  title(ylab = x, line = 2.5)
  
  legend("topright", legend = c("In Cellulo", "In Vitro", "Cellular"), col = c(brewer.pal(9,"Blues")[5], brewer.pal(9,"Oranges")[5], "grey"), pch = 15, bty = "n", cex = 0.95)
}

dev.off()
#X Interesting extra informative view. Second page of previous figure-set: Used for panel in figure after some adaptations in Adobe Illustrator such as orientation and color scheme to fit with the other panels

######## Figure 2
#### Fig.2a

#. Panel3A09122021.pdf (#)
#. HeatmapOfTheLTPLocationsWithoutLipidInformationWithoutDomainsAdded251120214hdrBlackTrianglesAddedlwd1ZOrderEnhancedAndLegendAdded2.pdf (#)

#. (HeatmapOfTheLTPLocationsWithLipidInformationWithDomainsAdded11072020c.pdf # Domain and subcellular locations added: more extensive version that is not relevant for article, so not included)
#. HeatmapOfTheLTPLocationsWithoutLipidInformationWithoutDomainsAdded251120214hdrBlackTrianglesAddedlwd2 #


#### HPTLC-data import as dataframe with factors in columns # Was in post-production of the figure updated by restricting some of these observations from the figures when uncertainties present.

HPTLCDataInVivoAndInVitro <- read.csv(file = "./InputData/KnownHPTLCResultsScraped31032020.tsv", header = FALSE, sep = "\t", as.is = FALSE)
colnames(HPTLCDataInVivoAndInVitro) <- c("LTPProtein", "Lipid", "Screen")


#### HPTLC-data: first clean-up and formatting # Was in post-production of the figure updated by restricting some of these observations from the figures when uncertainties present.

HPTLCDataInVivoAndInVitrob <- HPTLCDataInVivoAndInVitro[HPTLCDataInVivoAndInVitro[,1] != "SEC14L1",] # Removed because of a potential experimental issue
HPTLCDataInVivoAndInVitrob$LTPProtein <- factor(HPTLCDataInVivoAndInVitrob$LTPProtein, levels = levels(HPTLCDataInVivoAndInVitrob$LTPProtein)[levels(HPTLCDataInVivoAndInVitrob$LTPProtein) != "SEC14L1"])

# Updated with widen5 function
HPTLCSpecificitiesPerScreen <- lapply(c("A","E"), function(x){Col1ToRowNames(widen5(inputdf = HPTLCDataInVivoAndInVitrob[HPTLCDataInVivoAndInVitrob[,"Screen"] == x,], ColumnsLong = "LTPProtein", ColumnWide = "Lipid", ColumnValue = "Screen", AggregatingFunction = length, DropUnusedRowsLevels = FALSE, FunctionOutputValueType = integer(1)))})


HPTLCSpecificitiesPerScreen2 <- HPTLCSpecificitiesPerScreen

colnames(HPTLCSpecificitiesPerScreen2[[1]]) <- paste0(colnames(HPTLCSpecificitiesPerScreen2[[1]]), "*")
colnames(HPTLCSpecificitiesPerScreen2[[2]]) <- paste0(colnames(HPTLCSpecificitiesPerScreen2[[2]]), "*")

#### MS-data: clean-up and formatting
LTPProteins <- unique(c(unique(PureAntonella32b$LTPProtein), unique(PureEnric32$LTPProtein)))

HeadgroupOrderlnl <- c("d*Cer", "dCer", "DHCer", "DHOH*Cer", "tCer", "d*CerP", "d*HexCer", "t*HexCer", "t*Hex2Cer", "d*SHexCer", "d*SM", "DHSM", "t*SM", "FA", "FAL", "LPC",
                       "LPE", "LPE-O", "LPG", "DAG", "PA", "PC", "PC-O", "PE", "PE-O", "PI", "PS", "PG", "PG/BMP", "BMP", "TAG", "PGP", "CL", "VA")

# Updated with the widen5 function.
CovertToHeatmapMatrixlnl <- function(InputDataStartMatrixlnl, HeadgroupOrderlnl, LTPOrder, IonMode){ # CompactionMatrixForTotal, 
  
  reshapedstartofmatrixlnlx <- widen5(inputdf = InputDataStartMatrixlnl[InputDataStartMatrixlnl$IonMode == IonMode,], ColumnsLong = "LikelySubclass", ColumnWide = "LTPProtein", ColumnValue = "Intensity", AggregatingFunction = sum, FunctionOutputValueType = double(1), RowNamesColumnsLong = FALSE, FillerValue = 0)
  reshapedstartofmatrixlnl2x <- reshapedstartofmatrixlnlx[reshapedstartofmatrixlnlx$LikelySubclass != "P40", colnames(reshapedstartofmatrixlnlx) != "NaN"]
  
  reshapedstartofmatrixlnl4x <- Col1ToRowNames(reshapedstartofmatrixlnl2x)
  reshapedstartofmatrixlnl4x[,setdiff(LTPOrder, colnames(reshapedstartofmatrixlnl4x))] <- 0
  
  reshapedstartofmatrixlnl4x[setdiff(HeadgroupOrderlnl, rownames(reshapedstartofmatrixlnl4x)),] <- 0
  return(reshapedstartofmatrixlnl4x[HeadgroupOrderlnl, LTPOrder])
  
}


LTPMatrixPosInVivo <- CovertToHeatmapMatrixlnl(InputDataStartMatrixlnl = PureAntonella32b,
                                               HeadgroupOrderlnl = HeadgroupOrderlnl,
                                               
                                               LTPOrder = LTPProteins, 
                                               IonMode = "pos")

LTPMatrixNegInVivo <- CovertToHeatmapMatrixlnl(InputDataStartMatrixlnl = PureAntonella32b,
                                               HeadgroupOrderlnl = HeadgroupOrderlnl,
                                               
                                               LTPOrder = LTPProteins, 
                                               IonMode = "neg")

LTPMatrixPosInVitro <- CovertToHeatmapMatrixlnl(InputDataStartMatrixlnl = PureEnric32,
                                                HeadgroupOrderlnl = HeadgroupOrderlnl,
                                                
                                                LTPOrder = LTPProteins, 
                                                IonMode = "pos")

LTPMatrixNegInVitro <- CovertToHeatmapMatrixlnl(InputDataStartMatrixlnl = PureEnric32,
                                                HeadgroupOrderlnl = HeadgroupOrderlnl,
                                                
                                                LTPOrder = LTPProteins, 
                                                IonMode = "neg")

LTPMatrixPosInVivommn <- MinMaxNormMatrixFunc(inputdata = LTPMatrixPosInVivo)*9+1
LTPMatrixPosInVivommn[is.infinite(as.matrix(LTPMatrixPosInVivommn))] <- 0

LTPMatrixNegInVivommn <- MinMaxNormMatrixFunc(inputdata = LTPMatrixNegInVivo)*9+1
LTPMatrixNegInVivommn[is.infinite(as.matrix(LTPMatrixNegInVivommn))] <- 0

LTPMatrixPosInVitrommn <- MinMaxNormMatrixFunc(inputdata = LTPMatrixPosInVitro)*9+1
LTPMatrixPosInVitrommn[is.infinite(as.matrix(LTPMatrixPosInVitrommn))] <- 0

LTPMatrixNegInVitrommn <- MinMaxNormMatrixFunc(inputdata = LTPMatrixNegInVitro)*9+1 # Added in, because was missing.
LTPMatrixNegInVitrommn[is.infinite(as.matrix(LTPMatrixNegInVitrommn))] <- 0 # Added in, because was missing.

LTPMatrixTopInVivommn <- do.call("cbind", lapply(1:dim(LTPMatrixNegInVivommn)[2], function(x){pmax(LTPMatrixNegInVivommn[,x], LTPMatrixPosInVivommn[,x])}))
LTPMatrixTopInVitrommn <- do.call("cbind", lapply(1:dim(LTPMatrixNegInVitrommn)[2], function(x){pmax(LTPMatrixNegInVitrommn[,x], LTPMatrixPosInVitrommn[,x])}))

rownames(LTPMatrixTopInVivommn) <- rownames(LTPMatrixPosInVivommn)
colnames(LTPMatrixTopInVivommn) <- colnames(LTPMatrixPosInVivommn)

rownames(LTPMatrixTopInVitrommn) <- rownames(LTPMatrixPosInVitrommn)
colnames(LTPMatrixTopInVitrommn) <- colnames(LTPMatrixPosInVitrommn)

MeltedInVivoCombined <- rbind(meltmatrix(LTPMatrixTopInVivommn), meltmatrix(t(as.matrix(5*HPTLCSpecificitiesPerScreen2[[1]]))))
MeltedInVitroCombined <- rbind(meltmatrix(LTPMatrixTopInVitrommn), meltmatrix(t(as.matrix(5*HPTLCSpecificitiesPerScreen2[[2]]))))

CastInVivoCombined <- Col1ToRowNames(widen5(inputdf = MeltedInVivoCombined, ColumnsLong = "Var1", ColumnWide = "Var2", ColumnValue = "value", AggregatingFunction = sum, FunctionOutputValueType = double(1), SortRows = FALSE, FillerValue = 0))
CastInVitroCombined <- Col1ToRowNames(widen5(inputdf = MeltedInVitroCombined, ColumnsLong = "Var1", ColumnWide = "Var2", ColumnValue = "value", AggregatingFunction = sum, FunctionOutputValueType = double(1), SortRows = FALSE, FillerValue = 0))

# Make overview of LTP classes
MainDomainsOfTheLTPs <- unique(rbind(PureAntonella32b[,1:2], PureEnric32[,1:2])) 

MainDomainsOfTheLTPs2 <- cbind(LTPProtein = colnames(CastInVivoCombined), MainDomain = MainDomainsOfTheLTPs[match(colnames(CastInVivoCombined), MainDomainsOfTheLTPs[,1]),2])
MainDomainsOfTheLTPs2[is.na(MainDomainsOfTheLTPs2[,"MainDomain"]),2] <- "OSBP" # Updated to protect from potential future changes

MainDomainsOfTheLTPs4 <- MainDomainsOfTheLTPs2[order(MainDomainsOfTheLTPs2[,2], MainDomainsOfTheLTPs2[,1]),][c(1:29,32:37,30:31,38:40,43,41,42),]


# Presence of LTP-lipid pairs for different technologies
MainDomainsOfTheLTPs5 <- as.data.frame(list(MainDomainsOfTheLTPs4,
                                            
                                            HPTLC = ifelse(MainDomainsOfTheLTPs4[,"LTPProtein"] %in% rownames(HPTLCSpecificitiesPerScreen2[[1]]),"Present","Absent"),
                                            LCMS = ifelse(MainDomainsOfTheLTPs4[,"LTPProtein"] %in% MainDomainsOfTheLTPs[,"LTPProtein"],"Present","Absent")))

# Extraction of information from Uniprot
# The original code used to extract the information is quoted out here below and a fixed dataset is provided, to protect from future changes to the package or database

# BiocManager::install("UniprotR")
# library(UniprotR) # Version 1.2.4

MainDomainsOfTheLTPs7 <- cbind(MainDomainsOfTheLTPs5, 
                               UniprotNames = c("Q86WG3", "Q7Z465", "P12271", "O76054", "Q9UDX3", "O43304", "B5MCN3", "P49638", "Q9BTX7", # "Q92503" removed
                                                
                                                "Q9NZD2", "Q5TA50", "Q00169", "P48739", "Q9UKF7", "P17213", "Q8N4F0", "P29373", "P07148", "Q01469",
                                                "O15540", "P31025", "Q6UWW0", "P02689", "P09455", "P02753", "P82980", "P17900", "P22059", "Q969R2",
                                                
                                                "Q9BXW6", "Q9H1P3", "Q9H0X9", "Q9BZF2", "Q9BZF1", "Q96SU4", "Q9BXB5", "Q9BXB4", "Q6YN16", "P22307",
                                                "Q9UJQ7", "Q9UKL6", "Q9Y365", "Q9Y5P4"))


# MainDomainsOfTheLTPs8 <- cbind(MainDomainsOfTheLTPs7, GetFamily_Domains(MainDomainsOfTheLTPs7$UniprotNames))

# Make sure that we have an exact saved RDS-image of this file in time to avoid issues with Uniprot
# saveRDS(object = MainDomainsOfTheLTPs8, file = "./Output/BackUpOfMainDomainsOfTheLTPs8.rds")

# Load saved RDS-image of this file (if needed, otherwise this can be skipped)
MainDomainsOfTheLTPs8FromTheStorage <- readRDS("./InputData/BackUpOfMainDomainsOfTheLTPs8.rds")

# Check that nothing went wrong during the conversion (if needed, otherwise this can be skipped)
# identical(MainDomainsOfTheLTPs8FromTheStorage, MainDomainsOfTheLTPs8)

# Reassign reloaded variable to the original one
MainDomainsOfTheLTPs8 <- MainDomainsOfTheLTPs8FromTheStorage

EliminateThirds <- function(x){x[seq_along(x)%%3 != 0]}
MainDomainsOfTheLTPs8LongVersionDomainsx <- do.call("rbind", lapply(1:dim(MainDomainsOfTheLTPs8)[1], function(i){do.call("cbind", list(as.character(MainDomainsOfTheLTPs8$LTPProtein[i]), "Domain", matrix(EliminateThirds(unlist(strsplit(as.character(MainDomainsOfTheLTPs8[i,"Domain..FT."]), split = "\\; "))), ncol = 2, byrow = TRUE)))})) # Optimized to x version, to avoid the warnings for differing lengths.

MainDomainsOfTheLTPs8LongVersionDomains2x <- cbind(MainDomainsOfTheLTPs8LongVersionDomainsx, do.call("rbind",strsplit(gsub("DOMAIN ", "", MainDomainsOfTheLTPs8LongVersionDomainsx[,3]), "\\.\\.")))
MainDomainsOfTheLTPs8LongVersionDomains4x <- cbind(MainDomainsOfTheLTPs8LongVersionDomains2x, RegionName = gsub(" /note=", "", MainDomainsOfTheLTPs8LongVersionDomains2x[,4]))  

MainDomainsOfTheLTPs8LongVersionCoils <- do.call("rbind", lapply(1:dim(MainDomainsOfTheLTPs8)[1], function(i){do.call("cbind", list(as.character(MainDomainsOfTheLTPs8$LTPProtein[i]), "CC", matrix(unlist(strsplit(as.character(MainDomainsOfTheLTPs8[i,"Coiled.coil"]), split = "\\; ")),ncol = 2, byrow = TRUE)))}))
AmountOfCoiledCoilsKnown <- sapply(strsplit(as.character(MainDomainsOfTheLTPs8[,"Coiled.coil"]), split = "\\; "), length)/2

MainDomainsOfTheLTPs8LongVersionCoils2 <- cbind(MainDomainsOfTheLTPs8LongVersionCoils, RegionName = unlist(sapply(AmountOfCoiledCoilsKnown, function(x){if(x != 0.5){paste0("CC",1:x)}else{NA}})))
MainDomainsOfTheLTPs8LongVersionCoils4 <- cbind(MainDomainsOfTheLTPs8LongVersionCoils2, do.call("rbind", strsplit(gsub("COILED ", "", MainDomainsOfTheLTPs8LongVersionCoils2[,3]), split = "\\.\\.")))

ListCompositionalBiasElements <- strsplit(MainDomainsOfTheLTPs8$Compositional.bias, "\\;")
MainDomainsOfTheLTPs8LongVersionBias <- do.call("cbind", list(as.character(MainDomainsOfTheLTPs8[,1]), "CompBias", do.call("rbind", lapply(ListCompositionalBiasElements,function(x){if(is.na(x[1])){c(NA,NA,NA)}else{c(strsplit(gsub("COMPBIAS ", "", x[1]), split = "\\.\\.")[[1]], gsub("  /note=", "", x[2], fixed = TRUE))}})))) # NA handling was updated. 

ListOfTheMotifParts <- strsplit(MainDomainsOfTheLTPs8$Motif, "\\;")
MainDomainsOfTheLTPs8LongVersionMotifs <- do.call("cbind", list(as.character(MainDomainsOfTheLTPs8[,1]), "Motif", do.call("rbind", lapply(ListOfTheMotifParts, function(x){if(is.na(x[1])){NA}else{c(strsplit(gsub("MOTIF ", "", x[1]), split = "\\.\\.")[[1]], gsub("  /note=", "", x[2], fixed =TRUE))}})))) # Changed is.na(x) to is.na(x[1]) to solve warning because of different input-lengths.

SequenceAndDomainHighlightsLTPs <- as.data.frame(do.call("rbind", list(MainDomainsOfTheLTPs8LongVersionDomains4x[,c(1:2,5:7)], MainDomainsOfTheLTPs8LongVersionCoils4[,c(1,2,6,7,5)], MainDomainsOfTheLTPs8LongVersionBias, MainDomainsOfTheLTPs8LongVersionMotifs)))
colnames(SequenceAndDomainHighlightsLTPs) <- c("LTPProtein", "TypeRegion", "StartRegion", "StopRegion", "RegionName")

SequenceAndDomainHighlightsLTPs[,3] <- as.numeric(as.character(SequenceAndDomainHighlightsLTPs[,3]))
SequenceAndDomainHighlightsLTPs[,4] <- as.numeric(as.character(SequenceAndDomainHighlightsLTPs[,4]))

SequenceAndDomainHighlightsLTPs2 <- SequenceAndDomainHighlightsLTPs[!is.na(SequenceAndDomainHighlightsLTPs$StartRegion),]


CastProteinDomainsOfLTPs <- ZerosToNAsConverter(Col1ToRowNames(widen5(inputdf = SequenceAndDomainHighlightsLTPs2, ColumnsLong = "LTPProtein", ColumnWide = "RegionName", ColumnValue = "StartRegion", AggregatingFunction = sum, FunctionOutputValueType = double(1))))
CastProteinDomainsOfLTPs2 <- CastProteinDomainsOfLTPs

CastProteinDomainsOfLTPs2[setdiff(MainDomainsOfTheLTPs4[,"LTPProtein"], rownames(CastProteinDomainsOfLTPs2)),] <- NA
CastProteinDomainsOfLTPs4 <- CastProteinDomainsOfLTPs2[MainDomainsOfTheLTPs4[,"LTPProtein"],]


CastMainProteinDomainsLTPs <- Col1ToRowNames(widen5(inputdf = as.data.frame(MainDomainsOfTheLTPs4), ColumnsLong = "LTPProtein", ColumnWide = "MainDomain", ColumnValue = "MainDomain", AggregatingFunction = length, FunctionOutputValueType = integer(1), SortRows = FALSE))[MainDomainsOfTheLTPs4[,"LTPProtein"],]

CastMainProteinDomainsLTPswnas <- CastMainProteinDomainsLTPs
CastMainProteinDomainsLTPswnas[CastMainProteinDomainsLTPswnas == 0] <- NA

CastManyProteinDomainsLTPs <- cbind(CastProteinDomainsOfLTPs4,CastMainProteinDomainsLTPswnas*round(mean(as.matrix(CastProteinDomainsOfLTPs4), na.rm = TRUE)))


OrderDecreasingPerMainDomain <- unlist(lapply(unique(MainDomainsOfTheLTPs4[,2]), function(x){names(sort(rowSums(CastManyProteinDomainsLTPs[MainDomainsOfTheLTPs4[MainDomainsOfTheLTPs4[,2] == x,1],], na.rm = TRUE), decreasing = TRUE))}))
MainDomainsOfTheLTPs4x <- as.data.frame(cbind(OrderDecreasingPerMainDomain, MainDomainsOfTheLTPs4[,2]))

# Sort by type of domains and split in the types
# Explicitly defined in next code, saved and reloaded, to avoid possible future changes (cleans rownames too) (original code used for construction is quoted out below) 

# DomainTypes2 <- as.matrix(cbind(TypeRegion = c(rep("Domain", 6), rep("CC", 2), rep("CompBias", 3), rep("Motif", 3), rep("MainDomain", 9)), 
#                                 RegionName = c("GOLD", "CRAL-TRIO", "PRELI/MSF1", "START", "PH", "SCP2", "CC1", "CC2", "Ala/Gly-rich", "Poly-Leu", "Ser-rich", "FFAT", "Microbody targeting signal", "Nuclear localization signal", "CRAL-TRIO.1", "GLTP", "IP_trans", "START.1", "OSBP", "LBP_BPI_CETP", "lipocalin", "ML", "scp2")))

# Write this table away to ensure that always the same can be read in
# write.table(DomainTypes2, file="./Output/DomainTypes2.tsv", sep="\t", row.names = FALSE, quote = FALSE)

# Import the domain types from input file
DomainTypes2 <- as.matrix(read.csv(file = "./InputData/DomainTypes2.tsv", header = TRUE, sep = "\t", as.is = TRUE))

# Also here explicit definition in characters introduced instead of numeric subsetting of levels, to protect from possible future changes (original code is quoted out below)
ListDomainsOfTheLTPs <- lapply(c("CRAL-TRIO", "GLTP", "IP_trans", "START", "OSBP", "LBP_BPI_CETP", "lipocalin", "ML", "scp2"), function(x){MainDomainsOfTheLTPs4x[MainDomainsOfTheLTPs4x[,2] == x,]})

# Original code
# ListDomainsOfTheLTPs <- lapply(levels(MainDomainsOfTheLTPs4x[,2])[c(1:3, 9, 7, 4:6,8)], function(x){MainDomainsOfTheLTPs4x[MainDomainsOfTheLTPs4x[,2] == x,]})

ListDomainsOfTheLTPs[[4]] <- ListDomainsOfTheLTPs[[4]][3:1,]
ListDomainsOfTheLTPs[[7]][1:2,] <- ListDomainsOfTheLTPs[[7]][2:1,]

ListDomainsOfTheLTPs[[5]][1:2,] <- ListDomainsOfTheLTPs[[5]][2:1,]
MainDomainsOfTheLTPs4xx <- do.call("rbind",ListDomainsOfTheLTPs)

OrderDecreasingPerMainDomain <- unlist(lapply(unique(MainDomainsOfTheLTPs4[,2]), function(x){names(sort(rowSums(CastManyProteinDomainsLTPs[MainDomainsOfTheLTPs4[MainDomainsOfTheLTPs4[,2] == x,1],], na.rm = TRUE), decreasing = TRUE))}))
CastManyProteinDomainsLTPsDecreasinglyOrderDomainsPerMainDomain <- CastManyProteinDomainsLTPs[OrderDecreasingPerMainDomain, ]

DomainsByRowColumnAssociations <- CastManyProteinDomainsLTPsDecreasinglyOrderDomainsPerMainDomain[as.character(MainDomainsOfTheLTPs4xx[,1]),c(11,4,6,13,14,9,5,2,3,1,10,8,7,12,15:17,23,21,18:20,22)]
DomainsByRowColumnAssociationsReordered <- DomainsByRowColumnAssociations[,DomainTypes2[,2]]


LipidCompactionMatrixSphingolipids <- cbind(OriginalRow = rownames(LTPMatrixPosInVivo), NewRow = c(rep("Cer*",5), rownames(LTPMatrixPosInVivo)[6], rep("HexCer*",2), rownames(LTPMatrixPosInVivo)[9:10], rep("SM*",3), rownames(LTPMatrixPosInVivo)[14:34]))

LTPMatrixPosInVivoslc <- do.call("rbind", lapply(unique(LipidCompactionMatrixSphingolipids[,"NewRow"]), function(x){colSums(LTPMatrixPosInVivo[LipidCompactionMatrixSphingolipids[LipidCompactionMatrixSphingolipids[,"NewRow"] == x, "OriginalRow"],, drop = FALSE])}))
rownames(LTPMatrixPosInVivoslc) <- unique(LipidCompactionMatrixSphingolipids[,"NewRow"])

LTPMatrixPosInVivommnslc <- MinMaxNormMatrixFunc(inputdata = LTPMatrixPosInVivoslc)*9+1
LTPMatrixPosInVivommnslc[is.infinite(as.matrix(LTPMatrixPosInVivommnslc))] <- 0

LTPMatrixNegInVivoslc <- do.call("rbind", lapply(unique(LipidCompactionMatrixSphingolipids[,"NewRow"]), function(x){colSums(LTPMatrixNegInVivo[LipidCompactionMatrixSphingolipids[LipidCompactionMatrixSphingolipids[,"NewRow"] == x, "OriginalRow"],, drop = FALSE])}))
rownames(LTPMatrixNegInVivoslc) <- unique(LipidCompactionMatrixSphingolipids[,"NewRow"])

LTPMatrixNegInVivommnslc <- MinMaxNormMatrixFunc(inputdata = LTPMatrixNegInVivoslc)*9+1
LTPMatrixNegInVivommnslc[is.infinite(as.matrix(LTPMatrixNegInVivommnslc))] <- 0

LTPMatrixPosInVitroslc <- do.call("rbind", lapply(unique(LipidCompactionMatrixSphingolipids[,"NewRow"]), function(x){colSums(LTPMatrixPosInVitro[LipidCompactionMatrixSphingolipids[LipidCompactionMatrixSphingolipids[,"NewRow"] == x, "OriginalRow"],, drop = FALSE])}))
rownames(LTPMatrixPosInVitroslc) <- unique(LipidCompactionMatrixSphingolipids[,"NewRow"])

LTPMatrixPosInVitrommnslc <- MinMaxNormMatrixFunc(inputdata = LTPMatrixPosInVitroslc)*9+1
LTPMatrixPosInVitrommnslc[is.infinite(as.matrix(LTPMatrixPosInVitrommnslc))] <- 0

LTPMatrixNegInVitroslc <- do.call("rbind", lapply(unique(LipidCompactionMatrixSphingolipids[,"NewRow"]), function(x){colSums(LTPMatrixNegInVitro[LipidCompactionMatrixSphingolipids[LipidCompactionMatrixSphingolipids[,"NewRow"] == x, "OriginalRow"],, drop = FALSE])}))
rownames(LTPMatrixNegInVitroslc) <- unique(LipidCompactionMatrixSphingolipids[,"NewRow"])

LTPMatrixNegInVitrommnslc <- MinMaxNormMatrixFunc(inputdata = LTPMatrixNegInVitroslc)*9+1
LTPMatrixNegInVitrommnslc[is.infinite(as.matrix(LTPMatrixNegInVitrommnslc))] <- 0

LTPMatrixTopInVitrommnslc <- do.call("cbind", lapply(1:dim(LTPMatrixNegInVitrommnslc)[2], function(x){pmax(LTPMatrixNegInVitrommnslc[,x], LTPMatrixPosInVitrommnslc[,x])}))
LTPMatrixTopInVivommnslc <- do.call("cbind", lapply(1:dim(LTPMatrixNegInVivommnslc)[2], function(x){pmax(LTPMatrixNegInVivommnslc[,x], LTPMatrixPosInVivommnslc[,x])}))

rownames(LTPMatrixTopInVitrommnslc) <- rownames(LTPMatrixPosInVitrommnslc)
colnames(LTPMatrixTopInVitrommnslc) <- colnames(LTPMatrixPosInVitrommnslc)

rownames(LTPMatrixTopInVivommnslc) <- rownames(LTPMatrixPosInVivommnslc)
colnames(LTPMatrixTopInVivommnslc) <- colnames(LTPMatrixPosInVivommnslc)

#### HPTLC-data: further clean-up and formatting # Was in post-production of the figure updated by restricting some of these observations from the figures when uncertainties present.
HPTLCSpecificitiesPerScreen2hdr <- HPTLCSpecificitiesPerScreen2

colnames(HPTLCSpecificitiesPerScreen2hdr[[1]]) <- c("Cer*", "Sterol", "DAG", "PC", "PE", "PG", "PIPs", "PS")
colnames(HPTLCSpecificitiesPerScreen2hdr[[2]]) <- c("Cer*", "Sterol", "DAG", "PC", "PE", "PG", "PIPs", "PS")

HPTLCSpecificitiesPerScreen2hdrm <- lapply(HPTLCSpecificitiesPerScreen2hdr, function(x){meltmatrix(t(as.matrix(5*x)))})
HPTLCSpecificitiesPerScreen2hdrm2 <- lapply(HPTLCSpecificitiesPerScreen2hdrm, function(x){x[x$value != 0,]})

HPTLCSpecificitiesPerScreen2hdrm4 <- list()
HPTLCSpecificitiesPerScreen2hdrm4[[1]] <- HPTLCSpecificitiesPerScreen2hdrm2[[1]][sapply(1:dim(HPTLCSpecificitiesPerScreen2hdrm2[[1]])[1], function(x){
  
  !((as.character(HPTLCSpecificitiesPerScreen2hdrm2[[1]][x,1]) %in% rownames(LTPMatrixTopInVivommnslc)) &&
      (as.character(HPTLCSpecificitiesPerScreen2hdrm2[[1]][x,2]) %in% colnames(LTPMatrixTopInVivommnslc)) &&
      
      (LTPMatrixTopInVivommnslc[as.character(HPTLCSpecificitiesPerScreen2hdrm2[[1]][x,1]), as.character(HPTLCSpecificitiesPerScreen2hdrm2[[1]][x,2])] != 0))
}),] # Leaves PE_STARD10 in, but also Cer*_STARD11 out.


HPTLCSpecificitiesPerScreen2hdrm4[[2]] <- HPTLCSpecificitiesPerScreen2hdrm2[[2]][sapply(1:dim(HPTLCSpecificitiesPerScreen2hdrm2[[2]])[1], function(x){
  
  !((as.character(HPTLCSpecificitiesPerScreen2hdrm2[[2]][x,1]) %in% rownames(LTPMatrixTopInVitrommnslc)) &&
      (as.character(HPTLCSpecificitiesPerScreen2hdrm2[[2]][x,2]) %in% colnames(LTPMatrixTopInVitrommnslc)) &&
      
      (LTPMatrixTopInVitrommnslc[as.character(HPTLCSpecificitiesPerScreen2hdrm2[[2]][x,1]), as.character(HPTLCSpecificitiesPerScreen2hdrm2[[2]][x,2])] != 0))
}),] 

MeltedInVivoCombinedslchdr <- rbind(meltmatrix(LTPMatrixTopInVivommnslc), HPTLCSpecificitiesPerScreen2hdrm4[[1]])
MeltedInVitroCombinedslchdr <- rbind(meltmatrix(LTPMatrixTopInVitrommnslc), HPTLCSpecificitiesPerScreen2hdrm4[[2]])

CastInVivoCombinedslchdr <- as.matrix(Col1ToRowNames(widen5(inputdf = MeltedInVivoCombinedslchdr, ColumnsLong = "Var1", ColumnWide = "Var2", ColumnValue = "value", AggregatingFunction = function(x){if(length(x)>0){max(x,na.rm = TRUE)}else{0}}, FunctionOutputValueType = double(1), SortRows = FALSE, FillerValue = 0)))
CastInVitroCombinedslchdr <- as.matrix(Col1ToRowNames(widen5(inputdf = MeltedInVitroCombinedslchdr, ColumnsLong = "Var1", ColumnWide = "Var2", ColumnValue = "value", AggregatingFunction = function(x){if(length(x)>0){max(x,na.rm = TRUE)}else{0}}, FunctionOutputValueType = double(1), SortRows = FALSE, FillerValue = 0)))

# Internal conversions simplified/eliminated by introduction of more advanced widen5 function instead of the original script snippet, and thus no missing colnames anymore


# all(colnames(CastInVitroCombinedslchdr) %in% rownames(DomainsByRowColumnAssociationsReordered)) # TRUE
# all(rownames(DomainsByRowColumnAssociationsReordered) %in% colnames(CastInVitroCombinedslchdr)) # FALSE (in a stringAsFactors = FALSE environment)

# Make in vitro matrix matching with in cellulo matrix because it has some columns less in environments not following the initialized options (reintroduced code from the previous version)
CastInVitroCombinedslchdr2 <- cbind(CastInVitroCombinedslchdr, matrix(data = 0, 
                                                                      
                                                                      nrow = nrow(CastInVitroCombinedslchdr),
                                                                      ncol = length(rownames(DomainsByRowColumnAssociationsReordered)[!(rownames(DomainsByRowColumnAssociationsReordered) %in% colnames(CastInVitroCombinedslchdr))]),
                                                                      
                                                                      dimnames = list(rownames(CastInVitroCombinedslchdr),
                                                                                      rownames(DomainsByRowColumnAssociationsReordered)[!(rownames(DomainsByRowColumnAssociationsReordered) %in% colnames(CastInVitroCombinedslchdr))])))

InVivoDataSetTotalslchdr <- as.matrix(CastInVivoCombinedslchdr[, rownames(DomainsByRowColumnAssociationsReordered)])
InVitroDataSetTotalslchdr <- as.matrix(CastInVitroCombinedslchdr2[, rownames(DomainsByRowColumnAssociationsReordered)]) # Reintroduced back from previous code to get matching columns again.

c(rownames(InVivoDataSetTotalslchdr), rownames(InVitroDataSetTotalslchdr))[!unique(c(rownames(InVivoDataSetTotalslchdr), rownames(InVitroDataSetTotalslchdr))) %in% rownames(InVivoDataSetTotalslchdr)]
# character(0)

MissingLipidEntriesInVitro <- c(rownames(InVivoDataSetTotalslchdr), rownames(InVitroDataSetTotalslchdr))[!unique(c(rownames(InVivoDataSetTotalslchdr), rownames(InVitroDataSetTotalslchdr))) %in% rownames(InVitroDataSetTotalslchdr)]
# "Sterol" "PIPs"

MeltedInVivoCombinedslc <- rbind(meltmatrix(LTPMatrixTopInVivommnslc), meltmatrix(t(as.matrix(5*HPTLCSpecificitiesPerScreen2[[1]]))))
MeltedInVitroCombinedslc <- rbind(meltmatrix(LTPMatrixTopInVitrommnslc), meltmatrix(t(as.matrix(5*HPTLCSpecificitiesPerScreen2[[2]]))))

# Again reduction in script size by use of widen5 and internally defined function to handle infinites
# Should also handle possible warnings for introduction of infinites

CastInVivoCombinedslc <- as.matrix(Col1ToRowNames(widen5(inputdf = MeltedInVivoCombinedslc, ColumnsLong = "Var1", ColumnWide = "Var2", ColumnValue = "value", AggregatingFunction = function(x){if(length(x)>0){max(x,na.rm = TRUE)}else{0}}, FunctionOutputValueType = double(1), SortRows = FALSE)))
CastInVitroCombinedslc <- as.matrix(Col1ToRowNames(widen5(inputdf = MeltedInVitroCombinedslc, ColumnsLong = "Var1", ColumnWide = "Var2", ColumnValue = "value", AggregatingFunction = function(x){if(length(x)>0){max(x,na.rm = TRUE)}else{0}}, FunctionOutputValueType = double(1), SortRows = FALSE)))

InVivoDataSetTotalslc <- as.matrix(CastInVivoCombinedslc[, rownames(DomainsByRowColumnAssociationsReordered)])
InVitroDataSetTotalslc <- as.matrix(CastInVitroCombinedslc[, rownames(DomainsByRowColumnAssociationsReordered)])

InVivoDataSetWithoutRedundantStarsslc <- InVivoDataSetTotalslc[-match(c("DAG*", "PC*", "PE*", "PG*", "PS*"), rownames(InVivoDataSetTotalslc)),]
InVitroDataSetWithoutRedundantStarsslc <- InVitroDataSetTotalslc[-match(c("DAG*", "PC*", "PE*", "PG*", "PS"), rownames(InVitroDataSetTotalslc)),]

rownames(InVivoDataSetWithoutRedundantStarsslc)[rownames(InVivoDataSetWithoutRedundantStarsslc) == "PIP*"] <- "PIPs"
rownames(InVitroDataSetWithoutRedundantStarsslc)[rownames(InVitroDataSetWithoutRedundantStarsslc) == "PIP*"] <- "PIPs"

rownames(InVivoDataSetWithoutRedundantStarsslc)[rownames(InVivoDataSetWithoutRedundantStarsslc) == "CH*"] <- "CH"
rownames(InVitroDataSetWithoutRedundantStarsslc)[rownames(InVitroDataSetWithoutRedundantStarsslc) == "CH*"] <- "CH"

# PS for in vivo already corrected before
rownames(InVitroDataSetWithoutRedundantStarsslc)[rownames(InVitroDataSetWithoutRedundantStarsslc) == "PS*"] <- "PS"

InVivoDataSetslc <- InVivoDataSetWithoutRedundantStarsslc[c(1:12,14:19,29,20,25,21:23,26,13,24,28,27),]
InVitroDataSetslc <- InVitroDataSetWithoutRedundantStarsslc[c(1:12,14:19,28,29,24,20:22,25,13,23,27,26),]

InVivoDataSetslchdr <- InVivoDataSetTotalslchdr[c(1:12,14:19,29,20,25,21:23,26,13,24,28,27),]
InVitroDataSetslchdr <- rbind(InVitroDataSetTotalslchdr, matrix(data = 0, 
                                                                
                                                                nrow = length(MissingLipidEntriesInVitro),
                                                                ncol = ncol(InVitroDataSetTotalslchdr),
                                                                
                                                                dimnames = list(MissingLipidEntriesInVitro,
                                                                                colnames(InVitroDataSetTotalslchdr))))[rownames(InVivoDataSetslchdr),]

InVivoDataSetslc4hdr <- InVivoDataSetslchdr
InVivoDataSetslc4hdr[InVivoDataSetslc4hdr == 0] <- NA

InVitroDataSetslc4hdr <- InVitroDataSetslchdr
InVitroDataSetslc4hdr[InVitroDataSetslc4hdr == 0] <- NA

# Extract for comparison with the literature data
MainGroupsOfLipids <- cbind(Lipids = rownames(CastInVivoCombined)[c(35,1:19,37,20,21,38,22:23,39,24:26,41:42,27,40,28:30,32,33,31,36,34)], Group = c(rep("SL",14), rep("FA",2), rep("LL",4), rep("GPL",17), "PGP", "CL", "TAG", "CH", "VA"))

InVivoDataSetTotal <- as.matrix(CastInVivoCombined[MainGroupsOfLipids[,"Lipids"], rownames(DomainsByRowColumnAssociationsReordered)])
InVitroDataSetTotal <- as.matrix(CastInVitroCombined[MainGroupsOfLipids[,"Lipids"], rownames(DomainsByRowColumnAssociationsReordered)])

InVivoDataSetWithoutRedundantStars <- InVivoDataSetTotal[-match(c("Cer*", "DAG*", "PC*", "PE*", "PG*", "PS*"), rownames(InVivoDataSetTotal)),]
InVitroDataSetWithoutRedundantStars <- InVitroDataSetTotal[-match(c("Cer*", "DAG*", "PC*", "PE*", "PG*", "PS"), rownames(InVitroDataSetTotal)),]

rownames(InVivoDataSetWithoutRedundantStars)[rownames(InVivoDataSetWithoutRedundantStars) == "PIP*"] <- "PIPs"
rownames(InVitroDataSetWithoutRedundantStars)[rownames(InVitroDataSetWithoutRedundantStars) == "PIP*"] <- "PIPs"

rownames(InVivoDataSetWithoutRedundantStars)[rownames(InVivoDataSetWithoutRedundantStars) == "CH*"] <- "CH/CE"
rownames(InVitroDataSetWithoutRedundantStars)[rownames(InVitroDataSetWithoutRedundantStars) == "CH*"] <- "CH/CE"

# PS for in vivo already corrected before
rownames(InVitroDataSetWithoutRedundantStars)[rownames(InVitroDataSetWithoutRedundantStars) == "PS*"] <- "PS"

# Literature coverage # Simplified version used with only the data that are relevant here
LiteratureDataLTPsIntegratedSimplified <- read.csv(file = "./InputData/LiteratureDataLTPsIntegratedSimplified.txt", header = TRUE, sep = "\t", as.is = TRUE)

LiteratureDataLTPsIntegratedSimplified$InVivo <- as.numeric(gsub(",",".", LiteratureDataLTPsIntegratedSimplified$InVivo))
LiteratureDataLTPsIntegratedSimplified$InVitro <- as.numeric(gsub(",",".", LiteratureDataLTPsIntegratedSimplified$InVitro))

# 29.6% Literature consensus
sum(LiteratureDataLTPsIntegratedSimplified$LiteratureConsensus)/length(LiteratureDataLTPsIntegratedSimplified$LiteratureConsensus) # 0.296

# Conversion of script to widen5 gives slightly different column order because of the filling in of columns, but that is no issue considering the steps after.
LiteratureConsensusLinksWideVersion <- Col1ToRowNames(widen5(inputdf = LiteratureDataLTPsIntegratedSimplified, ColumnsLong = "Lipid", ColumnWide = "Protein", ColumnValue = "LiteratureConsensus", AggregatingFunction = sum, FunctionOutputValueType = integer(1), FillerValue = as.integer(0)))

LiteratureConsensusLinksWideVersion2 <- LiteratureConsensusLinksWideVersion[rownames(InVivoDataSetWithoutRedundantStars), colnames(InVivoDataSetWithoutRedundantStars)] # Simplified steps before this here to not necessitate input of two literature files.
LiteratureConsensusLinksWideVersion4 <- !((InVivoDataSetWithoutRedundantStars == 0) & (InVitroDataSetWithoutRedundantStars == 0)) & (LiteratureConsensusLinksWideVersion2 == 0) # InVivoDataSetWithoutRedundantStars instead of InVivoDataSet to avoid loop and same for in vitro data


LipidClassLiteratureDataSet <- as.matrix(LiteratureConsensusLinksWideVersion4)[c(1:19,21:28,32,29:31,33,20,34,35,36),]

LipidClassLiteratureDataSetslc <- LipidClassLiteratureDataSet[c(1,6,7,9:11,14:36),]
rownames(LipidClassLiteratureDataSetslc) <- rownames(InVivoDataSetslc)

# Correction on original heatmap: literature data for OSBPL2-PIPs known now & PC-O and PE-O should not be seen as known only because PC and PE are known
# Was in post-production of the figure updated by restricting some of these observations from the figures when uncertainties present.

LipidClassLiteratureDataSetslc4hdr <- LipidClassLiteratureDataSetslc
rownames(LipidClassLiteratureDataSetslc) <- rownames(InVivoDataSetslc4hdr)

LipidClassLiteratureDataSetslc4hdr["PIPs","OSBPL2"] <- FALSE
LipidClassLiteratureDataSetslc4hdr["PC-O","STARD2"] <- TRUE

LipidClassLiteratureDataSetslc4hdr["PC-O","STARD10"] <- TRUE
LipidClassLiteratureDataSetslc4hdr["PC-O","GM2A"] <- TRUE

LipidClassLiteratureDataSetslc4hdr["PC-O","LCN1"] <- TRUE
LipidClassLiteratureDataSetslc4hdr["PE-O","STARD10"] <- TRUE


ORPDomains <- read.csv(file = "./InputData/ORPDomains270520202.txt", header = TRUE, sep = "\t", as.is = TRUE)

ORPDomains2 <- do.call("cbind", list(ORPDomains[,1:3], ORPDomains[,4:7]-2, FFAT = ifelse(as.logical(ORPDomains[,8]), ORPDomains[,5]-24, NA)))
ORPDomains2[9,4] <- 148 # Longer version

ORPDomains2[10,4] <- 2 # Somewhere completely in front according to InterPro, and according to UniProt 2
colnames(ORPDomains2)[5] <- "ORD"

DomainsByRowColumnAssociationsReordered2 <- DomainsByRowColumnAssociationsReordered
for(i in ORPDomains2[ORPDomains2$HGNC %in% rownames(DomainsByRowColumnAssociationsReordered2),"HGNC"]){
  
  DomainsByRowColumnAssociationsReordered2[i, "PH"] <- ORPDomains2[ORPDomains2$HGNC == i, "PH"]
  DomainsByRowColumnAssociationsReordered2[i, "FFAT"] <- ORPDomains2[ORPDomains2$HGNC == i, "FFAT"]
  
}


DomainsByRowColumnAssociationsReordered4 <- Col1ToRowNames(merge(x = DomainsByRowColumnAssociationsReordered2, 
                                                                 y = ORPDomains2[ORPDomains2$HGNC %in% rownames(DomainsByRowColumnAssociationsReordered2),c("HGNC","TMpd","Ankyrin", "ORD")],
                                                                 
                                                                 by.x = 0,
                                                                 by.y = 1,
                                                                 
                                                                 all = TRUE))[rownames(DomainsByRowColumnAssociationsReordered2),]



# OSBPL1A: ANK 47,80,175 ==> PFAM: 7, 160; CC: 123,356,433,879 <-> 430, 877

DomainsByRowColumnAssociationsReordered4["OSBPL1A", "CC1"] <- 121 
DomainsByRowColumnAssociationsReordered4["OSBPL1A", "CC2"] <- 354

DomainsByRowColumnAssociationsReordered4["OSBPL1A", "CC3"] <- 430
DomainsByRowColumnAssociationsReordered4["OSBPL1A", "CC4"] <- 877

# 1 to 3 less --> 2; for TMpds: 18 & 33
# OSBPL5: PFAM: CC: 96, 666; TMpd: 860

DomainsByRowColumnAssociationsReordered4["OSBPL5", "CC1"] <- 93 # As before
DomainsByRowColumnAssociationsReordered4["OSBPL5", "CC2"] <- 664

# OSBP2: PFAM: CC: 389
DomainsByRowColumnAssociationsReordered4["OSBP2", "CC1"] <- 387


# OSBPL8: PFAM: CC: 118, 834; TMpd: 871

DomainsByRowColumnAssociationsReordered4["OSBPL8", "CC1"] <- 116 
DomainsByRowColumnAssociationsReordered4["OSBPL8", "CC2"] <- 832


# OSBPL7: PFAM: CC: 201, 776

DomainsByRowColumnAssociationsReordered4["OSBPL7", "CC1"] <- 199 
DomainsByRowColumnAssociationsReordered4["OSBPL7", "CC2"] <- 774


# OSBPL9: PFAM: CC: 130, 681

DomainsByRowColumnAssociationsReordered4["OSBPL9", "CC1"] <- 128 
DomainsByRowColumnAssociationsReordered4["OSBPL9", "CC2"] <- 679

# OSBPL2: PFAM: CC: 416
DomainsByRowColumnAssociationsReordered4["OSBPL2", "CC1"] <- 416 


# GLTP domain start sites

DomainsByRowColumnAssociationsReordered4["GLTPD1", "GLTP"] <- 28 #Interpro & PFAM
DomainsByRowColumnAssociationsReordered4["GLTP", "GLTP"] <- 18 #Interpro & PFAM


# IP_trans

DomainsByRowColumnAssociationsReordered4["PITPNA", "IP_trans"] <- 2 # PFAM And also CC 240 
DomainsByRowColumnAssociationsReordered4["PITPNB", "IP_trans"] <- 2 # PFAM And also CC 240 

DomainsByRowColumnAssociationsReordered4["PITPNA", "CC1"] <- 240
DomainsByRowColumnAssociationsReordered4["PITPNB", "CC1"] <- 240

DomainsByRowColumnAssociationsReordered4["PITPNC1", "IP_trans"] <- 1 # PFAM



# BPIFB2: PFAM LBP_BPI_CETP at 33 (with signal peptide of 20)

DomainsByRowColumnAssociationsReordered4["BPIFB2", "LBP_BPI_CETP"] <- 33
DomainsByRowColumnAssociationsReordered4["BPIFB2", "Signal peptide"] <- 1


# BPI: PFAM LBP_BPI_CETP at 43 (with signal peptide of 30)

DomainsByRowColumnAssociationsReordered4["BPI", "LBP_BPI_CETP"] <- 43
DomainsByRowColumnAssociationsReordered4["BPI", "Signal peptide"] <- 1

# lipocalin --> capital
DomainsByRowColumnAssociationsReordered4["CRABP2", "lipocalin"] <- 5 # RABR <- 133 uniprot

DomainsByRowColumnAssociationsReordered4["FABP5", "lipocalin"] <- 8 # FABR <- 129 uniprot
DomainsByRowColumnAssociationsReordered4["FABP1", "lipocalin"] <- 1 # No uniprot entry for FABR. No PFAM.

DomainsByRowColumnAssociationsReordered4["FABP7", "lipocalin"] <- 7
DomainsByRowColumnAssociationsReordered4["LCN1", "lipocalin"] <- 32 # PFAM: With 23 AA signal peptide 

DomainsByRowColumnAssociationsReordered4["LCN1", "Signal peptide"] <- 1
DomainsByRowColumnAssociationsReordered4["LCN15", "Signal peptide"] <- 1

DomainsByRowColumnAssociationsReordered4["LCN15", "lipocalin"] <- 34 # PFAM: With 20 AA signal peptide 
DomainsByRowColumnAssociationsReordered4["PMP2", "lipocalin"] <- 6 # FABR <- 127 uniprot

DomainsByRowColumnAssociationsReordered4["RBP1", "lipocalin"] <- 6 # STRA6 interaction region at 22
DomainsByRowColumnAssociationsReordered4["RBP4", "lipocalin"] <- 39 # PFAM: With 18 AA signal peptide  

DomainsByRowColumnAssociationsReordered4["RBP4", "Signal peptide"] <- 1
DomainsByRowColumnAssociationsReordered4["RBP5", "lipocalin"] <- 6


# ML

DomainsByRowColumnAssociationsReordered4["GM2A", "ML"] <- 33 # PFAM: With 23 AA signal peptide 
DomainsByRowColumnAssociationsReordered4["GM2A", "Signal peptide"] <- 1


# SCP2: PFAM: Thiolase: 14; SCP2: 437 <=> 433 UniProt

DomainsByRowColumnAssociationsReordered4["SCP2", "Thiolase"] <- 14 
# DomainsByRowColumnAssociationsReordered4["SCP2", "Signal peptide"] <- NA

# SCP2D1: 44 UniProt; PFAM: 48
# HSDL2: PFAM: adh_short: 11; SCP2: 310 <=> 306 UniProt

DomainsByRowColumnAssociationsReordered4["HSDL2", "adh_short"] <- 11 
# DomainsByRowColumnAssociationsReordered4["HSDL2", "Signal peptide"] <- NA

# Peroxisome targeting signal at C-terminus: ARL(415) vs AKL(544) for SCP2 (while SKL consensus) <==> AKF(153) for SCP2D1 (all at C-terminus)
DomainsByRowColumnAssociationsReordered4["HSDL2", "PTS1"] <- 415

DomainsByRowColumnAssociationsReordered4["SCP2", "PTS1"] <- 544
DomainsByRowColumnAssociationsReordered4["SCP2D1", "PTS1"] <- 153 # Not know yet necessarily: maybe indicate this!

# The START-domain and SEC14L proteins look fine


# BNIPL: PFAM: BNIP2: 83, CRAL_TRIO_2: 212 <=> 191 --> Keep original and put 83 anyway.
DomainsByRowColumnAssociationsReordered4["BNIPL", "BNIP2"] <- 83

# ATCAY: PFAM: BNIP2: 59; CRAL_TRIO_2: 188 <=> 171 --> change the following to 17 less: TMpd: 182, 214 --> 165, 199
DomainsByRowColumnAssociationsReordered4["ATCAY", "BNIP2"] <- 59

DomainsByRowColumnAssociationsReordered4["ATCAY", "TMpd"] <- 165
DomainsByRowColumnAssociationsReordered4["ATCAY", "TMpd2"] <- 199

# Check other CRAL_TRIO is also N
# RLBP1: PFAM: CRAL_TRIO: 139 <=> 136; But CRAL_TRIO_N: 15 --> keep; and CC1: 45 --> keep too

DomainsByRowColumnAssociationsReordered4["RLBP1", "CC1"] <- 45
DomainsByRowColumnAssociationsReordered4["RLBP1", "CRAL_TRIO_N"] <- 15

# TTPAL: PFAM: CRAL_TRIO_N: 38; CRAL_TRIO: 123 <=> 117
DomainsByRowColumnAssociationsReordered4["TTPAL", "CRAL_TRIO_N"] <- 38

# TTPA: PFAM: CRAL_TRIO_N: 38; CRAL_TRIO: 95 <=> 88
DomainsByRowColumnAssociationsReordered4["TTPA", "CRAL_TRIO_N"] <- 38


# Search types of CRAL_TRIO in previous CRAL_TRIOs

DomainsByRowColumnAssociationsReordered4["SEC14L5", "CRAL_TRIO_N"] <- 239 # PFAM & +-Interpro
DomainsByRowColumnAssociationsReordered4["SEC14L2", "CRAL_TRIO_N"] <- 2 # Interpro

DomainsByRowColumnAssociationsReordered4["SEC14L4", "CRAL_TRIO_N"] <- 2 # Interpro
DomainsByRowColumnAssociationsReordered4["SEC14L6", "CRAL_TRIO_N"] <- 2 # Interpro although at 34 also smaller domain called similar

# BNIPL and ATCAY do not have a CRAL_TRIO_N domain even in Interpro and PFAM calls the CRAL_TRIO CRAL_TRIO_2 instead
# Makes sense based on the detailed discussion in (Gupta et al., 2012, PLoS ONE).

DomainsByRowColumnAssociationsReordered4["ATCAY", "CRAL_TRIO_2"] <- DomainsByRowColumnAssociationsReordered4["ATCAY", "CRAL-TRIO"]
DomainsByRowColumnAssociationsReordered4["BNIPL", "CRAL_TRIO_2"] <- DomainsByRowColumnAssociationsReordered4["BNIPL", "CRAL-TRIO"]

DomainsByRowColumnAssociationsReordered4["ATCAY", "CRAL-TRIO"] <- NA
DomainsByRowColumnAssociationsReordered4["BNIPL", "CRAL-TRIO"] <- NA

DomainsByRowColumnAssociationsReordered5 <- DomainsByRowColumnAssociationsReordered4[,c("PRELI/MSF1","CRAL_TRIO_N","CRAL-TRIO", "GOLD", "BNIP2", "CRAL_TRIO_2", "GLTP", "IP_trans", "START", "PH", "ORD", "Ankyrin", "LBP_BPI_CETP", "lipocalin", "ML", "Thiolase", "adh_short","SCP2", "TMpd", "TMpd2", "CC1", "CC2", "CC3", "CC4", "Ser-rich", "Ala/Gly-rich", "Poly-Leu", "FFAT", "Signal peptide", "Nuclear localization signal", "PTS1")]
ReorderedLTPsByManualSeriation <- rownames(DomainsByRowColumnAssociationsReordered5)[c(1:4,9:5,10:17,22,18,25,19,27,26,21,23,20,24,28:29,40,34:35,38,37,39,36,32,33,31,30,43:41)]

# See reference to MainGroupsOfLipids earlier
MainGroupsOfLipidsWithoutRedundantStars <- MainGroupsOfLipids[-match(c("Cer*", "DAG*", "PC*", "PE*", "PG*", "PS*"), MainGroupsOfLipids[,"Lipids"]),]

MainGroupsOfLipidsWithoutRedundantStars2 <- MainGroupsOfLipidsWithoutRedundantStars[c(1:19,21:28,32,29:31,33,20,34,35,36),]
MainGroupsOfLipidsWithoutRedundantStars2[,"Group"] <- c(rep("SL", 13),rep("FA", 2), rep("LL", 4), rep("GPL", 13), rep("NL", 4)) # (D)GPL --> GPL

LegendName <- "Legend"
LegendColor <- col_fung

WidthAdaptor <- 160
HeightAdaptor <- 160

SplitByRowsVector <- factor(MainGroupsOfLipidsWithoutRedundantStars2[,"Group"], levels = unique(MainGroupsOfLipidsWithoutRedundantStars2[,"Group"]))
# SplitByColsVector <- factor(MainDomainsOfTheLTPs4xx[,2], levels = levels(MainDomainsOfTheLTPs4x[,2])[c(1:3, 9, 7, 4:6,8)]) # Original code

# Updated code to avoid numerical subsetting of levels to avoid future changes (see above for the original code-line)
SplitByColsVector <- factor(MainDomainsOfTheLTPs4xx[,2], levels = c("CRAL-TRIO", "GLTP", "IP_trans", "START", "OSBP", "LBP_BPI_CETP", "lipocalin", "ML", "scp2"))

GraphNameslc <- "White circles with black borders for novelty (vs. consensus literature knowledge)(sphingolipids aggregated)"
library(ComplexHeatmap)

pdf("./Output/HeatmapOfTheLTPLocationsWithoutLipidInformationWithoutDomainsAdded251120214hdrBlackTrianglesAddedlwd2.pdf",
    width = unit(20, "mm"), height = unit(20, "mm"))

Heatmap(InVivoDataSetslc4hdr[,ReorderedLTPsByManualSeriation], name = LegendName, col = LegendColor, rect_gp = gpar(type = "none"), column_title = GraphNameslc, 
        width = unit(WidthAdaptor*dim(InVivoDataSetslc4hdr)[2]/min(dim(InVivoDataSetslc4hdr)), "mm"), height = unit(HeightAdaptor*dim(InVivoDataSetslc4hdr)[1]/min(dim(InVivoDataSetslc4hdr)), "mm"),
        
        show_heatmap_legend = FALSE,
        cell_fun = function(j, i, x, y, width, height, fill) {
          
          grid.polygon(x = unit(as.numeric(gsub("[a-zA-Z ]", "", x)) + c(-as.numeric(gsub("[a-zA-Z ]", "", width))/2, as.numeric(gsub("[a-zA-Z ]", "", width))/2, -as.numeric(gsub("[a-zA-Z ]", "", width))/2), "npc"), 
                       y = unit(as.numeric(gsub("[a-zA-Z ]", "", y)) + c(-as.numeric(gsub("[a-zA-Z ]", "", height))/2, as.numeric(gsub("[a-zA-Z ]", "", height))/2, as.numeric(gsub("[a-zA-Z ]", "", height))/2), "npc"), 
                       
                       gp = gpar(fill = col_funb(InVivoDataSetslc4hdr[,ReorderedLTPsByManualSeriation][i, j]), 
                                 col = NA))
          
          grid.polygon(x = unit(as.numeric(gsub("[a-zA-Z ]", "", x)) + c(-as.numeric(gsub("[a-zA-Z ]", "", width))/2, as.numeric(gsub("[a-zA-Z ]", "", width))/2, as.numeric(gsub("[a-zA-Z ]", "", width))/2), "npc"), 
                       y = unit(as.numeric(gsub("[a-zA-Z ]", "", y)) + c(-as.numeric(gsub("[a-zA-Z ]", "", height))/2, as.numeric(gsub("[a-zA-Z ]", "", height))/2, -as.numeric(gsub("[a-zA-Z ]", "", height))/2), "npc"), 
                       
                       gp = gpar(fill = col_funo(InVitroDataSetslc4hdr[,ReorderedLTPsByManualSeriation][i, j]), 
                                 col = NA))
          
          grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "LightGrey", fill = NA, lwd = 1))
          
          
          
          if(!is.na(InVivoDataSetslc4hdr[,ReorderedLTPsByManualSeriation][i, j]) && (InVivoDataSetslc4hdr[,ReorderedLTPsByManualSeriation][i, j] == 5)){
            
            grid.polygon(x = unit(as.numeric(gsub("[a-zA-Z ]", "", x)) + c(-as.numeric(gsub("[a-zA-Z ]", "", width))/2, as.numeric(gsub("[a-zA-Z ]", "", width))/2, -as.numeric(gsub("[a-zA-Z ]", "", width))/2), "npc"), 
                         y = unit(as.numeric(gsub("[a-zA-Z ]", "", y)) + c(-as.numeric(gsub("[a-zA-Z ]", "", height))/2, as.numeric(gsub("[a-zA-Z ]", "", height))/2, as.numeric(gsub("[a-zA-Z ]", "", height))/2), "npc"), 
                         
                         gp = gpar(fill = NA, 
                                   col = "black", lwd = 2))}
          
          
          if(!is.na(InVitroDataSetslc4hdr[,ReorderedLTPsByManualSeriation][i, j]) && (InVitroDataSetslc4hdr[,ReorderedLTPsByManualSeriation][i, j] == 5)){
            
            grid.polygon(x = unit(as.numeric(gsub("[a-zA-Z ]", "", x)) + c(-as.numeric(gsub("[a-zA-Z ]", "", width))/2, as.numeric(gsub("[a-zA-Z ]", "", width))/2, as.numeric(gsub("[a-zA-Z ]", "", width))/2), "npc"), 
                         y = unit(as.numeric(gsub("[a-zA-Z ]", "", y)) + c(-as.numeric(gsub("[a-zA-Z ]", "", height))/2, as.numeric(gsub("[a-zA-Z ]", "", height))/2, -as.numeric(gsub("[a-zA-Z ]", "", height))/2), "npc"), 
                         
                         gp = gpar(fill = NA, 
                                   col = "black", lwd = 2))}
          
          grid.circle(x = x, y = y, r = unit(WidthAdaptor*LipidClassLiteratureDataSetslc4hdr[,ReorderedLTPsByManualSeriation][i,j]/(4*min(dim(InVivoDataSetslc[,ReorderedLTPsByManualSeriation]))), "mm"), gp = gpar(fill = "White", col = "Black"))
        }, cluster_rows = FALSE, cluster_columns = FALSE, 
        
        row_split = factor(SplitByRowsVector[8:36]), column_split = factor(SplitByColsVector[c(1:4,9:5,10:17,22,18,25,19,27,26,21,23,20,24,28:29,40,34:35,38,37,39,36,32,33,31,30,43:41)], levels = unique(SplitByColsVector[c(1:4,9:5,10:17,22,18,25,19,27,26,21,23,20,24,28:29,40,34:35,38,37,39,36,32,33,31,30,43:41)])),
        cluster_row_slices = FALSE, cluster_column_slices = FALSE,
        
        row_title = " ", row_title_gp = gpar(fontsize = 10))
dev.off()

# Previous figure: Basis for Figure 2a, with several stylistic updates in Adobe Illustrator and the inclusion of a legend.
# Was in post-production of the figure updated by restricting some of these observations from the figures when uncertainties present.


#### Protein-domain-based ordering of LTPs

#. SupplementaryInformationOverviewSeriationOfProteinDomains061120222.pdf (#)
#. HeatmapListDomainsReorderedTypesDomainsEnhanced06112022.pdf #

library(ComplexHeatmap)
LegendName <- "Legend"

WidthAdaptor <- 160
HeightAdaptor <- 160

# SplitByColsVector <- factor(MainDomainsOfTheLTPs4xx[,2], levels = levels(MainDomainsOfTheLTPs4x[,2])[c(1:3, 9, 7, 4:6,8)]) # Again updated to a version without levels subsetting to avoid future changes (original version on this line)
SplitByColsVector <- factor(MainDomainsOfTheLTPs4xx[,2], levels = c("CRAL-TRIO", "GLTP", "IP_trans", "START", "OSBP", "LBP_BPI_CETP", "lipocalin", "ML", "scp2"))


DomainsLTPsReorderedAfterSeriation <- colnames(DomainsByRowColumnAssociationsReordered5)[c(1:13,15,14,17,16,18:31)]

ReorderedDomainsInputMatrix <- t(as.matrix(DomainsByRowColumnAssociationsReordered5[ReorderedLTPsByManualSeriation, DomainsLTPsReorderedAfterSeriation]))
rownames(ReorderedDomainsInputMatrix)[rownames(ReorderedDomainsInputMatrix) %in% c("CRAL-TRIO", "Ankyrin", "lipocalin", "adh_short", "TMpd", "TMpd2")] <-  c("CRAL_TRIO", "Ankyrin repeats", "Lipocalin", "Adh_short", "TM1", "TM2")

ColFunctionForDomains <- colorRamp2(1:max(CastManyProteinDomainsLTPs, na.rm = TRUE), colorRampPalette(c("#fafa6e", "#2a4858"))(max(CastManyProteinDomainsLTPs, na.rm = TRUE)), space = "RGB")


DomainTypes5 <- cbind(TypeRegion = c(rep("Domain", 18), rep("Region", 9), rep("Motif", 4)),
                      RegionName = c("PRELI/MSF1", "CRAL_TRIO_N", "CRAL-TRIO", "GOLD", "BNIP2", "CRAL_TRIO_2", "GLTP", "IP_trans", "START", "PH", "ORD", "Ankyrin", "LBP_BPI_CETP", "lipocalin", "ML", "Thiolase", "adh_short","SCP2", "TMpd", "TMpd2", "CC1", "CC2", "CC3", "CC4", "Ser-rich", "Ala/Gly-rich", "Poly-Leu", "FFAT", "Signal peptide", "Nuclear localization signal", "PTS1"))

pdf("./Output/HeatmapListDomainsReorderedTypesDomainsEnhanced06112022.pdf",
    width = unit(14, "mm"), height = unit(14, "mm"))

Heatmap(matrix = ReorderedDomainsInputMatrix, 
        col = ColFunctionForDomains, name = "Domain Start" ,cluster_rows = FALSE, cluster_columns = FALSE, 
        
        
        width = unit(WidthAdaptor*dim(ReorderedDomainsInputMatrix)[2]/min(dim(ReorderedDomainsInputMatrix)), "mm"), height = unit(HeightAdaptor*dim(ReorderedDomainsInputMatrix)[1]/min(dim(ReorderedDomainsInputMatrix)), "mm"),
        
        border = "grey", show_column_names = TRUE, 
        column_labels = colnames(ReorderedDomainsInputMatrix), na_col = "white",
        
        column_split = factor(SplitByColsVector[c(1:4,9:5,10:17,22,18,25,19,27,26,21,23,20,24,28:29,40,34:35,38,37,39,36,32,33,31,30,43:41)], levels = unique(SplitByColsVector[c(1:4,9:5,10:17,22,18,25,19,27,26,21,23,20,24,28:29,40,34:35,38,37,39,36,32,33,31,30,43:41)])),
        row_split = factor(DomainTypes5[,1], levels = unique(DomainTypes5[,1])),
        
        cluster_row_slices = FALSE, row_title = " ",
        cluster_column_slices = FALSE, column_title = " ",
        
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "LightGrey", fill = NA, lwd = 1))
          
        })
dev.off()

# Previous figure is the basis for supplementary file on seriation of LTPs based on domains/regions/motifs while taking domain starts also in account
# Is only lightly further updated by adding a broader legend and a title.

#### Fig.6c & Extended Data Fig.6a
#. PITPFocussedPIPCPARemakeGraph25022022LipidsClusterdTogetherWithThinAndThickLinesOverlappingAndWithGuideLinesAddedToDifferentLevelsSeveralFigures28022022.pdf (#)

#. PITPFocussedPIPCPARemakeStripedVerticalLinesGraph140320223d.pdf (#)
#. BarplotHeatmapPITransportersCondensedRowEntriesOnlyScreenColumnsAt5ProcentCutoffAndAnnotationRows12102021.pdf #

LTPLipidConnectionsDataSet <- rbind(PureAntonella32b, PureEnric32)
LipidSpeciesLTPSpecificList <- setNames(lapply(unique(LTPLipidConnectionsDataSet[,"LikelySubclass"]),
                                               
                                               function(y){widen5(inputdf = aggregate(LTPLipidConnectionsDataSet[LTPLipidConnectionsDataSet[,"LikelySubclass"] == y, "Intensity"], 
                                                                                      by=list(paste0(LTPLipidConnectionsDataSet[LTPLipidConnectionsDataSet[,"LikelySubclass"] == y, "LTPProtein"], "(", LTPLipidConnectionsDataSet[LTPLipidConnectionsDataSet[,"LikelySubclass"] == y, "Screen"], ")"),
                                                                                              
                                                                                              LTPLipidConnectionsDataSet[LTPLipidConnectionsDataSet[,"LikelySubclass"] == y, "Lipid"]), 
                                                                                      FUN=function(x){sum(x, na.rm = TRUE)}), 
                                                                  
                                                                  ColumnsLong = "Group.2", ColumnWide = "Group.1", ColumnValue = "x", AggregatingFunction = sum)}),
                                        nm = unique(LTPLipidConnectionsDataSet[,"LikelySubclass"]))

# Again by switch to widen5 here above gives different column order in dataframes in the list, but ok further
LipidSpeciesLTPSpecificList <- lapply(LipidSpeciesLTPSpecificList, function(x){
  
  if(length(grep("vivo", colnames(x))) > 0){x2 <- cbind(x, 'in vivo' = rowSumsDFS(x[,grep("vivo", colnames(x))]))}else{x2 <- x}
  if(length(grep("vitro", colnames(x2))) > 0){x4 <- cbind(x2, 'in vitro' = rowSumsDFS(x2[,grep("vitro", colnames(x2))]))}else{x4 <- x2}
  
  return(x4)})


LipidSpeciesLTPFractionList <- lapply(LipidSpeciesLTPSpecificList, function(y){
  cbind(Lipid = y[,1,drop=F], setNames(as.data.frame(do.call("cbind", lapply(2:dim(y)[2], function(x){100*y[,x]/sum(y[,x], na.rm = TRUE)}))), nm = colnames(y)[2:dim(y)[2]]))})

LipidSpeciesLTPFractionOfTheMaxList <- lapply(LipidSpeciesLTPSpecificList, function(y){
  cbind(Lipid = y[,1,drop=F], setNames(as.data.frame(do.call("cbind", lapply(2:dim(y)[2], function(x){100*y[,x]/max(y[,x], na.rm = TRUE)}))), nm = colnames(y)[2:dim(y)[2]]))})

LipidSubclassesAddedToBackground170620214 <- setNames(lapply(unique(LipidSubclassesAddedToBackground17062021[,1]), function(x){
  if(x %in% names(LipidSpeciesLTPFractionOfTheMaxList)){
    
    IntermediateForEachSubclass <- Col1ToRowNames(merge(LipidSubclassesAddedToBackground170620212[LipidSubclassesAddedToBackground170620212[,2] == x,c(1,4)], LipidSpeciesLTPFractionOfTheMaxList[[x]], by = 1, all = TRUE))
    IntermediateForEachSubclass[,1] <- as.numeric(gsub(",",".", IntermediateForEachSubclass[,1]))*100/max(as.numeric(gsub(",",".", IntermediateForEachSubclass[,1])), na.rm=TRUE)
    
    IntermediateForEachSubclass <- IntermediateForEachSubclass[order(IntermediateForEachSubclass[,1], rowSums(as.matrix(IntermediateForEachSubclass), na.rm = TRUE), decreasing = TRUE),]
    colnames(IntermediateForEachSubclass)[1] <- "Cellular"
    
    return(IntermediateForEachSubclass)}}), nm = unique(LipidSubclassesAddedToBackground17062021[,1]))



LTPLipidConnectionsDataSubset <- cbind(LTPLipidConnectionsDataSet[,c("LTPProtein", "LikelySubclass", "Screen", "Intensity")], CarbonChain = paste(LTPLipidConnectionsDataSet[,"TotalCarbonChainLength"],LTPLipidConnectionsDataSet[,"TotalCarbonChainUnsaturations"], sep = ":"))

# Change to widen5 changes the order of the rows, so I set up a new series of variables here with x at the end. PC come e.g. below PC-O. The widen5 columns are more similarly ordered as before however.
LTPLipidConnectionsDataAggregatedx <- widen5(inputdf = LTPLipidConnectionsDataSubset, ColumnsLong = c("LTPProtein", "LikelySubclass", "Screen"), ColumnWide = "CarbonChain", ColumnValue = "Intensity", AggregatingFunction = function(x){sum(x, na.rm = TRUE)}, RowNamesColumnsLong = FALSE, FillerValue = 0)

PITransportersComobilizedPatternsx <- LTPLipidConnectionsDataAggregatedx[LTPLipidConnectionsDataAggregatedx[, "LTPProtein"] %in% unique(LTPLipidConnectionsDataAggregatedx[LTPLipidConnectionsDataAggregatedx$LikelySubclass == "PI", "LTPProtein"])[c(1,2,6,3,4,5)],]
PITransportersComobilizedPatterns2x <- cbind(PITransportersComobilizedPatternsx[,c("LTPProtein", "LikelySubclass", "Screen")], PITransportersComobilizedPatternsx[,!(colnames(PITransportersComobilizedPatternsx) %in% c("LTPProtein", "LikelySubclass", "Screen", "NaN:NaN"))][,colSums(PITransportersComobilizedPatternsx[,!(colnames(PITransportersComobilizedPatternsx) %in% c("LTPProtein", "LikelySubclass", "Screen", "NaN:NaN"))]) != 0])

PITransportersComobilizedPatterns4x <- PITransportersComobilizedPatterns2x
PITransportersComobilizedPatterns4x[,-(1:3)] <- PITransportersComobilizedPatterns4x[,-(1:3)]*100/do.call(pmax, PITransportersComobilizedPatterns4x[,-(1:3)])

PITransportersComobilizedPatterns4x <- cbind(LTP_Lipid = paste(PITransportersComobilizedPatterns4x[,1], PITransportersComobilizedPatterns4x[,2], sep = "_"), PITransportersComobilizedPatterns4x)
NewRownamesPITransportersx <- unique(PITransportersComobilizedPatterns4x$LTP_Lipid)

GroupingListPITransporters <- list(list(22, 21, 4 , c(1:3,5), 10, c(6:9), 13, 12, 11, 18, 15, 14, c(16:17), 20, 19),
                                   list(NA, NA, NA, c("BPI", NA, "in vivo"), NA, c("GM2A", NA, "in vivo"), NA, NA, NA, NA, NA, NA, c("PITPNB", "PG/BMP", "in vitro"), NA, NA))

PITransportersComobilizedPatterns2CondensedRowEntriesx <- as.matrix(do.call("rbind", lapply(1:length(GroupingListPITransporters[[1]]), function(x){
  if(length(GroupingListPITransporters[[1]][[x]]) == 1){
    
    PITransportersComobilizedPatterns2x[GroupingListPITransporters[[1]][[x]],]
  }else{
    
    c(GroupingListPITransporters[[2]][[x]], colSums(PITransportersComobilizedPatterns2x[GroupingListPITransporters[[1]][[x]],-(1:3)]))
    
    
  }
})))

PITransportersComobilizedPatterns2CondensedRowEntriesx[4,2] <- "Other"
PITransportersComobilizedPatterns2CondensedRowEntriesx[6,2] <- "Other"

PITransportersComobilizedPatterns2CondensedRowEntriesNumericSubsetDatax <- PITransportersComobilizedPatterns2CondensedRowEntriesx[,-(1:3)]
mode(PITransportersComobilizedPatterns2CondensedRowEntriesNumericSubsetDatax) <- "numeric" 

PITransportersComobilizedPatterns4CondensedRowEntriesx <- cbind(PITransportersComobilizedPatterns2CondensedRowEntriesx[,1:3], as.data.frame(PITransportersComobilizedPatterns2CondensedRowEntriesNumericSubsetDatax))
PITransportersComobilizedPatterns4CondensedRowEntriesx[,-(1:3)] <- PITransportersComobilizedPatterns4CondensedRowEntriesx[,-(1:3)]*100/do.call(pmax, PITransportersComobilizedPatterns4CondensedRowEntriesx[,-(1:3)])

PITransportersComobilizedPatterns4CondensedRowEntriesx <- cbind(LTP_Lipid = paste(PITransportersComobilizedPatterns4CondensedRowEntriesx[,1], PITransportersComobilizedPatterns4CondensedRowEntriesx[,2], sep = "_"), PITransportersComobilizedPatterns4CondensedRowEntriesx)
NewRownamesPITransportersCondensedRowEntriesx <- unique(PITransportersComobilizedPatterns4CondensedRowEntriesx$LTP_Lipid)

CellularListPITransportersx <- lapply(sapply(strsplit(sapply(strsplit(as.character(NewRownamesPITransportersx), "_"), "[[", 2), "/"), "[[", 1), function(x){rbind(rownames(LipidSubclassesAddedToBackground170620214[[x]]), 
                                                                                                                                                                  ZerosToNAsConverter(LipidSubclassesAddedToBackground170620214[[x]][,"Cellular"]))})

CellularListPITransportersVectorsWithoutTheNAsx <- lapply(1:length(CellularListPITransportersx), function(x){as.data.frame(t(setNames(as.numeric(CellularListPITransportersx[[x]][2,!is.na(CellularListPITransportersx[[x]][2,])]),
                                                                                                                                      gsub("O-","",sapply(strsplit(gsub(")","",CellularListPITransportersx[[x]][1,!is.na(CellularListPITransportersx[[x]][2,])]), "\\("), "[[", 2)))))})


CarbonChainColNamesx <- sort(unique(c(colnames(PITransportersComobilizedPatterns4x)[-(1:4)], unlist(lapply(CellularListPITransportersVectorsWithoutTheNAsx, names)))))[c(1:54,56:60,55)]

# Implemented consequences of earlier widen as well as introduced RbindMultipleDataframesFillNonmatches-Function to reduce package dependency
PITransportersComobilizedPatterns4WithCellularAddedx <- ZerosToNAsConverter(RbindMultipleDataframesFillNonmatches(list(PITransportersComobilizedPatterns4x, do.call("cbind", list(LTP_Lipid = "cellular",
                                                                                                                                                                                  
                                                                                                                                                                                  LTPProtein = "cellular",
                                                                                                                                                                                  LikelySubclass = sapply(strsplit(sapply(strsplit(as.character(NewRownamesPITransportersx), "_"), "[[", 2), "/"), "[[", 1),
                                                                                                                                                                                  
                                                                                                                                                                                  Screen = "cellular",
                                                                                                                                                                                  RbindMultipleDataframesFillNonmatches(CellularListPITransportersVectorsWithoutTheNAsx))))))[c("LTP_Lipid", "LTPProtein", "LikelySubclass", "Screen", CarbonChainColNamesx)]

# Implemented consequences of earlier widen as well as introduced RbindMultipleDataframesFillNonmatches-Function to reduce package dependency
PITransportersComobilizedPatterns4WithCellularAddedCondensedRowEntriesxx <- ZerosToNAsConverter(RbindMultipleDataframesFillNonmatches(list(PITransportersComobilizedPatterns4CondensedRowEntriesx, do.call("cbind", list(LTP_Lipid = "cellular",
                                                                                                                                                                                                                         
                                                                                                                                                                                                                         LTPProtein = "cellular",
                                                                                                                                                                                                                         LikelySubclass = sapply(strsplit(sapply(strsplit(as.character(NewRownamesPITransportersx), "_"), "[[", 2), "/"), "[[", 1),
                                                                                                                                                                                                                         
                                                                                                                                                                                                                         Screen = "cellular",
                                                                                                                                                                                                                         RbindMultipleDataframesFillNonmatches(CellularListPITransportersVectorsWithoutTheNAsx))))))[c("LTP_Lipid", "LTPProtein", "LikelySubclass", "Screen", CarbonChainColNamesx)]

rownames(PITransportersComobilizedPatterns4WithCellularAddedCondensedRowEntriesxx) <- 1:dim(PITransportersComobilizedPatterns4WithCellularAddedCondensedRowEntriesxx)[1]
PITransportersComobilizedPatterns4WithCellularAddedCondensedRowEntriesx <- PITransportersComobilizedPatterns4WithCellularAddedCondensedRowEntriesxx

PITransportersComobilizedPatterns4WithCellularAdded2x <- unique(PITransportersComobilizedPatterns4WithCellularAddedx)
PITransportersComobilizedPatterns4WithCellularAddedCondensedRowEntries2x <- unique(PITransportersComobilizedPatterns4WithCellularAddedCondensedRowEntriesx)


EmptyLayerPITransportersx <- matrix(NA, 
                                    
                                    nrow = length(NewRownamesPITransportersx), 
                                    ncol = length(CarbonChainColNamesx), 
                                    
                                    dimnames = list(NewRownamesPITransportersx,
                                                    CarbonChainColNamesx))

NewPITransLayersx <- list(InCellulo = EmptyLayerPITransportersx,
                          InVitro = EmptyLayerPITransportersx,
                          
                          Cellular = EmptyLayerPITransportersx)



EmptyLayerPITransportersCondensedRowEntriesx <- matrix(NA, 
                                                       
                                                       nrow = length(NewRownamesPITransportersCondensedRowEntriesx), 
                                                       ncol = length(CarbonChainColNamesx), 
                                                       
                                                       dimnames = list(NewRownamesPITransportersCondensedRowEntriesx,
                                                                       CarbonChainColNamesx))

NewPITransLayersCondensedRowEntriesx <- list(InCellulo = EmptyLayerPITransportersCondensedRowEntriesx,
                                             InVitro = EmptyLayerPITransportersCondensedRowEntriesx,
                                             
                                             Cellular = EmptyLayerPITransportersCondensedRowEntriesx)


for(i in rownames(NewPITransLayersx$InCellulox)){
  if(any((PITransportersComobilizedPatterns4WithCellularAdded2x$Screen == "in vivo") & (PITransportersComobilizedPatterns4WithCellularAdded2x$LTP_Lipid == i))){
    
    NewPITransLayersx$InCellulo[i,] <- as.numeric(ZerosToNAsConverter(PITransportersComobilizedPatterns4WithCellularAdded2x[(PITransportersComobilizedPatterns4WithCellularAdded2x$Screen == "in vivo") & (PITransportersComobilizedPatterns4WithCellularAdded2x$LTP_Lipid == i),-(1:4)]))
  }}

for(i in rownames(NewPITransLayersx$InVitro)){
  if(any((PITransportersComobilizedPatterns4WithCellularAdded2x$Screen == "in vitro") & (PITransportersComobilizedPatterns4WithCellularAdded2x$LTP_Lipid == i))){
    
    NewPITransLayersx$InVitro[i,] <- as.numeric(ZerosToNAsConverter(PITransportersComobilizedPatterns4WithCellularAdded2x[(PITransportersComobilizedPatterns4WithCellularAdded2x$Screen == "in vitro") & (PITransportersComobilizedPatterns4WithCellularAdded2x$LTP_Lipid == i),-(1:4)]))
  }}

for(i in rownames(NewPITransLayersx$Cellular)){
  NewPITransLayersx$Cellular[i,] <- as.numeric(ZerosToNAsConverter(PITransportersComobilizedPatterns4WithCellularAdded2x[(PITransportersComobilizedPatterns4WithCellularAdded2x$Screen == "cellular") & (PITransportersComobilizedPatterns4WithCellularAdded2x$LikelySubclass == strsplit(strsplit(i,"_")[[1]][2],"/")[[1]][1]),-(1:4)]))
  
}


for(i in rownames(NewPITransLayersCondensedRowEntriesx$InCellulo)){
  if(any((PITransportersComobilizedPatterns4WithCellularAddedCondensedRowEntries2x$Screen == "in vivo") & (PITransportersComobilizedPatterns4WithCellularAddedCondensedRowEntries2x$LTP_Lipid == i))){
    
    NewPITransLayersCondensedRowEntriesx$InCellulo[i,] <- as.numeric(ZerosToNAsConverter(PITransportersComobilizedPatterns4WithCellularAddedCondensedRowEntries2x[(PITransportersComobilizedPatterns4WithCellularAddedCondensedRowEntries2x$Screen == "in vivo") & (PITransportersComobilizedPatterns4WithCellularAddedCondensedRowEntries2x$LTP_Lipid == i),-(1:4)]))
  }}

for(i in rownames(NewPITransLayersCondensedRowEntriesx$InVitro)){
  if(any((PITransportersComobilizedPatterns4WithCellularAddedCondensedRowEntries2x$Screen == "in vitro") & (PITransportersComobilizedPatterns4WithCellularAddedCondensedRowEntries2x$LTP_Lipid == i))){
    
    NewPITransLayersCondensedRowEntriesx$InVitro[i,] <- as.numeric(ZerosToNAsConverter(PITransportersComobilizedPatterns4WithCellularAddedCondensedRowEntries2x[(PITransportersComobilizedPatterns4WithCellularAddedCondensedRowEntries2x$Screen == "in vitro") & (PITransportersComobilizedPatterns4WithCellularAddedCondensedRowEntries2x$LTP_Lipid == i),-(1:4)]))
  }}

GroupingListPITransportersAtLaterStagex <- list(20, 19, 4 , c(1:3,5), 10, c(6:9), 12, 11, 16, 13, c(14:15), 18, 17)
NewPITransLayersCondensedRowEntriesx$Cellular <- do.call("rbind", lapply(GroupingListPITransportersAtLaterStagex, function(x){
  
  if(length(unlist(x)) == 1){
    NewPITransLayersx$Cellular[unlist(x),]
    
  }else{
    colSums(NewPITransLayersx[["Cellular"]][unlist(x),], na.rm = TRUE)*100/max(colSums(NewPITransLayersx[["Cellular"]][unlist(x),], na.rm = TRUE))
    
  }
}))

rownames(NewPITransLayersCondensedRowEntriesx$Cellular) <- rownames(NewPITransLayersCondensedRowEntriesx$InCellulo)
NewPITransLayersCondensedRowEntriesx$Cellular <- ZerosToNAsConverter(NewPITransLayersCondensedRowEntriesx$Cellular)

library(RColorBrewer)
library(ComplexHeatmap)

LegendName <- "Legend"
LegendColor <- col_fung

WidthAdaptor <- 160
HeightAdaptor <- 160

TextDataframeWithLTPAnnotationx <- data.frame(LTP = factor(sapply(strsplit(rownames(NewPITransLayersCondensedRowEntriesx$InCellulo), "_"), "[[", 1), levels = unique(sapply(strsplit(rownames(NewPITransLayersCondensedRowEntriesx$InCellulo), "_"), "[[", 1))))
AnnotationDataframeWithLTPAnnotationx <- rowAnnotation(df = TextDataframeWithLTPAnnotationx, col = list(LTP = c("SEC14L2" = "#7E549F", "BPI" = "#FB836F", "GM2A" = "#C1549C", "PITPNA" = "#FFCB3E", "PITPNB" = "#E0B01C", "PITPNC1" = "#A47C00")))

CutOffPresenceProcent <- 5
CSPIDCREx <- (colSums(NewPITransLayersCondensedRowEntriesx$InCellulo, na.rm = TRUE) > CutOffPresenceProcent)|(colSums(NewPITransLayersCondensedRowEntriesx$InVitro, na.rm = TRUE) > CutOffPresenceProcent)

pdf(paste0("./Output/BarplotHeatmapPITransportersCondensedRowEntriesOnlyScreenColumnsAt", CutOffPresenceProcent,"ProcentCutoffAndAnnotationRows12102021.pdf"),
    width = unit(20, "mm"), height = unit(20, "mm"))

Heatmap(NewPITransLayersCondensedRowEntriesx$InCellulo[,CSPIDCREx], name = LegendName, col = LegendColor, rect_gp = gpar(type = "none"), column_title = "", 
        width = unit(WidthAdaptor*dim(NewPITransLayersCondensedRowEntriesx$InCellulo[,CSPIDCREx])[2]/min(dim(NewPITransLayersCondensedRowEntriesx$InCellulo[,CSPIDCREx])), "mm"), height = unit(HeightAdaptor*dim(NewPITransLayersCondensedRowEntriesx$InCellulo[,CSPIDCREx])[1]/min(dim(NewPITransLayersCondensedRowEntriesx$InCellulo[,CSPIDCREx])), "mm"),
        
        
        cell_fun = function(j, i, x, y, width, height, fill) {
          
          grid.rect(x = x, y = y-0.5*height, width = width, height = NewPITransLayersCondensedRowEntriesx$InCellulo[,CSPIDCREx][i,j]/100*height, gp = gpar(col = NA, fill = brewer.pal(9,"Blues")[5], lwd = 1), just = c("center","bottom"))
          grid.rect(x = x, y = y-0.5*height, width = 0.5*width, height = NewPITransLayersCondensedRowEntriesx$InVitro[,CSPIDCREx][i,j]/100*height, gp = gpar(col = NA, fill = brewer.pal(9,"Oranges")[5], lwd = 1), just = c("center","bottom"))
          
          grid.rect(x = x, y = y-0.5*height, width = width, height = NewPITransLayersCondensedRowEntriesx$Cellular[,CSPIDCREx][i,j]/100*height, gp = gpar(col = "LightGrey", fill = NA, lwd = 2), just = c("center","bottom"))
        }, cluster_rows = FALSE, cluster_columns = FALSE, 
        
        
        cluster_row_slices = FALSE, cluster_column_slices = FALSE,
        
        row_title = " ", row_title_gp = gpar(fontsize = 10),
        left_annotation = AnnotationDataframeWithLTPAnnotationx,
        
        row_labels = sapply(strsplit(rownames(NewPITransLayersCondensedRowEntriesx$InCellulo), "_"), "[[", 2),
        row_split = c(rep("A", 2), rep("B",2), rep("C",2), rep("D",2), rep("E",3), rep("F",2))
        
)
dev.off()

# Previous graph further in Adobe Illustrator used to subset for panel Figure 6c and for the Extended Data Figures 6a
# Note: only anything above 5% is visualized here to reduce the complexity of the figures

######## Figure 3
#### Fig.3a

#. CircosExtendedWithSpecies051120217WithTracksSwitchedDecreasedTextSize0412WhiteRingRemovedSmallerScreenTracksWithoutInnerHalfringByDiffHeight0ColorSchemeUpdatedWithoutSomeTLCEntriesWithBlackOutsides8.pdf #


PureAntonella32bs0 <- PureAntonella32b[, c("LTPProtein", "ProteinDomain", "Screen", "TotalCarbonChainLength", "TotalCarbonChainUnsaturations", "IonMode", "LikelySubclass", "Intensity")]
PureAntonella32bs0[is.na(PureAntonella32bs0)] <- 0

PureEnric32s0 <- PureEnric32[, c("LTPProtein", "ProteinDomain", "Screen", "TotalCarbonChainLength", "TotalCarbonChainUnsaturations", "IonMode", "LikelySubclass", "Intensity")]
PureEnric32s0[is.na(PureEnric32s0)] <- 0

AggregatedInCelluloMS <- aggregate(PureAntonella32bs0[, "Intensity"], 
                                   by = PureAntonella32bs0[, c("LTPProtein", "ProteinDomain", "Screen", "TotalCarbonChainLength", "TotalCarbonChainUnsaturations", "IonMode", "LikelySubclass")], 
                                   
                                   FUN = function(x){sum(x,na.rm=TRUE)})


AggregatedInVitroMS <- aggregate(PureEnric32s0[, "Intensity"], 
                                 by = PureEnric32s0[, c("LTPProtein", "ProteinDomain", "Screen", "TotalCarbonChainLength", "TotalCarbonChainUnsaturations", "IonMode", "LikelySubclass")], 
                                 
                                 FUN = function(x){sum(x,na.rm=TRUE)})


AggregatedInCelluloMS2 <- rbind(cbind(AggregatedInCelluloMS[AggregatedInCelluloMS[,"IonMode"] == "pos", ], NormInt = MinMaxNormMatrixFunc(inputdata = AggregatedInCelluloMS[AggregatedInCelluloMS[,"IonMode"] == "pos", "x"])*9+1),
                                cbind(AggregatedInCelluloMS[AggregatedInCelluloMS[,"IonMode"] == "neg", ], NormInt = MinMaxNormMatrixFunc(inputdata = AggregatedInCelluloMS[AggregatedInCelluloMS[,"IonMode"] == "neg", "x"])*9+1))

AggregatedInVitroMS2 <- rbind(cbind(AggregatedInVitroMS[AggregatedInVitroMS[,"IonMode"] == "pos", ], NormInt = MinMaxNormMatrixFunc(inputdata = AggregatedInVitroMS[AggregatedInVitroMS[,"IonMode"] == "pos", "x"])*9+1),
                              cbind(AggregatedInVitroMS[AggregatedInVitroMS[,"IonMode"] == "neg", ], NormInt = MinMaxNormMatrixFunc(inputdata = AggregatedInVitroMS[AggregatedInVitroMS[,"IonMode"] == "neg", "x"])*9+1))


# Aggregate pos & neg by max

AggregatedInCelluloMS4 <- aggregate(AggregatedInCelluloMS2[, "NormInt"], 
                                    by = AggregatedInCelluloMS2[, c("LTPProtein", "ProteinDomain", "Screen", "TotalCarbonChainLength", "TotalCarbonChainUnsaturations", "LikelySubclass")], 
                                    
                                    FUN = function(x){max(x,na.rm=TRUE)})


AggregatedInVitroMS4 <- aggregate(AggregatedInVitroMS2[, "NormInt"], 
                                  by = AggregatedInVitroMS2[, c("LTPProtein", "ProteinDomain", "Screen", "TotalCarbonChainLength", "TotalCarbonChainUnsaturations", "LikelySubclass")], 
                                  
                                  FUN = function(x){max(x,na.rm=TRUE)})


HPTLCaggregated <- rbind(cbind(meltmatrix(as.matrix(5*HPTLCSpecificitiesPerScreen2[[1]])), ProteinDomain = NA, Screen = "in vivo", TotalCarbonChainLength = NA, TotalCarbonChainUnsaturations = NA),
                         cbind(meltmatrix(as.matrix(5*HPTLCSpecificitiesPerScreen2[[2]])), ProteinDomain = NA, Screen = "in vitro", TotalCarbonChainLength = NA, TotalCarbonChainUnsaturations = NA))

HPTLCaggregated <- HPTLCaggregated[HPTLCaggregated[,"value"] != 0,]
HPTLCaggregated[,2] <- gsub("\\*","", HPTLCaggregated[,2])

colnames(HPTLCaggregated)[1:3] <- c("LTPProtein", "LikelySubclass", "x")
HPTLCaggregated <- HPTLCaggregated[,colnames(AggregatedInVitroMS4)]

HPTLCaggregated$ProteinDomain <- c("START",rep("OSBP",9), rep("START",3), "OSBP", "IP_trans", "IP_trans", "START", "IP_trans", "OSBP")
HPTLCaggregated[HPTLCaggregated$LikelySubclass == "PG","LikelySubclass"] <- "PG/BMP"

HPTLCaggregated[HPTLCaggregated$LikelySubclass == "Cer","LikelySubclass"] <- "Cer*"
HPTLCaggregated[HPTLCaggregated$LikelySubclass == "PIP","LikelySubclass"] <- "PIPs"

HPTLCaggregated[HPTLCaggregated$LikelySubclass == "CH","LikelySubclass"] <- "Sterol"
HPTLCaggregated[is.na(HPTLCaggregated)] <- 0

AggregatedAllScreenDataNormalized <- do.call("rbind", list(AggregatedInCelluloMS4, AggregatedInVitroMS4, HPTLCaggregated))
colnames(AggregatedAllScreenDataNormalized)[7] <- "NormInt"

AggregatedAllScreenDataNormalized2 <- aggregate(AggregatedAllScreenDataNormalized[, "NormInt"], 
                                                by = AggregatedAllScreenDataNormalized[, c("LTPProtein", "ProteinDomain", "Screen", "TotalCarbonChainLength", "TotalCarbonChainUnsaturations", "LikelySubclass")], 
                                                
                                                FUN = function(x){max(x, na.rm=TRUE)})
# widen5 gives slightly different orders below of rows and columns: were different indicated with x or xx added, until variable are similar again.

AggregatedAllScreenDataNormalized4x <- widen5(inputdf = AggregatedAllScreenDataNormalized2, ColumnsLong = c("LTPProtein", "ProteinDomain", "TotalCarbonChainLength", "TotalCarbonChainUnsaturations", "LikelySubclass"), ColumnWide = "Screen", ColumnValue = "x", AggregatingFunction = sum, FunctionOutputValueType = double(1), RowNamesColumnsLong = FALSE)
AggregatedAllScreenDataNormalized4x[,1] <- factor(AggregatedAllScreenDataNormalized4x[,1], levels = ReorderedLTPsByManualSeriation)

AggregatedAllScreenDataNormalized4x$LikelySubclass <- factor(AggregatedAllScreenDataNormalized4x$LikelySubclass, levels = c("Cer*", "d*Cer", "dCer", "DHCer", "DHOH*Cer", "tCer", "d*CerP", "d*HexCer", "t*HexCer", "t*Hex2Cer", "d*SHexCer", "d*SM", "DHSM", "t*SM", "FA",       
                                                                                                                            "FAL", "LPC", "LPE", "LPE-O", "LPG", "PA", "PC", "PC-O", "PE", "PE-O", "PI", "PIPs", "PS", "PGP", "PG", "PG/BMP", "BMP", "CL", "DAG", "TAG", "Sterol", "VA") )

AggregatedAllScreenDataNormalized4x$BandSize <- 1 


# Descrepancy in row sorting with row with unsaturations at 10 before 8 --> corrected
x <- AggregatedAllScreenDataNormalized4x[579,] 

AggregatedAllScreenDataNormalized4x[579,] <- AggregatedAllScreenDataNormalized4x[580,]
AggregatedAllScreenDataNormalized4x[580,] <- x

# Other differences changed to original order in xx version: in vivo after in vitro column, and for OSBPL2 PIPs after Sterol, so all is the same as original
AggregatedAllScreenDataNormalized4xx <- AggregatedAllScreenDataNormalized4x[c(1:213,215,214,216:nrow(AggregatedAllScreenDataNormalized4x)),c(1:5,7,6,8)]

rownames(AggregatedAllScreenDataNormalized4xx) <- 1:dim(AggregatedAllScreenDataNormalized4xx)[1]
AggregatedAllScreenDataNormalized4 <- AggregatedAllScreenDataNormalized4xx

HPTLCRowsInData <- which(AggregatedAllScreenDataNormalized4$TotalCarbonChainLength == 0)
AggregatedAllScreenDataNormalized4ifh <- AggregatedAllScreenDataNormalized4

HPTLCRowsInDataNAsToBeIntroduced <- do.call("cbind", list(AggregatedAllScreenDataNormalized4ifh[HPTLCRowsInData,],
                                                          RowToBeEvaluated = HPTLCRowsInData,
                                                          
                                                          InVitroRemovalData = sapply(HPTLCRowsInData, function(x){paste(AggregatedAllScreenDataNormalized4ifh[x,1], AggregatedAllScreenDataNormalized4ifh[x,5], sep = "_") %in% paste(AggregatedAllScreenDataNormalized4ifh[-x,][!is.na(AggregatedAllScreenDataNormalized4ifh[-x,][,"in vitro"]), ][,1], AggregatedAllScreenDataNormalized4ifh[-x,][!is.na(AggregatedAllScreenDataNormalized4ifh[-x,][,"in vitro"]), ][,5], sep = "_")}),
                                                          InVivoRemovalData = sapply(HPTLCRowsInData, function(x){paste(AggregatedAllScreenDataNormalized4ifh[x,1], AggregatedAllScreenDataNormalized4ifh[x,5], sep = "_") %in% paste(AggregatedAllScreenDataNormalized4ifh[-x,][!is.na(AggregatedAllScreenDataNormalized4ifh[-x,][,"in vivo"]), ][,1], AggregatedAllScreenDataNormalized4ifh[-x,][!is.na(AggregatedAllScreenDataNormalized4ifh[-x,][,"in vivo"]), ][,5], sep = "_")})))

AggregatedAllScreenDataNormalized4ifh[HPTLCRowsInDataNAsToBeIntroduced[HPTLCRowsInDataNAsToBeIntroduced[,"InVitroRemovalData"],"RowToBeEvaluated"], "in vitro"] <- NA
AggregatedAllScreenDataNormalized4ifh[HPTLCRowsInDataNAsToBeIntroduced[HPTLCRowsInDataNAsToBeIntroduced[,"InVivoRemovalData"],"RowToBeEvaluated"], "in vivo"] <- NA

# Has the corrected extra 2 entries for OSBPL9-PS and STARD10-PE: adapted to x version
AggregatedAllScreenDataNormalized4hdrx <- AggregatedAllScreenDataNormalized4ifh[!(is.na(AggregatedAllScreenDataNormalized4ifh[,"in vivo"]) & is.na(AggregatedAllScreenDataNormalized4ifh[,"in vitro"])),]

InCelluloInVitroPresence <- (InVivoDataSetslc != 0)|(InVitroDataSetslc != 0)
mode(InCelluloInVitroPresence) <- "numeric"


# Input of lipid colors

LipidColorsForCircos <- read.csv(file = "./InputData/LipidColorsUpdatesForCircos011120212.txt", header = TRUE, sep = "\t", as.is = TRUE)[,1:2]
LipidColorsForCircos[,2] <- paste0("#", LipidColorsForCircos[,2])

LipidColorsForCircos2 <- LipidColorsForCircos[,2]
names(LipidColorsForCircos2) <- LipidColorsForCircos[,1]

LipidColorsForCircos4 <- c(LipidColorsForCircos2, Sterol = "#100000")
LipidColorsForCircos4b <- LipidColorsForCircos4

LipidColorsForCircos4b["PI"] <- "#DFBF29"
LipidColorsForCircos4b["PIPs"] <- "#A0904D"

LipidColorsForCircos4b["PS"] <- "#F6E839"
all(unique(AggregatedAllScreenDataNormalized4[,5]) %in% names(LipidColorsForCircos4b)) # Check for missing color definitions # TRUE --> no missing definitions


# Actual circos visualization at species-level

library(circlize) # CircosExtendedWithSpecies051120217WithTracksSwitchedDecreasedTextSize0412WhiteRingRemovedSmallerScreenTracksWithoutInnerHalfringByDiffHeight0ColorSchemeUpdatedWithoutSomeTLCEntriesWithBlackOutsides7b.pdf
pdf("./Output/CEWS051120217WTSDTS0412WRRSSTWIHBDH0CSUWSTLCEWBO7b.pdf")

circos.clear()
circos.par(start.degree = -70, clock.wise = FALSE, cell.padding = c(0,0,0,0))

OSOW <- getOption("warn")
options(warn = -1)

testchord2x <- chordDiagram(AggregatedAllScreenDataNormalized4hdrx[,c(5,1,8)], annotationTrack = NULL, grid.col = LipidColorsForCircos4b, directional = -1, diffHeight = mm_h(0),
                            preAllocateTracks = list(list(track.height = 0.005),
                                                     
                                                     list(track.height = 0.16),
                                                     list(track.height = 0.05)),
                            
                            order = c(c("Cer*", "d*Cer", "dCer", "DHCer", "DHOH*Cer", "tCer", "d*CerP", "d*HexCer", "t*HexCer", "t*Hex2Cer", "d*SHexCer", "d*SM", "DHSM", "t*SM",
                                        "FA", "FAL", "LPC", "LPE", "LPE-O", "LPG", "PA", "PC", "PC-O", "PE", "PE-O", "PI", "PIPs", "PS", "PGP", "PG", "PG/BMP", "BMP", 
                                        
                                        "CL", "DAG", "TAG", "Sterol", "VA"),
                                      ReorderedLTPsByManualSeriation))

circos.track(track.index = 2, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.41)
}, bg.border = NA) 


for(i in list(1:9,10:11,12:14,15:17,18:27,28:29,30,31:40,41:43)){
  
  highlight.sector(colnames(InCelluloInVitroPresence[,ReorderedLTPsByManualSeriation])[i], track.index = 1, col = "grey", 
                   text = "", cex = 0.8, text.col = "white", niceFacing = TRUE)
  
}
ylim <- c(0, 1) 

y1 <- ylim[1]
y2 <- ylim[2]

cdmresaddedvaluesx <- cbind(testchord2x, AggregatedAllScreenDataNormalized4hdrx)
colnames(cdmresaddedvaluesx)[colnames(cdmresaddedvaluesx) == "in vivo"] <- "in cellulo"

suppressMessages(for(i in seq_len(nrow(cdmresaddedvaluesx))) {
  if(cdmresaddedvaluesx$value1[i] > 0) {
    
    circos.rect(cdmresaddedvaluesx[i, "x1"], y1, cdmresaddedvaluesx[i, "x1"] - abs(cdmresaddedvaluesx[i, "value1"]), y1 + (y2-y1)*0.45, 
                col = "#ededed", 
                
                border = "#ededed",
                sector.index = cdmresaddedvaluesx$rn[i], track.index = 3)
    
    circos.rect(cdmresaddedvaluesx[i, "x1"], y1 + (y2-y1)*0.55, cdmresaddedvaluesx[i, "x1"] - abs(cdmresaddedvaluesx[i, "value1"]), y2, 
                col = "#ededed", 
                
                border = "#ededed",
                sector.index = cdmresaddedvaluesx$rn[i], track.index = 3)
    
    circos.rect(cdmresaddedvaluesx[i, "x2"], y1, cdmresaddedvaluesx[i, "x2"] - abs(cdmresaddedvaluesx[i, "value1"]), y1 + (y2-y1)*0.45, 
                col = "#ededed", 
                
                border = "#ededed",
                sector.index = cdmresaddedvaluesx$cn[i], track.index = 3)
    
    circos.rect(cdmresaddedvaluesx[i, "x2"], y1 + (y2-y1)*0.55, cdmresaddedvaluesx[i, "x2"] - abs(cdmresaddedvaluesx[i, "value1"]), y2, 
                col = "#ededed", 
                
                border = "#ededed",
                sector.index = cdmresaddedvaluesx$cn[i], track.index = 3)
    
    circos.rect(cdmresaddedvaluesx[i, "x1"], y1, cdmresaddedvaluesx[i, "x1"] - abs(cdmresaddedvaluesx[i, "value1"]), y1 + (y2-y1)*0.45, 
                col = col_funo(cdmresaddedvaluesx[i, "in vitro"]), 
                
                border = col_funo(cdmresaddedvaluesx[i, "in vitro"]),
                sector.index = cdmresaddedvaluesx$rn[i], track.index = 3)
    
    circos.rect(cdmresaddedvaluesx[i, "x1"], y1 + (y2-y1)*0.55, cdmresaddedvaluesx[i, "x1"] - abs(cdmresaddedvaluesx[i, "value1"]), y2, 
                col = col_funb(cdmresaddedvaluesx[i, "in cellulo"]), 
                
                border = col_funb(cdmresaddedvaluesx[i, "in cellulo"]),
                sector.index = cdmresaddedvaluesx$rn[i], track.index = 3)
    
    circos.rect(cdmresaddedvaluesx[i, "x2"], y1, cdmresaddedvaluesx[i, "x2"] - abs(cdmresaddedvaluesx[i, "value1"]), y1 + (y2-y1)*0.45, 
                col = col_funo(cdmresaddedvaluesx[i, "in vitro"]), 
                
                border = col_funo(cdmresaddedvaluesx[i, "in vitro"]),
                sector.index = cdmresaddedvaluesx$cn[i], track.index = 3)
    
    circos.rect(cdmresaddedvaluesx[i, "x2"], y1 + (y2-y1)*0.55, cdmresaddedvaluesx[i, "x2"] - abs(cdmresaddedvaluesx[i, "value1"]), y2, 
                col = col_funb(cdmresaddedvaluesx[i, "in cellulo"]), 
                
                border = col_funb(cdmresaddedvaluesx[i, "in cellulo"]),
                sector.index = cdmresaddedvaluesx$cn[i], track.index = 3)
    
    
    if((cdmresaddedvaluesx[i, "TotalCarbonChainLength"] == 0) & !is.na(cdmresaddedvaluesx[i, "in vitro"]) & (cdmresaddedvaluesx[i, "rn"] != "VA")){
      
      circos.rect(cdmresaddedvaluesx[i, "x1"], y1-(y2-y1)*0.2, cdmresaddedvaluesx[i, "x1"] - abs(cdmresaddedvaluesx[i, "value1"]), y1, 
                  col = "black", 
                  
                  border = "black",
                  sector.index = cdmresaddedvaluesx$rn[i], track.index = 3)
      
      circos.rect(cdmresaddedvaluesx[i, "x2"], y1-(y2-y1)*0.2, cdmresaddedvaluesx[i, "x2"] - abs(cdmresaddedvaluesx[i, "value1"]), y1, 
                  col = "black", 
                  
                  border = "black",
                  sector.index = cdmresaddedvaluesx$cn[i], track.index = 3)
      
    }
    
    
    
    if((cdmresaddedvaluesx[i, "TotalCarbonChainLength"] == 0) & !is.na(cdmresaddedvaluesx[i, "in cellulo"]) & (cdmresaddedvaluesx[i, "rn"] != "VA")){
      
      circos.rect(cdmresaddedvaluesx[i, "x1"], (y2+0.2), cdmresaddedvaluesx[i, "x1"] - abs(cdmresaddedvaluesx[i, "value1"]), y2, 
                  col = "black", 
                  
                  border = "black",
                  sector.index = cdmresaddedvaluesx$rn[i], track.index = 3)
      
      circos.rect(cdmresaddedvaluesx[i, "x2"], (y2+0.2), cdmresaddedvaluesx[i, "x2"] - abs(cdmresaddedvaluesx[i, "value1"]), y2, 
                  col = "black", 
                  
                  border = "black",
                  sector.index = cdmresaddedvaluesx$cn[i], track.index = 3)
      
    }
  }  
  
}) # Messages suppressed because of on purpose plotting outside of sectors
options(warn = OSOW)

dev.off()
# Previous figure: Basis for figure panel 3a, after some stylistic enhancements in Adobe Illustrator, such as changes (of orientation of) some labels and optimization of some of the colors

#### (Subclasses version: not used in final figure.)
#. CircosExtended031120215SterolChangedAndTracksSwappedDifferentTrackColorsLinkColorsUpdatedRemovalWhiteSectionsWithoutSomeTLCEntries24112021BlackOuterBlocksWithLegendAdded.pdf (#)

#. CircosExtended031120215SterolChangedAndTracksSwappedDifferentTrackColorsLinkColorsUpdatedRemovalWhiteSectionsWithoutSomeTLCEntries24112021BlackOuterBlocks.pdf #


InCelluloInVitroPresence2 <- InCelluloInVitroPresence
rownames(InCelluloInVitroPresence2)[rownames(InCelluloInVitroPresence2) == "CH"] <- "Sterol"

InCelluloInVitroPresence2hdr <- InCelluloInVitroPresence2


# Updated version of sub-class - LTP connections circos with changed colors and removed white circles for TLC-reduced dataset, based on data that only has TLC if no other data there for the specific combo
pdf("./Output/CircosExtended031120215SterolChangedAndTracksSwappedDifferentTrackColorsLinkColorsUpdatedRemovalWhiteSectionsWithoutSomeTLCEntries24112021BlackOuterBlocks.pdf")

circos.clear()
circos.par(start.degree = -70, clock.wise = FALSE, cell.padding = c(0,0,0,0))

testchord <- chordDiagram(InCelluloInVitroPresence2hdr[,ReorderedLTPsByManualSeriation], annotationTrack = NULL, grid.col = LipidColorsForCircos4b, directional = -1, diffHeight = mm_h(0),
                          preAllocateTracks = list(list(track.height = 0.005),
                                                   
                                                   list(track.height = 0.16),
                                                   list(track.height = 0.05)))

circos.track(track.index = 2, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.64)
}, bg.border = NA) 


for(i in list(1:9,10:11,12:14,15:17,18:27,28:29,30,31:40,41:43)){
  
  highlight.sector(colnames(InCelluloInVitroPresence2hdr[,ReorderedLTPsByManualSeriation])[i], track.index = 1, col = "grey", 
                   text = "", cex = 0.8, text.col = "white", niceFacing = TRUE)
  
}


cdm_res <- testchord
ylim <- get.cell.meta.data("ylim", sector.index = rownames(InCelluloInVitroPresence2hdr[,ReorderedLTPsByManualSeriation])[1], track.index = 3)

y1 <- ylim[1]
y2 <- ylim[2]


i <- 10

suppressMessages(for(i in seq_len(nrow(cdm_res))) {
  if(cdm_res$value1[i] > 0) {
    
    circos.rect(cdm_res[i, "x1"], y1, cdm_res[i, "x1"] - abs(cdm_res[i, "value1"]), y1 + (y2-y1)*0.45, 
                col = "#ededed", 
                
                border = "#ededed",
                sector.index = cdm_res$rn[i], track.index = 3)
    
    circos.rect(cdm_res[i, "x1"], y1 + (y2-y1)*0.55, cdm_res[i, "x1"] - abs(cdm_res[i, "value1"]), y2, 
                col = "#ededed", 
                
                border = "#ededed",
                sector.index = cdm_res$rn[i], track.index = 3)
    
    circos.rect(cdm_res[i, "x2"], y1, cdm_res[i, "x2"] - abs(cdm_res[i, "value1"]), y1 + (y2-y1)*0.45, 
                col = "#ededed", 
                
                border = "#ededed",
                sector.index = cdm_res$cn[i], track.index = 3)
    
    circos.rect(cdm_res[i, "x2"], y1 + (y2-y1)*0.55, cdm_res[i, "x2"] - abs(cdm_res[i, "value1"]), y2, 
                col = "#ededed", 
                
                border = "#ededed",
                sector.index = cdm_res$cn[i], track.index = 3)
    
    circos.rect(cdm_res[i, "x1"], y1, cdm_res[i, "x1"] - abs(cdm_res[i, "value1"]), y1 + (y2-y1)*0.45, 
                col = col_funo(InVitroDataSetslc4hdr[,ReorderedLTPsByManualSeriation][cdm_res$rn[i], cdm_res$cn[i]]), 
                
                border = col_funo(InVitroDataSetslc4hdr[,ReorderedLTPsByManualSeriation][cdm_res$rn[i], cdm_res$cn[i]]),
                sector.index = cdm_res$rn[i], track.index = 3)
    
    circos.rect(cdm_res[i, "x1"], y1 + (y2-y1)*0.55, cdm_res[i, "x1"] - abs(cdm_res[i, "value1"]), y2, 
                col = col_funb(InVivoDataSetslc4hdr[,ReorderedLTPsByManualSeriation][cdm_res$rn[i], cdm_res$cn[i]]), 
                
                border =  col_funb(InVivoDataSetslc4hdr[,ReorderedLTPsByManualSeriation][cdm_res$rn[i], cdm_res$cn[i]]), 
                sector.index = cdm_res$rn[i], track.index = 3)
    
    circos.rect(cdm_res[i, "x2"], y1, cdm_res[i, "x2"] - abs(cdm_res[i, "value1"]), y1 + (y2-y1)*0.45, 
                col = col_funo(InVitroDataSetslc4hdr[,ReorderedLTPsByManualSeriation][cdm_res$rn[i], cdm_res$cn[i]]),
                
                border =  col_funo(InVitroDataSetslc4hdr[,ReorderedLTPsByManualSeriation][cdm_res$rn[i], cdm_res$cn[i]]),
                sector.index = cdm_res$cn[i], track.index = 3)
    
    circos.rect(cdm_res[i, "x2"], y1 + (y2-y1)*0.55, cdm_res[i, "x2"] - abs(cdm_res[i, "value1"]), y2, 
                col = col_funb(InVivoDataSetslc4hdr[,ReorderedLTPsByManualSeriation][cdm_res$rn[i], cdm_res$cn[i]]), 
                
                border =  col_funb(InVivoDataSetslc4hdr[,ReorderedLTPsByManualSeriation][cdm_res$rn[i], cdm_res$cn[i]]), 
                sector.index = cdm_res$cn[i], track.index = 3)
    
    if((!is.na(InVitroDataSetslc4hdr[,ReorderedLTPsByManualSeriation][cdm_res$rn[i], cdm_res$cn[i]])) && 
       (InVitroDataSetslc4hdr[,ReorderedLTPsByManualSeriation][cdm_res$rn[i], cdm_res$cn[i]] == 5)){
      
      circos.rect(cdm_res[i, "x1"], y1-(y2-y1)*0.2, cdm_res[i, "x1"] - abs(cdm_res[i, "value1"]), y1, 
                  col = "black", 
                  
                  border = "black",
                  sector.index = cdm_res$rn[i], track.index = 3)
      
      circos.rect(cdm_res[i, "x2"], y1-(y2-y1)*0.2, cdm_res[i, "x2"] - abs(cdm_res[i, "value1"]), y1, 
                  col = "black", 
                  
                  border = "black",
                  sector.index = cdm_res$cn[i], track.index = 3)
      
    }
    
    
    if((!is.na(InVivoDataSetslc4hdr[,ReorderedLTPsByManualSeriation][cdm_res$rn[i], cdm_res$cn[i]])) && 
       (InVivoDataSetslc4hdr[,ReorderedLTPsByManualSeriation][cdm_res$rn[i], cdm_res$cn[i]] == 5)){
      
      circos.rect(cdm_res[i, "x1"], (y2+0.2), cdm_res[i, "x1"] - abs(cdm_res[i, "value1"]), y2, 
                  col = "black", 
                  
                  border = "black",
                  sector.index = cdm_res$rn[i], track.index = 3)
      
      circos.rect(cdm_res[i, "x2"], (y2+0.2), cdm_res[i, "x2"] - abs(cdm_res[i, "value1"]), y2, 
                  col = "black", 
                  
                  border = "black",
                  sector.index = cdm_res$cn[i], track.index = 3)
      
    }
    
    
  }
}) # Messages suppressed because of on purpose plotting outside of sectors

dev.off()
# The previous figure is the alternative version of the circos focussed on only subclasses, instead of species, which was not used in the final figures of the article

#### Fig.3b
#. LTPSubclassDistributionMSAndTLCBothScreensInBargraphPlot170120222021AmountToNumberExtendedYAxisAndItalicLegend140320223g.pdf (#)

#. LTPSubclassDistributionMSAndTLCBothScreensInBargraphPlot170120222021AmountToNumber18012022.pdf #


UniqueMSConnectionsLTPsWithLipidSubclasses <- unique(rbind(PureAntonella32b, PureEnric32)[,c("LTPProtein", "LikelySubclass", "Screen")])


# With PG to PG/BMP
UniqueMSConnectionsLTPsWithLipidSubclassesPGConvertedToPGBMP <- UniqueMSConnectionsLTPsWithLipidSubclasses

UniqueMSConnectionsLTPsWithLipidSubclassesPGConvertedToPGBMP[UniqueMSConnectionsLTPsWithLipidSubclassesPGConvertedToPGBMP$LikelySubclass == "PG","LikelySubclass"] <- "PG/BMP"
UniqueMSConnectionsLTPsWithLipidSubclassesPGConvertedToPGBMP <- unique(UniqueMSConnectionsLTPsWithLipidSubclassesPGConvertedToPGBMP)


UniqueMSConnectionsLTPsWithLipidSubclassesPGConvertedToPGBMPUniqueHPTLCsAdded <- unique(do.call("rbind", list(UniqueMSConnectionsLTPsWithLipidSubclassesPGConvertedToPGBMP,
                                                                                                              
                                                                                                              setNames(cbind(HPTLCSpecificitiesPerScreen2hdrm4[[1]][,2:1], "in vivo"), colnames(UniqueMSConnectionsLTPsWithLipidSubclassesPGConvertedToPGBMP)),
                                                                                                              setNames(cbind(HPTLCSpecificitiesPerScreen2hdrm4[[2]][,2:1], "in vitro"), colnames(UniqueMSConnectionsLTPsWithLipidSubclassesPGConvertedToPGBMP)))))

LTPSubclassDistributionMSAndTLCInCellulo <- table(table(UniqueMSConnectionsLTPsWithLipidSubclassesPGConvertedToPGBMPUniqueHPTLCsAdded[UniqueMSConnectionsLTPsWithLipidSubclassesPGConvertedToPGBMPUniqueHPTLCsAdded[,3] == "in vivo",1]))
LTPSubclassDistributionMSAndTLCInVitro <- table(table(UniqueMSConnectionsLTPsWithLipidSubclassesPGConvertedToPGBMPUniqueHPTLCsAdded[UniqueMSConnectionsLTPsWithLipidSubclassesPGConvertedToPGBMPUniqueHPTLCsAdded[,3] == "in vitro",1]))

LTPSubclassDistributionMSAndTLCBothScreens <- merge(as.matrix(LTPSubclassDistributionMSAndTLCInCellulo), as.matrix(LTPSubclassDistributionMSAndTLCInVitro), by = 0, all = TRUE)
colnames(LTPSubclassDistributionMSAndTLCBothScreens) <- c("subclasses bound", "in cellulo", "in vitro") # All subclass amounts present: ok

LTPSubclassDistributionMSAndTLCBothScreens2 <- t(Col1ToRowNames(LTPSubclassDistributionMSAndTLCBothScreens))
LTPSubclassDistributionMSAndTLCBothScreens2[is.na(LTPSubclassDistributionMSAndTLCBothScreens2)] <- 0

par(mar = c(5.1, 4.1, 4.1, 2.1))
pdf("./Output/LTPSubclassDistributionMSAndTLCBothScreensInBargraphPlot170120222021AmountToNumber18012022.pdf")

barplot(LTPSubclassDistributionMSAndTLCBothScreens2, beside = TRUE, col = ColorMatrixTryOut["odd",], main = "For MS with HPTLC data", las = 1, xlab = "Number of lipid subclasses bound", ylab = "Number of LTPs") 
abline(h = 1:max(LTPSubclassDistributionMSAndTLCBothScreens2), col = "#FFFFFF33", lwd = 3.2)

legend("topright", legend = c("in cellulo", "in vitro"), col = ColorMatrixTryOut["odd",], pch = 15, bty = "n", cex = 1)
dev.off()

# Previous figure was used for panel 3b: only y-axis further extended


######## Figure 4
# Comparison of lipid co-mobilization with co-regulation and co-localization in tissues and subcellularly

# Note: first we describe the figures used for making the panels of figure 4, 
# afterwards we also describe an alternative version that is not part of the final figure 4 with Nrd0 in bandwidth selection for the Gaussian kernels for the densities.

# The final sets of figures that were used for the panels of figure 4 are: Panel5AWithTailoredScriptAndMoreDetailedVisualizationNowIncludesAllTooWithBinSize00117022022NoSegmentStacking18032022.pdf;
# Panel5BWithTailoredScriptAndMoreDetailedVisualizationNowIncludesAllTooWithBinSize001170220222NoSegmentStacking18032022.pdf; Panel5CWithTailoredScriptAndMoreDetailedVisualizationBinSize001NoSegmentStacking18032022.pdf


#### Panel 4A

#. Panel5AWithTheNrd0WdAdaptedWithScalesImplemented.pdf (#) (Non-used version of final file)
#. Panel5AWithoutTheNrd0WdAdaptedWithScalesImplemented.pdf (#) (Intermediate cleaned-up version of final file in Illustrator)

#. Panel5AWithTailoredScriptAndMoreDetailedVisualizationNowIncludesAllTooWithBinSize00117022022NoSegmentStacking18032022.pdf # (Basis for work to final version of article: see after part for panel 5C.)
#. (OverlayOfTheViolinPlotsWithBeanlineMediansAttemptToScaleInInverseOrderNrd0CorrectionDone11092020RemakeWithAllIncludedAndOverplottingWithReferenceUnity21032022.pdf # Not used version with Nrd0 and unity line.)


KoeberlinCorrelations <- read.csv(file = "./InputData/Koeberlin_Snijder_Cell_2015_lipid_lipid_correlations.txt", header = TRUE, sep = "\t", as.is = TRUE)

KoeberlinCorrelations2 <- KoeberlinCorrelations[KoeberlinCorrelations[,1] != "SM C20:0", colnames(KoeberlinCorrelations) != "SM.C20.0"]
rm(KoeberlinCorrelations)

SplitRownames <- strsplit(KoeberlinCorrelations2[,1], " ")
SplitRownamesmatrix <- do.call("rbind", lapply(SplitRownames, function(x){switch(length(x), "1" = c(x, "sl", "sl"), "2" = c(x[1], "sl", x[2]), "3" = x)}))

SplitRownamesmatrix21179 <- cbind(SplitRownamesmatrix[1:179,1:2], do.call("rbind", strsplit(do.call("rbind", strsplit(SplitRownamesmatrix[1:179,3], "C"))[,2], ":")))
SplitRownamesmatrix180244 <- do.call("rbind", strsplit(SplitRownamesmatrix[180:244,1], "-"))[,2:3]

frontsls <- do.call("rbind", lapply(strsplit(unlist(strsplit(sapply(strsplit(SplitRownamesmatrix180244[,1],"C"), "[[", 2),")")),"\\("), function(x){switch(length(x), "1" = c(x,""), "2" = x)}))
frontsls2 <- cbind(do.call("rbind",strsplit(frontsls[,1],":")), frontsls[,2])

slentriessplit <- cbind(frontsls2, do.call("rbind", lapply(strsplit(unlist(strsplit(SplitRownamesmatrix180244[,2],")")), "\\("), function(x){switch(length(x), "1" = c(x,""), "2" = x)})))


TypesOfLipidsKoeberlin <- rbind(do.call("cbind", list(SplitRownamesmatrix21179[,c(1,3,4,2)],rep("",179),rep("",179))),
                                do.call("cbind", list(slentriessplit[,c(4,1,2)], rep("sl",65), slentriessplit[,c(3,5)])))

TypesOfLipidsKoeberlin2 <- cbind(TypesOfLipidsKoeberlin, gsub("L", "Lyso", TypesOfLipidsKoeberlin[,1], fixed = T))


TypesOfLipidsKoeberlin4 <- cbind.data.frame(TypesOfLipidsKoeberlin2, V8 = ifelse(grepl("e",TypesOfLipidsKoeberlin2[,4], fixed = T), paste0(TypesOfLipidsKoeberlin2[,7],"-O"),TypesOfLipidsKoeberlin2[,7]))
TypesOfLipidsKoeberlin4x <- TypesOfLipidsKoeberlin4

TypesOfLipidsKoeberlin4x[,8] <- gsub("SM", "d*SM", gsub("Lyso", "L", TypesOfLipidsKoeberlin4x[,8], fixed = T), fixed = T)
TypesOfLipidsKoeberlin4x[,8] <- ifelse(TypesOfLipidsKoeberlin4x[,8] == "Cer", 
                                       
                                       ifelse(TypesOfLipidsKoeberlin4x[,5] != "", 
                                              ifelse(TypesOfLipidsKoeberlin4x[,6] != "", "DHOH*Cer", "dOHCer"), 
                                              
                                              ifelse(TypesOfLipidsKoeberlin4x[,6] != "", "DHCer", "dCer")), 
                                       TypesOfLipidsKoeberlin4x[,8])

TypesOfLipidsKoeberlin4x[,2] <- ifelse(TypesOfLipidsKoeberlin4x[,8] == "d*SM", as.numeric(as.character(TypesOfLipidsKoeberlin4x[,2]))+18, as.numeric(as.character(TypesOfLipidsKoeberlin4x[,2])))
TypesOfLipidsKoeberlin4x[,3] <- ifelse(TypesOfLipidsKoeberlin4x[,8] == "d*SM", as.numeric(as.character(TypesOfLipidsKoeberlin4x[,3]))+1, as.numeric(as.character(TypesOfLipidsKoeberlin4x[,3])))

TypesOfLipidsKoeberlin4x[,2] <- ifelse((TypesOfLipidsKoeberlin4x[,8] == "DHOH*Cer")|(TypesOfLipidsKoeberlin4x[,8] == "DHCer")|(TypesOfLipidsKoeberlin4x[,8] == "dOHCer")|(TypesOfLipidsKoeberlin4x[,8] == "dCer"), as.numeric(as.character(TypesOfLipidsKoeberlin4x[,2]))+18, as.numeric(as.character(TypesOfLipidsKoeberlin4x[,2])))
TypesOfLipidsKoeberlin4x[,3] <- ifelse((TypesOfLipidsKoeberlin4x[,8] == "dOHCer")|(TypesOfLipidsKoeberlin4x[,8] == "dCer"), as.numeric(as.character(TypesOfLipidsKoeberlin4x[,3]))+1, as.numeric(as.character(TypesOfLipidsKoeberlin4x[,3])))


KoeberlinCorrelationsConsensusNames <- KoeberlinCorrelations2[,-1]

rownames(KoeberlinCorrelationsConsensusNames) <- paste0(TypesOfLipidsKoeberlin4x[,8], "(", as.character(TypesOfLipidsKoeberlin4x[,2]), ":", as.character(TypesOfLipidsKoeberlin4x[,3]), ")")
colnames(KoeberlinCorrelationsConsensusNames) <- rownames(KoeberlinCorrelationsConsensusNames)

KoeberlinCorrelationsConsensusNames2 <- KoeberlinCorrelationsConsensusNames
KoeberlinCorrelationsConsensusNames2[upper.tri(KoeberlinCorrelationsConsensusNames2)] <- NA

KoeberlinCorrelationsConsensusNames2Long <- meltmatrix(as.matrix(Col1ToRowNames(cbind(rownames(KoeberlinCorrelationsConsensusNames2), KoeberlinCorrelationsConsensusNames2))))
KoeberlinCorrelationsConsensusNames2LongReducedVersion <- KoeberlinCorrelationsConsensusNames2Long[!is.na(KoeberlinCorrelationsConsensusNames2Long$value),]

colnames(KoeberlinCorrelationsConsensusNames2LongReducedVersion) <- c("lipidA", "lipidB", "correlation")


KoeberlinCorrelationsConsensusNames2LongReducedVersion2 <- data.frame(KoeberlinCorrelationsConsensusNames2LongReducedVersion,
                                                                      "HeadgroupA" = as.character(sapply(strsplit(as.character(KoeberlinCorrelationsConsensusNames2LongReducedVersion$lipidA), "\\("),"[[", 1)),
                                                                      
                                                                      "HeadgroupB" = as.character(sapply(strsplit(as.character(KoeberlinCorrelationsConsensusNames2LongReducedVersion$lipidB), "\\("),"[[", 1)),
                                                                      stringsAsFactors = FALSE)

KoeberlinCorrelationsConsensusNames2LongReducedVersion2 <- data.frame(KoeberlinCorrelationsConsensusNames2LongReducedVersion2,
                                                                      "MatchingHeadgroup" = (KoeberlinCorrelationsConsensusNames2LongReducedVersion2$HeadgroupA == KoeberlinCorrelationsConsensusNames2LongReducedVersion2$HeadgroupB), stringsAsFactors = FALSE)

ConversionMatrixForTotalKoeberlin_2 <- cbind(unique(KoeberlinCorrelationsConsensusNames2LongReducedVersion2$HeadgroupA), c("PC", "PC", "PC", "PC", "PE", "PE", "PE", "PE", "PG", "PG", "PG", "PS", "PS", "SM", "Cer", "Cer", "Cer", "Cer"))
KoeberlinCorrelationsConsensusNames2LongReducedVersion4 <- data.frame(KoeberlinCorrelationsConsensusNames2LongReducedVersion2, 
                                                                      
                                                                      "HeadgroupA2" = ConversionMatrixForTotalKoeberlin_2[match(KoeberlinCorrelationsConsensusNames2LongReducedVersion2$HeadgroupA, ConversionMatrixForTotalKoeberlin_2[,1]),2],
                                                                      "HeadgroupB2" = ConversionMatrixForTotalKoeberlin_2[match(KoeberlinCorrelationsConsensusNames2LongReducedVersion2$HeadgroupB, ConversionMatrixForTotalKoeberlin_2[,1]),2], stringsAsFactors = FALSE)

KoeberlinCorrelationsConsensusNames2LongReducedVersion5 <- data.frame(KoeberlinCorrelationsConsensusNames2LongReducedVersion4, "MatchingHeadgroup2" = (KoeberlinCorrelationsConsensusNames2LongReducedVersion4$HeadgroupA2 == KoeberlinCorrelationsConsensusNames2LongReducedVersion4$HeadgroupB2), stringsAsFactors = FALSE)
KoeberlinCorrelationsConsensusNames2LongReducedVersion7 <- data.frame(KoeberlinCorrelationsConsensusNames2LongReducedVersion5, "MatchingNumber" = (KoeberlinCorrelationsConsensusNames2LongReducedVersion5$MatchingHeadgroup + KoeberlinCorrelationsConsensusNames2LongReducedVersion5$MatchingHeadgroup2), stringsAsFactors = FALSE)

AggregatedLTPLipidPairsAntonella_2 <- aggregate(PureAntonella32b$Intensity, by=list(PureAntonella32b$LTPProtein, PureAntonella32b$Lipid, PureAntonella32b$LikelySubclass, PureAntonella32b$TotalCarbonChainLength, PureAntonella32b$TotalCarbonChainUnsaturations), FUN=sum)
AggregatedLTPLipidPairsEnric_2 <- aggregate(PureEnric32$Intensity, by=list(PureEnric32$LTPProtein, PureEnric32$Lipid, PureEnric32$LikelySubclass, PureEnric32$TotalCarbonChainLength, PureEnric32$TotalCarbonChainUnsaturations), FUN=sum)

colnames(AggregatedLTPLipidPairsAntonella_2) <- c("LTPProtein", "Lipid", "LikelySubclass", "TotalCarbonChainLength", "TotalCarbonChainUnsaturations", "Intensity")
colnames(AggregatedLTPLipidPairsEnric_2) <- c("LTPProtein", "Lipid", "LikelySubclass", "TotalCarbonChainLength", "TotalCarbonChainUnsaturations", "Intensity")

AggregatedLTPLipidPairsAntonella_4 <- cbind(AggregatedLTPLipidPairsAntonella_2, ConsensusWithKoeberlinSubclass = ifelse(as.character(AggregatedLTPLipidPairsAntonella_2$LikelySubclass) == "d*Cer", "dCer",
                                                                                                                        ifelse(as.character(AggregatedLTPLipidPairsAntonella_2$LikelySubclass) == "tCer", "DHOH*Cer", as.character(AggregatedLTPLipidPairsAntonella_2$LikelySubclass))))

AggregatedLTPLipidPairsEnric_4 <- cbind(AggregatedLTPLipidPairsEnric_2, ConsensusWithKoeberlinSubclass = ifelse(as.character(AggregatedLTPLipidPairsEnric_2$LikelySubclass) == "d*Cer", "dCer",
                                                                                                                ifelse(as.character(AggregatedLTPLipidPairsEnric_2$LikelySubclass) == "tCer", "DHOH*Cer", as.character(AggregatedLTPLipidPairsEnric_2$LikelySubclass))))

AggregatedLTPLipidPairsAntonella_5 <- cbind(AggregatedLTPLipidPairsAntonella_4, ConsensusWithKoeberlinSpecies = paste0(AggregatedLTPLipidPairsAntonella_4[,7], "(", as.character(AggregatedLTPLipidPairsAntonella_4[,4]), ":", as.character(AggregatedLTPLipidPairsAntonella_4[,5]), ")"))
AggregatedLTPLipidPairsEnric_5 <- cbind(AggregatedLTPLipidPairsEnric_4, ConsensusWithKoeberlinSpecies = paste0(AggregatedLTPLipidPairsEnric_4[,7], "(", as.character(AggregatedLTPLipidPairsEnric_4[,4]), ":", as.character(AggregatedLTPLipidPairsEnric_4[,5]), ")"))

AggregatedLTPLipidPairsAntonella_7 <- aggregate(AggregatedLTPLipidPairsAntonella_5$Intensity, by=list(AggregatedLTPLipidPairsAntonella_5$LTPProtein, AggregatedLTPLipidPairsAntonella_5$ConsensusWithKoeberlinSpecies, AggregatedLTPLipidPairsAntonella_5$ConsensusWithKoeberlinSubclass), FUN=sum)
AggregatedLTPLipidPairsEnric_7 <- aggregate(AggregatedLTPLipidPairsEnric_5$Intensity, by=list(AggregatedLTPLipidPairsEnric_5$LTPProtein, AggregatedLTPLipidPairsEnric_5$ConsensusWithKoeberlinSpecies, AggregatedLTPLipidPairsEnric_5$ConsensusWithKoeberlinSubclass), FUN=sum)

y <- AggregatedLTPLipidPairsAntonella_7 
LipidConnections_7 <- do.call("rbind", lapply(unique(y[,1]), function(x){format(if(sum(y[,1] == x) > 1){t(combn(y[y[,1] == x, 2],2))}, trim = TRUE)}))

LipidConnections2_7 <- LipidConnections_7[LipidConnections_7[,1] != "NULL",] 
LipidConnections2A_7 <- trimws(LipidConnections2_7)

y <- AggregatedLTPLipidPairsEnric_7
LipidConnections_7 <- do.call("rbind", lapply(unique(y[,1]), function(x){format(if(sum(y[,1] == x) > 1){t(combn(y[y[,1] == x, 2],2))}, trim = TRUE)}))

LipidConnections2_7 <- LipidConnections_7[LipidConnections_7[,1] != "NULL",]
LipidConnections2E_7 <- trimws(LipidConnections2_7)

LipidConnections4A_7 <- unique(LipidConnections2A_7)
LipidConnections4E_7 <- unique(LipidConnections2E_7)


CorrelationsLista_7 <- lapply(1:dim(LipidConnections4A_7)[1], function(x){KoeberlinCorrelationsConsensusNames[LipidConnections4A_7[x,1],LipidConnections4A_7[x,2]]})

CorrelationsMatrixa_7 <- cbind(LipidConnections4A_7[sapply(CorrelationsLista_7, function(x){!is.null(x)}),], do.call("rbind", CorrelationsLista_7))
CorrelationsMatrixa2_7 <- CorrelationsMatrixa_7[!is.na(CorrelationsMatrixa_7[,3]),]

CorrelationsMatrixa4_7 <- CorrelationsMatrixa2_7[order(as.numeric(CorrelationsMatrixa2_7[,3]), decreasing = TRUE),]


CorrelationsMatrixa5_7 <- do.call("cbind", list(CorrelationsMatrixa4_7,
                                                sapply(strsplit(as.character(CorrelationsMatrixa4_7[,1]), "\\("), "[[", 1),
                                                
                                                sapply(strsplit(as.character(CorrelationsMatrixa4_7[,2]), "\\("), "[[", 1)))
conversionmatrixa_7 <- cbind(unique(as.vector(CorrelationsMatrixa5_7[,4:5])), c("PC", "PE", "PC", "SM", "Cer", "PS", "Cer", "Cer", "PE"))


CorrelationsListe_7 <- lapply(1:dim(LipidConnections4E_7)[1], function(x){KoeberlinCorrelationsConsensusNames[LipidConnections4E_7[x,1],LipidConnections4E_7[x,2]]})

CorrelationsMatrixe_7 <- cbind(LipidConnections4E_7[sapply(CorrelationsListe_7, function(x){!is.null(x)}),], do.call("rbind", CorrelationsListe_7))
CorrelationsMatrixe2_7 <- CorrelationsMatrixe_7[!is.na(CorrelationsMatrixe_7[,3]),]

CorrelationsMatrixe4_7 <- CorrelationsMatrixe2_7[order(as.numeric(CorrelationsMatrixe2_7[,3]), decreasing = TRUE),]


CorrelationsMatrixe5_7 <- do.call("cbind", list(CorrelationsMatrixe4_7,
                                                sapply(strsplit(as.character(CorrelationsMatrixe4_7[,1]), "\\("), "[[", 1),
                                                
                                                sapply(strsplit(as.character(CorrelationsMatrixe4_7[,2]), "\\("), "[[", 1)))
conversionmatrixe_7 <- cbind(unique(as.vector(CorrelationsMatrixe5_7[,4:5])), c("PC", "PE", "PC", "PC", "SM", "PG", "PE"))

CorrelationsMatrixa7_7 <- do.call("cbind", list(CorrelationsMatrixa5_7, conversionmatrixa_7[match(CorrelationsMatrixa5_7[,4], conversionmatrixa_7[,1]),2],
                                                conversionmatrixa_7[match(CorrelationsMatrixa5_7[,5], conversionmatrixa_7[,1]),2]))

CorrelationsMatrixe7_7 <- do.call("cbind", list(CorrelationsMatrixe5_7, conversionmatrixe_7[match(CorrelationsMatrixe5_7[,4], conversionmatrixe_7[,1]),2],
                                                conversionmatrixe_7[match(CorrelationsMatrixe5_7[,5], conversionmatrixe_7[,1]),2]))

CorrelationsMatrixa8_7 <- cbind(CorrelationsMatrixa7_7, rowSums(cbind(as.numeric(CorrelationsMatrixa7_7[,4] == CorrelationsMatrixa7_7[,5]), as.numeric(CorrelationsMatrixa7_7[,6] == CorrelationsMatrixa7_7[,7]))))
CorrelationsMatrixe8_7 <- cbind(CorrelationsMatrixe7_7, rowSums(cbind(as.numeric(CorrelationsMatrixe7_7[,4] == CorrelationsMatrixe7_7[,5]), as.numeric(CorrelationsMatrixe7_7[,6] == CorrelationsMatrixe7_7[,7]))))

CorrelationsMatrixae8_8 <- rbind(cbind(CorrelationsMatrixa8_7, "A"), cbind(CorrelationsMatrixe8_7, "E"))
colnames(CorrelationsMatrixae8_8) <- c("From", "To", "Correlation", "From2", "To2", "From4", "To4", "MatchLevel", "Screen")

library(RColorBrewer)
pdf("./Output/Panel5AWithTailoredScriptAndMoreDetailedVisualizationNowIncludesAllTooWithBinSize00117022022NoSegmentStacking18032022.pdf")

BinSize <- 0.01
MaxList2 <- list()

StatList <- list()
StatList2 <- list()

ScreenColorsForFilling <- cbind(c("A","E"),c(brewer.pal(9,"Blues")[4], brewer.pal(9,"Oranges")[4]))
ScreenColorsForFillingPot <- cbind(c("A","E"),c(brewer.pal(9,"Blues")[7], brewer.pal(9,"Oranges")[7]))


SimilarityIndex <- "All"

x1a <- as.numeric(CorrelationsMatrixae8_8[(CorrelationsMatrixae8_8[,"Screen"] == "A"),"Correlation"])
x1e <- as.numeric(CorrelationsMatrixae8_8[(CorrelationsMatrixae8_8[,"Screen"] == "E"),"Correlation"])

x2 <- KoeberlinCorrelationsConsensusNames2LongReducedVersion7[, "correlation"]
MaxY <- max(c(unlist(density(x1a, from = -1, to = 1)["y"]), unlist(density(x1e, from = -1, to = 1)["y"]), unlist(density(x2, from = -1, to = 1)["y"])))

# Set stage
plot(unlist(density(x1a, from = -1, to = 1)["x"]), unlist(density(x1a, from = -1, to = 1)["y"]), xlab = "Co-regulation", ylab = "Normalized Density", ylim = c(-MaxY,MaxY), col = "white", main = paste0("Similarity: ", SimilarityIndex), type = "l", bty="n")

polygon(x = c(min(unlist(density(x1a, from = -1, to = 1)["x"])), unlist(density(x1a, from = -1, to = 1)["x"]), max(unlist(density(x1a, from = -1, to = 1)["x"]))),
        y = c(0, unlist(density(x1a, from = -1, to = 1)["y"]), 0), col = ScreenColorsForFilling[ScreenColorsForFilling[,1] == "A",2], border = "white")

polygon(x = c(min(unlist(density(x1e, from = -1, to = 1)["x"])), unlist(density(x1e, from = -1, to = 1)["x"]), max(unlist(density(x1e, from = -1, to = 1)["x"]))),
        y = c(0, -unlist(density(x1e, from = -1, to = 1)["y"]), 0), col = ScreenColorsForFilling[ScreenColorsForFilling[,1] == "E",2], border = "white")

segments(x1a, 0.16, x1 = x1a, y1 = 0, col = "black")
segments(x1e, -0.16, x1 = x1e, y1 = 0, col = "black")

lines(unlist(density(x2, from = -1, to = 1)["x"]), unlist(density(x2, from = -1, to = 1)["y"]), col = "DarkGrey", lwd = 2)
lines(unlist(density(x2, from = -1, to = 1)["x"]), -unlist(density(x2, from = -1, to = 1)["y"]), col = "DarkGrey", lwd = 2)

segments(x0 = median(x1a), x1 = median(x1a), y0 = MaxY, y1 = 0, lwd = 2)
segments(x0 = median(x2), x1 = median(x2), y0 = MaxY, y1 = 0, lwd = 2, col = "DarkGrey")

segments(x0 = median(x1e), x1 = median(x1e), y0 = -MaxY, y1 = 0, lwd = 2)
segments(x0 = median(x2), x1 = median(x2), y0 = -MaxY, y1 = 0, lwd = 2, col = "DarkGrey")

MaxList2[[SimilarityIndex]][["A"]] <- c(max(unlist(density(x1a, from = -1, to = 1)["y"])),max(unlist(density(x2, from = -1, to = 1)["y"])))
MaxList2[[SimilarityIndex]][["E"]] <- c(max(unlist(density(x1e, from = -1, to = 1)["y"])),max(unlist(density(x2, from = -1, to = 1)["y"])))

for(SimilarityIndex in as.character(0:2)){
  
  
  x1a <- as.numeric(CorrelationsMatrixae8_8[(CorrelationsMatrixae8_8[,"Screen"] == "A") & (CorrelationsMatrixae8_8[,"MatchLevel"] == SimilarityIndex),"Correlation"])
  x1e <- as.numeric(CorrelationsMatrixae8_8[(CorrelationsMatrixae8_8[,"Screen"] == "E") & (CorrelationsMatrixae8_8[,"MatchLevel"] == SimilarityIndex),"Correlation"])
  
  x2 <- KoeberlinCorrelationsConsensusNames2LongReducedVersion7[KoeberlinCorrelationsConsensusNames2LongReducedVersion7[,"MatchingNumber"] == SimilarityIndex, "correlation"]
  MaxY <- max(c(unlist(density(x1a, from = -1, to = 1)["y"]), unlist(density(x1e, from = -1, to = 1)["y"]), unlist(density(x2, from = -1, to = 1)["y"])))
  
  plot(unlist(density(x1a, from = -1, to = 1)["x"]), unlist(density(x1a, from = -1, to = 1)["y"]), xlab = "Co-regulation", ylab = "Normalized Density", ylim = c(-MaxY,MaxY), col = "white", main = paste0("Similarity: ", SimilarityIndex), type = "l", bty="n")
  
  
  polygon(x = c(min(unlist(density(x1a, from = -1, to = 1)["x"])), unlist(density(x1a, from = -1, to = 1)["x"]), max(unlist(density(x1a, from = -1, to = 1)["x"]))),
          y = c(0, unlist(density(x1a, from = -1, to = 1)["y"]), 0), col = ScreenColorsForFilling[ScreenColorsForFilling[,1] == "A",2], border = "white")
  
  polygon(x = c(min(unlist(density(x1e, from = -1, to = 1)["x"])), unlist(density(x1e, from = -1, to = 1)["x"]), max(unlist(density(x1e, from = -1, to = 1)["x"]))),
          y = c(0, -unlist(density(x1e, from = -1, to = 1)["y"]), 0), col = ScreenColorsForFilling[ScreenColorsForFilling[,1] == "E",2], border = "white")
  
  segments(x1a, 0.16, x1 = x1a, y1 = 0, col = "black")
  segments(x1e, -0.16, x1 = x1e, y1 = 0, col = "black")
  
  lines(unlist(density(x2, from = -1, to = 1)["x"]), unlist(density(x2, from = -1, to = 1)["y"]), col = "DarkGrey", lwd = 2)
  lines(unlist(density(x2, from = -1, to = 1)["x"]), -unlist(density(x2, from = -1, to = 1)["y"]), col = "DarkGrey", lwd = 2)
  
  segments(x0 = median(x1a), x1 = median(x1a), y0 = MaxY, y1 = 0, lwd = 2)
  segments(x0 = median(x2), x1 = median(x2), y0 = MaxY, y1 = 0, lwd = 2, col = "DarkGrey")
  
  segments(x0 = median(x1e), x1 = median(x1e), y0 = -MaxY, y1 = 0, lwd = 2)
  segments(x0 = median(x2), x1 = median(x2), y0 = -MaxY, y1 = 0, lwd = 2, col = "DarkGrey")
  
  MaxList2[[SimilarityIndex]][["A"]] <- c(max(unlist(density(x1a, from = -1, to = 1)["y"])),max(unlist(density(x2, from = -1, to = 1)["y"])))
  MaxList2[[SimilarityIndex]][["E"]] <- c(max(unlist(density(x1e, from = -1, to = 1)["y"])),max(unlist(density(x2, from = -1, to = 1)["y"])))
  
}
# Generalize the maximum for y: all figures on the same scale

SimilarityIndex <- "All"
MaxY <- max(unlist(MaxList2)) 

x1a <- as.numeric(CorrelationsMatrixae8_8[(CorrelationsMatrixae8_8[,"Screen"] == "A"),"Correlation"])
x1e <- as.numeric(CorrelationsMatrixae8_8[(CorrelationsMatrixae8_8[,"Screen"] == "E"),"Correlation"])

x2 <- KoeberlinCorrelationsConsensusNames2LongReducedVersion7[, "correlation"]
plot(unlist(density(x1a, from = -1, to = 1)["x"]), unlist(density(x1a, from = -1, to = 1)["y"]), xlab = "Co-regulation", ylab = "Normalized Density", ylim = c(-MaxY,MaxY), col = "white", main = paste0("Similarity: ", SimilarityIndex), type = "l", bty="n")

polygon(x = c(min(unlist(density(x1a, from = -1, to = 1)["x"])), unlist(density(x1a, from = -1, to = 1)["x"]), max(unlist(density(x1a, from = -1, to = 1)["x"]))),
        y = c(0, unlist(density(x1a, from = -1, to = 1)["y"]), 0), col = ScreenColorsForFilling[ScreenColorsForFilling[,1] == "A",2], border = "white")

polygon(x = c(min(unlist(density(x1e, from = -1, to = 1)["x"])), unlist(density(x1e, from = -1, to = 1)["x"]), max(unlist(density(x1e, from = -1, to = 1)["x"]))),
        y = c(0, -unlist(density(x1e, from = -1, to = 1)["y"]), 0), col = ScreenColorsForFilling[ScreenColorsForFilling[,1] == "E",2], border = "white")

segments(x1a, 0.16, x1 = x1a, y1 = 0, col = "black")
segments(x1e, -0.16, x1 = x1e, y1 = 0, col = "black")

lines(unlist(density(x2, from = -1, to = 1)["x"]), unlist(density(x2, from = -1, to = 1)["y"]), col = "DarkGrey", lwd = 2)
lines(unlist(density(x2, from = -1, to = 1)["x"]), -unlist(density(x2, from = -1, to = 1)["y"]), col = "DarkGrey", lwd = 2)

segments(x0 = median(x1a), x1 = median(x1a), y0 = MaxY, y1 = 0, lwd = 2)
segments(x0 = median(x2), x1 = median(x2), y0 = MaxY, y1 = 0, lwd = 2, col = "DarkGrey")

segments(x0 = median(x1e), x1 = median(x1e), y0 = -MaxY, y1 = 0, lwd = 2)
segments(x0 = median(x2), x1 = median(x2), y0 = -MaxY, y1 = 0, lwd = 2, col = "DarkGrey")


MaxList2 <- list()

MaxList2[[SimilarityIndex]][["A"]] <- c(max(unlist(density(x1a, from = -1, to = 1)["y"])),max(unlist(density(x2, from = -1, to = 1)["y"])))
MaxList2[[SimilarityIndex]][["E"]] <- c(max(unlist(density(x1e, from = -1, to = 1)["y"])),max(unlist(density(x2, from = -1, to = 1)["y"])))

StatList[[SimilarityIndex]] <- structure(rbind(FishersExactTestNegPos2(SubsetX = x1a, SupersetY = x2), FishersExactTestNegPos2(SubsetX = x1e, SupersetY = x2)), 
                                         dimnames = list(c("A","E"), c("NoNeg", "NoPos", "YesNeg", "YesPos", "LowConf", "HighConf", "OddsRatio", "pValue")))

StatList2[[SimilarityIndex]] <- structure(rbind(FishersExactTestNegPosEqualBackground2(SubsetX = x1a, SupersetY = x2), FishersExactTestNegPosEqualBackground2(SubsetX = x1e, SupersetY = x2)), 
                                          dimnames = list(c("A","E"), c("NoNeg", "NoPos", "YesNeg", "YesPos", "LowConf", "HighConf", "OddsRatio", "pValue")))


for(SimilarityIndex in as.character(0:2)){
  
  x1a <- as.numeric(CorrelationsMatrixae8_8[(CorrelationsMatrixae8_8[,"Screen"] == "A") & (CorrelationsMatrixae8_8[,"MatchLevel"] == SimilarityIndex),"Correlation"])
  x1e <- as.numeric(CorrelationsMatrixae8_8[(CorrelationsMatrixae8_8[,"Screen"] == "E") & (CorrelationsMatrixae8_8[,"MatchLevel"] == SimilarityIndex),"Correlation"])
  
  x2 <- KoeberlinCorrelationsConsensusNames2LongReducedVersion7[KoeberlinCorrelationsConsensusNames2LongReducedVersion7[,"MatchingNumber"] == SimilarityIndex, "correlation"]
  plot(unlist(density(x1a, from = -1, to = 1)["x"]), unlist(density(x1a, from = -1, to = 1)["y"]), xlab = "Co-regulation", ylab = "Normalized Density", ylim = c(-MaxY,MaxY), col = "white", main = paste0("Similarity: ", SimilarityIndex), type = "l", bty="n")
  
  polygon(x = c(min(unlist(density(x1a, from = -1, to = 1)["x"])), unlist(density(x1a, from = -1, to = 1)["x"]), max(unlist(density(x1a, from = -1, to = 1)["x"]))),
          y = c(0, unlist(density(x1a, from = -1, to = 1)["y"]), 0), col = ScreenColorsForFilling[ScreenColorsForFilling[,1] == "A",2], border = "white")
  
  polygon(x = c(min(unlist(density(x1e, from = -1, to = 1)["x"])), unlist(density(x1e, from = -1, to = 1)["x"]), max(unlist(density(x1e, from = -1, to = 1)["x"]))),
          y = c(0, -unlist(density(x1e, from = -1, to = 1)["y"]), 0), col = ScreenColorsForFilling[ScreenColorsForFilling[,1] == "E",2], border = "white")
  
  segments(x1a, 0.16, x1 = x1a, y1 = 0, col = "black")
  segments(x1e, -0.16, x1 = x1e, y1 = 0, col = "black")
  
  lines(unlist(density(x2, from = -1, to = 1)["x"]), unlist(density(x2, from = -1, to = 1)["y"]), col = "DarkGrey", lwd = 2)
  lines(unlist(density(x2, from = -1, to = 1)["x"]), -unlist(density(x2, from = -1, to = 1)["y"]), col = "DarkGrey", lwd = 2)
  
  segments(x0 = median(x1a), x1 = median(x1a), y0 = MaxY, y1 = 0, lwd = 2)
  segments(x0 = median(x2), x1 = median(x2), y0 = MaxY, y1 = 0, lwd = 2, col = "DarkGrey")
  
  segments(x0 = median(x1e), x1 = median(x1e), y0 = -MaxY, y1 = 0, lwd = 2)
  segments(x0 = median(x2), x1 = median(x2), y0 = -MaxY, y1 = 0, lwd = 2, col = "DarkGrey")
  
  MaxList2[[SimilarityIndex]][["A"]] <- c(max(unlist(density(x1a, from = -1, to = 1)["y"])),max(unlist(density(x2, from = -1, to = 1)["y"])))
  MaxList2[[SimilarityIndex]][["E"]] <- c(max(unlist(density(x1e, from = -1, to = 1)["y"])),max(unlist(density(x2, from = -1, to = 1)["y"])))
  
  StatList[[SimilarityIndex]] <- structure(rbind(FishersExactTestNegPos2(SubsetX = x1a, SupersetY = x2), FishersExactTestNegPos2(SubsetX = x1e, SupersetY = x2)), 
                                           dimnames = list(c("A","E"), c("NoNeg", "NoPos", "YesNeg", "YesPos", "LowConf", "HighConf", "OddsRatio", "pValue")))
  
  StatList2[[SimilarityIndex]] <- structure(rbind(FishersExactTestNegPosEqualBackground2(SubsetX = x1a, SupersetY = x2), FishersExactTestNegPosEqualBackground2(SubsetX = x1e, SupersetY = x2)), 
                                            dimnames = list(c("A","E"), c("NoNeg", "NoPos", "YesNeg", "YesPos", "LowConf", "HighConf", "OddsRatio", "pValue")))
  
}
dev.off()

# All of the figures in the second set of figure pages in the previous document were combined into Panel 4a and further esthetically enhanced and integrated with the other parts of Figure 4 in Adobe Illustrator.
# For the statistical significances the results of the Fisher exact tests on the difference in the medians was used (see further), and also added in Adobe Illustrator.


#### Panel 4B

#. Panel5BWithoutTheNrd0WdAdaptedWithScalesImplemented22032022.pdf (#)
#. Panel5BWithTheNrd0WdAdaptedWithScalesImplemented22032022.pdf (#)

#. Panel5BWithTailoredScriptAndMoreDetailedVisualizationNowIncludesAllTooWithBinSize001170220222NoSegmentStacking18032022.pdf # (Basis for used version (without the Nrd0))
#. Panel5BWithGreenUnityDistributionWithNrd0AndWd02322032022b # See after the 5C part. Not finally used version.

# Original input for the Manders' overlap coefficient calculations on METASPACE. Has been processed according to what is described below and restored in a different format for use here (see further).
# Coexp_Mander_NoTreshold_Rework21022019 <- read.csv(file = "./InputData/Coexp_Mander_NoTreshold_Rework21022019.csv", header = TRUE, sep = ";", as.is = TRUE)

# Which entries are lipids and lipid-like molecules according to ClassyFire software (see SuperClass column)
Lipid_Classes_Sergio_Subsets_210220194 <- read.csv(file = "./InputData/Lipid_Classes_Sergio_Subsets_210220194.txt", header = TRUE, sep = "\t", as.is = TRUE)

# table(Lipid_Classes_Sergio_Subsets_210220194$SuperClass)
# 1032 lipid-like molecules; 179 other molecules; 1 (224): wrong annotation: -hydroxyprednisolone"

# -hydroxyprednisolone" Lipids and lipid-like molecules                           Other 
#                               1                            1032                             179 

# Eliminate non-lipids from METASPACE subsetting
SubsetOfLipidClassesSergio <- Lipid_Classes_Sergio_Subsets_210220194[Lipid_Classes_Sergio_Subsets_210220194$SuperClass != "Other",]


# Previous reprocessing of the Coexp_Mander_NoTreshold_Rework210220192 before storage in another format for use in this repository (see below)

# Coexp_Mander_NoTreshold_Rework210220192 <- Coexp_Mander_NoTreshold_Rework21022019[,3:47]
# rownames(Coexp_Mander_NoTreshold_Rework210220192) <- paste(Coexp_Mander_NoTreshold_Rework21022019[,1], Coexp_Mander_NoTreshold_Rework21022019[,2], sep = "_")

# Save in smaller and easier to handle RDS format. GitHub issues with file size avoided in this way too.
# saveRDS(object = Coexp_Mander_NoTreshold_Rework210220192, file = "./InputData/Coexp_Mander_NoTreshold_Rework210220192.rds") 

# Read in new input-file in RDS-format for the Manders' overlap coefficient calculations of METASPACE.
Coexp_Mander_NoTreshold_Rework210220192 <- readRDS("./InputData/Coexp_Mander_NoTreshold_Rework210220192.rds")

MList <- lapply(1:dim(Coexp_Mander_NoTreshold_Rework210220192)[2],function(y){sapply(Coexp_Mander_NoTreshold_Rework210220192[,y],function(x){if(is.na(x)){NA}else{strsplit(GetStuffBetweenBrackets(x)[[1]], ",")[[1]][1]}})})
MMatrix <- do.call("cbind", MList)

mode(MMatrix) <- "numeric"
MMeansx <- cbind.data.frame(do.call("rbind", strsplit(rownames(Coexp_Mander_NoTreshold_Rework210220192), "_")), rowMeans(MMatrix, na.rm = TRUE)) # Updated to second version to keep script coherent & changed extraction accordingly and updated to x form to indicate the change.

MMeans2x <- cbind(MMeansx,cbind(sapply(MMeansx[,1], function(x){x %in% SubsetOfLipidClassesSergio[,1]}),sapply(MMeansx[,2], function(x){x %in% SubsetOfLipidClassesSergio[,1]})))
LipidSubsetMMeansx <- MMeans2x[MMeans2x[,4] & MMeans2x[,5],]

LipidSubsetMMeans_2x <- do.call("cbind", list(LipidSubsetMMeansx,
                                              "From" = gsub(pattern = "([A-Z](?![0-9]))", replacement = "\\11", x = LipidSubsetMMeansx[,1], perl = TRUE),
                                              
                                              "To" = gsub(pattern = "([A-Z](?![0-9]))", replacement = "\\11", x = LipidSubsetMMeansx[,2], perl = TRUE)))


LipidSubsetMMeans_4x <- data.frame(LipidSubsetMMeans_2x[,1:5], 
                                   FromTo = paste(LipidSubsetMMeans_2x[,"From"], LipidSubsetMMeans_2x[,"To"], sep = "_"), 
                                   
                                   LipidSubsetMMeans_2x[,6:7], 
                                   stringsAsFactors = FALSE)

HeadGroupConversionMatrix <- do.call("cbind", list("HeadGroup" = c("PC", "LPC", "PE", "LPE", "PG", "LPG", "PG/BMP", "BMP", "PS", "PI", "FA", "PA", "TAG", "DAG", "VE", "VA"),
                                                   "C" = c(8, 8, 5, 5, 6, 6, 6, 6, 6, 9, 0, 3, 3, 3, 26, 20),
                                                   
                                                   "H" = c(18, 18, 12, 12, 13, 13, 13, 13, 12, 16, 1, 7, 5, 6, 50, 30),
                                                   "N" = c(1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0),
                                                   
                                                   "O" = c(6, 6, 6, 6, 8, 8, 8, 8, 8, 11, 1, 6, 3, 3, 2, 1),
                                                   "P" = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0)))


LipidSubsetMMeans_4Subset <- data.frame(Cooccurrence = LipidSubsetMMeans_4x[,3], LipidPairMatchingType = "All") # Inflow in same previous system again: no x needed anymore.

AggregatedLTPLipidPairsCombined_7wvi <- unique(rbind(cbind(PureAntonella32b[,c("LTPProtein", "Lipid", "LikelySubclass")], Screen = "A"), cbind(PureEnric32[,c("LTPProtein", "Lipid", "LikelySubclass")], Screen = "E")))
colnames(AggregatedLTPLipidPairsCombined_7wvi) <- c("LTPProtein", "LipidSpecies", "LipidSubclass", "Screen")

# Convertion: a: acyl; e: ether; d: d sphingosine; h: DH sphingosine; t: t sphingosine)
# Alternative nomenclature: similar to Koeberlin but not exactly the same: lysolipids indicated by headgroup, so PC_a_C and not LPC_a_C

ConversionMatrixForLipidsForFullDataAlternativeKoeberlinNomenclaturewvi <- cbind(c("BMP", "CL", "d*Cer", "d*CerP", "d*HexCer", "d*SHexCer", "d*SM", "DAG", "dCer", "DHCer", "DHOH*Cer", "DHSM", "FA", "FAL", "LPC", "LPE", "LPE-O", "LPG", 
                                                                                   "PA", "PC", "PC-O", "PE", "PE-O", "PG", "PG/BMP", "PGP", "PI", "PS", "t*Hex2Cer", "t*HexCer", "t*SM", "TAG", "tCer", "VA"), 
                                                                                 
                                                                                 c("BMP_aa_C", "CL_aaaa_C", "Cer_da_C", "CerP_da_C", "HexCer_da_C", "SHexCer_da_C", "SM_da_C", "DAG_aa_C", "Cer_da_C", "Cer_ha_C", "Cer_ta_C", "SM_ha_C",
                                                                                   "FA_a_C", "FA_e_C", "PC_a_C", "PE_a_C", "PE_e_C", "PG_a_C", "PA_aa_C", "PC_aa_C", "PC_ae_C", "PE_aa_C", "PE_ae_C", "PG_aa_C", "PG/BMP_aa_C", "PGP_aa_C", 
                                                                                   
                                                                                   "PI_aa_C",  "PS_aa_C", "Hex2Cer_ta_C", "HexCer_ta_C", "SM_ta_C", "TAG_aaa_C", "Cer_ta_C", "VA_nn_C"))


AggregatedLTPLipidPairsCombined_8wvi <- cbind(AggregatedLTPLipidPairsCombined_7wvi, 
                                              LipidSpeciesAlternativeNomenclature = paste0(ConversionMatrixForLipidsForFullDataAlternativeKoeberlinNomenclaturewvi[match(AggregatedLTPLipidPairsCombined_7wvi$LipidSubclass, ConversionMatrixForLipidsForFullDataAlternativeKoeberlinNomenclaturewvi[,1]),2], 
                                                                                           
                                                                                           unlist(lapply(GetStuffBetweenBrackets(AggregatedLTPLipidPairsCombined_7wvi$LipidSpecies), function(x){if(length(x) != 0){mgsub::mgsub(x[length(x)], c("d", "DH", "t", "-OH", "\\*", "O-"), rep("",6))}else{"nn:nn"}}))))
write.table(AggregatedLTPLipidPairsCombined_8wvi, file="./Output/LTPLipidConnectionsWithSphingolipidsIncludedAndNewestDataKevinTiteca24092020.csv", sep="\t", row.names = FALSE, quote = FALSE)


# Corrected version of PI Headgroup

HeadGroupConversionMatrixc <- HeadGroupConversionMatrix
HeadGroupConversionMatrixc[HeadGroupConversionMatrixc[,"HeadGroup"] == "PI","H"] <- "17"


HeadGroupConversionMatrix2 <- rbind(cbind(HeadGroupConversionMatrixc, "S" = rep(0,16)),
                                    
                                    do.call("cbind", list("HeadGroup" = c("PGP", "CL", "Cer", "CerP", "SM", "HexCer", "Hex2Cer", "SHexCer"),
                                                          "C" = c(6,9,3,3,8,9,15,9), 
                                                          
                                                          "H" = c(14,18,7,8,19,17,27,17),
                                                          "N" = c(0,0,1,1,2,1,1,1),
                                                          
                                                          "O" = c(11,13,2,5,5,7,12,10),
                                                          "P" = c(2,2,0,1,1,0,0,0),
                                                          
                                                          "S" = c(0,0,0,0,0,0,0,1)))) 


ConnectorPieceConversionMatrix2 <-cbind("Connector" = c("a", "e", "d", "h", "t", "n"),
                                        "C" = c(1,1,2,2,2,0), 
                                        
                                        "H" = c(0,2,2,4,4,0),
                                        "N" = c(0,0,0,0,0,0),
                                        
                                        "O" = c(1,0,0,0,1,0),
                                        "P" = c(0,0,0,0,0,0),
                                        
                                        "S" = c(0,0,0,0,0,0))


HTips2 <- do.call("cbind", list("HeadGroup" = HeadGroupConversionMatrix2[,1],
                                "HAmount" = c(2,2,2,2,2,2,2,2,2,2,1,2,3,2,0,0,2,4,2,2,2,2,2,2))) # Similar to just counting the amount of linkages, which could be done as alternative if necessary.

AggregatedLTPLipidPairsCombined_8wvi2 <- cbind(AggregatedLTPLipidPairsCombined_8wvi, do.call("rbind",strsplit(as.character(AggregatedLTPLipidPairsCombined_8wvi$LipidSpeciesAlternativeNomenclature), "\\_C|\\:")))
AggregatedLTPLipidPairsCombined_8wvi2 <- cbind(AggregatedLTPLipidPairsCombined_8wvi2, do.call("rbind", strsplit(as.character(AggregatedLTPLipidPairsCombined_8wvi2[,6]),"\\_")))

colnames(AggregatedLTPLipidPairsCombined_8wvi2)[6:10] <- c("SubclassKoeberlinlike", "TotalCarb", "TotalUnsat", "Headgroup", "Linkages")


library(stringr)
AggregatedLTPLipidPairsCombined_8wvi4 <- cbind(AggregatedLTPLipidPairsCombined_8wvi2, do.call("cbind", structure(lapply(ConnectorPieceConversionMatrix2[,1], function(x){str_count(AggregatedLTPLipidPairsCombined_8wvi2$Linkages, x)}), names = paste0(ConnectorPieceConversionMatrix2[,1],"Counts"))))

Carb <- suppressWarnings(as.numeric(as.character(AggregatedLTPLipidPairsCombined_8wvi4[, "TotalCarb"]))) # NAs introduced on purpose (warnings therefor suppressed)
MissingCarbByLinkage <- colSums(t(AggregatedLTPLipidPairsCombined_8wvi4[,paste0(ConnectorPieceConversionMatrix2[,1],"Counts")])*c(1,1,5,5,5,0))

# C residuals fatty acids
AggregatedLTPLipidPairsCombined_8wvi4c <- cbind(AggregatedLTPLipidPairsCombined_8wvi4, CResidualFattyAcids = (ifelse(!is.na(Carb), Carb, 0) - MissingCarbByLinkage))

# H residuals fatty acids
AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms <- cbind(AggregatedLTPLipidPairsCombined_8wvi4c, HResidualFattyAcids = 2*AggregatedLTPLipidPairsCombined_8wvi4c$CResidualFattyAcids + as.numeric(HTips2[match(as.character(AggregatedLTPLipidPairsCombined_8wvi4c$Headgroup), HTips2[,"HeadGroup"]), "HAmount"]) - 2*suppressWarnings(as.numeric(as.character(AggregatedLTPLipidPairsCombined_8wvi4c$TotalUnsat)))) # NAs introduction wantes: warnings suppressed.

AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms[is.na(AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms$HResidualFattyAcids), "HResidualFattyAcids"] <- 0


# Correction for the da sphingolipids, because otherwise 2H to much subtracted 
AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms$HResidualFattyAcids <- ifelse(AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms$Linkages == "da", AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms[, "HResidualFattyAcids"] + 2, AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms[, "HResidualFattyAcids"])

AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms2 <- cbind(AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms, do.call("rbind", lapply(strsplit(as.character(AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms$Linkages),""), 
                                                                                                                                       function(x){if(length(x) > 1){colSums(apply(ConnectorPieceConversionMatrix2[match(x, ConnectorPieceConversionMatrix2[,1]),-1],2, as.numeric))
                                                                                                                                         
                                                                                                                                       }else{as.numeric(ConnectorPieceConversionMatrix2[match(x, ConnectorPieceConversionMatrix2[,1]),-1])}})) + # Linkages
                                                              do.call("cbind", list(AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms[,c("CResidualFattyAcids","HResidualFattyAcids")], 0, 0, 0, 0)) + # Residuals of fatty acids
                                                              
                                                              apply(HeadGroupConversionMatrix2[match(as.character(AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms$Headgroup),HeadGroupConversionMatrix2),-1], 2,  as.numeric))
colnames(AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms2)[19:24] <- c("C", "H", "N", "O", "P", "S")

AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms4 <- cbind(AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms2, ChemicalFormulaWithOnes = apply(AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms2[,c("C", "H", "N", "O", "P", "S")], MARGIN = 1, function(x){paste0(c("C", "H", "N", "O", "P", "S")[x != 0], x[x != 0], collapse = "")}))


LTPsScreens <- unique(AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms4[,c("LTPProtein","Screen")])
DataSet <- AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms4

ObservedPossibleChemicalCombinationsForLipidsOfLTPsx <- do.call("rbind", lapply(1:dim(LTPsScreens)[1], function(y){
  if(length(unique(as.character(DataSet[(DataSet$LTPProtein == LTPsScreens[y,]$LTPProtein) & (DataSet$Screen == LTPsScreens[y,]$Screen),"ChemicalFormulaWithOnes"]))) > 1){
    
    cbind(LTPsScreens[y,], t(combn(unique(as.character(DataSet[(DataSet$LTPProtein == LTPsScreens[y,]$LTPProtein) & (DataSet$Screen == LTPsScreens[y,]$Screen),"ChemicalFormulaWithOnes"])), 2)), row.names = NULL)}}))
# Simplified rownames: x version

AllPossibleChemicalCombinationsForLipidsOfLTPs <- t(combn(as.character(unique(DataSet$ChemicalFormulaWithOnes)),2)) 
#All possible pairs: 30381 <=> observed pairs: 10256

dim(ObservedPossibleChemicalCombinationsForLipidsOfLTPsx[ObservedPossibleChemicalCombinationsForLipidsOfLTPsx$Screen == "A",])[1] #2777
dim(ObservedPossibleChemicalCombinationsForLipidsOfLTPsx[ObservedPossibleChemicalCombinationsForLipidsOfLTPsx$Screen == "E",])[1] #7479

ObservedPossibleChemicalCombinationsForLipidsOfLTPs2x <- do.call("cbind", list(ObservedPossibleChemicalCombinationsForLipidsOfLTPsx, 
                                                                               ChemicalPairs1 = apply(ObservedPossibleChemicalCombinationsForLipidsOfLTPsx[,3:4], 1, function(x){paste0(x,collapse = "_")}),
                                                                               
                                                                               ChemicalPairs2 = apply(ObservedPossibleChemicalCombinationsForLipidsOfLTPsx[,4:3], 1, function(x){paste0(x,collapse = "_")})))


AllPossibleChemicalCombinationsForLipidsOfLTPs2 <- do.call("cbind", list(AllPossibleChemicalCombinationsForLipidsOfLTPs, 
                                                                         ChemicalPairs1 = apply(AllPossibleChemicalCombinationsForLipidsOfLTPs[,1:2], 1, function(x){paste0(x,collapse = "_")}),
                                                                         
                                                                         ChemicalPairs2 = apply(AllPossibleChemicalCombinationsForLipidsOfLTPs[,2:1], 1, function(x){paste0(x,collapse = "_")})))


AllPossibleChemicalCombinationsForLipidsOfLTPsa <- t(combn(as.character(unique(DataSet[DataSet$Screen == "A",]$ChemicalFormulaWithOnes)),2))
AllPossibleChemicalCombinationsForLipidsOfLTPse <- t(combn(as.character(unique(DataSet[DataSet$Screen == "E",]$ChemicalFormulaWithOnes)),2))

dim(AllPossibleChemicalCombinationsForLipidsOfLTPsa)[1] #8256
dim(AllPossibleChemicalCombinationsForLipidsOfLTPse)[1] #15931

AllPossibleChemicalCombinationsForLipidsOfLTPs2a <- do.call("cbind", list(AllPossibleChemicalCombinationsForLipidsOfLTPsa, 
                                                                          ChemicalPairs1 = apply(AllPossibleChemicalCombinationsForLipidsOfLTPsa[,1:2], 1, function(x){paste0(x,collapse = "_")}),
                                                                          
                                                                          ChemicalPairs2 = apply(AllPossibleChemicalCombinationsForLipidsOfLTPsa[,2:1], 1, function(x){paste0(x,collapse = "_")})))


AllPossibleChemicalCombinationsForLipidsOfLTPs2e <- do.call("cbind", list(AllPossibleChemicalCombinationsForLipidsOfLTPse, 
                                                                          ChemicalPairs1 = apply(AllPossibleChemicalCombinationsForLipidsOfLTPse[,1:2], 1, function(x){paste0(x,collapse = "_")}),
                                                                          
                                                                          ChemicalPairs2 = apply(AllPossibleChemicalCombinationsForLipidsOfLTPse[,2:1], 1, function(x){paste0(x,collapse = "_")})))


CooccurrenceOfChemicalPairsListx <- lapply(list(ObservedPossibleChemicalCombinationsForLipidsOfLTPs2x, 
                                                as.data.frame(AllPossibleChemicalCombinationsForLipidsOfLTPs2), 
                                                
                                                as.data.frame(AllPossibleChemicalCombinationsForLipidsOfLTPs2a), 
                                                as.data.frame(AllPossibleChemicalCombinationsForLipidsOfLTPs2e)), 
                                           
                                           function(y){cbind(y, Cooccurrence = LipidSubsetMMeans_4x[match(as.character(y[,"ChemicalPairs1"]), as.character(LipidSubsetMMeans_4x[,"FromTo"])),3])})


# Write to do in Excel # (Entered earlier to avoid repetition and looping)
write.table(t(t(unique(as.character(AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms4[,"LipidSubclass"])))), file="./Output/SubclassesForObservedLipids02102020.csv", sep="\t", row.names = FALSE, quote = FALSE)

# Input again with the consensus subclasses and classes, BMP & PG/BMP seen as a type of PG # (Entered earlier to avoid repetition and looping)
SubclassesMatchingToConsensusSubclassesAndClasses <- read.csv(file = "./InputData/SubclassesForObservedLipidsConsensusSubclassesAndClasses02102020.txt", header = TRUE, sep = "\t", as.is = TRUE, quote = "")


CooccurrenceOfChemicalPairsList2x <- list()

y <- CooccurrenceOfChemicalPairsListx[[1]]
CooccurrenceOfChemicalPairsList2x[[1]] <- do.call("cbind", list(y,
                                                                
                                                                Subclass1 = AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms4[match(apply(y[,1:3], 1, function(x){paste0(x,collapse = "_")}),
                                                                                                                                     apply(AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms4[,c("LTPProtein", "Screen","ChemicalFormulaWithOnes")], 1, function(x){paste0(x,collapse = "_")})), "LipidSubclass"],
                                                                
                                                                Subclass2 = AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms4[match(apply(y[,c(1,2,4)], 1, function(x){paste0(x,collapse = "_")}),
                                                                                                                                     apply(AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms4[,c("LTPProtein", "Screen","ChemicalFormulaWithOnes")], 1, function(x){paste0(x,collapse = "_")})), "LipidSubclass"]
                                                                
))



y <- CooccurrenceOfChemicalPairsListx[[2]]

CooccurrenceOfChemicalPairsList2x[[2]] <- do.call("cbind", list(y,
                                                                Subclass1 = AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms4[match(as.character(y$V1), as.character(AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms4$ChemicalFormulaWithOnes)), "LipidSubclass"],
                                                                
                                                                Subclass2 = AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms4[match(as.character(y$V2), as.character(AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms4$ChemicalFormulaWithOnes)), "LipidSubclass"],
                                                                SubclassesMatchingToConsensusSubclassesAndClasses[match(AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms4[match(as.character(y$V1), as.character(AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms4$ChemicalFormulaWithOnes)), "LipidSubclass"], SubclassesMatchingToConsensusSubclassesAndClasses[,1]),2:3],
                                                                
                                                                SubclassesMatchingToConsensusSubclassesAndClasses[match(AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms4[match(as.character(y$V2), as.character(AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms4$ChemicalFormulaWithOnes)), "LipidSubclass"], SubclassesMatchingToConsensusSubclassesAndClasses[,1]),2:3]))



y <- CooccurrenceOfChemicalPairsListx[[3]]

CooccurrenceOfChemicalPairsList2x[[3]] <- do.call("cbind", list(y,
                                                                Subclass1 = AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms4[match(as.character(y$V1), as.character(AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms4$ChemicalFormulaWithOnes)), "LipidSubclass"],
                                                                
                                                                Subclass2 = AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms4[match(as.character(y$V2), as.character(AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms4$ChemicalFormulaWithOnes)), "LipidSubclass"],
                                                                SubclassesMatchingToConsensusSubclassesAndClasses[match(AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms4[match(as.character(y$V1), as.character(AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms4$ChemicalFormulaWithOnes)), "LipidSubclass"], SubclassesMatchingToConsensusSubclassesAndClasses[,1]),2:3],
                                                                
                                                                SubclassesMatchingToConsensusSubclassesAndClasses[match(AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms4[match(as.character(y$V2), as.character(AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms4$ChemicalFormulaWithOnes)), "LipidSubclass"], SubclassesMatchingToConsensusSubclassesAndClasses[,1]),2:3]))



y <- CooccurrenceOfChemicalPairsListx[[4]]

CooccurrenceOfChemicalPairsList2x[[4]] <- do.call("cbind", list(y,
                                                                Subclass1 = AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms4[match(as.character(y$V1), as.character(AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms4$ChemicalFormulaWithOnes)), "LipidSubclass"],
                                                                
                                                                Subclass2 = AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms4[match(as.character(y$V2), as.character(AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms4$ChemicalFormulaWithOnes)), "LipidSubclass"],
                                                                SubclassesMatchingToConsensusSubclassesAndClasses[match(AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms4[match(as.character(y$V1), as.character(AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms4$ChemicalFormulaWithOnes)), "LipidSubclass"], SubclassesMatchingToConsensusSubclassesAndClasses[,1]),2:3],
                                                                
                                                                SubclassesMatchingToConsensusSubclassesAndClasses[match(AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms4[match(as.character(y$V2), as.character(AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms4$ChemicalFormulaWithOnes)), "LipidSubclass"], SubclassesMatchingToConsensusSubclassesAndClasses[,1]),2:3]))


for(i in 2:4){
  CooccurrenceOfChemicalPairsList2x[[i]] <- cbind(CooccurrenceOfChemicalPairsList2x[[i]], LipidPairMatchingType = (CooccurrenceOfChemicalPairsList2x[[i]][,8] == CooccurrenceOfChemicalPairsList2x[[i]][,10]) + (CooccurrenceOfChemicalPairsList2x[[i]][,9] == CooccurrenceOfChemicalPairsList2x[[i]][,11]))
  
}



CooccurrenceOfChemicalPairsList2x[[1]] <- do.call("cbind", list(CooccurrenceOfChemicalPairsList2x[[1]],
                                                                
                                                                SubclassesMatchingToConsensusSubclassesAndClasses[match(CooccurrenceOfChemicalPairsList2x[[1]]$Subclass1, SubclassesMatchingToConsensusSubclassesAndClasses[,1]),2:3],
                                                                SubclassesMatchingToConsensusSubclassesAndClasses[match(CooccurrenceOfChemicalPairsList2x[[1]]$Subclass2, SubclassesMatchingToConsensusSubclassesAndClasses[,1]),2:3]))

CooccurrenceOfChemicalPairsList2x[[1]] <- cbind(CooccurrenceOfChemicalPairsList2x[[1]], LipidPairMatchingType = (CooccurrenceOfChemicalPairsList2x[[1]][,10] == CooccurrenceOfChemicalPairsList2x[[1]][,12]) + (CooccurrenceOfChemicalPairsList2x[[1]][,11] == CooccurrenceOfChemicalPairsList2x[[1]][,13]))


CooccurrenceOfChemicalPairsListExpandedWithAllx <- CooccurrenceOfChemicalPairsList2x[[1]]
CooccurrenceOfChemicalPairsListExpandedWithAllx$LipidPairMatchingType <- "All"

CooccurrenceOfChemicalPairsListExpandedWithAll2x <- rbind(CooccurrenceOfChemicalPairsList2x[[1]], CooccurrenceOfChemicalPairsListExpandedWithAllx)



PossibleCooccurrenceOfChemicalPairsx <- rbind(cbind(CooccurrenceOfChemicalPairsList2x[[3]], Screen = "A"), cbind(CooccurrenceOfChemicalPairsList2x[[4]], Screen = "E"))

PossibleCooccurrenceOfChemicalPairs2x <- PossibleCooccurrenceOfChemicalPairsx
PossibleCooccurrenceOfChemicalPairs2x$LipidPairMatchingType <- "All"

PossibleCooccurrenceOfChemicalPairs2x <- rbind(PossibleCooccurrenceOfChemicalPairsx, PossibleCooccurrenceOfChemicalPairs2x)


library(RColorBrewer)
pdf("./Output/Panel5BWithTailoredScriptAndMoreDetailedVisualizationNowIncludesAllTooWithBinSize001170220222NoSegmentStacking18032022.pdf")

BinSize <- 0.01
MaxList2TissuesColoc <- list()

StatListTissuesColoc <- list()
StatList2TissuesColoc <- list()

StatListTissuesColocX2 <- list()
StatList2TissuesColocX2 <- list()

ScreenColorsForFilling <- cbind(c("A","E"),c(brewer.pal(9,"Blues")[4], brewer.pal(9,"Oranges")[4]))
ScreenColorsForFillingPot <- cbind(c("A","E"),c(brewer.pal(9,"Blues")[7], brewer.pal(9,"Oranges")[7]))


for(SimilarityIndex in c("All", as.character(0:2))){
  
  x1a <- na.omit(as.numeric(CooccurrenceOfChemicalPairsListExpandedWithAll2x[(CooccurrenceOfChemicalPairsListExpandedWithAll2x[,"Screen"] == "A") & (CooccurrenceOfChemicalPairsListExpandedWithAll2x[,"LipidPairMatchingType"] == SimilarityIndex),"Cooccurrence"]))
  x1e <- na.omit(as.numeric(CooccurrenceOfChemicalPairsListExpandedWithAll2x[(CooccurrenceOfChemicalPairsListExpandedWithAll2x[,"Screen"] == "E") & (CooccurrenceOfChemicalPairsListExpandedWithAll2x[,"LipidPairMatchingType"] == SimilarityIndex),"Cooccurrence"]))
  
  x1apot <- na.omit(as.numeric(PossibleCooccurrenceOfChemicalPairs2x[(PossibleCooccurrenceOfChemicalPairs2x[,"Screen"] == "A") & (PossibleCooccurrenceOfChemicalPairs2x[,"LipidPairMatchingType"] == SimilarityIndex),"Cooccurrence"]))
  x1epot <- na.omit(as.numeric(PossibleCooccurrenceOfChemicalPairs2x[(PossibleCooccurrenceOfChemicalPairs2x[,"Screen"] == "E") & (PossibleCooccurrenceOfChemicalPairs2x[,"LipidPairMatchingType"] == SimilarityIndex),"Cooccurrence"]))
  
  x2 <- LipidSubsetMMeans_4Subset[LipidSubsetMMeans_4Subset[,"LipidPairMatchingType"] == SimilarityIndex, "Cooccurrence"]
  if(length(x2) > 0){
    
    MaxY <- max(c(unlist(density(x1a, from = -1, to = 1)["y"]), unlist(density(x1e, from = -1, to = 1)["y"]), unlist(density(x1apot, from = -1, to = 1)["y"]), unlist(density(x1epot, from = -1, to = 1)["y"]), unlist(density(x2, from = -1, to = 1)["y"])))
    plot(x = unlist(density(x1a, from = 0, to = 1)["x"]), y = unlist(density(x1a, from = -1, to = 1)["y"]), xlab = "Co-occurrence (tissues)", ylab = "Normalized Density", ylim = c(-MaxY,MaxY), xlim = c(0,1), col = "white", main = paste0("Similarity: ", SimilarityIndex), type = "l", bty="n")
    
    polygon(x = c(min(unlist(density(x1a, from = -1, to = 1)["x"])), unlist(density(x1a, from = -1, to = 1)["x"]), max(unlist(density(x1a, from = -1, to = 1)["x"]))),
            y = c(0, unlist(density(x1a, from = -1, to = 1)["y"]), 0), col = ScreenColorsForFilling[ScreenColorsForFilling[,1] == "A",2], border = "white")
    
    polygon(x = c(min(unlist(density(x1e, from = -1, to = 1)["x"])), unlist(density(x1e, from = -1, to = 1)["x"]), max(unlist(density(x1e, from = -1, to = 1)["x"]))),
            y = c(0, -unlist(density(x1e, from = -1, to = 1)["y"]), 0), col = ScreenColorsForFilling[ScreenColorsForFilling[,1] == "E",2], border = "white")
    
    segments(x1a, 0.16, x1 = x1a, y1 = 0, col = "black")
    segments(x1e, -0.16, x1 = x1e, y1 = 0, col = "black")
    
    lines(c(0,unlist(density(x1apot, from = 0, to = 1)["x"]),1), c(0,unlist(density(x1apot, from = 0, to = 1)["y"]),0), col = ScreenColorsForFillingPot[ScreenColorsForFillingPot[,1] == "A",2], lwd = 2)
    lines(c(0,unlist(density(x1epot, from = 0, to = 1)["x"]),1), c(0,-unlist(density(x1epot, from = 0, to = 1)["y"]),0), col = ScreenColorsForFillingPot[ScreenColorsForFillingPot[,1] == "E",2], lwd = 2)
    
    lines(c(0,unlist(density(x2, from = 0, to = 1)["x"]),1), c(0,unlist(density(x2, from = 0, to = 1)["y"]),0), col = "DarkGrey", lwd = 2)
    lines(c(0,unlist(density(x2, from = 0, to = 1)["x"]),1), c(0,-unlist(density(x2, from = 0, to = 1)["y"]),0), col = "DarkGrey", lwd = 2)
    
    segments(x0 = median(x1a), x1 = median(x1a), y0 = MaxY, y1 = 0, lwd = 2)
    segments(x0 = median(x1e), x1 = median(x1e), y0 = -MaxY, y1 = 0, lwd = 2)
    
    segments(x0 = median(x1apot), x1 = median(x1apot), y0 = MaxY, y1 = 0, lwd = 2, col = ScreenColorsForFillingPot[ScreenColorsForFillingPot[,1] == "A",2])
    segments(x0 = median(x1epot), x1 = median(x1epot), y0 = -MaxY, y1 = 0, lwd = 2, col = ScreenColorsForFillingPot[ScreenColorsForFillingPot[,1] == "E",2])
    
    segments(x0 = median(x2), x1 = median(x2), y0 = MaxY, y1 = 0, lwd = 2, col = "DarkGrey")
    segments(x0 = median(x2), x1 = median(x2), y0 = -MaxY, y1 = 0, lwd = 2, col = "DarkGrey")
    
    MaxList2TissuesColoc[[SimilarityIndex]][["A"]] <- c(max(unlist(density(x1a, from = 0, to = 1)["y"])),max(unlist(density(x2, from = 0, to = 1)["y"])))
    MaxList2TissuesColoc[[SimilarityIndex]][["E"]] <- c(max(unlist(density(x1e, from = 0, to = 1)["y"])),max(unlist(density(x2, from = 0, to = 1)["y"])))
    
    
  }else{
    
    MaxY <- max(c(unlist(density(x1a, from = -1, to = 1)["y"]), unlist(density(x1e, from = -1, to = 1)["y"]), unlist(density(x1apot, from = -1, to = 1)["y"]), unlist(density(x1epot, from = -1, to = 1)["y"])))
    plot(x = unlist(density(x1a, from = 0, to = 1)["x"]), y = unlist(density(x1a, from = -1, to = 1)["y"]), xlab = "Co-occurrence (tissues)", ylab = "Normalized Density", ylim = c(-MaxY,MaxY), xlim = c(0,1), col = "white", main = paste0("Similarity: ", SimilarityIndex), type = "l", bty="n")
    
    polygon(x = c(min(unlist(density(x1a, from = -1, to = 1)["x"])), unlist(density(x1a, from = -1, to = 1)["x"]), max(unlist(density(x1a, from = -1, to = 1)["x"]))),
            y = c(0, unlist(density(x1a, from = -1, to = 1)["y"]), 0), col = ScreenColorsForFilling[ScreenColorsForFilling[,1] == "A",2], border = "white")
    
    polygon(x = c(min(unlist(density(x1e, from = -1, to = 1)["x"])), unlist(density(x1e, from = -1, to = 1)["x"]), max(unlist(density(x1e, from = -1, to = 1)["x"]))),
            y = c(0, -unlist(density(x1e, from = -1, to = 1)["y"]), 0), col = ScreenColorsForFilling[ScreenColorsForFilling[,1] == "E",2], border = "white")
    
    segments(x1a, 0.16, x1 = x1a, y1 = 0, col = "black")
    segments(x1e, -0.16, x1 = x1e, y1 = 0, col = "black")
    
    lines(c(0,unlist(density(x1apot, from = 0, to = 1)["x"]),1), c(0,unlist(density(x1apot, from = 0, to = 1)["y"]),0), col = ScreenColorsForFillingPot[ScreenColorsForFillingPot[,1] == "A",2], lwd = 2)
    lines(c(0,unlist(density(x1epot, from = 0, to = 1)["x"]),1), c(0,-unlist(density(x1epot, from = 0, to = 1)["y"]),0), col = ScreenColorsForFillingPot[ScreenColorsForFillingPot[,1] == "E",2], lwd = 2)
    
    segments(x0 = median(x1a), x1 = median(x1a), y0 = MaxY, y1 = 0, lwd = 2)
    segments(x0 = median(x1e), x1 = median(x1e), y0 = -MaxY, y1 = 0, lwd = 2)
    
    segments(x0 = median(x1apot), x1 = median(x1apot), y0 = MaxY, y1 = 0, lwd = 2, col = ScreenColorsForFillingPot[ScreenColorsForFillingPot[,1] == "A",2])
    segments(x0 = median(x1epot), x1 = median(x1epot), y0 = -MaxY, y1 = 0, lwd = 2, col = ScreenColorsForFillingPot[ScreenColorsForFillingPot[,1] == "E",2])
    
    MaxList2TissuesColoc[[SimilarityIndex]][["A"]] <- c(max(unlist(density(x1a, from = 0, to = 1)["y"])),0)
    MaxList2TissuesColoc[[SimilarityIndex]][["E"]] <- c(max(unlist(density(x1e, from = 0, to = 1)["y"])),0)    
    
  }
}


# Generalize the maximum for y: all figures on the same scale

SimilarityIndex <- "All"
MaxY <- max(unlist(MaxList2TissuesColoc)) 

MaxList2TissuesColoc <- list()



for(SimilarityIndex in c("All", as.character(0:2))){
  
  x1a <- na.omit(as.numeric(CooccurrenceOfChemicalPairsListExpandedWithAll2x[(CooccurrenceOfChemicalPairsListExpandedWithAll2x[,"Screen"] == "A") & (CooccurrenceOfChemicalPairsListExpandedWithAll2x[,"LipidPairMatchingType"] == SimilarityIndex),"Cooccurrence"]))
  x1e <- na.omit(as.numeric(CooccurrenceOfChemicalPairsListExpandedWithAll2x[(CooccurrenceOfChemicalPairsListExpandedWithAll2x[,"Screen"] == "E") & (CooccurrenceOfChemicalPairsListExpandedWithAll2x[,"LipidPairMatchingType"] == SimilarityIndex),"Cooccurrence"]))
  
  x1apot <- na.omit(as.numeric(PossibleCooccurrenceOfChemicalPairs2x[(PossibleCooccurrenceOfChemicalPairs2x[,"Screen"] == "A") & (PossibleCooccurrenceOfChemicalPairs2x[,"LipidPairMatchingType"] == SimilarityIndex),"Cooccurrence"]))
  x1epot <- na.omit(as.numeric(PossibleCooccurrenceOfChemicalPairs2x[(PossibleCooccurrenceOfChemicalPairs2x[,"Screen"] == "E") & (PossibleCooccurrenceOfChemicalPairs2x[,"LipidPairMatchingType"] == SimilarityIndex),"Cooccurrence"]))
  
  x2 <- LipidSubsetMMeans_4Subset[LipidSubsetMMeans_4Subset[,"LipidPairMatchingType"] == SimilarityIndex, "Cooccurrence"]
  if(length(x2) > 0){
    
    
    plot(x = unlist(density(x1a, from = 0, to = 1)["x"]), y = unlist(density(x1a, from = -1, to = 1)["y"]), xlab = "Co-occurrence (tissues)", ylab = "Normalized Density", ylim = c(-MaxY,MaxY), xlim = c(0,1), col = "white", main = paste0("Similarity: ", SimilarityIndex), type = "l", bty="n")
    
    polygon(x = c(min(unlist(density(x1a, from = -1, to = 1)["x"])), unlist(density(x1a, from = -1, to = 1)["x"]), max(unlist(density(x1a, from = -1, to = 1)["x"]))),
            y = c(0, unlist(density(x1a, from = -1, to = 1)["y"]), 0), col = ScreenColorsForFilling[ScreenColorsForFilling[,1] == "A",2], border = "white")
    
    polygon(x = c(min(unlist(density(x1e, from = -1, to = 1)["x"])), unlist(density(x1e, from = -1, to = 1)["x"]), max(unlist(density(x1e, from = -1, to = 1)["x"]))),
            y = c(0, -unlist(density(x1e, from = -1, to = 1)["y"]), 0), col = ScreenColorsForFilling[ScreenColorsForFilling[,1] == "E",2], border = "white")
    
    segments(x1a, 0.16, x1 = x1a, y1 = 0, col = "black")
    segments(x1e, -0.16, x1 = x1e, y1 = 0, col = "black")
    
    lines(c(0,unlist(density(x1apot, from = 0, to = 1)["x"]),1), c(0,unlist(density(x1apot, from = 0, to = 1)["y"]),0), col = ScreenColorsForFillingPot[ScreenColorsForFillingPot[,1] == "A",2], lwd = 2)
    lines(c(0,unlist(density(x1epot, from = 0, to = 1)["x"]),1), c(0,-unlist(density(x1epot, from = 0, to = 1)["y"]),0), col = ScreenColorsForFillingPot[ScreenColorsForFillingPot[,1] == "E",2], lwd = 2)
    
    lines(c(0,unlist(density(x2, from = 0, to = 1)["x"]),1), c(0,unlist(density(x2, from = 0, to = 1)["y"]),0), col = "DarkGrey", lwd = 2)
    lines(c(0,unlist(density(x2, from = 0, to = 1)["x"]),1), c(0,-unlist(density(x2, from = 0, to = 1)["y"]),0), col = "DarkGrey", lwd = 2)
    
    segments(x0 = median(x1a), x1 = median(x1a), y0 = MaxY, y1 = 0, lwd = 2)
    segments(x0 = median(x1e), x1 = median(x1e), y0 = -MaxY, y1 = 0, lwd = 2)
    
    segments(x0 = median(x1apot), x1 = median(x1apot), y0 = MaxY, y1 = 0, lwd = 2, col = ScreenColorsForFillingPot[ScreenColorsForFillingPot[,1] == "A",2])
    segments(x0 = median(x1epot), x1 = median(x1epot), y0 = -MaxY, y1 = 0, lwd = 2, col = ScreenColorsForFillingPot[ScreenColorsForFillingPot[,1] == "E",2])
    
    segments(x0 = median(x2), x1 = median(x2), y0 = MaxY, y1 = 0, lwd = 2, col = "DarkGrey")
    segments(x0 = median(x2), x1 = median(x2), y0 = -MaxY, y1 = 0, lwd = 2, col = "DarkGrey")
    
    MaxList2TissuesColoc[[SimilarityIndex]][["A"]] <- c(max(unlist(density(x1a, from = 0, to = 1)["y"])),max(unlist(density(x2, from = 0, to = 1)["y"])))
    MaxList2TissuesColoc[[SimilarityIndex]][["E"]] <- c(max(unlist(density(x1e, from = 0, to = 1)["y"])),max(unlist(density(x2, from = 0, to = 1)["y"])))
    
    StatListTissuesColoc[[SimilarityIndex]] <- structure(rbind(FishersExactTestBelowAbove05(SubsetX = x1a, SupersetY = c(x1apot,x1apot)), FishersExactTestBelowAbove05(SubsetX = x1e, SupersetY = c(x1epot,x1epot))), 
                                                         dimnames = list(c("A","E"), c("NoNeg", "NoPos", "YesNeg", "YesPos", "LowConf", "HighConf", "OddsRatio", "pValue")))
    
    StatList2TissuesColoc[[SimilarityIndex]] <- structure(rbind(FishersExactTestNegPosEqualBackground2(SubsetX = x1a, SupersetY = c(x1apot,x1apot)), FishersExactTestNegPosEqualBackground2(SubsetX = x1e, SupersetY = c(x1epot,x1epot))), 
                                                          dimnames = list(c("A","E"), c("NoNeg", "NoPos", "YesNeg", "YesPos", "LowConf", "HighConf", "OddsRatio", "pValue")))
    
    StatListTissuesColocX2[[SimilarityIndex]] <- structure(rbind(FishersExactTestBelowAbove05(SubsetX = x1a, SupersetY = c(x2,x2)), FishersExactTestBelowAbove05(SubsetX = x1e, SupersetY = c(x2,x2))), 
                                                           dimnames = list(c("A","E"), c("NoNeg", "NoPos", "YesNeg", "YesPos", "LowConf", "HighConf", "OddsRatio", "pValue")))
    
    StatList2TissuesColocX2[[SimilarityIndex]] <- structure(rbind(FishersExactTestNegPosEqualBackground2(SubsetX = x1a, SupersetY = c(x2,x2)), FishersExactTestNegPosEqualBackground2(SubsetX = x1e, SupersetY = c(x2,x2))), 
                                                            dimnames = list(c("A","E"), c("NoNeg", "NoPos", "YesNeg", "YesPos", "LowConf", "HighConf", "OddsRatio", "pValue")))
    
  }else{
    plot(x = unlist(density(x1a, from = 0, to = 1)["x"]), y = unlist(density(x1a, from = -1, to = 1)["y"]), xlab = "Co-occurrence (tissues)", ylab = "Normalized Density", ylim = c(-MaxY,MaxY), xlim = c(0,1), col = "white", main = paste0("Similarity: ", SimilarityIndex), type = "l", bty="n")
    
    polygon(x = c(min(unlist(density(x1a, from = -1, to = 1)["x"])), unlist(density(x1a, from = -1, to = 1)["x"]), max(unlist(density(x1a, from = -1, to = 1)["x"]))),
            y = c(0, unlist(density(x1a, from = -1, to = 1)["y"]), 0), col = ScreenColorsForFilling[ScreenColorsForFilling[,1] == "A",2], border = "white")
    
    polygon(x = c(min(unlist(density(x1e, from = -1, to = 1)["x"])), unlist(density(x1e, from = -1, to = 1)["x"]), max(unlist(density(x1e, from = -1, to = 1)["x"]))),
            y = c(0, -unlist(density(x1e, from = -1, to = 1)["y"]), 0), col = ScreenColorsForFilling[ScreenColorsForFilling[,1] == "E",2], border = "white")
    
    segments(x1a, 0.16, x1 = x1a, y1 = 0, col = "black")
    segments(x1e, -0.16, x1 = x1e, y1 = 0, col = "black")
    
    lines(c(0,unlist(density(x1apot, from = 0, to = 1)["x"]),1), c(0,unlist(density(x1apot, from = 0, to = 1)["y"]),0), col = ScreenColorsForFillingPot[ScreenColorsForFillingPot[,1] == "A",2], lwd = 2)
    lines(c(0,unlist(density(x1epot, from = 0, to = 1)["x"]),1), c(0,-unlist(density(x1epot, from = 0, to = 1)["y"]),0), col = ScreenColorsForFillingPot[ScreenColorsForFillingPot[,1] == "E",2], lwd = 2)
    
    segments(x0 = median(x1a), x1 = median(x1a), y0 = MaxY, y1 = 0, lwd = 2)
    segments(x0 = median(x1e), x1 = median(x1e), y0 = -MaxY, y1 = 0, lwd = 2)
    
    segments(x0 = median(x1apot), x1 = median(x1apot), y0 = MaxY, y1 = 0, lwd = 2, col = ScreenColorsForFillingPot[ScreenColorsForFillingPot[,1] == "A",2])
    segments(x0 = median(x1epot), x1 = median(x1epot), y0 = -MaxY, y1 = 0, lwd = 2, col = ScreenColorsForFillingPot[ScreenColorsForFillingPot[,1] == "E",2])
    
    MaxList2TissuesColoc[[SimilarityIndex]][["A"]] <- c(max(unlist(density(x1a, from = 0, to = 1)["y"])),0)
    MaxList2TissuesColoc[[SimilarityIndex]][["E"]] <- c(max(unlist(density(x1e, from = 0, to = 1)["y"])),0)    
    
    StatListTissuesColoc[[SimilarityIndex]] <- structure(rbind(FishersExactTestBelowAbove05(SubsetX = x1a, SupersetY = c(x1apot,x1apot)), FishersExactTestBelowAbove05(SubsetX = x1e, SupersetY = c(x1epot,x1epot))), 
                                                         dimnames = list(c("A","E"), c("NoNeg", "NoPos", "YesNeg", "YesPos", "LowConf", "HighConf", "OddsRatio", "pValue")))
    
    StatList2TissuesColoc[[SimilarityIndex]] <- structure(rbind(FishersExactTestNegPosEqualBackground2(SubsetX = x1a, SupersetY = c(x1apot,x1apot)), FishersExactTestNegPosEqualBackground2(SubsetX = x1e, SupersetY = c(x1epot,x1epot))), 
                                                          dimnames = list(c("A","E"), c("NoNeg", "NoPos", "YesNeg", "YesPos", "LowConf", "HighConf", "OddsRatio", "pValue")))
    
  }
}

dev.off()  


# All of the figures in the second set of figure pages in the previous document were combined into Panel 4b and further esthetically enhanced and integrated with the other parts of Figure 4 in Adobe Illustrator.
# For the statistical significances the results of the Fisher exact tests on the difference in the medians was used (see further), and also added in Adobe Illustrator.


#### Panel 4C

#. Panel5CWithoutTheNrd0WdAdaptedWithScalesImplemented22032022.pdf (#) (Name of the finally used version for inclusion in figure 5.)
#. Panel5CWithTheNrd0WdAdaptedWithScalesImplemented22032022.pdf (#) (Not the basis of the final version.)

#. Panel5CWithTailoredScriptAndMoreDetailedVisualizationBinSize001NoSegmentStacking18032022.pdf # (Finally used version for inclusion in figure 5 based on this.)
#. Panel5CWithGreenUnityDistributionAndBwWithNrd0Wd01422032022.pdf # (Not the basis of the final version.)

# How the document of the average subcellular localizations of lipid species was generated based on previous data 
# Now as input the further downstream document

# MolPctOfMads <- Col1ToRowNames(read.csv(file = "./InputData/processed.molpct_generatedbyMads.csv", header = TRUE, sep = ",", as.is = TRUE))[,1:38] # Get rid of empty space at end columns
# ExperimentTypesOfMads <- sapply(strsplit(colnames(MolPctOfMads), "\\."), "[",1)
# 
# 
# AggregateMPOfMads <- sapply(unique(ExperimentTypesOfMads), function(x){rowMeans(MolPctOfMads[,ExperimentTypesOfMads == x], na.rm = TRUE)})
# 
# AggregateMPOfMads0 <- AggregateMPOfMads
# AggregateMPOfMads0[is.na(AggregateMPOfMads0)] <- 0
# 
# 
# AggregateMPOfMads0RTL <- setNames(split(AggregateMPOfMads0, row(AggregateMPOfMads0)), rownames(AggregateMPOfMads0))
# 
# AggregateMPOfMads0RTLMOC <- do.call("cbind", setNames(lapply(AggregateMPOfMads0RTL, function(y){unlist(lapply(AggregateMPOfMads0RTL, function(x){MOCFunction(y, x)}))}), rownames(AggregateMPOfMads0)))


# Save the cleaned matrix of the averages of the subcellular lipid localization as a compact RDS-file
# saveRDS(object = AggregateMPOfMads0RTLMOC, file = "./InputData/SubcellularLocalizationAveragesLipids.rds")

# Export Extended Data Table 9A
# See further

# Load saved RDS-image of this file (if needed, otherwise this can be skipped)
AggregateMPOfMads0RTLMOCx <- readRDS("./InputData/SubcellularLocalizationAveragesLipids.rds")

# Check that nothing went wrong during the conversion (if needed, otherwise this can be skipped)
# identical(AggregateMPOfMads0RTLMOCx, AggregateMPOfMads0RTLMOC)

# Reassign reloaded cleaned-up document with averages of the subcellular localizations of the lipid species
AggregateMPOfMads0RTLMOC <- AggregateMPOfMads0RTLMOCx

# Remove the original large file
rm(AggregateMPOfMads0RTLMOCx)

MadsLipids <- do.call("rbind", strsplit(gsub(";2","", rownames(AggregateMPOfMads0RTLMOC)), " "))
MadsLipids <- cbind(MadsLipids, do.call("rbind", strsplit(MadsLipids[,2], ":")))

colnames(MadsLipids) <- c("LipidSubclassFromMads", "TotalLipidCarbonChainMads", "LipidCarbonChainLengthMads", "LipidCarbonChainUnsaturationMads")
AggregateMPOfMads0RTLMOC2 <- cbind(MadsLipids, AggregateMPOfMads0RTLMOC)


MadsLipidsSubclassesConverted <- MadsLipids[,1]

for(i in list(c("Cer","d*Cer"), c("HexCer","d*HexCer"), c("LPCO-","LPC-O"), c("LPEO-","LPE-O"), c("PCO-","PC-O"), c("PEO-","PE-O"), c("SM","d*SM"))){
  MadsLipidsSubclassesConverted[MadsLipidsSubclassesConverted == i[1]] <- i[2]
  
}
AggregateMPOfMads0RTLMOC4 <- cbind(CovertedSubclasses = MadsLipidsSubclassesConverted, AggregateMPOfMads0RTLMOC2)

rm(AggregateMPOfMads0RTLMOC)
rm(AggregateMPOfMads0RTLMOC2)

# Optional but advisable here because of size
# gc()

x <- AggregateMPOfMads0RTLMOC4[,-(1:5)]
y <- which(lower.tri(x), arr.ind = TRUE)

AggregateMPOfMads0RTLMOC5 <- do.call("cbind", list(AggregateMPOfMads0RTLMOC4[y[,1],1:5], AggregateMPOfMads0RTLMOC4[y[,2],1:5], SubcellularColocalization = sapply(1:dim(y)[1], function(z){x[y[z,1],y[z,2]]})))
SubclassToClassConversionMatrix <- cbind(Subclasses = unique(AggregateMPOfMads0RTLMOC5[,"CovertedSubclasses"]), Classes = c("CE","Cer","Chol","CL","DAG","FA","HexCer","PC","PC","PE","PE","PG","PI","PA","PC","PC","PE","PE","PG","PI","PS","SM"))

AggregateMPOfMads0RTLMOC7 <- do.call("cbind", list(ClassFrom = SubclassToClassConversionMatrix[match(AggregateMPOfMads0RTLMOC5[,1], SubclassToClassConversionMatrix[,1]),2],
                                                   ClassTo = SubclassToClassConversionMatrix[match(AggregateMPOfMads0RTLMOC5[,6], SubclassToClassConversionMatrix[,1]),2],
                                                   
                                                   AggregateMPOfMads0RTLMOC5))


AggregateMPOfMads0RTLMOC8 <- cbind(PairType = sapply(1:dim(AggregateMPOfMads0RTLMOC7)[1], function(x){if(AggregateMPOfMads0RTLMOC7[x,1] == AggregateMPOfMads0RTLMOC7[x,2]){
  if(AggregateMPOfMads0RTLMOC7[x,3] == AggregateMPOfMads0RTLMOC7[x,8]){"Species"}else{"Sub-Class"}
  
}else{"Class"}}), 
AggregateMPOfMads0RTLMOC7)

AggregateMPOfMads0RTLMOC10 <- cbind(Pairs = paste0(AggregateMPOfMads0RTLMOC8[,4], "(", AggregateMPOfMads0RTLMOC8[,6], ")_", AggregateMPOfMads0RTLMOC8[,9], "(", AggregateMPOfMads0RTLMOC8[,11], ")"),
                                    AggregateMPOfMads0RTLMOC8)


SubclassGeneralizationMatrixOurData <- do.call("rbind", list(c("BMP","PG"), c("PG/BMP","PG"), c("DHSM", "d*SM"), c("dCer", "d*Cer"), c("DHCer", "d*Cer")))

SubclassesForPureAntonella <- PureAntonella32b[,"LikelySubclass"]
SubclassesForPureEnric <- PureEnric32[,"LikelySubclass"]

for(i in 1:dim(SubclassGeneralizationMatrixOurData)[1]){
  SubclassesForPureAntonella[SubclassesForPureAntonella == SubclassGeneralizationMatrixOurData[i,1]] <- SubclassGeneralizationMatrixOurData[i,2]}

InCelluloReducedLinks <- unique(do.call("cbind", list(PureAntonella32b[,c("LTPProtein","Screen")], 
                                                      SpeciesName = paste0(SubclassesForPureAntonella,"(",PureAntonella32b[,"TotalCarbonChainLength"],":",PureAntonella32b[,"TotalCarbonChainUnsaturations"], ")")
                                                      
)))


for(i in 1:dim(SubclassGeneralizationMatrixOurData)[1]){
  SubclassesForPureEnric[SubclassesForPureEnric == SubclassGeneralizationMatrixOurData[i,1]] <- SubclassGeneralizationMatrixOurData[i,2]}

InVitroReducedLinks <- unique(do.call("cbind", list(PureEnric32[,c("LTPProtein","Screen")], 
                                                    SpeciesName = paste0(SubclassesForPureEnric,"(",PureEnric32[,"TotalCarbonChainLength"],":",PureEnric32[,"TotalCarbonChainUnsaturations"], ")")
                                                    
)))



InCelluloReducedLinksCouples <- do.call("rbind", lapply(unique(InCelluloReducedLinks[,1]), function(x){
  
  if(sum(InCelluloReducedLinks[,1] == x) > 1){
    cbind(unique(InCelluloReducedLinks[InCelluloReducedLinks[,1] == x, 1:2]), t(combn(as.character(InCelluloReducedLinks[InCelluloReducedLinks[,1] == x, 3]), m = 2)), row.names = NULL)
    
  }}))
# Cleaned rownames allocation up.


InVitroReducedLinksCouples <- do.call("rbind", lapply(unique(InVitroReducedLinks[,1]), function(x){
  
  if(sum(InVitroReducedLinks[,1] == x) > 1){
    cbind(unique(InVitroReducedLinks[InVitroReducedLinks[,1] == x, 1:2]), t(combn(as.character(InVitroReducedLinks[InVitroReducedLinks[,1] == x, 3]), m = 2)), row.names = NULL)
    
  }}))
# Cleaned rownames allocation up.

InCelluloReducedLinksCouples2 <- do.call("cbind", list(InCelluloReducedLinksCouples, paste(InCelluloReducedLinksCouples[,3], InCelluloReducedLinksCouples[,4], sep = "_"), paste(InCelluloReducedLinksCouples[,4], InCelluloReducedLinksCouples[,3], sep = "_")))
colnames(InCelluloReducedLinksCouples2) <- c("LTPProtein", "Screen", "FromLipid", "ToLipid", "FromTo", "ToFrom")

InVitroReducedLinksCouples2 <- do.call("cbind", list(InVitroReducedLinksCouples, paste(InVitroReducedLinksCouples[,3], InVitroReducedLinksCouples[,4], sep = "_"), paste(InVitroReducedLinksCouples[,4], InVitroReducedLinksCouples[,3], sep = "_")))
colnames(InVitroReducedLinksCouples2) <- c("LTPProtein", "Screen", "FromLipid", "ToLipid", "FromTo", "ToFrom")

AggregateMPOfMads0RTLMOC11x <- do.call("cbind", list(SubclassPairPresenceInData = as.character(AggregateMPOfMads0RTLMOC10[,5]) %in% unique(c(as.character(SubclassesForPureAntonella), as.character(SubclassesForPureEnric))) &
                                                       as.character(AggregateMPOfMads0RTLMOC10[,10]) %in% unique(c(as.character(SubclassesForPureAntonella), as.character(SubclassesForPureEnric))), 
                                                     
                                                     AggregateMPOfMads0RTLMOC10,
                                                     
                                                     
                                                     SubclassPairPresenceInCellulo = as.character(AggregateMPOfMads0RTLMOC10[,5]) %in% unique(as.character(SubclassesForPureAntonella)) &
                                                       as.character(AggregateMPOfMads0RTLMOC10[,10]) %in% unique(as.character(SubclassesForPureAntonella)), 
                                                     
                                                     SubclassPairPresenceInVitro = as.character(AggregateMPOfMads0RTLMOC10[,5]) %in% unique(as.character(SubclassesForPureEnric)) &
                                                       as.character(AggregateMPOfMads0RTLMOC10[,10]) %in% unique(as.character(SubclassesForPureEnric))
                                                     
))
AggregateMPOfMads0RTLMOC11xs <- AggregateMPOfMads0RTLMOC11x[as.logical(AggregateMPOfMads0RTLMOC11x[,"SubclassPairPresenceInData"]),]

AggregateMPOfMads0RTLMOC13xs <- cbind(MatchWithInCellulo = AggregateMPOfMads0RTLMOC11xs[,"Pairs"] %in% unique(c(as.character(InCelluloReducedLinksCouples2$FromTo), as.character(InCelluloReducedLinksCouples2$ToFrom))),
                                      AggregateMPOfMads0RTLMOC11xs)

AggregateMPOfMads0RTLMOC13xs2 <- cbind(MatchWithInVitro = AggregateMPOfMads0RTLMOC13xs[,"Pairs"] %in% unique(c(as.character(InVitroReducedLinksCouples2$FromTo), as.character(InVitroReducedLinksCouples2$ToFrom))),
                                       AggregateMPOfMads0RTLMOC13xs)

SubcellularLocalizationOverlapOurDatax <- do.call("rbind", list(cbind(AggregateMPOfMads0RTLMOC13xs2[as.logical(AggregateMPOfMads0RTLMOC13xs2[,"MatchWithInCellulo"]), c("SubcellularColocalization", "Pairs", "PairType")], Screen = "In Cellulo"),
                                                                cbind(AggregateMPOfMads0RTLMOC13xs2[as.logical(AggregateMPOfMads0RTLMOC13xs2[,"MatchWithInVitro"]), c("SubcellularColocalization", "Pairs","PairType")], Screen = "In Vitro"),
                                                                
                                                                do.call("cbind", list(AggregateMPOfMads0RTLMOC13xs2[as.logical(AggregateMPOfMads0RTLMOC13xs2[,"MatchWithInCellulo"]), c("SubcellularColocalization", "Pairs")], PairType = "All", Screen = "In Cellulo")),
                                                                do.call("cbind", list(AggregateMPOfMads0RTLMOC13xs2[as.logical(AggregateMPOfMads0RTLMOC13xs2[,"MatchWithInVitro"]), c("SubcellularColocalization", "Pairs")], PairType = "All", Screen = "In Vitro"))))

PossibleSubcellularLocalizationOverlapDatax <- do.call("rbind", list(cbind(AggregateMPOfMads0RTLMOC13xs2[as.logical(AggregateMPOfMads0RTLMOC13xs2[,"SubclassPairPresenceInCellulo"]), c("SubcellularColocalization", "Pairs", "PairType")], Screen = "In Cellulo"),
                                                                     cbind(AggregateMPOfMads0RTLMOC13xs2[as.logical(AggregateMPOfMads0RTLMOC13xs2[,"SubclassPairPresenceInVitro"]), c("SubcellularColocalization", "Pairs","PairType")], Screen = "In Vitro"),
                                                                     
                                                                     do.call("cbind", list(AggregateMPOfMads0RTLMOC13xs2[as.logical(AggregateMPOfMads0RTLMOC13xs2[,"SubclassPairPresenceInCellulo"]), c("SubcellularColocalization", "Pairs")], PairType = "All", Screen = "In Cellulo")),
                                                                     do.call("cbind", list(AggregateMPOfMads0RTLMOC13xs2[as.logical(AggregateMPOfMads0RTLMOC13xs2[,"SubclassPairPresenceInVitro"]), c("SubcellularColocalization", "Pairs")], PairType = "All", Screen = "In Vitro"))))

AllPairsInSubcellularData <- rbind(AggregateMPOfMads0RTLMOC11x[,c("PairType", "SubcellularColocalization", "Pairs")], cbind(PairType = "All", AggregateMPOfMads0RTLMOC11x[,c("SubcellularColocalization", "Pairs")]))


library(RColorBrewer)
pdf("./Output/Panel5CWithTailoredScriptAndMoreDetailedVisualizationBinSize001NoSegmentStacking18032022.pdf")

BinSize <- 0.01
MaxList2CellsColoc <- list()

StatListCellsColoc <- list()
StatList2CellsColoc <- list()

StatListCellsColocX2 <- list()
StatList2CellsColocX2 <- list()

ScreenColorsForFilling <- cbind(c("A","E"),c(brewer.pal(9,"Blues")[4], brewer.pal(9,"Oranges")[4]))
ScreenColorsForFillingPot <- cbind(c("A","E"),c(brewer.pal(9,"Blues")[7], brewer.pal(9,"Oranges")[7]))


for(SimilarityIndex in c("All", "Class", "Sub-Class", "Species")){
  
  x1a <- na.omit(as.numeric(SubcellularLocalizationOverlapOurDatax[(SubcellularLocalizationOverlapOurDatax[,"Screen"] == "In Cellulo") & (SubcellularLocalizationOverlapOurDatax[,"PairType"] == SimilarityIndex),"SubcellularColocalization"]))
  x1e <- na.omit(as.numeric(SubcellularLocalizationOverlapOurDatax[(SubcellularLocalizationOverlapOurDatax[,"Screen"] == "In Vitro") & (SubcellularLocalizationOverlapOurDatax[,"PairType"] == SimilarityIndex),"SubcellularColocalization"]))
  
  x1apot <- na.omit(as.numeric(PossibleSubcellularLocalizationOverlapDatax[(PossibleSubcellularLocalizationOverlapDatax[,"Screen"] == "In Cellulo") & (PossibleSubcellularLocalizationOverlapDatax[,"PairType"] == SimilarityIndex),"SubcellularColocalization"]))
  x1epot <- na.omit(as.numeric(PossibleSubcellularLocalizationOverlapDatax[(PossibleSubcellularLocalizationOverlapDatax[,"Screen"] == "In Vitro") & (PossibleSubcellularLocalizationOverlapDatax[,"PairType"] == SimilarityIndex),"SubcellularColocalization"]))
  
  x2 <- na.omit(as.numeric(AllPairsInSubcellularData[AllPairsInSubcellularData[,"PairType"] == SimilarityIndex, "SubcellularColocalization"]))
  
  
  MaxY <- max(c(unlist(density(x1a, from = -1, to = 1)["y"]), unlist(density(x1e, from = -1, to = 1)["y"]), unlist(density(x1apot, from = -1, to = 1)["y"]), unlist(density(x1epot, from = -1, to = 1)["y"]), unlist(density(x2, from = -1, to = 1)["y"])))
  plot(x = unlist(density(x1a, from = 0, to = 1)["x"]), y = unlist(density(x1a, from = -1, to = 1)["y"]), xlab = "Co-occurrence (sub-cellular)", ylab = "Normalized Density", ylim = c(-MaxY,MaxY), xlim = c(0,1), col = "white", main = paste0("Similarity: ", SimilarityIndex), type = "l", bty="n")
  
  polygon(x = c(min(unlist(density(x1a, from = -1, to = 1)["x"])), unlist(density(x1a, from = -1, to = 1)["x"]), max(unlist(density(x1a, from = -1, to = 1)["x"]))),
          y = c(0, unlist(density(x1a, from = -1, to = 1)["y"]), 0), col = ScreenColorsForFilling[ScreenColorsForFilling[,1] == "A",2], border = "white")
  
  polygon(x = c(min(unlist(density(x1e, from = -1, to = 1)["x"])), unlist(density(x1e, from = -1, to = 1)["x"]), max(unlist(density(x1e, from = -1, to = 1)["x"]))),
          y = c(0, -unlist(density(x1e, from = -1, to = 1)["y"]), 0), col = ScreenColorsForFilling[ScreenColorsForFilling[,1] == "E",2], border = "white")
  
  segments(x1a, 0.16, x1 = x1a, y1 = 0, col = "black")
  segments(x1e, -0.16, x1 = x1e, y1 = 0, col = "black")
  
  lines(c(0,unlist(density(x2, from = 0, to = 1)["x"]),1), c(0,unlist(density(x2, from = 0, to = 1)["y"]),0), col = "DarkGrey", lwd = 2)
  lines(c(0,unlist(density(x2, from = 0, to = 1)["x"]),1), c(0,-unlist(density(x2, from = 0, to = 1)["y"]),0), col = "DarkGrey", lwd = 2)
  
  lines(c(0,unlist(density(x1apot, from = 0, to = 1)["x"]),1), c(0,unlist(density(x1apot, from = 0, to = 1)["y"]),0), col = ScreenColorsForFillingPot[ScreenColorsForFillingPot[,1] == "A",2], lwd = 2)
  lines(c(0,unlist(density(x1epot, from = 0, to = 1)["x"]),1), c(0,-unlist(density(x1epot, from = 0, to = 1)["y"]),0), col = ScreenColorsForFillingPot[ScreenColorsForFillingPot[,1] == "E",2], lwd = 2)
  
  segments(x0 = median(x1a), x1 = median(x1a), y0 = MaxY, y1 = 0, lwd = 2)
  segments(x0 = median(x1e), x1 = median(x1e), y0 = -MaxY, y1 = 0, lwd = 2)
  
  segments(x0 = median(x1apot), x1 = median(x1apot), y0 = MaxY, y1 = 0, lwd = 2, col = ScreenColorsForFillingPot[ScreenColorsForFillingPot[,1] == "A",2])
  segments(x0 = median(x1epot), x1 = median(x1epot), y0 = -MaxY, y1 = 0, lwd = 2, col = ScreenColorsForFillingPot[ScreenColorsForFillingPot[,1] == "E",2])
  
  segments(x0 = median(x2), x1 = median(x2), y0 = MaxY, y1 = 0, lwd = 2, col = "DarkGrey")
  segments(x0 = median(x2), x1 = median(x2), y0 = -MaxY, y1 = 0, lwd = 2, col = "DarkGrey")
  
  MaxList2CellsColoc[[SimilarityIndex]][["In Cellulo"]] <- c(max(unlist(density(x1a, from = 0, to = 1)["y"])),max(unlist(density(x2, from = 0, to = 1)["y"])))
  MaxList2CellsColoc[[SimilarityIndex]][["In Vitro"]] <- c(max(unlist(density(x1e, from = 0, to = 1)["y"])),max(unlist(density(x2, from = 0, to = 1)["y"])))
  
}



# Generalize the maximum for y: all figures on the same scale

MaxY <- max(unlist(MaxList2CellsColoc)) 
MaxList2CellsColoc <- list()


for(SimilarityIndex in c("All", "Class", "Sub-Class", "Species")){
  
  x1a <- na.omit(as.numeric(SubcellularLocalizationOverlapOurDatax[(SubcellularLocalizationOverlapOurDatax[,"Screen"] == "In Cellulo") & (SubcellularLocalizationOverlapOurDatax[,"PairType"] == SimilarityIndex),"SubcellularColocalization"]))
  x1e <- na.omit(as.numeric(SubcellularLocalizationOverlapOurDatax[(SubcellularLocalizationOverlapOurDatax[,"Screen"] == "In Vitro") & (SubcellularLocalizationOverlapOurDatax[,"PairType"] == SimilarityIndex),"SubcellularColocalization"]))
  
  x1apot <- na.omit(as.numeric(PossibleSubcellularLocalizationOverlapDatax[(PossibleSubcellularLocalizationOverlapDatax[,"Screen"] == "In Cellulo") & (PossibleSubcellularLocalizationOverlapDatax[,"PairType"] == SimilarityIndex),"SubcellularColocalization"]))
  x1epot <- na.omit(as.numeric(PossibleSubcellularLocalizationOverlapDatax[(PossibleSubcellularLocalizationOverlapDatax[,"Screen"] == "In Vitro") & (PossibleSubcellularLocalizationOverlapDatax[,"PairType"] == SimilarityIndex),"SubcellularColocalization"]))
  
  x2 <- na.omit(as.numeric(AllPairsInSubcellularData[AllPairsInSubcellularData[,"PairType"] == SimilarityIndex, "SubcellularColocalization"]))
  plot(x = unlist(density(x1a, from = 0, to = 1)["x"]), y = unlist(density(x1a, from = -1, to = 1)["y"]), xlab = "Co-occurrence (sub-cellular)", ylab = "Normalized Density", ylim = c(-MaxY,MaxY), xlim = c(0,1), col = "white", main = paste0("Similarity: ", SimilarityIndex), type = "l", bty="n")
  
  polygon(x = c(min(unlist(density(x1a, from = -1, to = 1)["x"])), unlist(density(x1a, from = -1, to = 1)["x"]), max(unlist(density(x1a, from = -1, to = 1)["x"]))),
          y = c(0, unlist(density(x1a, from = -1, to = 1)["y"]), 0), col = ScreenColorsForFilling[ScreenColorsForFilling[,1] == "A",2], border = "white")
  
  polygon(x = c(min(unlist(density(x1e, from = -1, to = 1)["x"])), unlist(density(x1e, from = -1, to = 1)["x"]), max(unlist(density(x1e, from = -1, to = 1)["x"]))),
          y = c(0, -unlist(density(x1e, from = -1, to = 1)["y"]), 0), col = ScreenColorsForFilling[ScreenColorsForFilling[,1] == "E",2], border = "white")
  
  segments(x1a, 0.16, x1 = x1a, y1 = 0, col = "black")
  segments(x1e, -0.16, x1 = x1e, y1 = 0, col = "black")
  
  lines(c(0,unlist(density(x2, from = 0, to = 1)["x"]),1), c(0,unlist(density(x2, from = 0, to = 1)["y"]),0), col = "DarkGrey", lwd = 2)
  lines(c(0,unlist(density(x2, from = 0, to = 1)["x"]),1), c(0,-unlist(density(x2, from = 0, to = 1)["y"]),0), col = "DarkGrey", lwd = 2)
  
  lines(c(0,unlist(density(x1apot, from = 0, to = 1)["x"]),1), c(0,unlist(density(x1apot, from = 0, to = 1)["y"]),0), col = ScreenColorsForFillingPot[ScreenColorsForFillingPot[,1] == "A",2], lwd = 2)
  lines(c(0,unlist(density(x1epot, from = 0, to = 1)["x"]),1), c(0,-unlist(density(x1epot, from = 0, to = 1)["y"]),0), col = ScreenColorsForFillingPot[ScreenColorsForFillingPot[,1] == "E",2], lwd = 2)
  
  segments(x0 = median(x1a), x1 = median(x1a), y0 = MaxY, y1 = 0, lwd = 2)
  segments(x0 = median(x1e), x1 = median(x1e), y0 = -MaxY, y1 = 0, lwd = 2)
  
  segments(x0 = median(x1apot), x1 = median(x1apot), y0 = MaxY, y1 = 0, lwd = 2, col = ScreenColorsForFillingPot[ScreenColorsForFillingPot[,1] == "A",2])
  segments(x0 = median(x1epot), x1 = median(x1epot), y0 = -MaxY, y1 = 0, lwd = 2, col = ScreenColorsForFillingPot[ScreenColorsForFillingPot[,1] == "E",2])
  
  segments(x0 = median(x2), x1 = median(x2), y0 = MaxY, y1 = 0, lwd = 2, col = "DarkGrey")
  segments(x0 = median(x2), x1 = median(x2), y0 = -MaxY, y1 = 0, lwd = 2, col = "DarkGrey")
  
  MaxList2CellsColoc[[SimilarityIndex]][["In Cellulo"]] <- c(max(unlist(density(x1a, from = 0, to = 1)["y"])),max(unlist(density(x2, from = 0, to = 1)["y"])))
  MaxList2CellsColoc[[SimilarityIndex]][["E"]] <- c(max(unlist(density(x1e, from = 0, to = 1)["y"])),max(unlist(density(x2, from = 0, to = 1)["y"])))
  
  StatListCellsColoc[[SimilarityIndex]] <- structure(rbind(FishersExactTestBelowAbove05(SubsetX = x1a, SupersetY = c(x1apot,x1apot)), FishersExactTestBelowAbove05(SubsetX = x1e, SupersetY = c(x1epot,x1epot))), 
                                                     dimnames = list(c("In Cellulo","In Vitro"), c("NoNeg", "NoPos", "YesNeg", "YesPos", "LowConf", "HighConf", "OddsRatio", "pValue")))
  
  StatList2CellsColoc[[SimilarityIndex]] <- structure(rbind(FishersExactTestNegPosEqualBackground2(SubsetX = x1a, SupersetY = c(x1apot,x1apot)), FishersExactTestNegPosEqualBackground2(SubsetX = x1e, SupersetY = c(x1epot,x1epot))), 
                                                      dimnames = list(c("In Cellulo","In Vitro"), c("NoNeg", "NoPos", "YesNeg", "YesPos", "LowConf", "HighConf", "OddsRatio", "pValue")))
  
  StatListCellsColocX2[[SimilarityIndex]] <- structure(rbind(FishersExactTestBelowAbove05(SubsetX = x1a, SupersetY = c(x2,x2)), FishersExactTestBelowAbove05(SubsetX = x1e, SupersetY = c(x2,x2))), 
                                                       dimnames = list(c("In Cellulo","In Vitro"), c("NoNeg", "NoPos", "YesNeg", "YesPos", "LowConf", "HighConf", "OddsRatio", "pValue")))
  
  StatList2CellsColocX2[[SimilarityIndex]] <- structure(rbind(FishersExactTestNegPosEqualBackground2(SubsetX = x1a, SupersetY = c(x2,x2)), FishersExactTestNegPosEqualBackground2(SubsetX = x1e, SupersetY = c(x2,x2))), 
                                                        dimnames = list(c("In Cellulo","In Vitro"), c("NoNeg", "NoPos", "YesNeg", "YesPos", "LowConf", "HighConf", "OddsRatio", "pValue")))
  
}
dev.off()

# All of the figures in the second set of figure pages in the previous document were combined into Panel 4c and further esthetically enhanced and integrated with the other parts of Figure 4 in Adobe Illustrator.
# For the statistical significances the results of the Fisher exact tests on the difference in the medians was used (see further), and also added in Adobe Illustrator.

#### Similar figures for panels from figure 4, but with Nrd0 and with unity lines to enable internal comparisons, y-axes, and correct fitting with eachother during integration of figures in Adobe Illustrator.
#### Not the finally used version of these panels (see parts before, for the finally used versions).


#### Panel 4A Alternative With Nrd0

# After this "all" gets integrated in 4a in a similar set-up as for 4b and 4c too
CorrelationsMatrixae8_8g <- CorrelationsMatrixae8_8

CorrelationsMatrixae8_8g[,"MatchLevel"] <- "All"
CorrelationsMatrixae8_8combo <- rbind(CorrelationsMatrixae8_8g, CorrelationsMatrixae8_8)

KoeberlinCorrelationsConsensusNames2LongReducedVersion7g <- KoeberlinCorrelationsConsensusNames2LongReducedVersion7
KoeberlinCorrelationsConsensusNames2LongReducedVersion7g[,"MatchingNumber"] <- "All"

KoeberlinCorrelationsConsensusNames2LongReducedVersion7combo <- rbind(KoeberlinCorrelationsConsensusNames2LongReducedVersion7g, KoeberlinCorrelationsConsensusNames2LongReducedVersion7) 
library(beanplot)

pdf("./Output/OverlayOfTheViolinPlotsWithBeanlineMediansAttemptToScaleInInverseOrderNrd0CorrectionDone11092020RemakeWithAllIncludedAndOverplotting18032022.pdf")
TestBeanForSpecial218032022 <- beanplot(as.numeric(CorrelationsMatrixae8_8combo[,"Correlation"]) ~ factor(CorrelationsMatrixae8_8combo[,"Screen"], levels = c("E","A")) * factor(CorrelationsMatrixae8_8combo[,"MatchLevel"], levels = c(as.character(2:0), "All")), ll = 0.04,
                                        
                                        main = "Lipid hierarchy preference by correlation", side = "both", xlab="Correlation",
                                        col = list(brewer.pal(9,"Oranges")[4], c(brewer.pal(9,"Blues")[4], "black")), method = "overplot",
                                        
                                        axes=F, bw = "nrd0",
                                        horizontal = TRUE,  boxwex = 1,
                                        
                                        cutmin = -1,
                                        cutmax = 1, log = "",
                                        
                                        border = FALSE, 
                                        overalline = "median", beanlines = "median")


TestBeanForGeneral2standard18032022 <- beanplot(as.numeric(KoeberlinCorrelationsConsensusNames2LongReducedVersion7combo[,"correlation"]) ~ factor(KoeberlinCorrelationsConsensusNames2LongReducedVersion7combo[,"MatchingNumber"], levels = c(as.character(2:0),"All")), ll = 0,
                                                
                                                alpha = 0.5, bw = "nrd0",
                                                col = c(NA,"black", "black","darkgrey"),
                                                
                                                axes=F, log = "",
                                                horizontal = TRUE,  boxwex = 1,
                                                
                                                cutmin = -1,
                                                cutmax = 1,
                                                
                                                border = "DarkGrey", wd = TestBeanForSpecial218032022$wd,
                                                overalline = "median", beanlines = "median",
                                                
                                                add = TRUE)


axis(side = 2, at = 4:1, labels = c("all", "class", "subclass", "species"))
axis(1)

dev.off()
# Previous figure is Nrd0-version for panel 4A that was not finally used in the article, but without an unit line for the beanplots as reference for further steps

# Steps to include an unit line for all the beanplots
UnitBean <- do.call("rbind.data.frame", lapply(c(as.character(2:0),"All"), function(x){cbind.data.frame(as.numeric(seq(-1,1, by = 0.01)), x)}))

colnames(UnitBean) <- c("correlation", "MatchingNumber")


pdf("./Output/OverlayOfTheViolinPlotsWithBeanlineMediansAttemptToScaleInInverseOrderNrd0CorrectionDone11092020RemakeWithAllIncludedAndOverplottingWithReferenceUnity21032022.pdf")
TestBeanForSpecial218032022 <- beanplot(as.numeric(CorrelationsMatrixae8_8combo[,"Correlation"]) ~ factor(CorrelationsMatrixae8_8combo[,"Screen"], levels = c("E","A")) * factor(CorrelationsMatrixae8_8combo[,"MatchLevel"], levels = c(as.character(2:0), "All")), ll = 0.04,
                                        
                                        main = "Lipid hierarchy preference by correlation", side = "both", xlab="Correlation",
                                        col = list(brewer.pal(9,"Oranges")[4], c(brewer.pal(9,"Blues")[4], "black")), method = "overplot",
                                        
                                        axes=F, bw = "nrd0", 
                                        horizontal = TRUE,  boxwex = 1,
                                        
                                        cutmin = -1,
                                        cutmax = 1, log = "",
                                        
                                        border = FALSE, 
                                        overalline = "median", beanlines = "median")


TestBeanForGeneral2standard18032022 <- beanplot(as.numeric(KoeberlinCorrelationsConsensusNames2LongReducedVersion7combo[,"correlation"]) ~ factor(KoeberlinCorrelationsConsensusNames2LongReducedVersion7combo[,"MatchingNumber"], levels = c(as.character(2:0),"All")), ll = 0,
                                                
                                                alpha = 0.5, bw = "nrd0", 
                                                col = c(NA,"black", "black","darkgrey"),
                                                
                                                axes=F, log = "",
                                                horizontal = TRUE,  boxwex = 1,
                                                
                                                cutmin = -1,
                                                cutmax = 1,
                                                
                                                border = "DarkGrey", wd = TestBeanForSpecial218032022$wd,
                                                overalline = "median", beanlines = "median",
                                                
                                                add = TRUE)


axis(side = 2, at = 4:1, labels = c("all", "class", "subclass", "species"))
axis(1)


beanplot(as.numeric(UnitBean[,"correlation"]) ~ factor(UnitBean[,"MatchingNumber"], levels = c(as.character(2:0),"All")), ll = 0,
         
         alpha = 0.5, bw = "nrd0", 
         col = c(NA, "white", "white", NA), # updated from an NA-vector
         
         axes=F, log = "",
         horizontal = TRUE,  boxwex = 1,
         
         cutmin = -1,
         cutmax = 1,
         
         border = "green", wd = TestBeanForSpecial218032022$wd,
         overalline = "median", beanlines = "median",
         
         add = TRUE)
dev.off()

# Previous figure is the Nrd0 corrected version for panel 4A that was not finally used in the article, with an unit line for the beanplots as reference for further steps
rm(TestBeanForGeneral2standard18032022)


#### Panel 5C Alternative With Nrd0 (and unit line)

UnitBeanSubcellLoc <- do.call("rbind.data.frame", lapply(c("Species", "Sub-Class", "Class", "All"), function(x){cbind.data.frame(as.numeric(seq(0,1, by = 0.01)), x)}))
colnames(UnitBeanSubcellLoc) <- c("SubcellularColocalization", "PairType")

pdf("./Output/Panel5CWithGreenUnityDistributionAndBwWithNrd0Wd01422032022.pdf")
beanplot(as.numeric(as.character(SubcellularLocalizationOverlapOurDatax[,"SubcellularColocalization"])) ~ factor(SubcellularLocalizationOverlapOurDatax[,"Screen"], levels = c("In Vitro", "In Cellulo")) * factor(SubcellularLocalizationOverlapOurDatax[,"PairType"], levels = c("Species", "Sub-Class", "Class", "All")),
         
         main = "Co-transport: lipid hierarchy preference by co-occurrence in cells", side = "both", xlab="Sub-cellular co-occurrence", ll = 0.04, bw = "nrd0", wd = 0.14, # wd = 0.23,
         col = list(brewer.pal(9,"Oranges")[4], c(brewer.pal(9,"Blues")[4], "black")), method = "overplot", 
         
         axes=F,
         horizontal = TRUE,
         
         cutmin = 0,
         cutmax = 1,
         
         border = FALSE,
         overalline = "median")


beanplot(as.numeric(as.character(AllPairsInSubcellularData[,"SubcellularColocalization"])) ~ factor(AllPairsInSubcellularData[,"PairType"], levels = c("Species", "Sub-Class", "Class", "All")),
         
         alpha = 0.5, wd = 0.14, 
         col = c(NA, "black", "black", "DarkGrey"), side = "no",
         
         axes=F, ll = 0,
         horizontal = TRUE, bw = "nrd0",
         
         cutmin = 0,
         cutmax = 1,
         
         border = c("DarkGrey", "DarkGrey"),
         overalline = "median",
         
         add = TRUE,
         what = c(FALSE, TRUE, TRUE, TRUE))


beanplot(as.numeric(as.character(PossibleSubcellularLocalizationOverlapDatax[,"SubcellularColocalization"])) ~ factor(PossibleSubcellularLocalizationOverlapDatax[,"Screen"], levels = c("In Vitro", "In Cellulo")) * factor(PossibleSubcellularLocalizationOverlapDatax[,"PairType"], levels = c("Species", "Sub-Class", "Class", "All")),
         
         alpha = 0.5, wd = 0.14, 
         col = c(NA, "black", "black", "Pink"), side = "both",
         
         axes=F, ll = 0,
         horizontal = TRUE, bw = "nrd0",
         
         cutmin = 0,
         cutmax = 1,
         
         border = c(brewer.pal(9,"Oranges")[7], c(brewer.pal(9,"Blues")[7])),
         overalline = "median",
         
         add = TRUE,
         what = c(FALSE, TRUE, TRUE, TRUE))


beanplot(as.numeric(UnitBeanSubcellLoc[,"SubcellularColocalization"]) ~ factor(UnitBeanSubcellLoc[,"PairType"], levels = c("Species", "Sub-Class", "Class", "All")), ll = 0,
         
         bw = "nrd0", 
         col = c(NA, "white", "white", NA), # updated from an NA-vector
         
         axes=F, log = "",
         horizontal = TRUE,  boxwex = 1,
         
         cutmin = 0,
         cutmax = 1,
         
         border = "green", wd = 0.14, 
         overalline = "median", beanlines = "median",
         
         add = TRUE)


axis(side = 2, at = 1:4, labels = c("species", "subclass", "class", "all"))
axis(1)

dev.off()
# Previous figure is the Nrd0 corrected version for panel 4C that was not finally used in the article, with an unit line for the beanplots as reference for further steps


#### Panel 4B Alternative With Nrd0 (and unit line)

UnitBeanTissueLoc <- do.call("rbind.data.frame", lapply(c("Species", "Sub-Class", "Class", "All"), function(x){cbind.data.frame(as.numeric(seq(0,1, by = 0.01)), x)}))
colnames(UnitBeanTissueLoc) <- c("Cooccurrence", "LipidPairMatchingType")

pdf("./Output/Panel5BWithGreenUnityDistributionWithNrd0AndWd02322032022b.pdf")
beanplot(as.numeric(as.character(CooccurrenceOfChemicalPairsListExpandedWithAll2x[,"Cooccurrence"])) ~ factor(CooccurrenceOfChemicalPairsListExpandedWithAll2x[,"Screen"], levels = c("E", "A")) * factor(CooccurrenceOfChemicalPairsListExpandedWithAll2x[,"LipidPairMatchingType"], levels = c(as.character(2:0), "All")), 
         
         main = "Co-transport: lipid hierarchy preference by co-occurrence in tissues", side = "both", xlab="Co-occurrence", ll = 0.04, wd = 0.23,
         col = list(brewer.pal(9,"Oranges")[4], c(brewer.pal(9,"Blues")[4], "black")), method = "overplot",
         
         axes=F, bw = "nrd0", 
         horizontal = TRUE,
         
         cutmin = 0,
         cutmax = 1,
         
         border = FALSE,
         overalline = "median")


beanplot(as.numeric(as.character(PossibleCooccurrenceOfChemicalPairs2x[,"Cooccurrence"])) ~ factor(PossibleCooccurrenceOfChemicalPairs2x[,"Screen"], levels = c("E", "A")) * factor(PossibleCooccurrenceOfChemicalPairs2x[,"LipidPairMatchingType"], levels = c(as.character(2:0), "All")), 
         
         alpha = 0.5, wd = 0.23,
         col = c(NA, "black", "black", "DarkGrey"), side = "both",
         
         axes=F, ll = 0,
         horizontal = TRUE, bw = "nrd0", 
         
         cutmin = 0,
         cutmax = 1,
         
         border = c(brewer.pal(9,"Oranges")[7], c(brewer.pal(9,"Blues")[7])),
         overalline = "median",
         
         add = TRUE,
         what = c(FALSE, TRUE, TRUE, TRUE))


beanplot(as.numeric(as.character(LipidSubsetMMeans_4Subset[,"Cooccurrence"])) ~ factor(LipidSubsetMMeans_4Subset[,"LipidPairMatchingType"], levels = c(as.character(2:0), "All")), 
         
         alpha = 0.5, wd = 0.23,
         col = c(NA, "black", "black", "DarkGrey"), side = "no",
         
         axes=F, ll = 0,
         horizontal = TRUE, bw = "nrd0", 
         
         cutmin = 0,
         cutmax = 1,
         
         border = c("DarkGrey", "DarkGrey"),
         overalline = "median",
         
         add = TRUE,
         what = c(FALSE, TRUE, TRUE, TRUE))


beanplot(as.numeric(UnitBeanSubcellLoc[,"SubcellularColocalization"]) ~ factor(UnitBeanSubcellLoc[,"PairType"], levels = c("Species", "Sub-Class", "Class", "All")), ll = 0,
         
         
         col = c(NA, "white", "white", NA), # updated from an NA-vector
         
         axes=F, log = "",
         horizontal = TRUE,  boxwex = 1,
         
         cutmin = 0, bw = "nrd0", 
         cutmax = 1,
         
         border = "green", wd = 0.23, 
         overalline = "median", beanlines = "median",
         
         add = TRUE)


axis(side = 2, at = 1:4, labels = c("species", "subclass", "class", "all"))
axis(1)

dev.off()
# Previous figure is the Nrd0 corrected version for panel 4C that was not finally used in the article, with an unit line for the beanplots as reference for further steps


#### Results of Fisher exact tests

#. FisherExactTestResultsStandardCorrectedWithAll18082022.csv # Not finally used version because not all-round applicable & least strict
#. FisherExactTestResultsWithEqualizedBackgroundCorrectedWithAll18082022.csv # Finally used version of comparisons: applicable to all & even more strict

write.table(structure(do.call("rbind", StatList), dimnames = list(paste0(rep(c("A","E"),3), unlist(lapply(c("All",as.character(0:2)), function(x){rep(x, 2)}))), c("NoNeg", "NoPos", "YesNeg", "YesPos", "LowConf", "HighConf", "OddsRatio", "pValue"))),
            file="./Output/FisherExactTestResultsStandardCorrectedWithAll18082022.csv", sep="\t", row.names = TRUE, quote = FALSE)

write.table(structure(do.call("rbind", StatList2), dimnames = list(paste0(rep(c("A","E"),3), unlist(lapply(c("All",as.character(0:2)), function(x){rep(x, 2)}))), c("NoNeg", "NoPos", "YesNeg", "YesPos", "LowConf", "HighConf", "OddsRatio", "pValue"))),
            file="./Output/FisherExactTestResultsWithEqualizedBackgroundCorrectedWithAll18082022.csv", sep="\t", row.names = TRUE, quote = FALSE)

#. FisherExactTestResultsStandardForTheCooccurrences09102020.csv # Not finally used version because not all-round applicable & least strict
#. FisherExactTestResultsWithEqualizedBackgroundForTheCooccurrencesCorrected18082022.csv # Finally used version of comparisons: applicable to all & even more strict

StatListForCooccurrences <- list()
StatList2ForCooccurrences <- list()


SimilarityIndex <- "All"

x1a <- RemoveNAsFromVector(as.numeric(CooccurrenceOfChemicalPairsListExpandedWithAll2x[(CooccurrenceOfChemicalPairsListExpandedWithAll2x[,"Screen"] == "A") & (CooccurrenceOfChemicalPairsListExpandedWithAll2x[,"LipidPairMatchingType"] == SimilarityIndex),"Cooccurrence"]))
x1e <- RemoveNAsFromVector(as.numeric(CooccurrenceOfChemicalPairsListExpandedWithAll2x[(CooccurrenceOfChemicalPairsListExpandedWithAll2x[,"Screen"] == "E") & (CooccurrenceOfChemicalPairsListExpandedWithAll2x[,"LipidPairMatchingType"] == SimilarityIndex),"Cooccurrence"]))

pdForx1a <- RemoveNAsFromVector(as.numeric(as.character(PossibleCooccurrenceOfChemicalPairs2x[(PossibleCooccurrenceOfChemicalPairs2x[,"Screen"] == "A") & (PossibleCooccurrenceOfChemicalPairs2x[,"LipidPairMatchingType"] == SimilarityIndex),"Cooccurrence"])))
pdForx1e <- RemoveNAsFromVector(as.numeric(as.character(PossibleCooccurrenceOfChemicalPairs2x[(PossibleCooccurrenceOfChemicalPairs2x[,"Screen"] == "E") & (PossibleCooccurrenceOfChemicalPairs2x[,"LipidPairMatchingType"] == SimilarityIndex),"Cooccurrence"])))

StatListForCooccurrences[[SimilarityIndex]] <- structure(rbind(FishersExactTestBelowAbove05(SubsetX = x1a, SupersetY = c(pdForx1a,pdForx1a)), FishersExactTestBelowAbove05(SubsetX = x1e, SupersetY = c(pdForx1e,pdForx1e))), 
                                                         dimnames = list(c("A","E"), c("NoNeg", "NoPos", "YesNeg", "YesPos", "LowConf", "HighConf", "OddsRatio", "pValue")))

StatList2ForCooccurrences[[SimilarityIndex]] <- structure(rbind(FishersExactTestNegPosEqualBackground2(SubsetX = x1a, SupersetY = c(pdForx1a,pdForx1a)), FishersExactTestNegPosEqualBackground2(SubsetX = x1e, SupersetY = c(pdForx1e,pdForx1e))), 
                                                          dimnames = list(c("A","E"), c("NoNeg", "NoPos", "YesNeg", "YesPos", "LowConf", "HighConf", "OddsRatio", "pValue")))

Stat1FullAllForCooccurrences <- structure(rbind(FishersExactTestBelowAbove05(SubsetX = x1a, SupersetY = c(as.numeric(as.character(LipidSubsetMMeans_4Subset[,1])),as.numeric(as.character(LipidSubsetMMeans_4Subset[,1])))), FishersExactTestBelowAbove05(SubsetX = x1e, SupersetY = c(as.numeric(as.character(LipidSubsetMMeans_4Subset[,1])),as.numeric(as.character(LipidSubsetMMeans_4Subset[,1]))))), 
                                          dimnames = list(c("A","E"), c("NoNeg", "NoPos", "YesNeg", "YesPos", "LowConf", "HighConf", "OddsRatio", "pValue")))

Stat2FullAllForCooccurrencesCorrectedVersion <- structure(rbind(FishersExactTestNegPosEqualBackground2(SubsetX = x1a, SupersetY = c(as.numeric(as.character(LipidSubsetMMeans_4Subset[,1])),as.numeric(as.character(LipidSubsetMMeans_4Subset[,1])))), FishersExactTestNegPosEqualBackground2(SubsetX = x1e, SupersetY = c(as.numeric(as.character(LipidSubsetMMeans_4Subset[,1])),as.numeric(as.character(LipidSubsetMMeans_4Subset[,1]))))), 
                                                          dimnames = list(c("A","E"), c("NoNeg", "NoPos", "YesNeg", "YesPos", "LowConf", "HighConf", "OddsRatio", "pValue")))


for(SimilarityIndex in as.character(0:2)){
  
  x1a <- RemoveNAsFromVector(as.numeric(CooccurrenceOfChemicalPairsListExpandedWithAll2x[(CooccurrenceOfChemicalPairsListExpandedWithAll2x[,"Screen"] == "A") & (CooccurrenceOfChemicalPairsListExpandedWithAll2x[,"LipidPairMatchingType"] == SimilarityIndex),"Cooccurrence"]))
  x1e <- RemoveNAsFromVector(as.numeric(CooccurrenceOfChemicalPairsListExpandedWithAll2x[(CooccurrenceOfChemicalPairsListExpandedWithAll2x[,"Screen"] == "E") & (CooccurrenceOfChemicalPairsListExpandedWithAll2x[,"LipidPairMatchingType"] == SimilarityIndex),"Cooccurrence"]))
  
  pdForx1a <- RemoveNAsFromVector(as.numeric(as.character(PossibleCooccurrenceOfChemicalPairs2x[(PossibleCooccurrenceOfChemicalPairs2x[,"Screen"] == "A") & (PossibleCooccurrenceOfChemicalPairs2x[,"LipidPairMatchingType"] == SimilarityIndex),"Cooccurrence"])))
  pdForx1e <- RemoveNAsFromVector(as.numeric(as.character(PossibleCooccurrenceOfChemicalPairs2x[(PossibleCooccurrenceOfChemicalPairs2x[,"Screen"] == "E") & (PossibleCooccurrenceOfChemicalPairs2x[,"LipidPairMatchingType"] == SimilarityIndex),"Cooccurrence"])))
  
  StatListForCooccurrences[[SimilarityIndex]] <- structure(rbind(FishersExactTestBelowAbove05(SubsetX = x1a, SupersetY = c(pdForx1a,pdForx1a)), FishersExactTestBelowAbove05(SubsetX = x1e, SupersetY = c(pdForx1e,pdForx1e))), 
                                                           dimnames = list(c("A","E"), c("NoNeg", "NoPos", "YesNeg", "YesPos", "LowConf", "HighConf", "OddsRatio", "pValue")))
  
  StatList2ForCooccurrences[[SimilarityIndex]] <- structure(rbind(FishersExactTestNegPosEqualBackground2(SubsetX = x1a, SupersetY = c(pdForx1a,pdForx1a)), FishersExactTestNegPosEqualBackground2(SubsetX = x1e, SupersetY = c(pdForx1e,pdForx1e))), 
                                                            dimnames = list(c("A","E"), c("NoNeg", "NoPos", "YesNeg", "YesPos", "LowConf", "HighConf", "OddsRatio", "pValue")))
  
}


StatListForCooccurrences[["AllAll"]] <- Stat1FullAllForCooccurrences


StatList2ForCooccurrencesCorrected <- StatList2ForCooccurrences
StatList2ForCooccurrencesCorrected[["AllAll"]] <- Stat2FullAllForCooccurrencesCorrectedVersion

write.table(structure(do.call("rbind", StatListForCooccurrences[c("AllAll", "All", as.character(0:2))]), dimnames = list(paste0(rep(c("A","E"),5), unlist(lapply(c("AllAll", "All", as.character(0:2)), function(x){rep(x, 2)}))), c("NoNeg", "NoPos", "YesNeg", "YesPos", "LowConf", "HighConf", "OddsRatio", "pValue"))),
            file="./Output/FisherExactTestResultsStandardForTheCooccurrences09102020.csv", sep="\t", row.names = TRUE, quote = FALSE)

write.table(structure(do.call("rbind", StatList2ForCooccurrencesCorrected[c("AllAll", "All", as.character(0:2))]), dimnames = list(paste0(rep(c("A","E"),5), unlist(lapply(c("AllAll", "All", as.character(0:2)), function(x){rep(x, 2)}))), c("NoNeg", "NoPos", "YesNeg", "YesPos", "LowConf", "HighConf", "OddsRatio", "pValue"))),
            file="./Output/FisherExactTestResultsWithEqualizedBackgroundForTheCooccurrencesCorrected18082022.csv", sep="\t", row.names = TRUE, quote = FALSE)

# Subcellular Fisher Exact Tests # Finally used version of comparisons: only the median versions, for full reference (grey lines in graphs) and for what possible (colored lines in graphs) and not the others because applicable to all & even more strict, and double references used for comparible size in comparisons.
# They are all named ...21022022.csv . 

StatListForSubcellularCooccurrences <- list()
StatList2ForSubcellularCooccurrences <- list()

StatListForSubcellularCooccurrencesWithSingleRef <- list()
StatList2ForSubcellularCooccurrencesWithSingleRef <- list()

StatListForSubcellularCooccurrencesVersusAllShotgun <- list()
StatList2ForSubcellularCooccurrencesVersusAllShotgun <- list()

StatListForSubcellularCooccurrencesWithSingleRefVersusAllShotgun <- list()
StatList2ForSubcellularCooccurrencesWithSingleRefVersusAllShotgun <- list()

SimilarityIndex <- unique(SubcellularLocalizationOverlapOurDatax[,"PairType"])[1]
for(SimilarityIndex in as.character(unique(SubcellularLocalizationOverlapOurDatax[,"PairType"]))){
  
  x1a <- RemoveNAsFromVector(as.numeric(SubcellularLocalizationOverlapOurDatax[(SubcellularLocalizationOverlapOurDatax[,"Screen"] == "In Cellulo") & (SubcellularLocalizationOverlapOurDatax[,"PairType"] == SimilarityIndex),"SubcellularColocalization"]))
  x1e <- RemoveNAsFromVector(as.numeric(SubcellularLocalizationOverlapOurDatax[(SubcellularLocalizationOverlapOurDatax[,"Screen"] == "In Vitro") & (SubcellularLocalizationOverlapOurDatax[,"PairType"] == SimilarityIndex),"SubcellularColocalization"]))
  
  pdForx1a <- RemoveNAsFromVector(as.numeric(PossibleSubcellularLocalizationOverlapDatax[(PossibleSubcellularLocalizationOverlapDatax[,"Screen"] == "In Cellulo") & (PossibleSubcellularLocalizationOverlapDatax[,"PairType"] == SimilarityIndex),"SubcellularColocalization"]))
  pdForx1e <- RemoveNAsFromVector(as.numeric(PossibleSubcellularLocalizationOverlapDatax[(PossibleSubcellularLocalizationOverlapDatax[,"Screen"] == "In Vitro") & (PossibleSubcellularLocalizationOverlapDatax[,"PairType"] == SimilarityIndex),"SubcellularColocalization"]))
  
  allx1 <- RemoveNAsFromVector(as.numeric(AllPairsInSubcellularData[(AllPairsInSubcellularData[,"PairType"] == SimilarityIndex),"SubcellularColocalization"]))
  
  
  StatListForSubcellularCooccurrences[[SimilarityIndex]] <- structure(rbind(FishersExactTestBelowAbove05(SubsetX = x1a, SupersetY = c(pdForx1a,pdForx1a)), FishersExactTestBelowAbove05(SubsetX = x1e, SupersetY = c(pdForx1e,pdForx1e))), 
                                                                      dimnames = list(c("A","E"), c("NoNeg", "NoPos", "YesNeg", "YesPos", "LowConf", "HighConf", "OddsRatio", "pValue")))
  
  StatList2ForSubcellularCooccurrences[[SimilarityIndex]] <- structure(rbind(FishersExactTestNegPosEqualBackground2(SubsetX = x1a, SupersetY = c(pdForx1a,pdForx1a)), FishersExactTestNegPosEqualBackground2(SubsetX = x1e, SupersetY = c(pdForx1e,pdForx1e))), 
                                                                       dimnames = list(c("A","E"), c("NoNeg", "NoPos", "YesNeg", "YesPos", "LowConf", "HighConf", "OddsRatio", "pValue")))
  
  StatListForSubcellularCooccurrencesWithSingleRef[[SimilarityIndex]] <- structure(rbind(FishersExactTestBelowAbove05(SubsetX = x1a, SupersetY = pdForx1a), FishersExactTestBelowAbove05(SubsetX = x1e, SupersetY = pdForx1e)), 
                                                                                   dimnames = list(c("A","E"), c("NoNeg", "NoPos", "YesNeg", "YesPos", "LowConf", "HighConf", "OddsRatio", "pValue")))
  
  StatList2ForSubcellularCooccurrencesWithSingleRef[[SimilarityIndex]] <- structure(rbind(FishersExactTestNegPosEqualBackground2(SubsetX = x1a, SupersetY = pdForx1a), FishersExactTestNegPosEqualBackground2(SubsetX = x1e, SupersetY = pdForx1e)), 
                                                                                    dimnames = list(c("A","E"), c("NoNeg", "NoPos", "YesNeg", "YesPos", "LowConf", "HighConf", "OddsRatio", "pValue")))
  
  StatListForSubcellularCooccurrencesVersusAllShotgun[[SimilarityIndex]] <- structure(rbind(FishersExactTestBelowAbove05(SubsetX = x1a, SupersetY = c(allx1,allx1)), FishersExactTestBelowAbove05(SubsetX = x1e, SupersetY = c(allx1,allx1))), 
                                                                                      dimnames = list(c("A","E"), c("NoNeg", "NoPos", "YesNeg", "YesPos", "LowConf", "HighConf", "OddsRatio", "pValue")))
  
  StatList2ForSubcellularCooccurrencesVersusAllShotgun[[SimilarityIndex]] <- structure(rbind(FishersExactTestNegPosEqualBackground2(SubsetX = x1a, SupersetY = c(allx1,allx1)), FishersExactTestNegPosEqualBackground2(SubsetX = x1e, SupersetY = c(allx1,allx1))), 
                                                                                       dimnames = list(c("A","E"), c("NoNeg", "NoPos", "YesNeg", "YesPos", "LowConf", "HighConf", "OddsRatio", "pValue")))
  
  StatListForSubcellularCooccurrencesWithSingleRefVersusAllShotgun[[SimilarityIndex]] <- structure(rbind(FishersExactTestBelowAbove05(SubsetX = x1a, SupersetY = allx1), FishersExactTestBelowAbove05(SubsetX = x1e, SupersetY = allx1)), 
                                                                                                   dimnames = list(c("A","E"), c("NoNeg", "NoPos", "YesNeg", "YesPos", "LowConf", "HighConf", "OddsRatio", "pValue")))
  
  StatList2ForSubcellularCooccurrencesWithSingleRefVersusAllShotgun[[SimilarityIndex]] <- structure(rbind(FishersExactTestNegPosEqualBackground2(SubsetX = x1a, SupersetY = allx1), FishersExactTestNegPosEqualBackground2(SubsetX = x1e, SupersetY = allx1)), 
                                                                                                    dimnames = list(c("A","E"), c("NoNeg", "NoPos", "YesNeg", "YesPos", "LowConf", "HighConf", "OddsRatio", "pValue")))
  
}


FisherTestsForSubcellularColocalizations <- list(FisherExactTestForSubcellularColocalizationsAboveBelow05WithDoubleRefForWhatPossible = StatListForSubcellularCooccurrences,
                                                 FisherExactTestForSubcellularColocalizationsAboveBelowRefMedianWithDoubleRefForWhatPossible = StatList2ForSubcellularCooccurrences, 
                                                 
                                                 FisherExactTestForSubcellularColocalizationsAboveBelow05WithSingleRefForWhatPossible = StatListForSubcellularCooccurrencesWithSingleRef, 
                                                 FisherExactTestForSubcellularColocalizationsAboveBelowRefMedianWithSingleRefForWhatPossible = StatList2ForSubcellularCooccurrencesWithSingleRef,
                                                 
                                                 FisherExactTestForSubcellularColocalizationsAboveBelow05WithDoubleRefForFullRef = StatListForSubcellularCooccurrencesVersusAllShotgun, 
                                                 FisherExactTestForSubcellularColocalizationsAboveBelowRefMedianWithDoubleRefForFullRef = StatList2ForSubcellularCooccurrencesVersusAllShotgun, 
                                                 
                                                 FisherExactTestForSubcellularColocalizationsAboveBelow05WithSingleRefForFullRef = StatListForSubcellularCooccurrencesWithSingleRefVersusAllShotgun, 
                                                 FisherExactTestForSubcellularColocalizationsAboveBelowRefMedianWithSingleRefForFullRef = StatList2ForSubcellularCooccurrencesWithSingleRefVersusAllShotgun)


for(i in 1:length(FisherTestsForSubcellularColocalizations)){
  
  write.table(structure(do.call("rbind", FisherTestsForSubcellularColocalizations[[i]]), dimnames = list(paste0(rep(c("A","E"),4),rep(names(FisherTestsForSubcellularColocalizations[[i]]), each = 2)), colnames(FisherTestsForSubcellularColocalizations[[i]][[1]]))),
              file= paste0("./Output/",names(FisherTestsForSubcellularColocalizations)[i],"21022022.csv"), quote = FALSE)
  
}
# The previous table forms the basis for Extended Data Table 9B


######## Figure 6

#### Fig.6a
#. BarplotHeatmapTestSphingolipidTransporters06092021WithSidesAdded2CorrecteddSMZeroToNA.pdf # (R-basis of final figure.)

#. BarplotHeatmapTestSphingolipidTransporters06092021WithSidesAdded2CorrecteddSMZeroToNA10CorrectMetabolism5MoreAdvancedVersions4WithHighlighLinesWithStrongerCrosses.pdf (#) 
# (Last page is further intermediate figure after integration in Adobe Illustrator: Final figure has LCN1 information removed, legends simplified, CERS information added to the top, and highlights of most important part focusses.)

# LTPLipidConnectionsDataSubset and LTPLipidConnectionsDataAggregatedx creation does not need repeating: eliminated here, but updated for x-versionLTPLipidConnectionsDataSubset <- cbind(LTPLipidConnectionsDataSet[,c("LTPProtein", "LikelySubclass", "Screen", "Intensity")], CarbonChain = paste(LTPLipidConnectionsDataSet[,"TotalCarbonChainLength"],LTPLipidConnectionsDataSet[,"TotalCarbonChainUnsaturations"], sep = ":"))
# Change to widen5 changes the order of the rows, so I set up a new series of variables here with x at the end. PC come e.g. below PC-O.


SphingolipidSubclasses <- c("d*Cer", "dCer", "DHCer", "DHOH*Cer", "tCer", "d*CerP", "d*SM", "DHSM", "t*SM", "d*HexCer", "t*HexCer", "d*SHexCer", "t*Hex2Cer")

# Not exactly identical but has everything similar inside.
LTPLipidConnectionsDataSphingolipidsx <- LTPLipidConnectionsDataAggregatedx[LTPLipidConnectionsDataAggregatedx$LikelySubclass %in% SphingolipidSubclasses,]

LTPLipidConnectionsDataSphingolipids2x <- LTPLipidConnectionsDataSphingolipidsx[,c(rep(TRUE,3), colSums(LTPLipidConnectionsDataSphingolipidsx[,4:101]) != 0)]


# Converts some NAs to Zeros, but is equivalent otherwise after the widen5 introduction
CellularSphingolipidDistributionInfox <- Col1ToRowNames(widen5(inputdf = do.call("rbind", lapply(SphingolipidSubclasses[SphingolipidSubclasses %in% names(LipidSubclassesAddedToBackground170620214)], function(x){
  
  data.frame(CarbonChain = mgsub::mgsub(sapply(strsplit(rownames(LipidSubclassesAddedToBackground170620214[[x]]), "\\("), "[[", 2), c(")", "d", "\\*"), rep("",3)),
             Cellular = LipidSubclassesAddedToBackground170620214[[x]][,"Cellular"],
             
             Lipid = x)})),
  ColumnsLong = "Lipid", ColumnWide = "CarbonChain", ColumnValue = "Cellular", AggregatingFunction = sum, FunctionOutputValueType = double(1), SortRows = FALSE))

CellularSphingolipidDistributionInfo2x <- CellularSphingolipidDistributionInfox[rowSums(CellularSphingolipidDistributionInfox, na.rm = TRUE) > 0,
                                                                                colSums(CellularSphingolipidDistributionInfox, na.rm = TRUE) > 0]

CellularSphingolipidDistributionInfo4x <- do.call("cbind", list(LTPProtein = "Cellular", 
                                                                LikelySubclass = rownames(CellularSphingolipidDistributionInfo2x),
                                                                
                                                                Screen = "Cellular",
                                                                CellularSphingolipidDistributionInfo2x))

LTPLipidConnectionsDataSphingolipidsWithCellularAddedxx <- as.data.frame(t(Col1ToRowNames(merge(t(LTPLipidConnectionsDataSphingolipids2x), t(CellularSphingolipidDistributionInfo4x), by = 0, all = TRUE))))
LTPLipidConnectionsDataSphingolipidsWithCellularMatrixxx <- as.matrix(LTPLipidConnectionsDataSphingolipidsWithCellularAddedxx[,1:22])

mode(LTPLipidConnectionsDataSphingolipidsWithCellularMatrixxx) <- "numeric"
LTPLipidConnectionsDataSphingolipidsWithCellularMatrixxx[is.na(LTPLipidConnectionsDataSphingolipidsWithCellularMatrixxx)] <- 0

LTPLipidConnectionsDataSphingolipidsWithCellularMatrixToRowMaxxx <- do.call("rbind", lapply(1:dim(LTPLipidConnectionsDataSphingolipidsWithCellularMatrixxx)[1], 
                                                                                            function(x){LTPLipidConnectionsDataSphingolipidsWithCellularMatrixxx[x,]*100/max(LTPLipidConnectionsDataSphingolipidsWithCellularMatrixxx[x,], na.rm = TRUE)}))

rownames(LTPLipidConnectionsDataSphingolipidsWithCellularMatrixToRowMaxxx) <- apply(LTPLipidConnectionsDataSphingolipidsWithCellularAddedxx[,23:25], 1 , paste , collapse = "_" )
LTPLipidConnectionsDataSphingolipidsWithCellularMatrixToRowMax2xx <- LTPLipidConnectionsDataSphingolipidsWithCellularMatrixToRowMaxxx[c(19,14:18,20,7,21,9,8,11,10,12,13,22,2,1,6,5,3,4),]

OriginalRowsSphingoTransxx <- do.call("rbind", strsplit(rownames(LTPLipidConnectionsDataSphingolipidsWithCellularMatrixToRowMax2xx), split = "_"))
NewRowsSphingoTransxx <- unique(OriginalRowsSphingoTransxx[OriginalRowsSphingoTransxx[,3] != "Cellular", 1:2])


EmptyLayerxx <- matrix(NA, 
                       
                       nrow = dim(NewRowsSphingoTransxx)[1], 
                       ncol = dim(LTPLipidConnectionsDataSphingolipidsWithCellularMatrixToRowMax2xx)[2], 
                       
                       dimnames = list(paste(NewRowsSphingoTransxx[,1],NewRowsSphingoTransxx[,2],sep = "_"),
                                       colnames(LTPLipidConnectionsDataSphingolipidsWithCellularMatrixToRowMax2xx)))

NewSphingoTransLayersxx <- list(InCellulo = EmptyLayerxx,
                                InVitro = EmptyLayerxx,
                                
                                Cellular = EmptyLayerxx)


for(j in which(OriginalRowsSphingoTransxx[,3] == "in vivo")){
  for(i in which(paste0(NewRowsSphingoTransxx[,1],NewRowsSphingoTransxx[,2]) == paste0(OriginalRowsSphingoTransxx[j,1],OriginalRowsSphingoTransxx[j,2]))){
    
    NewSphingoTransLayersxx$InCellulo[i,] <- ZerosToNAsConverter(LTPLipidConnectionsDataSphingolipidsWithCellularMatrixToRowMax2xx[j,])
  }}

for(j in which(OriginalRowsSphingoTransxx[,3] == "in vitro")){
  for(i in which(paste0(NewRowsSphingoTransxx[,1],NewRowsSphingoTransxx[,2]) == paste0(OriginalRowsSphingoTransxx[j,1],OriginalRowsSphingoTransxx[j,2]))){
    
    NewSphingoTransLayersxx$InVitro[i,] <- ZerosToNAsConverter(LTPLipidConnectionsDataSphingolipidsWithCellularMatrixToRowMax2xx[j,])
  }}

for(j in which(OriginalRowsSphingoTransxx[,3] == "Cellular")){
  for(i in which(NewRowsSphingoTransxx[,1] == OriginalRowsSphingoTransxx[j,1])){
    
    NewSphingoTransLayersxx$Cellular[i,] <- ZerosToNAsConverter(LTPLipidConnectionsDataSphingolipidsWithCellularMatrixToRowMax2xx[j,])
  }}

library(ComplexHeatmap)
library(RColorBrewer)

TextDataframeLeftLipidsLTPs2 <- data.frame(Protein = c(rep("CERT",5), rep("GLTPD1",4), "LCN1", rep("GLTP",4)),
                                           Lipid = c(rep("Cer", 5),rep("CerP",1),rep("SM",4), rep("xHexCer",4)))

AnnotationDataframeLeftLipidsLTPs2 <- rowAnnotation(df = TextDataframeLeftLipidsLTPs2, col = list(Protein = c("CERT" = "#FFCB3E", "LCN1" = "#FB836F", "GLTPD1" = "#C1549C", "GLTP" = "#7E549F"),
                                                                                                  Lipid = c("Cer" = "#ffc09f", "CerP" = "#fcf5c7", "SM" = "#a0ced9", "xHexCer" = "#b8e0d2")))

# LegendName <- "Legend"
# LegendColor <- col_fung
# 
# WidthAdaptor <- 160
# HeightAdaptor <- 160

pdf("./Output/BarplotHeatmapTestSphingolipidTransporters06092021WithSidesAdded2CorrecteddSMZeroToNA.pdf",
    width = unit(20, "mm"), height = unit(25, "mm"))

Heatmap(NewSphingoTransLayersxx$InCellulo, name = LegendName, col = LegendColor, rect_gp = gpar(type = "none"), column_title = "", 
        width = unit(WidthAdaptor*dim(NewSphingoTransLayersxx$InCellulo)[2]/min(dim(NewSphingoTransLayersxx$InCellulo)), "mm"), height = unit(HeightAdaptor*dim(NewSphingoTransLayersxx$InCellulo)[1]/min(dim(NewSphingoTransLayersxx$InCellulo)), "mm"),
        
        
        cell_fun = function(j, i, x, y, width, height, fill) {
          
          
          grid.rect(x = x, y = y, width = width, height = NewSphingoTransLayersxx$InCellulo[i,j]/100*height, gp = gpar(col = NA, fill = brewer.pal(9,"Blues")[5], lwd = 1), just = c("center","bottom"))
          grid.rect(x = x, y = y, width = 0.5*width, height = NewSphingoTransLayersxx$InVitro[i,j]/100*height, gp = gpar(col = NA, fill = brewer.pal(9,"Oranges")[5], lwd = 1), just = c("center","bottom"))
          
          grid.rect(x = x, y = y, width = width, height = NewSphingoTransLayersxx$Cellular[i,j]/100*height, gp = gpar(col = "LightGrey", fill = NA, lwd = 2), just = c("center","bottom"))
        }, cluster_rows = FALSE, cluster_columns = FALSE, 
        
        
        cluster_row_slices = FALSE, cluster_column_slices = FALSE,
        
        row_title = " ", row_title_gp = gpar(fontsize = 10),
        left_annotation = AnnotationDataframeLeftLipidsLTPs2,
        
        row_labels = sapply(strsplit(rownames(NewSphingoTransLayersxx$InCellulo), "_"), "[[", 1),
        row_split = c(rep("A", 5), rep("B",1), rep("C",3), "D", rep("E",4))
        
)
dev.off() 

# Previous figure is R-basis of final figure. 
# Final figure has LCN1 information removed, legends added, CERS information added to the top, and highlights of most important part focusses by putting black rectagles in Adobe Illustrator.

#### Extra information of CERT changes not included in final article because it just validates already existing observations
#. CERTAndSEC14L1LipidChangesHEK293WelchsTTest120820192.pdf #

min(OverexpressionHEK5[,2:ncol(OverexpressionHEK5)][OverexpressionHEK5[,2:ncol(OverexpressionHEK5)] != 0])
# 0.000100165 --> 0.00001 should be fine as background (to avoid 0 issues in the calculations); corrected the x+1 in lines below 

OverexpressionHEKRatiosMatrixAll <- do.call("cbind", lapply(1:9, function(x){(OverexpressionHEK5[,(2*x+1)] + 0.00001)/(OverexpressionHEK5[,(2*x)] + 0.00001)}))
colnames(OverexpressionHEKRatiosMatrixAll) <- paste0(colnames(OverexpressionHEK5)[sapply(1:9,function(x){(2*x+1)})], "_ratio")

rownames(OverexpressionHEKRatiosMatrixAll) <- OverexpressionHEK5[,1]
OverexpressionHEKLogRatiosMatrix <- log10(OverexpressionHEKRatiosMatrixAll)

OverexpressionHEKLogRatiosGeneral2 <- t(OverexpressionHEKLogRatiosMatrix[c("Cer", "CerP", "SM", "HexCer", "SHexCer", "diHexCer", "GM3", "DAG", "PA", "PA O-", "PC", "PC O-", "PE", "PE O-", "PS", "PI", "PI O-", "PG", "CL", "LPA", "LPC", "LPC O-", "LPE", "LPE O-", "LPS", "LPS O-", "LPI", "LPI O-", "LPG", "LPG O-", "CE", "Chol :"),c(1,4,7,2,5,8,3,6,9)])


OverexpressionHEKLogRatiosMelted <- meltmatrix(OverexpressionHEKLogRatiosGeneral2)
colnames(OverexpressionHEKLogRatiosMelted) <- c("Sample", "Lipid", "Ratio")


OverexpressionHEKLogRatiosMelted2 <- cbind(OverexpressionHEKLogRatiosMelted, LTPType = sapply(strsplit(as.character(OverexpressionHEKLogRatiosMelted[,"Sample"]), split = "_"), "[[", 3))

StatHEKLogRatiosMelted2 <- do.call("rbind",lapply(unique(as.character(OverexpressionHEKLogRatiosMelted2[,"Lipid"])), function(y){c(y,sapply(c("CERT","Sec14L1"), function(x){unlist(t.test(OverexpressionHEKLogRatiosMelted2[(OverexpressionHEKLogRatiosMelted2[,"LTPType"] == "control") & (OverexpressionHEKLogRatiosMelted2[,"Lipid"] == y),"Ratio"], OverexpressionHEKLogRatiosMelted2[(OverexpressionHEKLogRatiosMelted2[,"LTPType"] == x) & (OverexpressionHEKLogRatiosMelted2[,"Lipid"] == y),"Ratio"])[c(3,5,4)])}))}))
colnames(StatHEKLogRatiosMelted2) <- c("Lipid", "CERTp.value", "CERTControlMeanEstimate", "CERTProteinMeanEstimate", "CERTDiffConfidenceLow", "CERTDiffConfidenceHigh", "SEC14L1p.value", "SEC14L1ControlMeanEstimate", "SEC14L1ProteinMeanEstimate", "SEC14L1DiffConfidenceLow", "SEC14L1DiffConfidenceHigh")

StatHEKLogRatiosMelted4 <- do.call("cbind", list(StatHEKLogRatiosMelted2, CERTPercentDifference = 100*10^(as.numeric(StatHEKLogRatiosMelted2[,4]) - as.numeric(StatHEKLogRatiosMelted2[,3]))-100, SEC14L1PercentDifference = 100*10^(as.numeric(StatHEKLogRatiosMelted2[,9]) - as.numeric(StatHEKLogRatiosMelted2[,8]))-100))  


pdf("./Output/CERTAndSEC14L1LipidChangesHEK293WelchsTTest120820192.pdf")

plot(as.numeric(StatHEKLogRatiosMelted4[,"CERTPercentDifference"]), -log10(as.numeric(StatHEKLogRatiosMelted4[,"CERTp.value"])), xlab = "Percent difference lipid abundance: protein vs. control", ylab = "-log10(p-value) [higher is better] (Welch's 2 sample t-test) ", pch = 16, col = "#FC4E07", main = "CERT (HEK293)", xlim = c(-100,550))
abline(h = -log10(0.05), col = "lightgrey")

points(as.numeric(StatHEKLogRatiosMelted4[as.numeric(StatHEKLogRatiosMelted4[,"CERTp.value"]) <= 0.075,"CERTPercentDifference"]), -log10(as.numeric(StatHEKLogRatiosMelted4[as.numeric(StatHEKLogRatiosMelted4[,"CERTp.value"]) <= 0.075,"CERTp.value"])), pch = 16, col = "#E7B800")
points(as.numeric(StatHEKLogRatiosMelted4[as.numeric(StatHEKLogRatiosMelted4[,"CERTp.value"]) <= 0.05,"CERTPercentDifference"]), -log10(as.numeric(StatHEKLogRatiosMelted4[as.numeric(StatHEKLogRatiosMelted4[,"CERTp.value"]) <= 0.05,"CERTp.value"])), pch = 16, col = "#00AFBB")

text(as.numeric(StatHEKLogRatiosMelted4[as.numeric(StatHEKLogRatiosMelted4[,"CERTp.value"]) <= 0.075,"CERTPercentDifference"]), -log10(as.numeric(StatHEKLogRatiosMelted4[as.numeric(StatHEKLogRatiosMelted4[,"CERTp.value"]) <= 0.075,"CERTp.value"])), labels = paste0(StatHEKLogRatiosMelted4[as.numeric(StatHEKLogRatiosMelted4[,"CERTp.value"]) <= 0.075,"Lipid"], " (", round(as.numeric(StatHEKLogRatiosMelted4[as.numeric(StatHEKLogRatiosMelted4[,"CERTp.value"]) <= 0.075,"CERTPercentDifference"])), "%)"), pos = 4, cex = 0.8)
text(500, -log10(0.05), labels = "p-value 0.05", pos = 3, cex = 0.8, col = "darkgrey")

dev.off()
# The above figure is extra information on HEK293 lipid subclass changes upon CERT overexpresion

#### Lipid subclass changes upon overexpression of CERT in HeLa cells
# See further

######## Expanded part after integration of more data and analyses upon the reviewer questions
########

#### Analysis of full cell lipidomics data after overexpression of different LTPs
# Preparation of the data

MolPFileFull <- read.csv(file = "./InputData/A4_Lipotype_Report_Gavin (QS-21-531)_molp_namechange2-1.csv", header = TRUE, sep = ",", as.is = TRUE)
colnames(MolPFileFull)[colnames(MolPFileFull) == "PITPNC1.noninduced.biorep.C"] <- "PITPNC1.noninduced.biorepC" # Correct wrong name conversion


ExperimentalMetadataForFull <- t(Col1ToRowNames(MolPFileFull[1990:1991, 8:266]))

MolPFileFull2 <- MolPFileFull[1:1989,]
rm(MolPFileFull)

LipidMetadataOfFull <- MolPFileFull2[,1:7]
colnames(LipidMetadataOfFull) <- c("Feature", "LipidSubclass", "StructuralCategory", "FunctionalCategory", "TotalChainlength", "TotalUnsat", "TotalHydroxyls")


LipidDataOfFull <- as.matrix(MolPFileFull2[,9:266])

rownames(LipidDataOfFull) <- MolPFileFull2[,1]
mode(LipidDataOfFull) <- "numeric"

range(colSums(LipidDataOfFull, na.rm = T) - 100)
# -2.563428e-08  2.667633e-08 # Fluctuations from a perfect 100%

rm(MolPFileFull2)
SampleMetaDataFull <- do.call("rbind", strsplit(colnames(LipidDataOfFull), split = "\\."))

min(LipidDataOfFull, na.rm = TRUE)
# 4.627692e-06 ==> best to use 10^-8 for possible background corrections

LipidDataOfFullx <- LipidDataOfFull
LipidDataOfFullx[is.na(LipidDataOfFullx)] <- 0

# Extended Data Table X1 with the abundances of lipid species upon overexpression of each LTP (3 induced and 3 non-induced for each LTP)
write.table(LipidDataOfFullx, file = "./Output/SupplentaryTable6CBasis21052025.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

# 10^-8 as imputation value (see above)
LipidDataOfFullxi <- LipidDataOfFullx + 10^-8

# After imputation above the below corresponds to log10((Induced+10^-8)/(Noninduced+10^-8))
LipidDataFullInducedNoninducedImputedLogRatio <- do.call("cbind", lapply(unique(SampleMetaDataFull[,1]), function(x){do.call("cbind", setNames(lapply(unique(SampleMetaDataFull[,3]), function(y){log10(LipidDataOfFullxi[, paste(x, "induced", y, sep = ".")] / LipidDataOfFullxi[, paste(x, "noninduced", y, sep = ".")])}), paste(x, unique(SampleMetaDataFull[,3]), sep = "_")))}))

# Make a lipid subclass matrix for LipidDataOfFullxi
# unique(LipidMetadataOfFull$LipidSubclass) # Intermediate test of subclasses

# [1] "Cer"    "DAG"    "HexCer" "LPA"    "LPC"    "LPC O-" "LPE"    "LPE O-" "LPG"    "LPI"    "LPS"    "PA"     "PC"     "PC O-"  "PE"     "PE O-"  "PG"     "PI"     "PS"     "SM"     "TAG"   
# [22] "CL"     "CE"     "Chol" 

OrderOfLipidSubclasses <- c("Cer","SM", "HexCer", "LPA", "PA", "LPC", "LPC O-", "PC", "PC O-", "LPE", "LPE O-", "PE", "PE O-", "LPS", "PS", "LPI", "PI", "LPG", "PG", "CL", "DAG", "TAG", "CE", "Chol")


LipidSubclassesDataFullxi <- do.call("rbind", lapply(OrderOfLipidSubclasses, function(x){if(is.matrix(LipidDataOfFullx[LipidMetadataOfFull$LipidSubclass == x,])){
  colSums(LipidDataOfFullx[LipidMetadataOfFull$LipidSubclass == x,]) + 10^-8
  
}else{LipidDataOfFullx[LipidMetadataOfFull$LipidSubclass == x,] + 10^-8}}))
rownames(LipidSubclassesDataFullxi) <- OrderOfLipidSubclasses

# # The below corresponds to log10((Induced+10^-8)/(Noninduced+10^-8)) for subclass data
SubclassDataFullInducedNoninducedImputedLogRatio <- do.call("cbind", lapply(unique(SampleMetaDataFull[,1]), function(x){do.call("cbind", setNames(lapply(unique(SampleMetaDataFull[,3]), function(y){log10(LipidSubclassesDataFullxi[, paste(x, "induced", y, sep = ".")] / LipidSubclassesDataFullxi[, paste(x, "noninduced", y, sep = ".")])}), paste(x, unique(SampleMetaDataFull[,3]), sep = "_")))}))

# All LTPs analyzed in the overexpression data
LTPColumns <- sapply(strsplit(colnames(SubclassDataFullInducedNoninducedImputedLogRatio), "_"), "[[", 1)

# Paired and non-paired t-test on data itself
# Each time for both lipid subclasses and species

NonPairedTTestResultsForSubclasses <- do.call("rbind", lapply(unique(SampleMetaDataFull[,1]), function(x){do.call("cbind.data.frame", list(do.call("rbind", lapply(rownames(LipidSubclassesDataFullxi), function(z){if(all(mean(LipidSubclassesDataFullxi[z,(SampleMetaDataFull[,1] == x)]) == LipidSubclassesDataFullxi[z,(SampleMetaDataFull[,1] == x)])){
  setNames(c(1,0,0,0,0), c("p.value", "estimate.mean of x", "estimate.mean of y", "conf.int1", "conf.int2")) 
  
}else{
  unlist(t.test(LipidSubclassesDataFullxi[z,(SampleMetaDataFull[,1] == x) & (SampleMetaDataFull[,2] == "noninduced")], LipidSubclassesDataFullxi[z,(SampleMetaDataFull[,1] == x) & (SampleMetaDataFull[,2] == "induced")])[c(3,5,4)])  
  
}})), 
Lipid = rownames(LipidSubclassesDataFullxi), LTPProtein = x))}))

PairedTTestResultsForSubclasses <- do.call("rbind", lapply(unique(SampleMetaDataFull[,1]), function(x){do.call("cbind.data.frame", list(do.call("rbind", lapply(rownames(LipidSubclassesDataFullxi), function(z){if(all(mean(LipidSubclassesDataFullxi[z,(SampleMetaDataFull[,1] == x)]) == LipidSubclassesDataFullxi[z,(SampleMetaDataFull[,1] == x)])){
  setNames(c(1,0,0,0,0), c("p.value", "estimate.mean of the differences", "conf.int1", "conf.int2","PercentDifference")) 
  
}else{
  c(unlist(t.test(LipidSubclassesDataFullxi[z,(SampleMetaDataFull[,1] == x) & (SampleMetaDataFull[,2] == "induced")], LipidSubclassesDataFullxi[z,(SampleMetaDataFull[,1] == x) & (SampleMetaDataFull[,2] == "noninduced")], paired = TRUE)[c(3,5,4)]), PercentDifference = mean(LipidSubclassesDataFullxi[z,(SampleMetaDataFull[,1] == x) & (SampleMetaDataFull[,2] == "induced")])*100/mean(LipidSubclassesDataFullxi[z,(SampleMetaDataFull[,1] == x) & (SampleMetaDataFull[,2] == "noninduced")])-100)  
  
}})), 
Lipid = rownames(LipidSubclassesDataFullxi), LTPProtein = x))}))[,c(1:4,6:7,5)] # Updated so differences are calculated as x-y and with difference.

NonPairedTTestResultsForSpecies <- do.call("rbind", lapply(unique(SampleMetaDataFull[,1]), function(x){do.call("cbind.data.frame", list(do.call("rbind", lapply(rownames(LipidDataOfFullxi), function(z){if(all(mean(LipidDataOfFullxi[z,(SampleMetaDataFull[,1] == x)]) == LipidDataOfFullxi[z,(SampleMetaDataFull[,1] == x)])){
  setNames(c(1,0,0,0,0), c("p.value", "estimate.mean of x", "estimate.mean of y", "conf.int1", "conf.int2")) 
  
}else{
  unlist(t.test(LipidDataOfFullxi[z,(SampleMetaDataFull[,1] == x) & (SampleMetaDataFull[,2] == "noninduced")], LipidDataOfFullxi[z,(SampleMetaDataFull[,1] == x) & (SampleMetaDataFull[,2] == "induced")])[c(3,5,4)])  
  
}})), 
Lipid = rownames(LipidDataOfFullxi), LTPProtein = x))}))

PairedTTestResultsForSpecies <- do.call("rbind", lapply(unique(SampleMetaDataFull[,1]), function(x){do.call("cbind.data.frame", list(do.call("rbind", lapply(rownames(LipidDataOfFullxi), function(z){if(all(mean(LipidDataOfFullxi[z,(SampleMetaDataFull[,1] == x)]) == LipidDataOfFullxi[z,(SampleMetaDataFull[,1] == x)])){
  setNames(c(1,0,0,0,0), c("p.value", "estimate.mean of the differences", "conf.int1", "conf.int2", "PercentDifference")) 
  
}else{
  c(unlist(t.test(LipidDataOfFullxi[z,(SampleMetaDataFull[,1] == x) & (SampleMetaDataFull[,2] == "induced")], LipidDataOfFullxi[z,(SampleMetaDataFull[,1] == x) & (SampleMetaDataFull[,2] == "noninduced")], paired = TRUE)[c(3,5,4)]), PercentDifference = mean(LipidDataOfFullxi[z,(SampleMetaDataFull[,1] == x) & (SampleMetaDataFull[,2] == "induced")])*100/mean(LipidDataOfFullxi[z,(SampleMetaDataFull[,1] == x) & (SampleMetaDataFull[,2] == "noninduced")])-100)  
  
}})), 
Lipid = rownames(LipidDataOfFullxi), LTPProtein = x))}))[,c(1:4,6:7,5)] # Updated so differences are calculated as x-y and the difference.

# Write t-tests for to a file
write.table(PairedTTestResultsForSpecies, file = "./Output/PairedTTestResultsSpecies18062024.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

NonPairedTTestResultsForSubclasses <- cbind(NonPairedTTestResultsForSubclasses, PercentDifference = ifelse(NonPairedTTestResultsForSubclasses[,"estimate.mean of x"] != 0, 100*(as.numeric(NonPairedTTestResultsForSubclasses[,"estimate.mean of y"]) / as.numeric(NonPairedTTestResultsForSubclasses[,"estimate.mean of x"]))-100, 0))
NonPairedTTestResultsForSpecies <- cbind(NonPairedTTestResultsForSpecies, PercentDifference = ifelse(NonPairedTTestResultsForSpecies[,"estimate.mean of x"] != 0, 100*(as.numeric(NonPairedTTestResultsForSpecies[,"estimate.mean of y"]) / as.numeric(NonPairedTTestResultsForSpecies[,"estimate.mean of x"]))-100, 0))

LogRatioDifferenceTTestsForLipidSubclasses <- do.call("rbind", lapply(unique(LTPColumns), function(x){do.call("cbind.data.frame", list(do.call("rbind", lapply(rownames(SubclassDataFullInducedNoninducedImputedLogRatio), function(z){unlist(t.test(SubclassDataFullInducedNoninducedImputedLogRatio[z, LTPColumns != x], SubclassDataFullInducedNoninducedImputedLogRatio[z, LTPColumns == x])[c(3,5,4)])})),
                                                                                                                                       Lipid = rownames(SubclassDataFullInducedNoninducedImputedLogRatio), LTPProtein = x))}))

LogRatioDifferenceTTestsForLipidSpecies <- do.call("rbind", lapply(unique(LTPColumns), function(x){do.call("cbind.data.frame", list(do.call("rbind", lapply(rownames(LipidDataFullInducedNoninducedImputedLogRatio), function(z){unlist(t.test(LipidDataFullInducedNoninducedImputedLogRatio[z, LTPColumns != x], LipidDataFullInducedNoninducedImputedLogRatio[z, LTPColumns == x])[c(3,5,4)])})),
                                                                                                                                    Lipid = rownames(LipidDataFullInducedNoninducedImputedLogRatio), LTPProtein = x))}))

LogRatioDifferenceTTestsForLipidSubclasses <- cbind(LogRatioDifferenceTTestsForLipidSubclasses, PercentDifference = 100*10^(as.numeric(LogRatioDifferenceTTestsForLipidSubclasses[,"estimate.mean of y"]) - as.numeric(LogRatioDifferenceTTestsForLipidSubclasses[,"estimate.mean of x"]))-100)
LogRatioDifferenceTTestsForLipidSpecies <- cbind(LogRatioDifferenceTTestsForLipidSpecies, PercentDifference = 100*10^(as.numeric(LogRatioDifferenceTTestsForLipidSpecies[,"estimate.mean of y"]) - as.numeric(LogRatioDifferenceTTestsForLipidSpecies[,"estimate.mean of x"]))-100)

LogRatioDifferenceTTestsForLipidSpeciesb <- cbind(LogRatioDifferenceTTestsForLipidSpecies, LogRatio = log10((10^as.numeric(LogRatioDifferenceTTestsForLipidSpecies[,"estimate.mean of y"]))/(10^as.numeric(LogRatioDifferenceTTestsForLipidSpecies[,"estimate.mean of x"]))))


# Export file in reduced format
write.table(LogRatioDifferenceTTestsForLipidSpeciesb[,c("LTPProtein", "Lipid", "PercentDifference", "LogRatio", "p.value")], file = "./Output/469777_SupplementaryTable6B_v1.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# Export the file again to form the basis of Ext. Data Table X2 but now expanded with the condensed lipid subclass names and updated names for LTPs (STARD11 will be switched to CERT naming also afterwards) # Maybe place later: is ok.
write.table(LogRatioDifferenceTTestsForLipidSpeciesbwls[,c("NewLTPs", "Lipid", "ConvertedOELipidName", "LogRatio", "p.value")], file = "./Output/469777_SupplementaryTable7C_v4.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE) 


pdf("./Output/Panel2bBeforeEstheticOptimizationOn3.2By3.2Dim24042025bWithoutAnyDingbats.pdf", width = 3.2, height = 3.2, useDingbats = FALSE)

LipidHighlightsInputList <- setNames(list(cbind(c("^DAG", "^TAG"),c("orange","brown")),
                                          do.call("rbind", list(c("^DAG", "orange"), c("^PC\\s\\d", "#72a11e"), c("^TAG", "brown"), c("^SM", "blue"))), # Maybe check without PC-O too.
                                          
                                          cbind(c("^PC O-", "^PC\\s\\d"),c("#a3e827", "#72a11e")),
                                          cbind(c("^PE O-", "^PE\\s\\d", "^PC O-", "^PC\\s\\d"),c("green", "darkgreen", "#a3e827", "#72a11e")),
                                          
                                          cbind(c("^DAG"),c("orange"))
), nm = c("HSDL2", "STARD11","STARD2", "STARD10", "SEC14L2"))


for(LTPFocus in names(LipidHighlightsInputList)){
  
  AdaptiveSubsetForLogRatioDifferenceTTestsForLipidSpeciesb <- LogRatioDifferenceTTestsForLipidSpeciesb[LogRatioDifferenceTTestsForLipidSpeciesb$LTPProtein == LTPFocus,]
  plot(AdaptiveSubsetForLogRatioDifferenceTTestsForLipidSpeciesb$LogRatio, -log10(AdaptiveSubsetForLogRatioDifferenceTTestsForLipidSpeciesb$p.value), pch = 16, col = rgb(0,0,0, 0.20), xlab = "Log10 of ratio of induced vs. non-induced", ylab = "-Log10 p.value", main = LTPFocus, xlim = c(-7.07,7.07), ylim = c(0,50))
  
  for(i in 1:dim(LipidHighlightsInputList[[LTPFocus]])[1]){points(AdaptiveSubsetForLogRatioDifferenceTTestsForLipidSpeciesb[grepl(pattern = LipidHighlightsInputList[[LTPFocus]][i,1], x = AdaptiveSubsetForLogRatioDifferenceTTestsForLipidSpeciesb$Lipid), "LogRatio"], -log10(AdaptiveSubsetForLogRatioDifferenceTTestsForLipidSpeciesb[grepl(pattern = LipidHighlightsInputList[[LTPFocus]][i,1], x = AdaptiveSubsetForLogRatioDifferenceTTestsForLipidSpeciesb$Lipid), "p.value"]), col = LipidHighlightsInputList[[LTPFocus]][i,2], pch = 16)
  }
  
  abline(h = -log10(0.05/dim(AdaptiveSubsetForLogRatioDifferenceTTestsForLipidSpeciesb)[1]), col = "grey")
  abline(v = 0, col = "grey")
  
  plot(AdaptiveSubsetForLogRatioDifferenceTTestsForLipidSpeciesb$LogRatio, jitter(rep(1,length(AdaptiveSubsetForLogRatioDifferenceTTestsForLipidSpeciesb$LogRatio)), amount = 0.14), pch = 16, col = rgb(0,0,0, 0.20), xlab = "Log10 of ratio of induced vs. non-induced", ylab = "", ylim = c(0,dim(LipidHighlightsInputList[[LTPFocus]])[1]+1.4), main = LTPFocus, xlim = c(-7.07,7.07), yaxt = "n")
  for(i in 1:dim(LipidHighlightsInputList[[LTPFocus]])[1]){points(AdaptiveSubsetForLogRatioDifferenceTTestsForLipidSpeciesb[grepl(pattern = LipidHighlightsInputList[[LTPFocus]][i,1], x = AdaptiveSubsetForLogRatioDifferenceTTestsForLipidSpeciesb$Lipid), "LogRatio"], jitter(rep(i+1,sum(grepl(pattern = LipidHighlightsInputList[[LTPFocus]][i,1], x = AdaptiveSubsetForLogRatioDifferenceTTestsForLipidSpeciesb$Lipid))), amount = 0.14), col = LipidHighlightsInputList[[LTPFocus]][i,2], pch = 16)}
  
  abline(v = 0, col = "grey")
  axis(2, at = 1:(dim(LipidHighlightsInputList[[LTPFocus]])[1]+1), label = c("All lipids", gsub(pattern = "\\^", replacement = "", LipidHighlightsInputList[[LTPFocus]][,1])), las = 2)
  
}


plot.new()
legend("center", c("SM", "PC", "PC-O", "PE", "PE-O", "DAG", "TAG", "Other"), fill=c("blue", "#72a11e", "#a3e827", "darkgreen", "green", "orange", "brown", "lightgrey"), border = FALSE, box.col = "lightgrey")

dev.off()
# Volcano plot outputs from previous figures formed the basis for Fig. 2b, later combined and esthetically refined in Adobe Illustrator.


ListOfLTPTTests <- list(NonPairedTTestResultsForSpecies, NonPairedTTestResultsForSubclasses, PairedTTestResultsForSpecies, PairedTTestResultsForSubclasses, LogRatioDifferenceTTestsForLipidSpecies, LogRatioDifferenceTTestsForLipidSubclasses)

ListOfLTPTTests <- lapply(ListOfLTPTTests, function(x){do.call("cbind", list(x, Star = ifelse(x$p.value <= 0.05, 
                                                                                              ifelse(x$p.value <= 0.01,
                                                                                                     
                                                                                                     "**","*"),""),
                                                                             PercentDifferenceCutOff = ifelse(x$PercentDifference >= 100, 100, x$PercentDifference)))})

write.table(ListOfLTPTTests[[1]], file = "./Output/NonPairedTTestResultsForSpecies.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(ListOfLTPTTests[[2]], file = "./Output/NonPairedTTestResultsForSubclasses.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

write.table(ListOfLTPTTests[[3]], file = "./Output/PairedTTestResultsForSpecies.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(ListOfLTPTTests[[4]], file = "./Output/PairedTTestResultsForSubclasses.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

write.table(ListOfLTPTTests[[5]], file = "./Output/LogRatioDifferenceTTestsForLipidSpecies.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(ListOfLTPTTests[[6]], file = "./Output/LogRatioDifferenceTTestsForLipidSubclasses.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


ListOfLTPTTests <- list()

ListOfLTPTTests <- lapply(c("NonPairedTTestResultsForSpecies", "NonPairedTTestResultsForSubclasses", "PairedTTestResultsForSpecies", "PairedTTestResultsForSubclasses", "LogRatioDifferenceTTestsForLipidSpecies", "LogRatioDifferenceTTestsForLipidSubclasses"), function(x){
  read.csv(file = paste0("./Output/", x, ".tsv"), header = TRUE, sep = "\t", as.is = TRUE)})

ListOfLTPTTestsfs <- ListOfLTPTTests[c(1,3,5)]
ListOfLTPTTestsfs <- lapply(ListOfLTPTTestsfs, function(x){cbind(x, LipidMetadataOfFull[match(x[,"Lipid"], LipidMetadataOfFull$Feature), ])})

# For subset of only true GPLs with 2 fatty acid chains
ListOfLTPTTestsfsSubsetoes <- lapply(ListOfLTPTTestsfs, function(x){x[x$LipidSubclass %in% c("PA", "PC", "PE", "PG", "PI", "PS"),]})

# Write intermediate file to save it for overview excel file after subsetting for true GPLs with 2 legs # Basis for the Extended Data Table 7D.
write.table(ListOfLTPTTestsfsSubsetoes[[2]], file = "./Output/PairedTTestResultsSpeciesOnlySubsetFocus18062024.tsv.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

ChainLengthsOfSignificancesSubsetoes <- lapply(ListOfLTPTTestsfsSubsetoes, function(x){setNames(merge(aggregate(x[x[,"Star"] == "","TotalChainlength"], by = list(x[x[,"Star"] == "","TotalChainlength"]), FUN=length),
                                                                                                      aggregate(x[x[,"Star"] != "","TotalChainlength"], by = list(x[x[,"Star"] != "","TotalChainlength"]), FUN=length), by = "Group.1", all = TRUE), c("TotalChainlength", "Non-significant", "Significant"))})

ChainLengthsOfSignificancesSubsetoes <- lapply(ChainLengthsOfSignificancesSubsetoes, function(x){x[is.na(x)] <- 0; return(x)})
ChainLengthsOfSignificancesSubsetoes <- lapply(ChainLengthsOfSignificancesSubsetoes, function(x){cbind(x, SignificantFraction = x[,"Significant"]*100/(x[,"Non-significant"]+x[,"Significant"]))})


# Input of the data from the mobilization to make the comparisons (see earlier)

LTPMobilizedGPL2sInCelluloSumIntensities <- colSums(HeatmapMatrixBothAntonellan[c("PA", "PC", "PE", "PI", "PS", "PGP", "PG", "PG/BMP", "BMP"),])
LTPMobilizedGPL2sInVitroSumIntensities <- colSums(HeatmapMatrixBothEnricn[c("PA", "PC", "PE", "PI", "PS", "PGP", "PG", "PG/BMP", "BMP"),])

LTPMobilizedGPL2sInCelluloMaxNorm <- LTPMobilizedGPL2sInCelluloSumIntensities*100/max(LTPMobilizedGPL2sInCelluloSumIntensities)
LTPMobilizedGPL2sInVitroMaxNorm <- LTPMobilizedGPL2sInVitroSumIntensities*100/max(LTPMobilizedGPL2sInVitroSumIntensities)


# Similar treatment of the overexpression data for comparisons

LTPOEResultsGPL2sInVitroMaxNormoes <- setNames(ChainLengthsOfSignificancesSubsetoes[[2]][,4]*100/max(ChainLengthsOfSignificancesSubsetoes[[2]][,4]), ChainLengthsOfSignificancesSubsetoes[[2]][,1])
LTPOEResultsGPL2sInVitroMaxNormoes <- setNames(LTPOEResultsGPL2sInVitroMaxNormoes[as.character(16:70)], as.character(16:70))

LTPOEResultsGPL2sInVitroMaxNormoes[is.na(LTPOEResultsGPL2sInVitroMaxNormoes)] <- 0



# From the cellular lipidome of the HEK293 cells (control background for mobilization data)

LTPBackgroundLipidsGPL2s <- colSums(DetailsOfOverexpressionHEK2Split2bl2WideVersion2[c("PA", "PC", "PE", "PI", "PS", "PGP", "PG", "PG/BMP", "BMP"),])[as.character(16:70)]
LTPBackgroundLipidsGPL2sWithMaxNorm <- LTPBackgroundLipidsGPL2s*100/max(LTPBackgroundLipidsGPL2s)


# Combination of these datasets

LTPMobilizedGPL2sMaxNormCombinedVersionoes <- do.call("rbind", list(LTPMobilizedGPL2sInCelluloMaxNorm, LTPMobilizedGPL2sInVitroMaxNorm, LTPBackgroundLipidsGPL2sWithMaxNorm, LTPOEResultsGPL2sInVitroMaxNormoes))
rownames(LTPMobilizedGPL2sMaxNormCombinedVersionoes) <- c("in cellulo mobilized lipids", "in vitro mobilized lipids", "general lipidome", "changed lipids after overexpression")


# t-tests of differences between chainlengths 28:38 vs. 39:44

t.test(LTPMobilizedGPL2sMaxNormCombinedVersionoes[4,as.character(28:38)], LTPMobilizedGPL2sMaxNormCombinedVersionoes[4,as.character(39:44)])
# Welch Two Sample t-test

# data:  LTPMobilizedGPL2sMaxNormCombinedVersionoes[4, as.character(28:38)] and LTPMobilizedGPL2sMaxNormCombinedVersionoes[4, as.character(39:44)]
# t = 6.7094, df = 12.756, p-value = 1.595e-05

# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:

#   35.57803 69.46737
# sample estimates:

#   mean of x mean of y 
# 66.92891  14.40621 


# Rework for clean-up and also use for unsaturations

ListOfLTPTTestsfsSubsetoes <- lapply(ListOfLTPTTestsfs, function(x){x[x$LipidSubclass %in% c("PA", "PC", "PE", "PG", "PI", "PS"),]})
# Only 2th element needed below and above


ListInputsChainlengthUnsat <- list(list(ListOfLTPTTestsfsSubsetoes[[2]], "TotalChainlength"), list(ListOfLTPTTestsfs[[2]], "TotalUnsat"), list(ListOfLTPTTestsfsSubsetoes[[2]], "TotalUnsat"))

VariousOfSignificancesSubsetoes <- lapply(ListInputsChainlengthUnsat, function(x){setNames(merge(aggregate(x[[1]][x[[1]][,"Star"] == "",x[[2]]], by = list(x[[1]][x[[1]][,"Star"] == "",x[[2]]]), FUN=length),
                                                                                                 aggregate(x[[1]][x[[1]][,"Star"] != "",x[[2]]], by = list(x[[1]][x[[1]][,"Star"] != "",x[[2]]]), FUN=length), by = "Group.1", all = TRUE), c(x[[2]], "Non-significant", "Significant"))})

VariousOfSignificancesSubsetoes <- lapply(VariousOfSignificancesSubsetoes, function(x){x[is.na(x)] <- 0; return(x)})
VariousOfSignificancesSubsetoes <- lapply(VariousOfSignificancesSubsetoes, function(x){cbind(x, SignificantFraction = x[,"Significant"]*100/(x[,"Non-significant"]+x[,"Significant"]))})

LTPOEResultsGPL2sOfVariousMaxNormoes <-lapply(VariousOfSignificancesSubsetoes, function(x){setNames(x[,4]*100/max(x[,4]), x[,1])})
LTPOEResultsGPL2sOfVariousMaxNormoes <- lapply(LTPOEResultsGPL2sOfVariousMaxNormoes, function(x){NAsToZerosConverter(setNames(x[as.character(min(as.numeric(names(x))):max(as.numeric(names(x))))], as.character(min(as.numeric(names(x))):max(as.numeric(names(x))))))})

# t-tests of differences between chainlengths 28:38 vs. 39:44 (alternative but fundamentally the same)
t.test(LTPOEResultsGPL2sOfVariousMaxNormoes[[1]][as.character(28:38)], LTPOEResultsGPL2sOfVariousMaxNormoes[[1]][as.character(39:44)])

# Result
# Welch Two Sample t-test

# data:  LTPOEResultsGPL2sOfVariousMaxNormoes[[1]][as.character(28:38)] and LTPOEResultsGPL2sOfVariousMaxNormoes[[1]][as.character(39:44)]
# t = 6.7094, df = 12.756, p-value = 1.595e-05

# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:

#   35.57803 69.46737
# sample estimates:

#   mean of x mean of y 
# 66.92891  14.40621 

# t-test for total unsaturations for all the lipids, comparison 1-2 vs. others (0,3-12)
t.test(LTPOEResultsGPL2sOfVariousMaxNormoes[[2]][as.character(1:2)], LTPOEResultsGPL2sOfVariousMaxNormoes[[2]][as.character(c(0,3:12))])

# Result
# Welch Two Sample t-test

# data:  LTPOEResultsGPL2sOfVariousMaxNormoes[[2]][as.character(1:2)] and LTPOEResultsGPL2sOfVariousMaxNormoes[[2]][as.character(c(0, 3:12))]
# t = 2.9146, df = 8.4571, p-value = 0.01833

# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:

#   9.126593 75.299253
# sample estimates:

#   mean of x mean of y 
# 92.77953  50.56661 

# t-test for total unsaturations for all the lipids, comparison 1-2 vs. others (0,3-12): subset of lipids under focus
t.test(LTPOEResultsGPL2sOfVariousMaxNormoes[[3]][as.character(1:2)], LTPOEResultsGPL2sOfVariousMaxNormoes[[3]][as.character(c(0,3:12))])

# Result 
# Welch Two Sample t-test

# data:  LTPOEResultsGPL2sOfVariousMaxNormoes[[3]][as.character(1:2)] and LTPOEResultsGPL2sOfVariousMaxNormoes[[3]][as.character(c(0, 3:12))]
# t = 4.2507, df = 4.9166, p-value = 0.008393

# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:

#   22.79129 93.44278
# sample estimates:

#   mean of x mean of y 
# 91.23815  33.12111 

########
# Analysis of vesicle contents of in vitro screen based on lipidomes of the input-vesicles

InVitroVesicleLipidomes <- list()
InVitroVesicleLipidomes2 <- list()

InVitroVesicleLipidomes[["LiverPos"]] <- read.csv(file = "./InputData/LIVER_pos_filtered_annotated.csv", header = TRUE, sep = ",", as.is = TRUE)
InVitroVesicleLipidomes[["LiverNeg"]] <- read.csv(file = "./InputData/LIVER_neg_filtered_annotated.csv", header = TRUE, sep = ",", as.is = TRUE)

InVitroVesicleLipidomes[["BrainPos"]] <- read.csv(file = "./InputData/BRAIN_pos_filtered_annotated.csv", header = TRUE, sep = ",", as.is = TRUE)
InVitroVesicleLipidomes[["BrainNeg"]] <- read.csv(file = "./InputData/BRAIN_neg_filtered_annotated.csv", header = TRUE, sep = ",", as.is = TRUE)

InVitroVesicleLipidomes2[["LiverPos"]] <- InVitroVesicleLipidomes[["LiverPos"]][rowSums(InVitroVesicleLipidomes[["LiverPos"]][,c("X.M.H..", "X.M.Na..", "X.M.NH4..")] != "[]") > 0,]
InVitroVesicleLipidomes2[["LiverNeg"]] <- InVitroVesicleLipidomes[["LiverNeg"]][rowSums(InVitroVesicleLipidomes[["LiverNeg"]][,c("X.M.H..", "X.M.HCOO..")] != "[]") > 0,]

InVitroVesicleLipidomes2[["BrainPos"]] <- InVitroVesicleLipidomes[["BrainPos"]][rowSums(InVitroVesicleLipidomes[["BrainPos"]][,c("X.M.H..", "X.M.Na..", "X.M.NH4..")] != "[]") > 0,]
InVitroVesicleLipidomes2[["BrainNeg"]] <- InVitroVesicleLipidomes[["BrainNeg"]][rowSums(InVitroVesicleLipidomes[["BrainNeg"]][,c("X.M.H..", "X.M.HCOO..")] != "[]") > 0,]

rm(InVitroVesicleLipidomes)


AllLipidSpeciesInData <- unique(unlist(lapply(list(InVitroVesicleLipidomes2[["LiverPos"]][,c("X.M.H..")],
                                                   InVitroVesicleLipidomes2[["LiverPos"]][,c("X.M.Na..")],
                                                   
                                                   InVitroVesicleLipidomes2[["LiverPos"]][,c("X.M.NH4..")],
                                                   InVitroVesicleLipidomes2[["LiverNeg"]][,c("X.M.H..")],
                                                   
                                                   InVitroVesicleLipidomes2[["LiverNeg"]][,c("X.M.HCOO..")],
                                                   
                                                   
                                                   InVitroVesicleLipidomes2[["BrainPos"]][,c("X.M.H..")],
                                                   InVitroVesicleLipidomes2[["BrainPos"]][,c("X.M.Na..")],
                                                   
                                                   InVitroVesicleLipidomes2[["BrainPos"]][,c("X.M.NH4..")],
                                                   InVitroVesicleLipidomes2[["BrainNeg"]][,c("X.M.H..")],
                                                   
                                                   InVitroVesicleLipidomes2[["BrainNeg"]][,c("X.M.HCOO..")]), function(y){
                                                     unique(unlist(strsplit(gsub("'\n '", "__", gsub("' '", "__", gsub("\\[\\]", "", gsub("'\\]","",gsub("\\['","", y))))), "__")))})))

HeadgroupConv <- cbind(FullName = unique(sapply(strsplit(AllLipidSpeciesInData , split = " \\("), "[[", 1)),
                       ShortName = c("Cer", "LPE", "DAG", "GM4", "PS", "PC", "PIP", "PE", "HexCer", "TAG", "DAG", "TAG", "SM", "CerP", "LPG", "FA(18:3)", "LPA", "PA", "LPS", "LPC", "PI", "Hex2DAG", "PG", "BMP", "PIP3", "PGP", "Gb3", "PECer", "LPI", "FA(16:0)", "HexDAG", "SHex2Cer", "GM3", "GA1", "PIP2", "Hex2Cer", "FA(20:3)", "SHexCer", "FA(20:5)", "MAG", "SE", "FA(16:1)", "GD3", "FA(18:0)", "FA(38:5)", "FA(17:0)", "GA2", "FA(20:4)", "FA(18:2)", "FA(22:4)", "FA(18:1)", "GM2", "FA(30:6)", "FA(24:6)", "GM1", "FA(22:6)"))

SpeciesMatches <- do.call("cbind", list(FullSpecies = AllLipidSpeciesInData, 
                                        FullHeadgroup = sapply(strsplit(AllLipidSpeciesInData , split = " \\("), "[[", 1),
                                        
                                        ShortSpecies = HeadgroupConv[match(sapply(strsplit(AllLipidSpeciesInData , split = " \\("), "[[", 1), HeadgroupConv[,"FullName"]), "ShortName"]
))

SpeciesMatches[grep("\\(O-", SpeciesMatches[,"FullSpecies"]), "ShortSpecies"] <- paste0(SpeciesMatches[grep("\\(O-", SpeciesMatches[,"FullSpecies"]), "ShortSpecies"], "-O")
SpeciesMatches[grep("\\(d", SpeciesMatches[,"FullSpecies"]), "ShortSpecies"] <- paste0("d", SpeciesMatches[grep("\\(d", SpeciesMatches[,"FullSpecies"]), "ShortSpecies"])

SpeciesMatches[grep("\\(t", SpeciesMatches[,"FullSpecies"]), "ShortSpecies"] <- paste0("t", SpeciesMatches[grep("\\(t", SpeciesMatches[,"FullSpecies"]), "ShortSpecies"])
SpeciesMatches <- cbind(SpeciesMatches, SubClass = sapply(strsplit(SpeciesMatches[,"ShortSpecies"], "\\("), "[[", 1))

SpeciesMatches <- cbind(SpeciesMatches, SpeciesDetails = sapply(strsplit(sapply(strsplit(paste0(SpeciesMatches[,"FullSpecies"], SpeciesMatches[,"ShortSpecies"]), "\\)"), "[[", 1), "\\("), "[[", 2))
SpeciesMatches <- do.call("cbind", list(SpeciesMatches, 
                                        
                                        OnlyChains = sapply(strsplit(SpeciesMatches[,"SpeciesDetails"], ":"), "[[", 1),
                                        TotalUnsat = sapply(strsplit(SpeciesMatches[,"SpeciesDetails"], ":"), "[[", 2)))

SpeciesMatches <- cbind(SpeciesMatches, TotalChainlengths = gsub("O-|d|t","", SpeciesMatches[,"OnlyChains"]))
SpeciesMatches <- cbind(SpeciesMatches, ShortFull = paste0(SpeciesMatches[,"SubClass"], "(", SpeciesMatches[,"TotalChainlengths"], ":", SpeciesMatches[,"TotalUnsat"], ")") )

# Information on the adduct needed
library(reshape2) #! Package dependency that is removable

InVitroVesicleLipidomes4 <- list()
InVitroVesicleLipidomes5 <- list()

InVitroVesicleLipidomes4[["LiverPos"]] <- melt(InVitroVesicleLipidomes2[["LiverPos"]][, c("index", "m.z", "RT.range", "LIVER_N50_pos.m.z", "LIVER_N50_pos.RT.mean", "LIVER_N50_pos.Normalized.Area", "RT.mean", "X.M.H..", "X.M.Na..", "X.M.NH4..")],
                                               id.vars=c("index", "m.z", "RT.range", "LIVER_N50_pos.m.z", "LIVER_N50_pos.RT.mean", "LIVER_N50_pos.Normalized.Area", "RT.mean"),
                                               
                                               variable.name="Adduct",
                                               value.name="Lipids")

InVitroVesicleLipidomes4[["BrainPos"]] <- melt(InVitroVesicleLipidomes2[["BrainPos"]][, c("index", "m.z", "RT.range", "BRAIN_N50_pos.m.z", "BRAIN_N50_pos.RT.mean", "BRAIN_N50_pos.Normalized.Area", "RT.mean", "X.M.H..", "X.M.Na..", "X.M.NH4..")],
                                               id.vars=c("index", "m.z", "RT.range", "BRAIN_N50_pos.m.z", "BRAIN_N50_pos.RT.mean", "BRAIN_N50_pos.Normalized.Area", "RT.mean"),
                                               
                                               variable.name="Adduct",
                                               value.name="Lipids")

InVitroVesicleLipidomes4[["LiverNeg"]] <- melt(InVitroVesicleLipidomes2[["LiverNeg"]][, c("index", "m.z", "RT.range", "LIVER_N50_neg.m.z", "LIVER_N50_neg.RT.mean", "LIVER_N50_neg.Normalized.Area", "RT.mean", "X.M.H..", "X.M.HCOO..")],
                                               id.vars=c("index", "m.z", "RT.range", "LIVER_N50_neg.m.z", "LIVER_N50_neg.RT.mean", "LIVER_N50_neg.Normalized.Area", "RT.mean"),
                                               
                                               variable.name="Adduct",
                                               value.name="Lipids")

InVitroVesicleLipidomes4[["BrainNeg"]] <- melt(InVitroVesicleLipidomes2[["BrainNeg"]][, c("index", "m.z", "RT.range", "BRAIN_N50_neg.m.z", "BRAIN_N50_neg.RT.mean", "BRAIN_N50_neg.Normalized.Area", "RT.mean", "X.M.H..", "X.M.HCOO..")],
                                               id.vars=c("index", "m.z", "RT.range", "BRAIN_N50_neg.m.z", "BRAIN_N50_neg.RT.mean", "BRAIN_N50_neg.Normalized.Area", "RT.mean"),
                                               
                                               variable.name="Adduct",
                                               value.name="Lipids")

InVitroVesicleLipidomes4[["LiverPos"]] <- cbind(InVitroVesicleLipidomes4[["LiverPos"]], CleanedLipids = gsub("'\n '", "__", gsub("' '", "__", gsub("\\[\\]", "", gsub("'\\]","",gsub("\\['","", InVitroVesicleLipidomes4[["LiverPos"]][,"Lipids"]))))))
InVitroVesicleLipidomes4[["LiverPos"]] <- InVitroVesicleLipidomes4[["LiverPos"]][InVitroVesicleLipidomes4[["LiverPos"]][,"CleanedLipids"] != "",]

InVitroVesicleLipidomes4[["BrainPos"]] <- cbind(InVitroVesicleLipidomes4[["BrainPos"]], CleanedLipids = gsub("'\n '", "__", gsub("' '", "__", gsub("\\[\\]", "", gsub("'\\]","",gsub("\\['","", InVitroVesicleLipidomes4[["BrainPos"]][,"Lipids"]))))))
InVitroVesicleLipidomes4[["BrainPos"]] <- InVitroVesicleLipidomes4[["BrainPos"]][InVitroVesicleLipidomes4[["BrainPos"]][,"CleanedLipids"] != "",]

InVitroVesicleLipidomes4[["LiverNeg"]] <- cbind(InVitroVesicleLipidomes4[["LiverNeg"]], CleanedLipids = gsub("'\n '", "__", gsub("' '", "__", gsub("\\[\\]", "", gsub("'\\]","",gsub("\\['","", InVitroVesicleLipidomes4[["LiverNeg"]][,"Lipids"]))))))
InVitroVesicleLipidomes4[["LiverNeg"]] <- InVitroVesicleLipidomes4[["LiverNeg"]][InVitroVesicleLipidomes4[["LiverNeg"]][,"CleanedLipids"] != "",]

InVitroVesicleLipidomes4[["BrainNeg"]] <- cbind(InVitroVesicleLipidomes4[["BrainNeg"]], CleanedLipids = gsub("'\n '", "__", gsub("' '", "__", gsub("\\[\\]", "", gsub("'\\]","",gsub("\\['","", InVitroVesicleLipidomes4[["BrainNeg"]][,"Lipids"]))))))
InVitroVesicleLipidomes4[["BrainNeg"]] <- InVitroVesicleLipidomes4[["BrainNeg"]][InVitroVesicleLipidomes4[["BrainNeg"]][,"CleanedLipids"] != "",]

InVitroVesicleLipidomes4[["LiverPos"]][,"Adduct"] <- as.character(InVitroVesicleLipidomes4[["LiverPos"]][,"Adduct"])
InVitroVesicleLipidomes4[["BrainPos"]][,"Adduct"] <- as.character(InVitroVesicleLipidomes4[["BrainPos"]][,"Adduct"])

InVitroVesicleLipidomes4[["LiverNeg"]][,"Adduct"] <- as.character(InVitroVesicleLipidomes4[["LiverNeg"]][,"Adduct"])
InVitroVesicleLipidomes4[["BrainNeg"]][,"Adduct"] <- as.character(InVitroVesicleLipidomes4[["BrainNeg"]][,"Adduct"])

InVitroVesicleLipidomes5[["LiverPos"]] <- do.call("rbind.data.frame", lapply(1:dim(InVitroVesicleLipidomes4[["LiverPos"]])[1], function(y){
  do.call("rbind", lapply(unlist(strsplit(InVitroVesicleLipidomes4[["LiverPos"]][y,"CleanedLipids"], "__")), function(x){unlist(c(InVitroVesicleLipidomes4[["LiverPos"]][y,], CleanedFurther = x))}))}))

InVitroVesicleLipidomes5[["BrainPos"]] <- do.call("rbind.data.frame", lapply(1:dim(InVitroVesicleLipidomes4[["BrainPos"]])[1], function(y){
  do.call("rbind", lapply(unlist(strsplit(InVitroVesicleLipidomes4[["BrainPos"]][y,"CleanedLipids"], "__")), function(x){unlist(c(InVitroVesicleLipidomes4[["BrainPos"]][y,], CleanedFurther = x))}))}))

InVitroVesicleLipidomes5[["LiverNeg"]] <- do.call("rbind.data.frame", lapply(1:dim(InVitroVesicleLipidomes4[["LiverNeg"]])[1], function(y){
  do.call("rbind", lapply(unlist(strsplit(InVitroVesicleLipidomes4[["LiverNeg"]][y,"CleanedLipids"], "__")), function(x){unlist(c(InVitroVesicleLipidomes4[["LiverNeg"]][y,], CleanedFurther = x))}))}))

InVitroVesicleLipidomes5[["BrainNeg"]] <- do.call("rbind.data.frame", lapply(1:dim(InVitroVesicleLipidomes4[["BrainNeg"]])[1], function(y){
  do.call("rbind", lapply(unlist(strsplit(InVitroVesicleLipidomes4[["BrainNeg"]][y,"CleanedLipids"], "__")), function(x){unlist(c(InVitroVesicleLipidomes4[["BrainNeg"]][y,], CleanedFurther = x))}))}))

InVitroVesicleLipidomes5[["LiverPos"]] <- cbind.data.frame(InVitroVesicleLipidomes5[["LiverPos"]], SpeciesMatches[match(InVitroVesicleLipidomes5[["LiverPos"]][,"CleanedFurther"],
                                                                                                                        SpeciesMatches[,"FullSpecies"]), c("ShortFull", "SubClass", "TotalChainlengths", "TotalUnsat")])

InVitroVesicleLipidomes5[["BrainPos"]] <- cbind.data.frame(InVitroVesicleLipidomes5[["BrainPos"]], SpeciesMatches[match(InVitroVesicleLipidomes5[["BrainPos"]][,"CleanedFurther"],
                                                                                                                        SpeciesMatches[,"FullSpecies"]), c("ShortFull", "SubClass", "TotalChainlengths", "TotalUnsat")])

InVitroVesicleLipidomes5[["LiverNeg"]] <- cbind.data.frame(InVitroVesicleLipidomes5[["LiverNeg"]], SpeciesMatches[match(InVitroVesicleLipidomes5[["LiverNeg"]][,"CleanedFurther"],
                                                                                                                        SpeciesMatches[,"FullSpecies"]), c("ShortFull", "SubClass", "TotalChainlengths", "TotalUnsat")])

InVitroVesicleLipidomes5[["BrainNeg"]] <- cbind.data.frame(InVitroVesicleLipidomes5[["BrainNeg"]], SpeciesMatches[match(InVitroVesicleLipidomes5[["BrainNeg"]][,"CleanedFurther"],
                                                                                                                        SpeciesMatches[,"FullSpecies"]), c("ShortFull", "SubClass", "TotalChainlengths", "TotalUnsat")])

InVitroVesicleLipidomes5[["LiverPos"]][, c("index", "m.z", "LIVER_N50_pos.m.z", "LIVER_N50_pos.RT.mean", "LIVER_N50_pos.Normalized.Area", "RT.mean", "TotalChainlengths", "TotalUnsat")] <- apply(InVitroVesicleLipidomes5[["LiverPos"]][, c("index", "m.z", "LIVER_N50_pos.m.z", "LIVER_N50_pos.RT.mean", "LIVER_N50_pos.Normalized.Area", "RT.mean", "TotalChainlengths", "TotalUnsat")], 2, function(x) as.numeric(as.character(x)))
InVitroVesicleLipidomes5[["BrainPos"]][, c("index", "m.z", "BRAIN_N50_pos.m.z", "BRAIN_N50_pos.RT.mean", "BRAIN_N50_pos.Normalized.Area", "RT.mean", "TotalChainlengths", "TotalUnsat")] <- apply(InVitroVesicleLipidomes5[["BrainPos"]][, c("index", "m.z", "BRAIN_N50_pos.m.z", "BRAIN_N50_pos.RT.mean", "BRAIN_N50_pos.Normalized.Area", "RT.mean", "TotalChainlengths", "TotalUnsat")], 2, function(x) as.numeric(as.character(x)))

InVitroVesicleLipidomes5[["LiverNeg"]][, c("index", "m.z", "LIVER_N50_neg.m.z", "LIVER_N50_neg.RT.mean", "LIVER_N50_neg.Normalized.Area", "RT.mean", "TotalChainlengths", "TotalUnsat")] <- apply(InVitroVesicleLipidomes5[["LiverNeg"]][, c("index", "m.z", "LIVER_N50_neg.m.z", "LIVER_N50_neg.RT.mean", "LIVER_N50_neg.Normalized.Area", "RT.mean", "TotalChainlengths", "TotalUnsat")], 2, function(x) as.numeric(as.character(x)))
InVitroVesicleLipidomes5[["BrainNeg"]][, c("index", "m.z", "BRAIN_N50_neg.m.z", "BRAIN_N50_neg.RT.mean", "BRAIN_N50_neg.Normalized.Area", "RT.mean", "TotalChainlengths", "TotalUnsat")] <- apply(InVitroVesicleLipidomes5[["BrainNeg"]][, c("index", "m.z", "BRAIN_N50_neg.m.z", "BRAIN_N50_neg.RT.mean", "BRAIN_N50_neg.Normalized.Area", "RT.mean", "TotalChainlengths", "TotalUnsat")], 2, function(x) as.numeric(as.character(x)))

# Determine retention time minima and maxima per subclass in the mobilization-data (for data filtering)
PureMobilizationDataCombined <- rbind(PureAntonella32b, PureEnric32)

OverviewRTsMobilizedLipids <- merge(merge(aggregate(PureMobilizationDataCombined$MinimumOfRetentionTime, by = list(PureMobilizationDataCombined$LikelySubclass), FUN = function(x){min(x,na.rm = TRUE)}),
                                          aggregate(PureMobilizationDataCombined$MaximumOfRetentionTime, by = list(PureMobilizationDataCombined$LikelySubclass), FUN = function(x){max(x,na.rm = TRUE)}), by = "Group.1", all = TRUE),
                                    
                                    merge(merge(aggregate(PureMobilizationDataCombined$TotalCarbonChainLength, by = list(PureMobilizationDataCombined$LikelySubclass), FUN = function(x){min(x,na.rm = TRUE)}),
                                                aggregate(PureMobilizationDataCombined$TotalCarbonChainLength, by = list(PureMobilizationDataCombined$LikelySubclass), FUN = function(x){max(x,na.rm = TRUE)}), by = "Group.1", all = TRUE),
                                          
                                          merge(aggregate(PureMobilizationDataCombined$TotalCarbonChainUnsaturations, by = list(PureMobilizationDataCombined$LikelySubclass), FUN = function(x){min(x,na.rm = TRUE)}),
                                                aggregate(PureMobilizationDataCombined$TotalCarbonChainUnsaturations, by = list(PureMobilizationDataCombined$LikelySubclass), FUN = function(x){max(x,na.rm = TRUE)}), by = "Group.1", all = TRUE), by = "Group.1", all = TRUE), by = "Group.1", all = TRUE)

colnames(OverviewRTsMobilizedLipids) <- c("SubClass", "MinOfRT", "MaxOfRT", "MinOfChainLength", "MaxOfChainLength", "MinUnsat", "MaxUnsat")
OverviewRTsMobilizedLipids <- cbind(RTAmbiguousClasses = c("PG/BMP", "CL", "d*Cer", "d*CerP", "d*HexCer", "d*SHexCer", "d*SM", "DAG", "d*Cer", "d*Cer", "t*Cer", "d*SM", "FA", "FA", "LPC", "LPE", "LPE", "LPG", "PA", "PC", "PC", "PE", "PE", "PG/BMP", "PG/BMP", "PGP", "PI", "PS", "t*Hex2Cer", "t*HexCer", "t*SM", "TAG", "t*Cer", "VA"), OverviewRTsMobilizedLipids)

# Just a back-up of these intermediate data
write.table(OverviewRTsMobilizedLipids, file = "./Output/OverviewRTsMobilizedLipids28062024.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

OverviewRTsMobilizedLipidsasc <- merge(merge(aggregate(OverviewRTsMobilizedLipids$MinOfRT, by = list(OverviewRTsMobilizedLipids$RTAmbiguousClasses), FUN = function(x){min(x, na.rm = TRUE)}),
                                             aggregate(OverviewRTsMobilizedLipids$MaxOfRT, by = list(OverviewRTsMobilizedLipids$RTAmbiguousClasses), FUN = function(x){max(x, na.rm = TRUE)}), by = "Group.1", all = TRUE),
                                       
                                       merge(merge(aggregate(OverviewRTsMobilizedLipids$MinOfChainLength, by = list(OverviewRTsMobilizedLipids$RTAmbiguousClasses), FUN = function(x){min(x, na.rm = TRUE)}),
                                                   aggregate(OverviewRTsMobilizedLipids$MaxOfChainLength, by = list(OverviewRTsMobilizedLipids$RTAmbiguousClasses), FUN = function(x){max(x, na.rm = TRUE)}), by = "Group.1", all = TRUE),
                                             
                                             merge(aggregate(OverviewRTsMobilizedLipids$MinUnsat, by = list(OverviewRTsMobilizedLipids$RTAmbiguousClasses), FUN = function(x){min(x, na.rm = TRUE)}),
                                                   aggregate(OverviewRTsMobilizedLipids$MaxUnsat, by = list(OverviewRTsMobilizedLipids$RTAmbiguousClasses), FUN = function(x){max(x, na.rm = TRUE)}), by = "Group.1", all = TRUE), by = "Group.1", all = TRUE), by = "Group.1", all = TRUE)

colnames(OverviewRTsMobilizedLipidsasc) <- c("RTAmbiguousClasses", "MinOfRT", "MaxOfRT", "MinOfChainLength", "MaxOfChainLength", "MinUnsat", "MaxUnsat")
write.table(OverviewRTsMobilizedLipidsasc, file = "./Output/OverviewRTsMobilizedLipidsasc28062024.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

OverviewRTsMobilizedLipids <- read.csv(file = "./Output/OverviewRTsMobilizedLipids28062024.tsv", header = TRUE, sep = "\t", as.is = TRUE)
OverviewRTsMobilizedLipidsasc <- read.csv(file = "./Output/OverviewRTsMobilizedLipidsasc28062024.tsv", header = TRUE, sep = "\t", as.is = TRUE)

# Assign data to ambiguous retention time classes
AmbiguityConverterRTClasses <- cbind(SubClass = unique(unlist(lapply(InVitroVesicleLipidomes5, function(x){x[,"SubClass"]}))),
                                     
                                     RTAmbiguousClasses = c("d*Cer", "LPE", "DAG", "d*GM4", "PS", "PS", "PC", "t*Cer", "PIP", "PC", 
                                                            "PE", "d*HexCer", "TAG", "PE", "DAG", "TAG", "d*SM", "d*CerP", "LPG", "FA",
                                                            
                                                            "LPA", "PA", "LPS", "LPC", "LPC", "PI", "Hex2DAG", "t*HexCer", "PG/BMP", "PG/BMP",
                                                            "PIP3", "PGP", "d*Gb3", "d*PECer", "LPI", "PI", "t*SM", "t*PECer", "LPE", "HexDAG", 
                                                            
                                                            "d*SHex2Cer", "PA", "d*GM3", "d*GA1", "PIP2", "LPA", "d*Hex2Cer", "d*GM2", "LPS", "SE",
                                                            "d*GA2", "d*SHexCer", "MAG", "LPI", "d*GD3", "d*GM1"))

InVitroVesicleLipidomes7 <- lapply(InVitroVesicleLipidomes5, function(x){cbind(x, RTAmbiguousClasses = AmbiguityConverterRTClasses[match(x[,"SubClass"], AmbiguityConverterRTClasses[,"SubClass"]), "RTAmbiguousClasses"])
})

# Determine if exists in mobilization data and if in correct retention time range
InVitroVesicleLipidomes7 <- lapply(InVitroVesicleLipidomes7, function(x){cbind(x, OverviewRTsMobilizedLipidsasc[match(x[,"RTAmbiguousClasses"], OverviewRTsMobilizedLipidsasc[,"RTAmbiguousClasses"]), 2:7])})

InVitroVesicleLipidomes7[["LiverPos"]] <- cbind(InVitroVesicleLipidomes7[["LiverPos"]], RTIsFit = (InVitroVesicleLipidomes7[["LiverPos"]][,"LIVER_N50_pos.RT.mean"] > InVitroVesicleLipidomes7[["LiverPos"]][,"MinOfRT"]) &
                                                  (InVitroVesicleLipidomes7[["LiverPos"]][,"LIVER_N50_pos.RT.mean"] < InVitroVesicleLipidomes7[["LiverPos"]][,"MaxOfRT"]))

InVitroVesicleLipidomes7[["BrainPos"]] <- cbind(InVitroVesicleLipidomes7[["BrainPos"]], RTIsFit = (InVitroVesicleLipidomes7[["BrainPos"]][,"BRAIN_N50_pos.RT.mean"] > InVitroVesicleLipidomes7[["BrainPos"]][,"MinOfRT"]) &
                                                  (InVitroVesicleLipidomes7[["BrainPos"]][,"BRAIN_N50_pos.RT.mean"] < InVitroVesicleLipidomes7[["BrainPos"]][,"MaxOfRT"]))

InVitroVesicleLipidomes7[["LiverNeg"]] <- cbind(InVitroVesicleLipidomes7[["LiverNeg"]], RTIsFit = (InVitroVesicleLipidomes7[["LiverNeg"]][,"LIVER_N50_neg.RT.mean"] > InVitroVesicleLipidomes7[["LiverNeg"]][,"MinOfRT"]) &
                                                  (InVitroVesicleLipidomes7[["LiverNeg"]][,"LIVER_N50_neg.RT.mean"] < InVitroVesicleLipidomes7[["LiverNeg"]][,"MaxOfRT"]))

InVitroVesicleLipidomes7[["BrainNeg"]] <- cbind(InVitroVesicleLipidomes7[["BrainNeg"]], RTIsFit = (InVitroVesicleLipidomes7[["BrainNeg"]][,"BRAIN_N50_neg.RT.mean"] > InVitroVesicleLipidomes7[["BrainNeg"]][,"MinOfRT"]) &
                                                  (InVitroVesicleLipidomes7[["BrainNeg"]][,"BRAIN_N50_neg.RT.mean"] < InVitroVesicleLipidomes7[["BrainNeg"]][,"MaxOfRT"]))


InVitroVesicleLipidomes7 <- lapply(InVitroVesicleLipidomes7, function(x){do.call("cbind", list(x, 
                                                                                               
                                                                                               FitsInChainLengthRange = (x[,"TotalChainlengths"] >= x[,"MinOfChainLength"]) &
                                                                                                 (x[,"TotalChainlengths"] <= x[,"MaxOfChainLength"]),
                                                                                               
                                                                                               FitsInUnsatRange = (x[,"TotalUnsat"] >= x[,"MinUnsat"]) &
                                                                                                 (x[,"TotalUnsat"] <= x[,"MaxUnsat"])
                                                                                               
))})


InVitroVesicleLipidomes7 <- lapply(InVitroVesicleLipidomes7, function(x){cbind(x, Match = ifelse(x$RTIsFit, "yes", 
                                                                                                 ifelse(!(x$FitsInChainLengthRange & x$FitsInUnsatRange), "maybe", "no")))})

InVitroVesicleLipidomes7[["LiverPos"]] <- do.call("cbind", list(InVitroVesicleLipidomes7[["LiverPos"]], Mode = "Pos", Source = "Liver"))
InVitroVesicleLipidomes7[["BrainPos"]] <- do.call("cbind", list(InVitroVesicleLipidomes7[["BrainPos"]], Mode = "Pos", Source = "Brain"))

InVitroVesicleLipidomes7[["LiverNeg"]] <- do.call("cbind", list(InVitroVesicleLipidomes7[["LiverNeg"]], Mode = "Neg", Source = "Liver"))
InVitroVesicleLipidomes7[["BrainNeg"]] <- do.call("cbind", list(InVitroVesicleLipidomes7[["BrainNeg"]], Mode = "Neg", Source = "Brain"))

for(i in 1:4){colnames(InVitroVesicleLipidomes7[[i]])[4:6] <- c("N50_m.z", "N50_RTMean", "N50_Area")}
InVitroVesicleLipidomes8 <- do.call("rbind", InVitroVesicleLipidomes7)

InVitroVesicleLipidomes8 <- cbind(InVitroVesicleLipidomes8, RTChainlengthSameDirectionOutsideKnownRange = ((InVitroVesicleLipidomes8$N50_RTMean < InVitroVesicleLipidomes8$MinOfRT) & 
                                                                                                             (InVitroVesicleLipidomes8$TotalChainlengths < InVitroVesicleLipidomes8$MinOfChainLength))|
                                    
                                    ((InVitroVesicleLipidomes8$N50_RTMean > InVitroVesicleLipidomes8$MaxOfRT) & 
                                       (InVitroVesicleLipidomes8$TotalChainlengths > InVitroVesicleLipidomes8$MaxOfChainLength)))  

InVitroVesicleLipidomes8 <- cbind(InVitroVesicleLipidomes8, StricterMatch = ifelse(InVitroVesicleLipidomes8$Match != "maybe", InVitroVesicleLipidomes8$Match,
                                                                                   ifelse(InVitroVesicleLipidomes8$RTChainlengthSameDirectionOutsideKnownRange, "maybe", "no")))

# Overview Of Entries: Yes: 1553; Maybe: 659; No: 427; NA: 374


InVitroVesicleLipidomes8 <- cbind(TotalIndex = paste0(InVitroVesicleLipidomes8$Source, InVitroVesicleLipidomes8$Mode, "_", InVitroVesicleLipidomes8$index), InVitroVesicleLipidomes8)
write.table(InVitroVesicleLipidomes8, file = "./Output/InVitroVesicleLipidomes828062024.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

InVitroVesicleLipidomesCleaner <- InVitroVesicleLipidomes8[!is.na(InVitroVesicleLipidomes8$StricterMatch) & InVitroVesicleLipidomes8$StricterMatch != "no",]
write.table(InVitroVesicleLipidomesCleaner, file = "./Output/InVitroVesicleLipidomesCleaner28062024.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

AdductCommonness <- aggregate(InVitroVesicleLipidomesCleaner$SubClass, by = list(InVitroVesicleLipidomesCleaner$Mode, InVitroVesicleLipidomesCleaner$Adduct, InVitroVesicleLipidomesCleaner$SubClass), FUN = length)
colnames(AdductCommonness) <- c("Mode", "Adduct", "SubClass", "Counts")

AdductCommonness <- cbind(AdductCommonness, AdductScore = unlist(lapply(unique(AdductCommonness$SubClass), function(x){AdductCommonness[AdductCommonness$SubClass == x, "Counts"]*100/max(AdductCommonness[AdductCommonness$SubClass == x, "Counts"])})))


InVitroVesicleLipidomesCleaner2 <- cbind(InVitroVesicleLipidomesCleaner, AdductScore = AdductCommonness[match(paste0(InVitroVesicleLipidomesCleaner$Mode, InVitroVesicleLipidomesCleaner$Adduct, InVitroVesicleLipidomesCleaner$SubClass), paste0(AdductCommonness$Mode, AdductCommonness$Adduct, AdductCommonness$SubClass)), "AdductScore"])
InVitroVesicleLipidomesCleaner2 <- cbind(InVitroVesicleLipidomesCleaner2, CombinedScore = InVitroVesicleLipidomesCleaner2$AdductScore + (InVitroVesicleLipidomesCleaner2$StricterMatch == "yes")*50 + (InVitroVesicleLipidomesCleaner2$TotalChainlengths%%2 == 0)*50)

MaxCombinedScoreExtractor <- function(y){y[y[,"CombinedScore"] == max(y[,"CombinedScore"]),]}
ListOfHighScoringLines <- lapply(unique(InVitroVesicleLipidomesCleaner2$TotalIndex), function(x){MaxCombinedScoreExtractor(InVitroVesicleLipidomesCleaner2[InVitroVesicleLipidomesCleaner2$TotalIndex == x,])})


lapply(ListOfHighScoringLines[sapply(ListOfHighScoringLines, function(x){dim(x)[1]}) > 1], function(x){sort(x$SubClass)}) 

InVitroVesicleLipidomesCleaner4 <- do.call("rbind", ListOfHighScoringLines)
write.table(InVitroVesicleLipidomesCleaner4, file = "./Output/InVitroVesicleLipidomesCleaner429062024.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# Make an overview of only ambiguities to help with selection
RemainingAmbiguousLipidIdentifications <- do.call("rbind", ListOfHighScoringLines[sapply(ListOfHighScoringLines, function(x){dim(x)[1]}) > 1])

write.table(RemainingAmbiguousLipidIdentifications, file = "./Output/RemainingAmbiguousLipidIdentifications01072024b.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
RemainingAmbiguousLipidIdentificationsCorrections <- read.csv(file = "./InputData/RemainingAmbiguousLipidIdentifications01072024bCorrected.tsv", header = TRUE, sep = "\t", as.is = TRUE) # This contains manual corrections.

RemainingAmbiguousLipidIdentificationsCorrections[RemainingAmbiguousLipidIdentificationsCorrections$manual.confirmed.annotation == "PC(18:2/19:0)", "manual.confirmed.annotation"] <- "PC(37:2)"
RemainingAmbiguousLipidIdentificationsCorrections2 <- RemainingAmbiguousLipidIdentificationsCorrections[RemainingAmbiguousLipidIdentificationsCorrections$manual.confirmed.annotation != "",]

# Currently the most important columns to update: c("Subclass", "TotalChainlengths", "N50_Area", "Source", "TotalUnsat", "Mode")
# Watch out: most other columns will not be updated further, because not used!

RemainingAmbiguousLipidIdentificationsCorrections2$SubClass <- sapply(strsplit(RemainingAmbiguousLipidIdentificationsCorrections2$manual.confirmed.annotation, "\\("), "[[", 1)
RemainingAmbiguousLipidIdentificationsCorrections2$RTAmbiguousClasses <- sapply(strsplit(RemainingAmbiguousLipidIdentificationsCorrections2$manual.confirmed.annotation, "\\("), "[[", 1)

RemainingAmbiguousLipidIdentificationsCorrections2$TotalChainlengths <- as.numeric(sapply(strsplit(gsub("\\)", "", sapply(strsplit(RemainingAmbiguousLipidIdentificationsCorrections2$manual.confirmed.annotation, "\\("), "[[", 2)), ":"), "[[", 1))
RemainingAmbiguousLipidIdentificationsCorrections2$TotalUnsat <- sapply(strsplit(gsub("\\)", "", sapply(strsplit(RemainingAmbiguousLipidIdentificationsCorrections2$manual.confirmed.annotation, "\\("), "[[", 2)), ":"), "[[", 2)

RemainingAmbiguousLipidIdentificationsCorrections2$Adduct <- paste0("X.", gsub("\\+|-",".", RemainingAmbiguousLipidIdentificationsCorrections2$ion), "..")
RemainingAmbiguousLipidIdentificationsCorrections2$ShortFull <- RemainingAmbiguousLipidIdentificationsCorrections2$manual.confirmed.annotation


# Make an overview of only non-ambiguities to work further with

NonAmbiguousaLipidIdentifications <- do.call("rbind", ListOfHighScoringLines[sapply(ListOfHighScoringLines, function(x){dim(x)[1]}) == 1])
write.table(NonAmbiguousaLipidIdentifications, file = "./Output/NonAmbiguousaLipidIdentifications03072024.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# Extension of non-ambiguous set with the corrected entries
NonAmbiguousaLipidIdentificationsx <- rbind(NonAmbiguousaLipidIdentifications, RemainingAmbiguousLipidIdentificationsCorrections2[,colnames(NonAmbiguousaLipidIdentifications)])

# Limits of chain length ranges
LimitsOnChainLengthRanges <- cbind.data.frame(unique(NonAmbiguousaLipidIdentifications$SubClass), do.call("rbind.data.frame", list(c(26,50), c(38,76), c(26,50), c(26,50), c(26,50), c(26,50), c(26,50), c(26,50),c(38,76), c(26,50), c(26,50), c(26,50), c(26,50), c(26,50), c(26,50), c(0,25), c(0,25), c(0,25), c(26,50), c(26,50), c(26,50), c(26,50), c(0,25), c(26,50), c(26,50), c(26,50), c(0,25), c(26,50), c(0,25))))

colnames(LimitsOnChainLengthRanges) <- c("SubClass", "ChainlengthMin", "ChainlengthMax")
LimitsOnChainLengthRanges <- rbind(LimitsOnChainLengthRanges, c("PG/BMP",26,50))

NonAmbiguousaLipidIdentifications2x <- cbind(NonAmbiguousaLipidIdentificationsx, CommonRange = (NonAmbiguousaLipidIdentificationsx$TotalChainlengths >= LimitsOnChainLengthRanges[match(NonAmbiguousaLipidIdentificationsx$SubClass, LimitsOnChainLengthRanges$SubClass), "ChainlengthMin"]) &
                                               (NonAmbiguousaLipidIdentificationsx$TotalChainlengths <= LimitsOnChainLengthRanges[match(NonAmbiguousaLipidIdentificationsx$SubClass, LimitsOnChainLengthRanges$SubClass), "ChainlengthMax"]))

NonAmbiguousaLipidIdentifications4x <- NonAmbiguousaLipidIdentifications2x[NonAmbiguousaLipidIdentifications2x$CommonRange,]


# Write out to work further on for expansion of figurepanel 2c # Basis for the Extended Data Table 10A
write.table(NonAmbiguousaLipidIdentifications4x, file = "./Output/NonAmbiguousaLipidIdentifications4x09092024.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# Reload data for expansion of figurepanel 2c
NonAmbiguousaLipidIdentifications4x <- read.csv(file = "./Output/NonAmbiguousaLipidIdentifications4x09092024.tsv", header = TRUE, sep = "\t", as.is = TRUE)

# For all lipid subclasses: aggregation of negative and positive modes
AggregatedVesicleLipidomesSubclassesx <- setNames(lapply(unique(NonAmbiguousaLipidIdentifications4x$SubClass), function(x){
  
  aggregate(NonAmbiguousaLipidIdentifications4x[NonAmbiguousaLipidIdentifications4x$SubClass == x,"N50_Area"], 
            by = NonAmbiguousaLipidIdentifications4x[NonAmbiguousaLipidIdentifications4x$SubClass == x,c("Source", "TotalChainlengths", "TotalUnsat")],
            
            FUN = sum)
}), unique(NonAmbiguousaLipidIdentifications4x$SubClass))


AggregatedVesicleLipidomesSubclassesx <- lapply(AggregatedVesicleLipidomesSubclassesx, function(x){cbind(x, LipidChains = paste0(x$TotalChainlengths,":", x$TotalUnsat))})

saveRDS(AggregatedVesicleLipidomesSubclassesx, file="./Output/AggregatedVesicleLipidomesSubclassesx.RData")
AggregatedVesicleLipidomesSubclassesx <- readRDS("./Output/AggregatedVesicleLipidomesSubclassesx.RData")

AggregatedVesicleLipidomesSubclassesMatricesx <- lapply(AggregatedVesicleLipidomesSubclassesx, function(y){if(dim(y)[1] == 1){matrix(y$x, dimnames = list(y$Source,y$LipidChains))}else{as.matrix(NAsToZerosConverter(Col1ToRowNames(dcast(data = y[,c("Source","LipidChains", "x")], formula = y[,"Source"] ~ y[,"LipidChains"], value.var= "x"))))}}) # Done again on other computer
PercentagesVesicleLipidomesSubclassesMatricesOfCombinedTissuesx <- lapply(AggregatedVesicleLipidomesSubclassesMatricesx, function(y){colSums(y)*100/max(colSums(y))})

# Only columns above 5% for mobilized data relevant here
SpeciesColumnsSelectedAbove0.05 <- c("30:0", "30:1", "32:1", "33:1", "34:1", "34:2", "35:2", "36:2", "38:4")

SubselectedListOfPIPCPAInMobilizedRangex <- setNames(lapply(match(c("PI", "PC", "PA"), names(PercentagesVesicleLipidomesSubclassesMatricesOfCombinedTissuesx)), function(y){unlist(lapply(list(NAsToZerosConverter(setNames(PercentagesVesicleLipidomesSubclassesMatricesOfCombinedTissuesx[[y]][SpeciesColumnsSelectedAbove0.05], SpeciesColumnsSelectedAbove0.05))),
                                                                                                                                                                                          function(x){x*100/max(x)}))}), c("PI", "PC", "PA"))

pdf("./Output/VesicleLipidomesTissuesCombinedMobilizationSubsetFocusx04092024.pdf")
for(y in c("PI", "PC", "PA")){
  
  barplot(SubselectedListOfPIPCPAInMobilizedRangex[[y]], las = 2, main = paste0("Vesicle lipidomes: tissues combined & subset focus: ", names(SubselectedListOfPIPCPAInMobilizedRangex[y])), col = "NA", ylab = "Fraction of max. intensity", xlab = "Total carbon chain length and unsaturation of lipids")

  
  }
dev.off() # The above about is used in the construction of the vesicle layer for Fig.6c

# Similar visualization for sphingolipids for Fig.6a
# Later we removed tSM: Is empty after selection and not needed in final visualization.

SpeciesSelectionForSphingolipids <- c("32:1", "32:2", "33:1", "34:0", "34:1", "34:2", "35:1", "36:1", "36:2", "38:1", "38:2", "39:1", "40:1", "40:2", "41:1", "41:2", "42:1", "42:2", "42:3", "44:2", "46:0", "48:0")
HeadgroupSelectionForSphingolipids <- c("dCer", "tCer", "dCerP", "dSM", "dHexCer", "tHexCer", "dSHexCer")

SubselectedListOfSphingolipidsResultsMobilizedRangex <- setNames(lapply(match(HeadgroupSelectionForSphingolipids, names(PercentagesVesicleLipidomesSubclassesMatricesOfCombinedTissuesx)), function(y){unlist(lapply(list(NAsToZerosConverter(setNames(PercentagesVesicleLipidomesSubclassesMatricesOfCombinedTissuesx[[y]][SpeciesSelectionForSphingolipids], SpeciesSelectionForSphingolipids))),
                                                                                                                                                                                                                     function(x){x*100/max(x)}))}), HeadgroupSelectionForSphingolipids)

pdf("./Output/VesicleLipidomesTissuesCombinedMobilizationSubsetFocusOnSphingolipidTransportersx24092024.pdf")
for(y in HeadgroupSelectionForSphingolipids){
  
  barplot(SubselectedListOfSphingolipidsResultsMobilizedRangex[[y]], las = 2, main = paste0("Vesicle lipidomes: tissues combined & subset focus: ", names(SubselectedListOfSphingolipidsResultsMobilizedRangex[y])), col = "NA", ylab = "Fraction of max. intensity", xlab = "Total carbon chain length and unsaturation of lipids")
}

dev.off() # The above about is used in the construction of the vesicle layer for Fig.6a



MatrixOfVesiclesLipidomesx <- do.call("rbind", setNames(lapply(1:length(AggregatedVesicleLipidomesSubclassesx), function(y){
  
  z <- aggregate(x = AggregatedVesicleLipidomesSubclassesx[[y]][,4], list(AggregatedVesicleLipidomesSubclassesx[[y]][,2]), FUN = sum)
  return(NAsToZerosConverter(setNames(setNames(z[,2], z[,1])[as.character(14:74)], as.character(14:74))))}), names(AggregatedVesicleLipidomesSubclassesx)))


ListOfSubclassPartsx <- list()

for(j in 1:dim(MatrixOfVesiclesLipidomesx)[2]){
  z <- aggregate(MatrixOfVesiclesLipidomesx[,j], by = list(c("Cer", "TAG", "PI", "Cer", "PS-O", "PC-O", "TAG-O", "PC", "HexCer", "HexCer", "PE-O", "PA-O", "PI-O", "PE", "DAG", "SM", "DAG-O", "PE", "FA", "PA", "CerP", "SM", "PGP", "SHexCer", "PC", "PE-O", "PG/BMP", "PS", "PG/BMP")), FUN = sum)
  
  ListOfSubclassPartsx[[j]] <- setNames(z[,2], z[,1])}
AggregatedMatrixOfVesiclesLipidomesx <- do.call("cbind", setNames(ListOfSubclassPartsx, colnames(MatrixOfVesiclesLipidomesx)))

MobilizedLipidHeadgroupsLines <- c("Cer", "CerP", "HexCer", "Hex2Cer", "SHexCer", "SM", "FA", "PA", "PC", "PE", "PI", "PS", "PGP", "PG", "PG/BMP", "BMP", "CL", "DAG", "TAG")


AggregatedMatrixOfVesiclesLipidomes2x <- AggregatedMatrixOfVesiclesLipidomesx[rownames(AggregatedMatrixOfVesiclesLipidomesx) %in% MobilizedLipidHeadgroupsLines, as.character(14:72)]
AggregatedMatrixOfVesiclesLipidomes2x <- AggregatedMatrixOfVesiclesLipidomes2x[MobilizedLipidHeadgroupsLines[MobilizedLipidHeadgroupsLines %in% rownames(AggregatedMatrixOfVesiclesLipidomesx)],]

AggregatedMatrixOfVesiclesLipidomes4x <- MinMaxNormMatrixFunc(inputdata = AggregatedMatrixOfVesiclesLipidomes2x)*9+1
AggregatedMatrixOfVesiclesLipidomes4x[is.infinite(as.matrix(AggregatedMatrixOfVesiclesLipidomes4x))] <- 0

# Visualization of only vesicles lipidomes
library(ComplexHeatmap)

library(RColorBrewer)
library(circlize)


col_fung <- colorRamp2(0:10, c("white", colorRampPalette(brewer.pal(9,"Greys")[2:7])(10)))

LegendName <- "Legend"
LegendColor <- col_fung

pdf("./Output/FullVesicleLipidomeWithTissuesAndModesCombinedx04092024.pdf",
    width = unit(23, "mm"), height = unit(14, "mm"))

Heatmap(AggregatedMatrixOfVesiclesLipidomes4x, name = LegendName, col = LegendColor, rect_gp = gpar(type = "none"),  
        width = unit(100*dim(AggregatedMatrixOfVesiclesLipidomes4x)[2]/min(dim(AggregatedMatrixOfVesiclesLipidomes4x)), "mm"), height = unit(100*dim(AggregatedMatrixOfVesiclesLipidomes4x)[1]/min(dim(AggregatedMatrixOfVesiclesLipidomes4x)), "mm"), show_heatmap_legend = FALSE, column_title = "Empty circles: vesicle lipidome", 
     
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.circle(x = x, y = y, r = abs(AggregatedMatrixOfVesiclesLipidomes4x[i, j])/8.5 * min(unit.c(width, height)), gp = gpar(fill = NA, col = "Black", lwd = 2))
          
          grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "LightGrey", fill = NA, lwd = 1))
        }, cluster_rows = FALSE, cluster_columns = FALSE)

dev.off()
# The above is a layer of Fig.5a

# RemainingAmbiguousLipidIdentificationsCorrections <- read.csv(file = "D:/Rest/Rest5/ForUpdatesToFigure2/InVitroVesicleBackgroundFiles/RemainingAmbiguousLipidIdentifications01072024bCorrected.tsv", header = TRUE, sep = "\t", as.is = TRUE)
# Duplicated entry # Can be removed if necessary

# Write subcellular localizations file # This file is the basis for Extended Data Table 9A.
write.table(AggregateMPOfMads0RTLMOCx, file = "./Output/SubcellularLocalizationAveragesLipids30102024.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

ComparisonOfEffectsSetbackWithInterAverage <- apply(do.call("cbind", lapply(list(TestMatrixa2[c(12:19,21),], TestMatrixb[c(12:19,21),], TestMatrixe2[c(12:19,21),], AggregatedMatrixOfVesiclesLipidomes4[7:13,]), 
                                                                            function(x){log10(colSums(10^x)/dim(x)[1])})), 2, DivMaxToPercent)

ComparisonOfEffectsSetback <- apply(do.call("cbind", lapply(list(TestMatrixa2[c(12:19,21),], TestMatrixb[c(12:19,21),], TestMatrixe2[c(12:19,21),], AggregatedMatrixOfVesiclesLipidomes4[7:13,]), 
                                                            function(x){log10(colSums(10^x))})), 2, DivMaxToPercent)

ComparisonOfEffectsColsumsWithInterAverage <- apply(do.call("cbind", lapply(list(TestMatrixa2[c(12:19,21),], TestMatrixb[c(12:19,21),], TestMatrixe2[c(12:19,21),], AggregatedMatrixOfVesiclesLipidomes4[7:13,]), 
                                                                            function(x){colSums(x)/dim(x)[1]})), 2, DivMaxToPercent)

ComparisonOfEffectsColsums <- apply(do.call("cbind", lapply(list(TestMatrixa2[c(12:19,21),], TestMatrixb[c(12:19,21),], TestMatrixe2[c(12:19,21),], AggregatedMatrixOfVesiclesLipidomes4[7:13,]), 
                                                            function(x){colSums(x)})), 2, DivMaxToPercent)

CombinedGPL <- do.call("cbind", list(InCellulo = DivMaxToPercent(colSums(TestMatrixa2[c(12:19,21),])/length(c(12:19,21))),
                                     Cellular = DivMaxToPercent(colSums(TestMatrixb[c(12:19,21),])/length(c(12:19,21))),
                                     
                                     InVitro = DivMaxToPercent(colSums(TestMatrixe2[c(12:19,21),])/length(c(12:19,21))),
                                     Liposomes = DivMaxToPercent(colSums(AggregatedMatrixOfVesiclesLipidomes4[7:13,])/length(7:13))))

InputFocusChainLengths <- cbind(CombinedGPL[as.character(27:51),], t(t(LTPMobilizedGPL2sMaxNormCombinedVersionoes[4,as.character(27:51)])))
colnames(InputFocusChainLengths) <- c("InCellulo", "Cellular", "InVitro", "Liposomes", "SignificantChangesUponOE")


SpeciesOverviewLipidBodyConnectionsFocusOnGPLs <- SpeciesOverviewLipidBodyConnections[as.character(SpeciesOverviewLipidBodyConnections$Head) %in% c("PA", "PC", "PE", "PI", "PS", "PGP", "PG", "PG/BMP", "BMP") & !grepl(pattern = "O-", x = SpeciesOverviewLipidBodyConnections$Body),]

IntensitiesOfObservedVsExpectedUnsatLipidsGPLs <- setNames(Col1ToRowNames(merge(aggregate(PureAntonella32b[as.character(PureAntonella32b$LikelySubclass) %in% c("PA", "PC", "PE", "PI", "PS", "PGP", "PG", "PG/BMP", "BMP"), "Intensity"], by = list(PureAntonella32b[as.character(PureAntonella32b$LikelySubclass) %in% c("PA", "PC", "PE", "PI", "PS", "PGP", "PG", "PG/BMP", "BMP"), "TotalCarbonChainUnsaturations"]), FUN = function(x){sum(x,na.rm=TRUE)}), 
                                                                                aggregate(SpeciesOverviewLipidBodyConnectionsFocusOnGPLs$Amount, by = list(as.numeric(sapply(strsplit(as.character(SpeciesOverviewLipidBodyConnectionsFocusOnGPLs$Body), ":"), "[[", 2))), FUN = function(x){sum(x,na.rm=TRUE)}), 
                                                                                
                                                                                by = 1, all = TRUE)),
                                                           nm = c("Observed", "Expected"))

InVitroSummedIntensitiesUnsaturationsForTheGPLs <- aggregate(PureEnric32[as.character(PureEnric32$LikelySubclass) %in% c("PA", "PC", "PE", "PI", "PS", "PGP", "PG", "PG/BMP", "BMP"),"Intensity"], by = list(PureEnric32[as.character(PureEnric32$LikelySubclass) %in% c("PA", "PC", "PE", "PI", "PS", "PGP", "PG", "PG/BMP", "BMP"),"TotalCarbonChainUnsaturations"]), FUN = function(x){sum(x, na.rm = TRUE)})


IntensitiesUnsatLevelsScreensComparisonForTheGPLs <- setNames(Col1ToRowNames(rbind(merge(IntensitiesOfObservedVsExpectedUnsatLipidsGPLs, InVitroSummedIntensitiesUnsaturationsForTheGPLs, by.x = 0, by.y = 1, all = TRUE), c(11,0,0,0)))[as.character(0:12),c(1,3,2)],
                                                              nm = c("In Cellulo", "In Vitro", "Cellular"))

IntensitiesUnsatLevelsScreensComparisonForTheGPLs[is.na(IntensitiesUnsatLevelsScreensComparisonForTheGPLs)] <- 0
VesiclesSummedIntensitiesUnsaturationsGPLs <- aggregate(NonAmbiguousaLipidIdentifications4x[NonAmbiguousaLipidIdentifications4x$SubClass %in% c("PI", "PC", "PE", "PA", "PGP", "PS", "PG/BMP"), "N50_Area"], by = list(NonAmbiguousaLipidIdentifications4x[NonAmbiguousaLipidIdentifications4x$SubClass %in% c("PI", "PC", "PE", "PA", "PGP", "PS", "PG/BMP"), "TotalUnsat"]), FUN = function(x){sum(x, na.rm = TRUE)})

IntensitiesUnsatLevelsScreensComparisonGPLsx <- NAsToZerosConverter(Col1ToRowNames(merge(IntensitiesUnsatLevelsScreensComparisonForTheGPLs, VesiclesSummedIntensitiesUnsaturationsGPLs, by.x = 0, by.y = 1, all = TRUE))[as.character(0:12), c(1,3,2,4)])     
colnames(IntensitiesUnsatLevelsScreensComparisonGPLsx) <- c("In Cellulo", "Cellular", "In Vitro", "Liposomes")   

IntensitiesUnsatLevelsScreensComparisonGPLsx2 <- cbind(IntensitiesUnsatLevelsScreensComparisonGPLsx, SignificantChangesUponOE = t(t(LTPOEResultsGPL2sOfVariousMaxNormoes[[3]]))) # Fine to use this as it is just normalized to the maximum number, an earlier step could also be used.
IntensitiesUnsatLevelsScreensComparisonNormForTheGPLsByMaxx <- apply(IntensitiesUnsatLevelsScreensComparisonGPLsx2, 2, function(x){x*100/max(x,na.rm=TRUE)})


# Visualization with circles for chainlengths

library(ComplexHeatmap)
CirclesAveragesChainLengthsGLTPs <- t(InputFocusChainLengths/10)

pdf("./Output/AveragesChainLengthsGLTPs15122024.pdf",
    width = unit(23, "mm"), height = unit(10, "mm"))

Heatmap(CirclesAveragesChainLengthsGLTPs, name = LegendName, col = LegendColor, rect_gp = gpar(type = "none"),  
        width = unit(100*dim(CirclesAveragesChainLengthsGLTPs)[2]/min(dim(CirclesAveragesChainLengthsGLTPs)), "mm"), height = unit(100*dim(CirclesAveragesChainLengthsGLTPs)[1]/min(dim(CirclesAveragesChainLengthsGLTPs)), "mm"), show_heatmap_legend = FALSE, column_title = "Averages fatty acid chainlengths of GPLs (min-max normalized)", 
        
        
        cell_fun = function(j, i, x, y, width, height, fill) {
          
          grid.circle(x = x, y = y, r = abs(CirclesAveragesChainLengthsGLTPs[i, j])/4 * min(unit.c(width, height)), gp = gpar(fill = NA, col = "black", lwd = 2.3))
         grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "LightGrey", fill = NA, lwd = 1))
         
        }, cluster_rows = FALSE, cluster_columns = FALSE)
dev.off()


# Visualization with circles for unsaturations

CirclesAveragesUnsaturationsGFLsScreens <- t(IntensitiesUnsatLevelsScreensComparisonNormForTheGPLsByMaxx[,c(1,3)]/10)
CirclesAveragesUnsaturationsGFLsBackgrounds <- t(IntensitiesUnsatLevelsScreensComparisonNormForTheGPLsByMaxx[,c(2,4)]/10)

SignificantChangesOEUnsaturationsGFLs <- t(IntensitiesUnsatLevelsScreensComparisonNormForTheGPLsByMaxx[,5]/10)
library(ComplexHeatmap)

pdf("./Output/AveragesUnsaturationsGPLsScreenvsBackground090120252024.pdf", useDingbats = FALSE,
    width = unit(28, "mm"), height = unit(8, "mm"))

Heatmap(CirclesAveragesUnsaturationsGFLsScreens, name = LegendName, col = LegendColor, rect_gp = gpar(type = "none"),  
        width = unit(100*dim(CirclesAveragesUnsaturationsGFLsScreens)[2]/min(dim(CirclesAveragesUnsaturationsGFLsScreens)), "mm"), height = unit(100*dim(CirclesAveragesUnsaturationsGFLsScreens)[1]/min(dim(CirclesAveragesUnsaturationsGFLsScreens)), "mm"), show_heatmap_legend = FALSE, column_title = "Averages fatty acid chainlengths of GPLs (min-max normalized)", 
        
        
        cell_fun = function(j, i, x, y, width, height, fill) {
          
          grid.circle(x = x, y = y, r = abs(CirclesAveragesUnsaturationsGFLsScreens[i, j])/3.1 * min(unit.c(width, height)), gp = gpar(fill = "lightgrey", col = "lightgrey", lwd = 2.3))
          grid.circle(x = x, y = y, r = abs(CirclesAveragesUnsaturationsGFLsBackgrounds[i, j])/3.1 * min(unit.c(width, height)), gp = gpar(fill = NA, col = "black", lwd = 2.3))
          
          grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "LightGrey", fill = NA, lwd = 1))
        }, cluster_rows = FALSE, cluster_columns = FALSE)

dev.off()
# The previous figure was used as basis for the upper part of Fig.5b_b, after integration with the rest and making the colors congruent in Adobe Illustrator.

pdf("./Output/SignificantChangesOEUnsaturationsGFLs19122024.pdf",
    width = unit(14, "mm"), height = unit(4, "mm"))

heatmap(rbind(SignificantChangesOEUnsaturationsGFLs, c(0:10,10,0)), Colv = NA, Rowv = NA, scale = "none", col = colorRampPalette(c("white","#EFEFEF", "#E9E9E9", "lightgrey", "black"))(256))
heatmap(rbind(SignificantChangesOEUnsaturationsGFLs, c(0:10,10,0)), Colv = NA, Rowv = NA, scale = "none", col = colorRampPalette(c("#440154", "#414487", "#2A788E", "#22A884", "#7AD151", "#DCE230", "#FDE725"))(256))

heatmap(rbind(SignificantChangesOEUnsaturationsGFLs, c(0:10,10,0)), Colv = NA, Rowv = NA, scale = "none", col = colorRampPalette(c("white", "grey", "black"))(256))
dev.off()

# The previous second figure was used as basis for the lower part of Fig.5b_b, after integration with the rest in Adobe Illustrator.
# The figure legend for Fig.5 is also derivable from this second figure, as from other related figures.


# Redone of visualization with circles for the chainlengths

CirclesAveragesChainlengthsGFLsScreens <- CirclesAveragesChainLengthsGLTPs[c(1,3),]
CirclesAveragesChainlengthsGFLsBackgrounds <- CirclesAveragesChainLengthsGLTPs[c(2,4),]

SignificantChangesOEChainlengthsGFLs <- CirclesAveragesChainLengthsGLTPs[5,]
library(ComplexHeatmap)

pdf("./Output/AveragesChainlengthsGPLsScreenvsBackground18122024.pdf",
    width = unit(64, "mm"), height = unit(10, "mm"))

Heatmap(CirclesAveragesChainlengthsGFLsScreens, name = LegendName, col = LegendColor, rect_gp = gpar(type = "none"),  
        width = unit(100*dim(CirclesAveragesChainlengthsGFLsScreens)[2]/min(dim(CirclesAveragesChainlengthsGFLsScreens)), "mm"), height = unit(100*dim(CirclesAveragesChainlengthsGFLsScreens)[1]/min(dim(CirclesAveragesChainlengthsGFLsScreens)), "mm"), show_heatmap_legend = FALSE, column_title = "Averages fatty acid chainlengths of GPLs (min-max normalized)", 
        
        
        cell_fun = function(j, i, x, y, width, height, fill) {
          
          grid.circle(x = x, y = y, r = abs(CirclesAveragesChainlengthsGFLsScreens[i, j])/1.64 * min(unit.c(width, height)), gp = gpar(fill = "lightgrey", col = "lightgrey", lwd = 2.3))
          grid.circle(x = x, y = y, r = abs(CirclesAveragesChainlengthsGFLsBackgrounds[i, j])/1.64 * min(unit.c(width, height)), gp = gpar(fill = NA, col = "black", lwd = 2.3))
          
          grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "LightGrey", fill = NA, lwd = 1))
        }, cluster_rows = FALSE, cluster_columns = FALSE)

dev.off()
# The previous figure was used as basis for the upper part of Fig.5b_a, after integration with the rest and making the colors congruent in Adobe Illustrator.

pdf("./Output/SignificantChangesOEChainlengthsGFLs19122024.pdf",
    width = unit(64, "mm"), height = unit(4, "mm"))

heatmap(rbind(SignificantChangesOEChainlengthsGFLs, NA), Colv = NA, Rowv = NA, scale = "none", col = colorRampPalette(c("white","#EFEFEF", "#E9E9E9", "lightgrey", "black"))(256))
heatmap(rbind(SignificantChangesOEChainlengthsGFLs, NA), Colv = NA, Rowv = NA, scale = "none", col = colorRampPalette(c("#440154", "#414487", "#2A788E", "#22A884", "#7AD151", "#DCE230", "#FDE725"))(256))

heatmap(rbind(SignificantChangesOEChainlengthsGFLs, NA), Colv = NA, Rowv = NA, scale = "none", col = colorRampPalette(c("white", "grey", "black"))(256))
dev.off()

# The previous second figure was used as basis for the lower part of Fig.5b_a, after integration with the rest in Adobe Illustrator.
# The figure legend for Fig.5 is also derivable from this second figure, as from other related figures.

library(reshape2) #! Package dependency that is removable
NoveltyHits <- melt(LipidClassLiteratureDataSetslc4hdr)

NoveltyHits2 <- NoveltyHits[NoveltyHits[,3],1:2] # Double-checked manually too: all good.
colnames(NoveltyHits2) <- c("NovelLipid", "LTP")

write.table(NoveltyHits2, file = "./Output/NoveltyHits20122024.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
ExtendedDataTable1Basis <- read.csv(file = "./InputData/Titeca_ExtendedDataTable1_v1(1).csv", header = TRUE, sep = "\t", as.is = TRUE)

# Clean up the table
CleanedExtendedDataTable1Basis <- do.call("rbind", strsplit(x = gsub(";;;;", "", ExtendedDataTable1Basis[,1]), split = "\t"))

colnames(CleanedExtendedDataTable1Basis) <- unlist(strsplit(gsub("\\.\\.\\.\\.","", colnames(ExtendedDataTable1Basis)), "\\."))
CleanedExtendedDataTable1Basis[,"LikelySubclass"] <- gsub(";", "", CleanedExtendedDataTable1Basis[,"LikelySubclass"])

LipidNameConverterHigherSubclasses <- cbind(StandardLevel = sort(unique(as.character(CleanedExtendedDataTable1Basis[,"LikelySubclass"]))),
                                            HigherLevel = c("BMP", "CL", "Cer*", "d*CerP", "HexCer*", "d*SHexCer",
                                                            
                                                            "SM*", "DAG", "Cer*", "Cer*", "Cer*", "SM*",
                                                            "FA", "FAL", "LPC", "LPE", "LPE-O", "LPG",
                                                            
                                                            "PA", "PC", "PC-O", "PE", "PE-O", "PG",
                                                            "PG/BMP", "PGP", "PI", "PS", "t*Hex2Cer", "HexCer*",
                                                            
                                                            "SM*", "TAG", "Cer*", "VA"))
# Matching of sphingolipids corrected

CleanedExtendedDataTable1Basis2 <- cbind(CleanedExtendedDataTable1Basis, NovelLinkLTPLipidSubclass = rowSums(do.call("cbind", lapply(1:dim(NoveltyHits2)[1], function(x){(CleanedExtendedDataTable1Basis[,"LTPProtein"] == as.character(NoveltyHits2[x,2])) & (LipidNameConverterHigherSubclasses[match(CleanedExtendedDataTable1Basis[,"LikelySubclass"], LipidNameConverterHigherSubclasses[,"StandardLevel"]), "HigherLevel"] == as.character(NoveltyHits2[x,1]))
}))))

write.table(CleanedExtendedDataTable1Basis2, file = "./Output/ExtendedDataTable1Basis28012025Corrected.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
# The above table is the basis for further analyses that resulted in parts of Extended Data Table 7 and Fig.1b&c

####
# Find and visualize which fraction of novel and known lipid subclass - LTP pairs is also changed in the overexpression data

# Converter for the overexpression data
x <- cbind(as.character(LogRatioDifferenceTTestsForLipidSpeciesb$Lipid), do.call("rbind", strsplit(as.character(LogRatioDifferenceTTestsForLipidSpeciesb$Lipid), split = " ")))

# O- ==> -O
x[,2] <- paste0(x[,2], ifelse(grepl(pattern = "O-", x[,3]), "-O", ""))

# ;2 ==> d*
x[,2] <- paste0(ifelse(grepl(pattern = ";2", x[,3]), "d*", ""), x[,2])

# ;3 ==> t*
x[,2] <- paste0(ifelse(grepl(pattern = ";3", x[,3]), "t*", ""), x[,2])

# ;4 ==> tOH*
x[,2] <- paste0(ifelse(grepl(pattern = ";4", x[,3]), "tOH*", ""), x[,2])

LogRatioDifferenceTTestsForLipidSpeciesbwls <- cbind(LogRatioDifferenceTTestsForLipidSpeciesb, LipidSubclass = as.character(x[,2]))


LTPNamesConvertingMatrix <- cbind(OriginalLTPs = sort(unique(as.character(LogRatioDifferenceTTestsForLipidSpeciesbwls[,"LTPProtein"]))),
                                  NewLTPs = c("ATCAY", "BNIPL", "BPI", "BPIFB2", "CRABP2", "FABP1", "FABP5", "FABP7", "GLTP", "GLTPD1", "GM2A", "HSDL2", "LCN1", "LCN15",
                                              
                                              "OSBPL1A", "OSBPL10", "OSBPL11", "OSBPL2", "OSBP2", "OSBPL5", "OSBPL7", "OSBPL8", "OSBPL9", "OSBP", "PITPNA", "PITPNB", "PITPNC1", "PMP2", "RBP1", "RBP4", "RBP5", "RLBP1", "SCP2", "SCP2D1", "SEC14L2", "SEC14L4",
                                              "SEC14L5", "SEC14L6", "STARD10", "STARD11", "STARD2", "TTPA", "TTPAL"))

LogRatioDifferenceTTestsForLipidSpeciesbwls <- cbind(LogRatioDifferenceTTestsForLipidSpeciesbwls, 
                                                     NewLTPs = LTPNamesConvertingMatrix[match(LogRatioDifferenceTTestsForLipidSpeciesbwls[,"LTPProtein"], LTPNamesConvertingMatrix[,"OriginalLTPs"]), "NewLTPs"])

LipidNamesOEConverter <- cbind(OENamesOfLipids = sort(unique(as.character(LogRatioDifferenceTTestsForLipidSpeciesbwls[,"LipidSubclass"]))),
                               ConsensusNamesOfLipids = c("CE", "ST", "CL", "Cer*", "HexCer*", "SM*", "DAG", "LPA", "LPC", "LPC-O", "LPE", "LPE-O", "LPG", "LPI", "LPS",
                                                          
                                                          "PA", "PC", "PC-O", "PE", "PE-O", "PG", "PI", "PS", "Cer*", "HexCer*", "SM*", "TAG", "Cer*",  "SM*"))
LogRatioDifferenceTTestsForLipidSpeciesbwls <- cbind(LogRatioDifferenceTTestsForLipidSpeciesbwls,
                                                     
                                                     ConvertedOELipidName = LipidNamesOEConverter[match(LogRatioDifferenceTTestsForLipidSpeciesbwls[,"LipidSubclass"], LipidNamesOEConverter[,"OENamesOfLipids"]), "ConsensusNamesOfLipids"])


BroadSubsetPairsUponOverexpression2 <- unique(LogRatioDifferenceTTestsForLipidSpeciesbwls[LogRatioDifferenceTTestsForLipidSpeciesbwls$p.value <= 0.05/1989, c("NewLTPs", "ConvertedOELipidName")])
NarrowSubsetPairsUponOverexpression2 <- unique(LogRatioDifferenceTTestsForLipidSpeciesbwls[LogRatioDifferenceTTestsForLipidSpeciesbwls$p.value <= 0.05/dim(LogRatioDifferenceTTestsForLipidSpeciesbwls)[1], c("NewLTPs", "ConvertedOELipidName")])

LipidNamesScreenConverter2 <- cbind(ScreenNamesOfLipids = sort(rownames(InVivoDataSetslc4hdr)),
                                    ConsensusNamesOfLipids = c("PG", "Cer*", "CL", "d*CerP", "d*SHexCer", "DAG", "FA", "FAL", "HexCer*", "LPC", "LPE", "LPE-O", "LPG",
                                                               
                                                               "PA", "PC", "PC-O", "PE", "PE-O", "PG", "PG", "PGP", "PI", "PIPs", "PS", "SM*", "ST",
                                                               "t*Hex2Cer", "TAG", "VA"))

library(reshape2) #! Package dependency that is removable


LipidClassLiteratureDataSetslc4hdrCHConvertedToSterolName <- LipidClassLiteratureDataSetslc4hdr
rownames(LipidClassLiteratureDataSetslc4hdrCHConvertedToSterolName)[rownames(LipidClassLiteratureDataSetslc4hdrCHConvertedToSterolName) == "CH"] <- "Sterol"

CombinedViewOfScreens <- merge(merge(melt(InVivoDataSetslc4hdr[,!(colnames(InVivoDataSetslc4hdr) %in% c("OSBPL7", "OSBPL8", "OSBPL10", "OSBPL11"))]), 
                                     melt(InVitroDataSetslc4hdr[,!(colnames(InVitroDataSetslc4hdr) %in% c("OSBPL7", "OSBPL8", "OSBPL10", "OSBPL11"))]),
                                     
                                     by = c("Var1","Var2")), 
                               melt(LipidClassLiteratureDataSetslc4hdrCHConvertedToSterolName[,!(colnames(LipidClassLiteratureDataSetslc4hdrCHConvertedToSterolName) %in% c("OSBPL7", "OSBPL8", "OSBPL10", "OSBPL11"))]), 
                               
                               by = c("Var1","Var2"))
colnames(CombinedViewOfScreens) <- c("Lipid","LTPProtein", "InCellulo", "InVitro", "Novelty")

CombinedViewOfScreens <- cbind(CombinedViewOfScreens, 
                               ConsensusNamesOfLipids = LipidNamesScreenConverter2[match(CombinedViewOfScreens[,"Lipid"], LipidNamesScreenConverter2[,"ScreenNamesOfLipids"]), "ConsensusNamesOfLipids"])

CombinedViewOfScreens <- cbind(CombinedViewOfScreens, 
                               LipidInOEData = CombinedViewOfScreens[,"ConsensusNamesOfLipids"] %in% LipidNamesOEConverter[,"ConsensusNamesOfLipids"])

# Combination
CombinedViewOfScreens2 <- do.call("cbind", list(CombinedViewOfScreens, 
                                                
                                                MatchesNarrowSignificantOEs = paste(CombinedViewOfScreens[,"LTPProtein"], CombinedViewOfScreens[,"ConsensusNamesOfLipids"], sep = "_") %in%
                                                  paste(as.character(NarrowSubsetPairsUponOverexpression2[,"NewLTPs"]), as.character(NarrowSubsetPairsUponOverexpression2[,"ConvertedOELipidName"]), sep = "_"),
                                                
                                                MatchesBroadSignificantOEs = paste(CombinedViewOfScreens[,"LTPProtein"], CombinedViewOfScreens[,"ConsensusNamesOfLipids"], sep = "_") %in%
                                                  paste(as.character(BroadSubsetPairsUponOverexpression2[,"NewLTPs"]), as.character(BroadSubsetPairsUponOverexpression2[,"ConvertedOELipidName"]), sep = "_")))

CombinedViewOfScreens2 <- cbind(CombinedViewOfScreens2, LTPProteinInOEData = CombinedViewOfScreens2[,"LTPProtein"] %in% LogRatioDifferenceTTestsForLipidSpeciesbwls[,"NewLTPs"])
# Mainly several ORPs not matching --> correct the matching. # Done: Corrected!

# Correct the PS-OSBPL9 that was kicked out in vitro, because we went for high certainty identifications only and it previously could be interpreted in several ways
CombinedViewOfScreens2[(CombinedViewOfScreens2[,"Lipid"] == "PS") & (CombinedViewOfScreens2[,"LTPProtein"] == "OSBPL9"), "InVitro"] <- NA

# How many LTP-lipid pairs are observed in the screens in cellulo
ScreenMatchingOEOverview <- do.call("cbind", list(rbind(c(ScreenPresence = sum(!is.na(CombinedViewOfScreens2$InCellulo)), # 43
                                                          
                                                          RestScreen = sum(is.na(CombinedViewOfScreens2$InCellulo)), # 1161
                                                          PercentPresenceInScreen = sum(!is.na(CombinedViewOfScreens2$InCellulo))*100/length(CombinedViewOfScreens2$InCellulo)), # 3.571429
                                                        
                                                        # How many LTP-lipid pairs are observed in the screens in vitro
                                                        c(ScreenPresence = sum(!is.na(CombinedViewOfScreens2$InVitro)), # 75
                                                          
                                                          RestScreen = sum(is.na(CombinedViewOfScreens2$InVitro)), # 1129
                                                          PercentPresenceInScreen = sum(!is.na(CombinedViewOfScreens2$InVitro))*100/length(CombinedViewOfScreens2$InVitro))), # 6.229236
                                                  
                                                  # How many of the observed LTP-lipid pairs are novel in the screens in cellulo
                                                  rbind(c(NovelPairsInScreen = sum(!is.na(CombinedViewOfScreens2$InCellulo) & CombinedViewOfScreens2$Novelty), # 29
                                                          
                                                          KnownPairsInScreen = sum(!is.na(CombinedViewOfScreens2$InCellulo) & !CombinedViewOfScreens2$Novelty), # 14
                                                          PercentNovelPairsInScreen = sum(!is.na(CombinedViewOfScreens2$InCellulo) & CombinedViewOfScreens2$Novelty)*100/sum(!is.na(CombinedViewOfScreens2$InCellulo))), # 67.44186
                                                        
                                                        # How many of the observed LTP-lipid pairs are novel in the screens in vitro
                                                        c(NovelPairsInScreen = sum(!is.na(CombinedViewOfScreens2$InVitro) & CombinedViewOfScreens2$Novelty), # 61
                                                          
                                                          KnownPairsInScreen = sum(!is.na(CombinedViewOfScreens2$InVitro) & !CombinedViewOfScreens2$Novelty), # 14
                                                          PercentNovelPairsInScreen = sum(!is.na(CombinedViewOfScreens2$InVitro) & CombinedViewOfScreens2$Novelty)*100/sum(!is.na(CombinedViewOfScreens2$InVitro)))), # 81.33333
                                                  
                                                  ##### How many of the observed LTP-lipid pairs had lipids that were identifiable in the overexpression data
                                                  # For novelty:
                                                  
                                                  rbind(c(NovelMatchableInOE = sum(!is.na(CombinedViewOfScreens2$InCellulo) & CombinedViewOfScreens2$Novelty & CombinedViewOfScreens2$LipidInOEData), # 28 OK
                                                          NovelNonmatchableInOE = sum(!is.na(CombinedViewOfScreens2$InCellulo) & CombinedViewOfScreens2$Novelty) - sum(!is.na(CombinedViewOfScreens2$InCellulo) & CombinedViewOfScreens2$Novelty & CombinedViewOfScreens2$LipidInOEData), # 1
                                                          
                                                          PercentageNovelMatchableOE = sum(!is.na(CombinedViewOfScreens2$InCellulo) & CombinedViewOfScreens2$Novelty & CombinedViewOfScreens2$LipidInOEData)*100/sum(!is.na(CombinedViewOfScreens2$InCellulo) & CombinedViewOfScreens2$Novelty)), # 96.55172
                                                        
                                                        
                                                        c(NovelMatchableInOE = sum(!is.na(CombinedViewOfScreens2$InVitro) & CombinedViewOfScreens2$Novelty & CombinedViewOfScreens2$LipidInOEData), # 50 OK
                                                          NovelNonmatchableInOE = sum(!is.na(CombinedViewOfScreens2$InVitro) & CombinedViewOfScreens2$Novelty) - sum(!is.na(CombinedViewOfScreens2$InVitro) & CombinedViewOfScreens2$Novelty & CombinedViewOfScreens2$LipidInOEData), # 11
                                                          
                                                          PercentageNovelMatchableOE = sum(!is.na(CombinedViewOfScreens2$InVitro) & CombinedViewOfScreens2$Novelty & CombinedViewOfScreens2$LipidInOEData)*100/sum(!is.na(CombinedViewOfScreens2$InVitro) & CombinedViewOfScreens2$Novelty))), # 81.96721
                                                  # For non-novelty:
                                                  
                                                  rbind(c(KnownMatchableInOE = sum(!is.na(CombinedViewOfScreens2$InCellulo) & !CombinedViewOfScreens2$Novelty & CombinedViewOfScreens2$LipidInOEData), # 10 OK
                                                          KnownNonmatchableInOE = sum(!is.na(CombinedViewOfScreens2$InCellulo) & !CombinedViewOfScreens2$Novelty) - sum(!is.na(CombinedViewOfScreens2$InCellulo) & !CombinedViewOfScreens2$Novelty & CombinedViewOfScreens2$LipidInOEData), # 4
                                                          
                                                          PercentageKnownMatchableOE = sum(!is.na(CombinedViewOfScreens2$InCellulo) & !CombinedViewOfScreens2$Novelty & CombinedViewOfScreens2$LipidInOEData)*100/sum(!is.na(CombinedViewOfScreens2$InCellulo) & !CombinedViewOfScreens2$Novelty)), # 71.42857
                                                        
                                                        
                                                        c(KnownMatchableInOE = sum(!is.na(CombinedViewOfScreens2$InVitro) & !CombinedViewOfScreens2$Novelty & CombinedViewOfScreens2$LipidInOEData), # 10 OK
                                                          KnownNonmatchableInOE = sum(!is.na(CombinedViewOfScreens2$InVitro) & !CombinedViewOfScreens2$Novelty) - sum(!is.na(CombinedViewOfScreens2$InVitro) & !CombinedViewOfScreens2$Novelty & CombinedViewOfScreens2$LipidInOEData), # 4
                                                          
                                                          PercentageKnownMatchableOE = sum(!is.na(CombinedViewOfScreens2$InVitro) & !CombinedViewOfScreens2$Novelty & CombinedViewOfScreens2$LipidInOEData)*100/sum(!is.na(CombinedViewOfScreens2$InVitro) & !CombinedViewOfScreens2$Novelty))), # 71.42857
                                                  
                                                  
                                                  # How many of the observed LTP-lipid pairs are observed to have significantly changed in the overexpression data
                                                  # For novelty:
                                                  
                                                  rbind(c(NovelMatchedInOE = sum(!is.na(CombinedViewOfScreens2$InCellulo) & CombinedViewOfScreens2$Novelty & CombinedViewOfScreens2$LipidInOEData & CombinedViewOfScreens2$MatchesBroadSignificantOEs), # 19
                                                          NovelUnmatchedInOE = sum(!is.na(CombinedViewOfScreens2$InCellulo) & CombinedViewOfScreens2$Novelty & CombinedViewOfScreens2$LipidInOEData) - sum(!is.na(CombinedViewOfScreens2$InCellulo) & CombinedViewOfScreens2$Novelty & CombinedViewOfScreens2$LipidInOEData & CombinedViewOfScreens2$MatchesBroadSignificantOEs), # 9
                                                          
                                                          PercentageNovelMatchedOE = sum(!is.na(CombinedViewOfScreens2$InCellulo) & CombinedViewOfScreens2$Novelty & CombinedViewOfScreens2$LipidInOEData & CombinedViewOfScreens2$MatchesBroadSignificantOEs)*100/sum(!is.na(CombinedViewOfScreens2$InCellulo) & CombinedViewOfScreens2$Novelty & CombinedViewOfScreens2$LipidInOEData)), # 67.85714
                                                        
                                                        
                                                        c(NovelMatchedInOE = sum(!is.na(CombinedViewOfScreens2$InVitro) & CombinedViewOfScreens2$Novelty & CombinedViewOfScreens2$LipidInOEData & CombinedViewOfScreens2$MatchesBroadSignificantOEs), # 17
                                                          NovelUnmatchedInOE = sum(!is.na(CombinedViewOfScreens2$InVitro) & CombinedViewOfScreens2$Novelty & CombinedViewOfScreens2$LipidInOEData) - sum(!is.na(CombinedViewOfScreens2$InVitro) & CombinedViewOfScreens2$Novelty & CombinedViewOfScreens2$LipidInOEData & CombinedViewOfScreens2$MatchesBroadSignificantOEs),
                                                          
                                                          PercentageNovelMatchedOE = sum(!is.na(CombinedViewOfScreens2$InVitro) & CombinedViewOfScreens2$Novelty & CombinedViewOfScreens2$LipidInOEData & CombinedViewOfScreens2$MatchesBroadSignificantOEs)*100/sum(!is.na(CombinedViewOfScreens2$InVitro) & CombinedViewOfScreens2$Novelty & CombinedViewOfScreens2$LipidInOEData))), # 34
                                                  # For non-novelty:
                                                  
                                                  rbind(c(KnownMatchedInOE = sum(!is.na(CombinedViewOfScreens2$InCellulo) & !CombinedViewOfScreens2$Novelty & CombinedViewOfScreens2$LipidInOEData & CombinedViewOfScreens2$MatchesBroadSignificantOEs), # 6
                                                          KnownUnmatchedInOE = sum(!is.na(CombinedViewOfScreens2$InCellulo) & !CombinedViewOfScreens2$Novelty & CombinedViewOfScreens2$LipidInOEData) - sum(!is.na(CombinedViewOfScreens2$InCellulo) & !CombinedViewOfScreens2$Novelty & CombinedViewOfScreens2$LipidInOEData & CombinedViewOfScreens2$MatchesBroadSignificantOEs), # 4
                                                          
                                                          PercentageKnownMatchedOE = sum(!is.na(CombinedViewOfScreens2$InCellulo) & !CombinedViewOfScreens2$Novelty & CombinedViewOfScreens2$LipidInOEData & CombinedViewOfScreens2$MatchesBroadSignificantOEs)*100/sum(!is.na(CombinedViewOfScreens2$InCellulo) & !CombinedViewOfScreens2$Novelty & CombinedViewOfScreens2$LipidInOEData)), # 60
                                                        
                                                        
                                                        c(KnownMatchedInOE = sum(!is.na(CombinedViewOfScreens2$InVitro) & !CombinedViewOfScreens2$Novelty & CombinedViewOfScreens2$LipidInOEData & CombinedViewOfScreens2$MatchesBroadSignificantOEs), # 4
                                                          KnownUnmatchedInOE = sum(!is.na(CombinedViewOfScreens2$InVitro) & !CombinedViewOfScreens2$Novelty & CombinedViewOfScreens2$LipidInOEData) - sum(!is.na(CombinedViewOfScreens2$InVitro) & !CombinedViewOfScreens2$Novelty & CombinedViewOfScreens2$LipidInOEData & CombinedViewOfScreens2$MatchesBroadSignificantOEs), #6
                                                          
                                                          PercentageKnownMatchedOE = sum(!is.na(CombinedViewOfScreens2$InVitro) & !CombinedViewOfScreens2$Novelty & CombinedViewOfScreens2$LipidInOEData & CombinedViewOfScreens2$MatchesBroadSignificantOEs)*100/sum(!is.na(CombinedViewOfScreens2$InVitro) & !CombinedViewOfScreens2$Novelty & CombinedViewOfScreens2$LipidInOEData))))) # 40
rownames(ScreenMatchingOEOverview) <- c("in cellulo", "in vitro")

ScreenMatchingOEFocus <- ScreenMatchingOEOverview[,c("NovelMatchedInOE", "NovelUnmatchedInOE", "NovelNonmatchableInOE", "KnownMatchedInOE", "KnownUnmatchedInOE", "KnownNonmatchableInOE")]
rownames(ScreenMatchingOEFocus) <- c("in cellulo", "in vitro")

write.table(ScreenMatchingOEOverview, file = "./Output/ScreenMatchingOEOverview28012025.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE) # This file is the basis for part of the data used in figure 1d
write.table(ScreenMatchingOEFocus, file = "./Output/ScreenMatchingOEFocus28012025.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# The output above forms the basis for several parts of Fig.1d. Most of the visualization and integration of these numbers was done in Illustrator.


####
SplitVersionOfOEs <- strsplit(as.character(LogRatioDifferenceTTestsForLipidSpeciesbwls[,"Lipid"]), split = " ")

SplitVersionOfOEsp <- strsplit(gsub("O-", "", lapply(SplitVersionOfOEs, function(x){if(length(x)==2){x[2]}else{"0:0;0_0:0;0"}})), split = "_|/")
y <- SplitVersionOfOEsp[[901]]

LogRatioDifferenceTTestsForLipidSpeciesbwlschu <- cbind(LogRatioDifferenceTTestsForLipidSpeciesbwls, do.call("rbind", lapply(SplitVersionOfOEsp, function(y){c(sum(as.numeric(sapply(y, function(x){strsplit(strsplit(x, split = ";")[[1]][1], split = ":")[[1]][1]}))),
                                                                                                                                                               sum(as.numeric(sapply(y, function(x){strsplit(strsplit(x, split = ";")[[1]][1], split = ":")[[1]][2]}))))})))

colnames(LogRatioDifferenceTTestsForLipidSpeciesbwlschu)[13:14] <- c("TotalChainlength", "TotalUnsaturation")
apply(PureMobilizationDataCombined[,c("LTPProtein", "LikelySubclass", "TotalCarbonChainLength", "TotalCarbonChainUnsaturations")] , 1 , paste0 , collapse = "_", sep = "")

y <- paste(PureMobilizationDataCombined[,"LTPProtein"], PureMobilizationDataCombined[,"LikelySubclass"], as.numeric(PureMobilizationDataCombined[,"TotalCarbonChainLength"]), as.numeric(PureMobilizationDataCombined[,"TotalCarbonChainUnsaturations"]), sep = "_")
x <- paste(as.character(LogRatioDifferenceTTestsForLipidSpeciesbwlschu[,"NewLTPs"]), as.character(LogRatioDifferenceTTestsForLipidSpeciesbwlschu[,"LipidSubclass"]), as.numeric(LogRatioDifferenceTTestsForLipidSpeciesbwlschu[,"TotalChainlength"]), as.numeric(LogRatioDifferenceTTestsForLipidSpeciesbwlschu[,"TotalUnsaturation"]), sep = "_")


MappingSpeciesScreensOE <- lapply(1:length(y), function(t){if(any(x %in% y[t])){cbind(PureMobilizationDataCombined[t,], LogRatioDifferenceTTestsForLipidSpeciesbwlschu[x %in% y[t],])}})

sum(sapply(MappingSpeciesScreensOE, is.null)) # 270
length(MappingSpeciesScreensOE) # 1137

sum(sapply(MappingSpeciesScreensOE, is.null))/length(MappingSpeciesScreensOE) # 0.237467
SubclassMatcher <- sort(as.character(unique(PureMobilizationDataCombined[,"LikelySubclass"])))

SubclassMatcher <- cbind(SubclassMatcher, SubclassMatcher)
colnames(SubclassMatcher) <- c("Screen", "OE")

SubclassMatcher[c(1, 9:12, 25, 33),2] <- c("PG", "d*Cer", "d*Cer", "t*Cer", "d*SM", "PG", "t*Cer")
y2 <- paste(PureMobilizationDataCombined[,"LTPProtein"], SubclassMatcher[match(PureMobilizationDataCombined[,"LikelySubclass"], SubclassMatcher[,1]),2], as.numeric(PureMobilizationDataCombined[,"TotalCarbonChainLength"]), as.numeric(PureMobilizationDataCombined[,"TotalCarbonChainUnsaturations"]), sep = "_")

# Redo previous with y2 with corrected sub-classes
MappingSpeciesScreensOE2 <- lapply(1:length(y2), function(t){if(any(x %in% y2[t])){cbind(PureMobilizationDataCombined[t,], LogRatioDifferenceTTestsForLipidSpeciesbwlschu[x %in% y2[t],])}})

sum(sapply(MappingSpeciesScreensOE2, is.null)) # 226
length(MappingSpeciesScreensOE2) # 1137

sum(sapply(MappingSpeciesScreensOE2, is.null))/length(MappingSpeciesScreensOE2) # 0.1987687
1 - sum(sapply(MappingSpeciesScreensOE2, is.null))/length(MappingSpeciesScreensOE2) # 0.8012313

MappingSpeciesScreensOEm <- do.call("rbind", MappingSpeciesScreensOE2)



# Fraction of overlap between screens and overexpression data being significant based on overlap-size

sum(as.numeric(MappingSpeciesScreensOEm$p.value) <= 0.05/dim(MappingSpeciesScreensOEm)[1]) # 28
dim(MappingSpeciesScreensOEm)[1] # 4164

sum(as.numeric(MappingSpeciesScreensOEm$p.value) <= 0.05/dim(MappingSpeciesScreensOEm)[1])/dim(MappingSpeciesScreensOEm)[1] # 0.006724304



# Fraction of overexpression data being significant based on size of the overexpression data

sum(as.numeric(LogRatioDifferenceTTestsForLipidSpeciesbwlschu$p.value) <= 0.05/dim(LogRatioDifferenceTTestsForLipidSpeciesbwlschu)[1]) # 422
dim(LogRatioDifferenceTTestsForLipidSpeciesbwlschu)[1] # 85527

sum(as.numeric(LogRatioDifferenceTTestsForLipidSpeciesbwlschu$p.value) <= 0.05/dim(LogRatioDifferenceTTestsForLipidSpeciesbwlschu)[1])/dim(LogRatioDifferenceTTestsForLipidSpeciesbwlschu)[1] # 0.004934114



# Fraction of overexpression data being significant based on size of the overlap data, to keep p-value cut-off similar for both

sum(as.numeric(LogRatioDifferenceTTestsForLipidSpeciesbwlschu$p.value) <= 0.05/dim(MappingSpeciesScreensOEm)[1]) # 463
sum(as.numeric(LogRatioDifferenceTTestsForLipidSpeciesbwlschu$p.value) <= 0.05/dim(MappingSpeciesScreensOEm)[1])/dim(LogRatioDifferenceTTestsForLipidSpeciesbwlschu)[1] # 0.005413495


# Fraction of overlap data being significant based on size of the overexpression data, to keep p-value cut-off similar for both

sum(as.numeric(MappingSpeciesScreensOEm$p.value) <= 0.05/dim(LogRatioDifferenceTTestsForLipidSpeciesbwlschu)[1]) # 26
sum(as.numeric(MappingSpeciesScreensOEm$p.value) <= 0.05/dim(LogRatioDifferenceTTestsForLipidSpeciesbwlschu)[1])/dim(MappingSpeciesScreensOEm)[1] # 0.006243996

sum(LogRatioDifferenceTTestsForLipidSpeciesbwlschu$NewLTPs == "ATCAY") # Size of single LTP overexpression results section: 1989
# Fraction of overlap data being significant based on size of amount of lipids observed in individual LTP overexpression, to keep p-value cut-off similar for both

sum(as.numeric(MappingSpeciesScreensOEm$p.value) <= 0.05/1989) # 28
sum(as.numeric(MappingSpeciesScreensOEm$p.value) <= 0.05/1989)/dim(MappingSpeciesScreensOEm)[1] # 0.006724304


# Fraction of overexpression data being significant based on size of amount of lipids observed in individual LTP overexpression, to keep p-value cut-off similar for both

sum(as.numeric(LogRatioDifferenceTTestsForLipidSpeciesbwlschu$p.value) <= 0.05/1989) # 475
sum(as.numeric(LogRatioDifferenceTTestsForLipidSpeciesbwlschu$p.value) <= 0.05/1989)/dim(LogRatioDifferenceTTestsForLipidSpeciesbwlschu)[1] # 0.005553802


# Fraction of overlap data being significant without correction for multiple testing

sum(as.numeric(MappingSpeciesScreensOEm$p.value) <= 0.05) # 174
sum(as.numeric(MappingSpeciesScreensOEm$p.value) <= 0.05)/dim(MappingSpeciesScreensOEm)[1] # 0.04178674


# Fraction of overexpression data being significant without correction for multiple testing

sum(as.numeric(LogRatioDifferenceTTestsForLipidSpeciesbwlschu$p.value) <= 0.05) # 3064
sum(as.numeric(LogRatioDifferenceTTestsForLipidSpeciesbwlschu$p.value) <= 0.05)/dim(LogRatioDifferenceTTestsForLipidSpeciesbwlschu)[1] # 0.03582494

InputMatrixOverlapVsOverexpression <- cbind(OE = c(No = dim(LogRatioDifferenceTTestsForLipidSpeciesbwlschu)[1] - sum(as.numeric(LogRatioDifferenceTTestsForLipidSpeciesbwlschu$p.value) <= 0.05/dim(LogRatioDifferenceTTestsForLipidSpeciesbwlschu)[1]),
                                                   Yes = sum(as.numeric(LogRatioDifferenceTTestsForLipidSpeciesbwlschu$p.value) <= 0.05/dim(LogRatioDifferenceTTestsForLipidSpeciesbwlschu)[1])),
                                            
                                            Overlap = c(No = dim(MappingSpeciesScreensOEm)[1] - sum(as.numeric(MappingSpeciesScreensOEm$p.value) <= 0.05/dim(MappingSpeciesScreensOEm)[1]), 
                                                        Yes = sum(as.numeric(MappingSpeciesScreensOEm$p.value) <= 0.05/dim(MappingSpeciesScreensOEm)[1])))

ListOfInputMatricesOverlapVsOverexpression <- list(CorrectionOnSampleSize = cbind(OE = c(No = dim(LogRatioDifferenceTTestsForLipidSpeciesbwlschu)[1] - sum(as.numeric(LogRatioDifferenceTTestsForLipidSpeciesbwlschu$p.value) <= 0.05/dim(LogRatioDifferenceTTestsForLipidSpeciesbwlschu)[1]),
                                                                                         Yes = sum(as.numeric(LogRatioDifferenceTTestsForLipidSpeciesbwlschu$p.value) <= 0.05/dim(LogRatioDifferenceTTestsForLipidSpeciesbwlschu)[1])),
                                                                                  
                                                                                  Overlap = c(No = dim(MappingSpeciesScreensOEm)[1] - sum(as.numeric(MappingSpeciesScreensOEm$p.value) <= 0.05/dim(MappingSpeciesScreensOEm)[1]), 
                                                                                              Yes = sum(as.numeric(MappingSpeciesScreensOEm$p.value) <= 0.05/dim(MappingSpeciesScreensOEm)[1]))),
                                                   
                                                   CorrectionOnLTPGroupSize = cbind(OE = c(No = dim(LogRatioDifferenceTTestsForLipidSpeciesbwlschu)[1] - sum(as.numeric(LogRatioDifferenceTTestsForLipidSpeciesbwlschu$p.value) <= 0.05/1989),
                                                                                           Yes = sum(as.numeric(LogRatioDifferenceTTestsForLipidSpeciesbwlschu$p.value) <= 0.05/1989)),
                                                                                    
                                                                                    Overlap = c(No = dim(MappingSpeciesScreensOEm)[1] - sum(as.numeric(MappingSpeciesScreensOEm$p.value) <= 0.05/1989), 
                                                                                                Yes = sum(as.numeric(MappingSpeciesScreensOEm$p.value) <= 0.05/1989))),
                                                   
                                                   CorrectionOnOESize = cbind(OE = c(No = dim(LogRatioDifferenceTTestsForLipidSpeciesbwlschu)[1] - sum(as.numeric(LogRatioDifferenceTTestsForLipidSpeciesbwlschu$p.value) <= 0.05/dim(LogRatioDifferenceTTestsForLipidSpeciesbwlschu)[1]),
                                                                                     Yes = sum(as.numeric(LogRatioDifferenceTTestsForLipidSpeciesbwlschu$p.value) <= 0.05/dim(LogRatioDifferenceTTestsForLipidSpeciesbwlschu)[1])),
                                                                              
                                                                              Overlap = c(No = dim(MappingSpeciesScreensOEm)[1] - sum(as.numeric(MappingSpeciesScreensOEm$p.value) <= 0.05/dim(LogRatioDifferenceTTestsForLipidSpeciesbwlschu)[1]), 
                                                                                          Yes = sum(as.numeric(MappingSpeciesScreensOEm$p.value) <= 0.05/dim(LogRatioDifferenceTTestsForLipidSpeciesbwlschu)[1]))),
                                                   
                                                   CorrectionOnOverlapSize = cbind(OE = c(No = dim(LogRatioDifferenceTTestsForLipidSpeciesbwlschu)[1] - sum(as.numeric(LogRatioDifferenceTTestsForLipidSpeciesbwlschu$p.value) <= 0.05/dim(MappingSpeciesScreensOEm)[1]),
                                                                                          Yes = sum(as.numeric(LogRatioDifferenceTTestsForLipidSpeciesbwlschu$p.value) <= 0.05/dim(MappingSpeciesScreensOEm)[1])),
                                                                                   
                                                                                   Overlap = c(No = dim(MappingSpeciesScreensOEm)[1] - sum(as.numeric(MappingSpeciesScreensOEm$p.value) <= 0.05/dim(MappingSpeciesScreensOEm)[1]), 
                                                                                               Yes = sum(as.numeric(MappingSpeciesScreensOEm$p.value) <= 0.05/dim(MappingSpeciesScreensOEm)[1]))),
                                                   
                                                   NoCorrection = cbind(OE = c(No = dim(LogRatioDifferenceTTestsForLipidSpeciesbwlschu)[1] - sum(as.numeric(LogRatioDifferenceTTestsForLipidSpeciesbwlschu$p.value) <= 0.05),
                                                                               Yes = sum(as.numeric(LogRatioDifferenceTTestsForLipidSpeciesbwlschu$p.value) <= 0.05)),
                                                                        
                                                                        Overlap = c(No = dim(MappingSpeciesScreensOEm)[1] - sum(as.numeric(MappingSpeciesScreensOEm$p.value) <= 0.05), 
                                                                                    Yes = sum(as.numeric(MappingSpeciesScreensOEm$p.value) <= 0.05))))


ResultsFisherExactTestsOverlapVsOverexpression <- lapply(ListOfInputMatricesOverlapVsOverexpression, fisher.test)

# Look at LTP - lipid subclass
sum(!(is.na(CombinedViewOfScreens2$InCellulo) & is.na(CombinedViewOfScreens2$InVitro)) &
      
      CombinedViewOfScreens2$LipidInOEData &
      CombinedViewOfScreens2$LTPProteinInOEData &
      
      CombinedViewOfScreens2$MatchesBroadSignificantOEs)
# Broad significance matches of overlap: 41 (in cellulo: 25, and in vitro: 21, both: 5)


sum(!(is.na(CombinedViewOfScreens2$InCellulo) & is.na(CombinedViewOfScreens2$InVitro)) &
      
      CombinedViewOfScreens2$LipidInOEData &
      CombinedViewOfScreens2$LTPProteinInOEData &
      
      !CombinedViewOfScreens2$MatchesBroadSignificantOEs)
# Broad non-significance matches of overlap: 53 (in cellulo: 18, and in vitro: 38, both: 3)


sum(!(is.na(CombinedViewOfScreens2$InCellulo) & is.na(CombinedViewOfScreens2$InVitro)) &
      
      !(CombinedViewOfScreens2$LipidInOEData & CombinedViewOfScreens2$LTPProteinInOEData))
# Not possible to match because of lipid/LTP absence in overexpression data: 20

dim(BroadSubsetPairsUponOverexpression2)[1] # How many significant LTP - lipid subclass matches present? # 225
dim(unique(LogRatioDifferenceTTestsForLipidSpeciesbwls[, c("NewLTPs", "ConvertedOELipidName")]))[1] # How many LTP - lipid subclass entries in total? # 1032

# Percentage of significant matches for the screen overlap
41*100/(53+41) # 43.61702

# Percentage of significant matches for the overexpression data
225*100/(1032) # 21.80233


# Test which LTPs in overexpression data but not in overlap data

unique(LogRatioDifferenceTTestsForLipidSpeciesbwls[,"NewLTPs"])[!(unique(LogRatioDifferenceTTestsForLipidSpeciesbwls[,"NewLTPs"]) %in% unique(as.character(CombinedViewOfScreens2$LTPProtein)))] 
# OSBPL10 OSBPL11 OSBPL7  OSBPL8 

unique(BroadSubsetPairsUponOverexpression2[,"NewLTPs"])[!(unique(BroadSubsetPairsUponOverexpression2[,"NewLTPs"]) %in% unique(as.character(CombinedViewOfScreens2$LTPProtein)))] 
# OSBPL10 OSBPL11 OSBPL7  OSBPL8 

dim(BroadSubsetPairsUponOverexpression2[BroadSubsetPairsUponOverexpression2[,"NewLTPs"] %in% unique(as.character(CombinedViewOfScreens2$LTPProtein)),])[1] # How many significant LTP - lipid subclass matches present in only LTPs present everywhere? # 207
dim(unique(LogRatioDifferenceTTestsForLipidSpeciesbwls[LogRatioDifferenceTTestsForLipidSpeciesbwls[,"NewLTPs"] %in% unique(as.character(CombinedViewOfScreens2$LTPProtein)), c("NewLTPs", "ConvertedOELipidName")]))[1] # How many LTP - lipid subclass entries in total for overexpression dataset with focus on only LTPs present in screens? # 936

ListOfInputMatricesOverlapVsOverexpressionForSubClasses <- list(CorrectionOnSampleSizeFocusOnKnown = cbind(OE = c(No = dim(unique(LogRatioDifferenceTTestsForLipidSpeciesbwls[, c("NewLTPs", "ConvertedOELipidName")]))[1] - dim(BroadSubsetPairsUponOverexpression2)[1],
                                                                                                                  Yes = dim(BroadSubsetPairsUponOverexpression2)[1]),
                                                                                                           
                                                                                                           Overlap = c(No = sum(ScreenMatchingOEOverview[,"KnownUnmatchedInOE"]), 
                                                                                                                       Yes = sum(ScreenMatchingOEOverview[,"KnownMatchedInOE"]))),
                                                                
                                                                CorrectionOnSampleSizeFocusOnScreenLTPSFocusOnKnown = cbind(OE = c(No = dim(unique(LogRatioDifferenceTTestsForLipidSpeciesbwls[LogRatioDifferenceTTestsForLipidSpeciesbwls[,"NewLTPs"] %in% unique(as.character(CombinedViewOfScreens2$LTPProtein)), c("NewLTPs", "ConvertedOELipidName")]))[1] - dim(BroadSubsetPairsUponOverexpression2[BroadSubsetPairsUponOverexpression2[,"NewLTPs"] %in% unique(as.character(CombinedViewOfScreens2$LTPProtein)),])[1],
                                                                                                                                   Yes = dim(BroadSubsetPairsUponOverexpression2[BroadSubsetPairsUponOverexpression2[,"NewLTPs"] %in% unique(as.character(CombinedViewOfScreens2$LTPProtein)),])[1]),
                                                                                                                            
                                                                                                                            Overlap = c(No = sum(ScreenMatchingOEOverview[,"KnownUnmatchedInOE"]), 
                                                                                                                                        Yes = sum(ScreenMatchingOEOverview[,"KnownMatchedInOE"]))),
                                                                
                                                                CorrectionOnSampleSizeFocusOnNovel = cbind(OE = c(No = dim(unique(LogRatioDifferenceTTestsForLipidSpeciesbwls[, c("NewLTPs", "ConvertedOELipidName")]))[1] - dim(BroadSubsetPairsUponOverexpression2)[1],
                                                                                                                  Yes = dim(BroadSubsetPairsUponOverexpression2)[1]),
                                                                                                           
                                                                                                           Overlap = c(No = sum(ScreenMatchingOEOverview[,"NovelUnmatchedInOE"]), 
                                                                                                                       Yes = sum(ScreenMatchingOEOverview[,"NovelMatchedInOE"]))),
                                                                
                                                                CorrectionOnSampleSizeFocusOnScreenLTPSFocusOnNovel = cbind(OE = c(No = dim(unique(LogRatioDifferenceTTestsForLipidSpeciesbwls[LogRatioDifferenceTTestsForLipidSpeciesbwls[,"NewLTPs"] %in% unique(as.character(CombinedViewOfScreens2$LTPProtein)), c("NewLTPs", "ConvertedOELipidName")]))[1] - dim(BroadSubsetPairsUponOverexpression2[BroadSubsetPairsUponOverexpression2[,"NewLTPs"] %in% unique(as.character(CombinedViewOfScreens2$LTPProtein)),])[1],
                                                                                                                                   Yes = dim(BroadSubsetPairsUponOverexpression2[BroadSubsetPairsUponOverexpression2[,"NewLTPs"] %in% unique(as.character(CombinedViewOfScreens2$LTPProtein)),])[1]),
                                                                                                                            
                                                                                                                            Overlap = c(No = sum(ScreenMatchingOEOverview[,"NovelUnmatchedInOE"]), 
                                                                                                                                        Yes = sum(ScreenMatchingOEOverview[,"NovelMatchedInOE"]))),
                                                                
                                                                CorrectionOnSampleSize = cbind(OE = c(No = dim(unique(LogRatioDifferenceTTestsForLipidSpeciesbwls[, c("NewLTPs", "ConvertedOELipidName")]))[1] - dim(BroadSubsetPairsUponOverexpression2)[1],
                                                                                                      Yes = dim(BroadSubsetPairsUponOverexpression2)[1]),
                                                                                               
                                                                                               
                                                                                               Overlap = c(No = sum(!(is.na(CombinedViewOfScreens2$InCellulo) & is.na(CombinedViewOfScreens2$InVitro)) &
                                                                                                                      
                                                                                                                      CombinedViewOfScreens2$LipidInOEData &
                                                                                                                      CombinedViewOfScreens2$LTPProteinInOEData &
                                                                                                                      
                                                                                                                      !CombinedViewOfScreens2$MatchesBroadSignificantOEs), 
                                                                                                           Yes = sum(!(is.na(CombinedViewOfScreens2$InCellulo) & is.na(CombinedViewOfScreens2$InVitro)) &
                                                                                                                       
                                                                                                                       CombinedViewOfScreens2$LipidInOEData &
                                                                                                                       CombinedViewOfScreens2$LTPProteinInOEData &
                                                                                                                       
                                                                                                                       CombinedViewOfScreens2$MatchesBroadSignificantOEs))),
                                                                
                                                                
                                                                CorrectionOnSampleSizeFocusOnScreenLTPS = cbind(OE = c(No = dim(unique(LogRatioDifferenceTTestsForLipidSpeciesbwls[LogRatioDifferenceTTestsForLipidSpeciesbwls[,"NewLTPs"] %in% unique(as.character(CombinedViewOfScreens2$LTPProtein)), c("NewLTPs", "ConvertedOELipidName")]))[1] - dim(BroadSubsetPairsUponOverexpression2[BroadSubsetPairsUponOverexpression2[,"NewLTPs"] %in% unique(as.character(CombinedViewOfScreens2$LTPProtein)),])[1],
                                                                                                                       Yes = dim(BroadSubsetPairsUponOverexpression2[BroadSubsetPairsUponOverexpression2[,"NewLTPs"] %in% unique(as.character(CombinedViewOfScreens2$LTPProtein)),])[1]),
                                                                                                                
                                                                                                                
                                                                                                                Overlap = c(No = sum(!(is.na(CombinedViewOfScreens2$InCellulo) & is.na(CombinedViewOfScreens2$InVitro)) &
                                                                                                                                       
                                                                                                                                       CombinedViewOfScreens2$LipidInOEData &
                                                                                                                                       CombinedViewOfScreens2$LTPProteinInOEData &
                                                                                                                                       
                                                                                                                                       !CombinedViewOfScreens2$MatchesBroadSignificantOEs), 
                                                                                                                            Yes = sum(!(is.na(CombinedViewOfScreens2$InCellulo) & is.na(CombinedViewOfScreens2$InVitro)) &
                                                                                                                                        
                                                                                                                                        CombinedViewOfScreens2$LipidInOEData &
                                                                                                                                        CombinedViewOfScreens2$LTPProteinInOEData &
                                                                                                                                        
                                                                                                                                        CombinedViewOfScreens2$MatchesBroadSignificantOEs))),
                                                                
                                                                
                                                                CorrectionOnSampleSizeOnlyInCellulo = cbind(OE = c(No = dim(unique(LogRatioDifferenceTTestsForLipidSpeciesbwls[, c("NewLTPs", "ConvertedOELipidName")]))[1] - dim(BroadSubsetPairsUponOverexpression2)[1],
                                                                                                                   Yes = dim(BroadSubsetPairsUponOverexpression2)[1]),
                                                                                                            
                                                                                                            
                                                                                                            Overlap = c(No = sum(!(is.na(CombinedViewOfScreens2$InCellulo)) &
                                                                                                                                   
                                                                                                                                   CombinedViewOfScreens2$LipidInOEData &
                                                                                                                                   CombinedViewOfScreens2$LTPProteinInOEData &
                                                                                                                                   
                                                                                                                                   !CombinedViewOfScreens2$MatchesBroadSignificantOEs), 
                                                                                                                        Yes = sum(!(is.na(CombinedViewOfScreens2$InCellulo)) &
                                                                                                                                    
                                                                                                                                    CombinedViewOfScreens2$LipidInOEData &
                                                                                                                                    CombinedViewOfScreens2$LTPProteinInOEData &
                                                                                                                                    
                                                                                                                                    CombinedViewOfScreens2$MatchesBroadSignificantOEs))),
                                                                
                                                                
                                                                CorrectionOnSampleSizeFocusOnScreenLTPSOnlyInCellulo = cbind(OE = c(No = dim(unique(LogRatioDifferenceTTestsForLipidSpeciesbwls[LogRatioDifferenceTTestsForLipidSpeciesbwls[,"NewLTPs"] %in% unique(as.character(CombinedViewOfScreens2$LTPProtein)), c("NewLTPs", "ConvertedOELipidName")]))[1] - dim(BroadSubsetPairsUponOverexpression2[BroadSubsetPairsUponOverexpression2[,"NewLTPs"] %in% unique(as.character(CombinedViewOfScreens2$LTPProtein)),])[1],
                                                                                                                                    Yes = dim(BroadSubsetPairsUponOverexpression2[BroadSubsetPairsUponOverexpression2[,"NewLTPs"] %in% unique(as.character(CombinedViewOfScreens2$LTPProtein)),])[1]),
                                                                                                                             
                                                                                                                             
                                                                                                                             Overlap = c(No = sum(!(is.na(CombinedViewOfScreens2$InCellulo)) &
                                                                                                                                                    
                                                                                                                                                    CombinedViewOfScreens2$LipidInOEData &
                                                                                                                                                    CombinedViewOfScreens2$LTPProteinInOEData &
                                                                                                                                                    
                                                                                                                                                    !CombinedViewOfScreens2$MatchesBroadSignificantOEs), 
                                                                                                                                         Yes = sum(!(is.na(CombinedViewOfScreens2$InCellulo)) &
                                                                                                                                                     
                                                                                                                                                     CombinedViewOfScreens2$LipidInOEData &
                                                                                                                                                     CombinedViewOfScreens2$LTPProteinInOEData &
                                                                                                                                                     
                                                                                                                                                     CombinedViewOfScreens2$MatchesBroadSignificantOEs))),
                                                                
                                                                
                                                                CorrectionOnSampleSizeOnlyInVitro = cbind(OE = c(No = dim(unique(LogRatioDifferenceTTestsForLipidSpeciesbwls[, c("NewLTPs", "ConvertedOELipidName")]))[1] - dim(BroadSubsetPairsUponOverexpression2)[1],
                                                                                                                 Yes = dim(BroadSubsetPairsUponOverexpression2)[1]),
                                                                                                          
                                                                                                          
                                                                                                          Overlap = c(No = sum(!(is.na(CombinedViewOfScreens2$InVitro)) &
                                                                                                                                 
                                                                                                                                 CombinedViewOfScreens2$LipidInOEData &
                                                                                                                                 CombinedViewOfScreens2$LTPProteinInOEData &
                                                                                                                                 
                                                                                                                                 !CombinedViewOfScreens2$MatchesBroadSignificantOEs), 
                                                                                                                      Yes = sum(!(is.na(CombinedViewOfScreens2$InVitro)) &
                                                                                                                                  
                                                                                                                                  CombinedViewOfScreens2$LipidInOEData &
                                                                                                                                  CombinedViewOfScreens2$LTPProteinInOEData &
                                                                                                                                  
                                                                                                                                  CombinedViewOfScreens2$MatchesBroadSignificantOEs))),
                                                                
                                                                
                                                                CorrectionOnSampleSizeFocusOnScreenLTPSOnlyInVitro = cbind(OE = c(No = dim(unique(LogRatioDifferenceTTestsForLipidSpeciesbwls[LogRatioDifferenceTTestsForLipidSpeciesbwls[,"NewLTPs"] %in% unique(as.character(CombinedViewOfScreens2$LTPProtein)), c("NewLTPs", "ConvertedOELipidName")]))[1] - dim(BroadSubsetPairsUponOverexpression2[BroadSubsetPairsUponOverexpression2[,"NewLTPs"] %in% unique(as.character(CombinedViewOfScreens2$LTPProtein)),])[1],
                                                                                                                                  Yes = dim(BroadSubsetPairsUponOverexpression2[BroadSubsetPairsUponOverexpression2[,"NewLTPs"] %in% unique(as.character(CombinedViewOfScreens2$LTPProtein)),])[1]),
                                                                                                                           
                                                                                                                           
                                                                                                                           Overlap = c(No = sum(!(is.na(CombinedViewOfScreens2$InVitro)) &
                                                                                                                                                  
                                                                                                                                                  CombinedViewOfScreens2$LipidInOEData &
                                                                                                                                                  CombinedViewOfScreens2$LTPProteinInOEData &
                                                                                                                                                  
                                                                                                                                                  !CombinedViewOfScreens2$MatchesBroadSignificantOEs), 
                                                                                                                                       Yes = sum(!(is.na(CombinedViewOfScreens2$InVitro)) &
                                                                                                                                                   
                                                                                                                                                   CombinedViewOfScreens2$LipidInOEData &
                                                                                                                                                   CombinedViewOfScreens2$LTPProteinInOEData &
                                                                                                                                                   
                                                                                                                                                   CombinedViewOfScreens2$MatchesBroadSignificantOEs))),
                                                                
                                                                
                                                                CorrectionOnSampleSizeBothInCelluloAndInVitro = cbind(OE = c(No = dim(unique(LogRatioDifferenceTTestsForLipidSpeciesbwls[, c("NewLTPs", "ConvertedOELipidName")]))[1] - dim(BroadSubsetPairsUponOverexpression2)[1],
                                                                                                                             Yes = dim(BroadSubsetPairsUponOverexpression2)[1]),
                                                                                                                      
                                                                                                                      
                                                                                                                      Overlap = c(No = sum(!(is.na(CombinedViewOfScreens2$InCellulo)) & !(is.na(CombinedViewOfScreens2$InVitro)) &
                                                                                                                                             
                                                                                                                                             CombinedViewOfScreens2$LipidInOEData &
                                                                                                                                             CombinedViewOfScreens2$LTPProteinInOEData &
                                                                                                                                             
                                                                                                                                             !CombinedViewOfScreens2$MatchesBroadSignificantOEs), 
                                                                                                                                  Yes = sum(!(is.na(CombinedViewOfScreens2$InCellulo)) & !(is.na(CombinedViewOfScreens2$InVitro)) &
                                                                                                                                              
                                                                                                                                              CombinedViewOfScreens2$LipidInOEData &
                                                                                                                                              CombinedViewOfScreens2$LTPProteinInOEData &
                                                                                                                                              
                                                                                                                                              CombinedViewOfScreens2$MatchesBroadSignificantOEs))),
                                                                
                                                                
                                                                CorrectionOnSampleSizeFocusOnScreenLTPSBothInCelluloAndInVitro = cbind(OE = c(No = dim(unique(LogRatioDifferenceTTestsForLipidSpeciesbwls[LogRatioDifferenceTTestsForLipidSpeciesbwls[,"NewLTPs"] %in% unique(as.character(CombinedViewOfScreens2$LTPProtein)), c("NewLTPs", "ConvertedOELipidName")]))[1] - dim(BroadSubsetPairsUponOverexpression2[BroadSubsetPairsUponOverexpression2[,"NewLTPs"] %in% unique(as.character(CombinedViewOfScreens2$LTPProtein)),])[1],
                                                                                                                                              Yes = dim(BroadSubsetPairsUponOverexpression2[BroadSubsetPairsUponOverexpression2[,"NewLTPs"] %in% unique(as.character(CombinedViewOfScreens2$LTPProtein)),])[1]),
                                                                                                                                       
                                                                                                                                       
                                                                                                                                       Overlap = c(No = sum(!(is.na(CombinedViewOfScreens2$InCellulo)) & !(is.na(CombinedViewOfScreens2$InVitro)) &
                                                                                                                                                              
                                                                                                                                                              CombinedViewOfScreens2$LipidInOEData &
                                                                                                                                                              CombinedViewOfScreens2$LTPProteinInOEData &
                                                                                                                                                              
                                                                                                                                                              !CombinedViewOfScreens2$MatchesBroadSignificantOEs), 
                                                                                                                                                   Yes = sum(!(is.na(CombinedViewOfScreens2$InCellulo)) & !(is.na(CombinedViewOfScreens2$InVitro)) &
                                                                                                                                                               
                                                                                                                                                               CombinedViewOfScreens2$LipidInOEData &
                                                                                                                                                               CombinedViewOfScreens2$LTPProteinInOEData &
                                                                                                                                                               
                                                                                                                                                               CombinedViewOfScreens2$MatchesBroadSignificantOEs))))
colSums(ScreenMatchingOEOverview[,c("NovelMatchedInOE", "NovelUnmatchedInOE", "KnownMatchedInOE", "KnownUnmatchedInOE")])

ListOfInputMatricesOverlapVsOverexpressionForSubClasses2 <- list(CorrectionOnSampleSize = cbind(OE = c(No = dim(unique(LogRatioDifferenceTTestsForLipidSpeciesbwls[, c("NewLTPs", "ConvertedOELipidName")]))[1] - dim(BroadSubsetPairsUponOverexpression2)[1],
                                                                                                       Yes = dim(BroadSubsetPairsUponOverexpression2)[1]),
                                                                                                
                                                                                                
                                                                                                Overlap = c(No = sum(!(is.na(CombinedViewOfScreens2$InCellulo) & is.na(CombinedViewOfScreens2$InVitro)) &
                                                                                                                       
                                                                                                                       CombinedViewOfScreens2$LipidInOEData &
                                                                                                                       CombinedViewOfScreens2$LTPProteinInOEData &
                                                                                                                       
                                                                                                                       !CombinedViewOfScreens2$MatchesBroadSignificantOEs), 
                                                                                                            Yes = sum(!(is.na(CombinedViewOfScreens2$InCellulo) & is.na(CombinedViewOfScreens2$InVitro)) &
                                                                                                                        
                                                                                                                        CombinedViewOfScreens2$LipidInOEData &
                                                                                                                        CombinedViewOfScreens2$LTPProteinInOEData &
                                                                                                                        
                                                                                                                        CombinedViewOfScreens2$MatchesBroadSignificantOEs))),
                                                                 
                                                                 
                                                                 CorrectionOnSampleSizeFocusOnScreenLTPS = cbind(OE = c(No = dim(unique(LogRatioDifferenceTTestsForLipidSpeciesbwls[LogRatioDifferenceTTestsForLipidSpeciesbwls[,"NewLTPs"] %in% unique(as.character(CombinedViewOfScreens2$LTPProtein)), c("NewLTPs", "ConvertedOELipidName")]))[1] - dim(BroadSubsetPairsUponOverexpression2[BroadSubsetPairsUponOverexpression2[,"NewLTPs"] %in% unique(as.character(CombinedViewOfScreens2$LTPProtein)),])[1],
                                                                                                                        Yes = dim(BroadSubsetPairsUponOverexpression2[BroadSubsetPairsUponOverexpression2[,"NewLTPs"] %in% unique(as.character(CombinedViewOfScreens2$LTPProtein)),])[1]),
                                                                                                                 
                                                                                                                 
                                                                                                                 Overlap = c(No = sum(!(is.na(CombinedViewOfScreens2$InCellulo) & is.na(CombinedViewOfScreens2$InVitro)) &
                                                                                                                                        
                                                                                                                                        CombinedViewOfScreens2$LipidInOEData &
                                                                                                                                        CombinedViewOfScreens2$LTPProteinInOEData &
                                                                                                                                        
                                                                                                                                        !CombinedViewOfScreens2$MatchesBroadSignificantOEs), 
                                                                                                                             Yes = sum(!(is.na(CombinedViewOfScreens2$InCellulo) & is.na(CombinedViewOfScreens2$InVitro)) &
                                                                                                                                         
                                                                                                                                         CombinedViewOfScreens2$LipidInOEData &
                                                                                                                                         CombinedViewOfScreens2$LTPProteinInOEData &
                                                                                                                                         
                                                                                                                                         CombinedViewOfScreens2$MatchesBroadSignificantOEs))),
                                                                 
                                                                 
                                                                 CorrectionOnSampleSizeOnlyInCellulo = cbind(OE = c(No = dim(unique(LogRatioDifferenceTTestsForLipidSpeciesbwls[, c("NewLTPs", "ConvertedOELipidName")]))[1] - dim(BroadSubsetPairsUponOverexpression2)[1],
                                                                                                                    Yes = dim(BroadSubsetPairsUponOverexpression2)[1]),
                                                                                                             
                                                                                                             
                                                                                                             Overlap = c(No = sum(!(is.na(CombinedViewOfScreens2$InCellulo)) &
                                                                                                                                    
                                                                                                                                    CombinedViewOfScreens2$LipidInOEData &
                                                                                                                                    CombinedViewOfScreens2$LTPProteinInOEData &
                                                                                                                                    
                                                                                                                                    !CombinedViewOfScreens2$MatchesBroadSignificantOEs), 
                                                                                                                         Yes = sum(!(is.na(CombinedViewOfScreens2$InCellulo)) &
                                                                                                                                     
                                                                                                                                     CombinedViewOfScreens2$LipidInOEData &
                                                                                                                                     CombinedViewOfScreens2$LTPProteinInOEData &
                                                                                                                                     
                                                                                                                                     CombinedViewOfScreens2$MatchesBroadSignificantOEs))),
                                                                 
                                                                 
                                                                 CorrectionOnSampleSizeFocusOnScreenLTPSOnlyInCellulo = cbind(OE = c(No = dim(unique(LogRatioDifferenceTTestsForLipidSpeciesbwls[LogRatioDifferenceTTestsForLipidSpeciesbwls[,"NewLTPs"] %in% unique(as.character(CombinedViewOfScreens2$LTPProtein)), c("NewLTPs", "ConvertedOELipidName")]))[1] - dim(BroadSubsetPairsUponOverexpression2[BroadSubsetPairsUponOverexpression2[,"NewLTPs"] %in% unique(as.character(CombinedViewOfScreens2$LTPProtein)),])[1],
                                                                                                                                     Yes = dim(BroadSubsetPairsUponOverexpression2[BroadSubsetPairsUponOverexpression2[,"NewLTPs"] %in% unique(as.character(CombinedViewOfScreens2$LTPProtein)),])[1]),
                                                                                                                              
                                                                                                                              
                                                                                                                              Overlap = c(No = sum(!(is.na(CombinedViewOfScreens2$InCellulo)) &
                                                                                                                                                     
                                                                                                                                                     CombinedViewOfScreens2$LipidInOEData &
                                                                                                                                                     CombinedViewOfScreens2$LTPProteinInOEData &
                                                                                                                                                     
                                                                                                                                                     !CombinedViewOfScreens2$MatchesBroadSignificantOEs), 
                                                                                                                                          Yes = sum(!(is.na(CombinedViewOfScreens2$InCellulo)) &
                                                                                                                                                      
                                                                                                                                                      CombinedViewOfScreens2$LipidInOEData &
                                                                                                                                                      CombinedViewOfScreens2$LTPProteinInOEData &
                                                                                                                                                      
                                                                                                                                                      CombinedViewOfScreens2$MatchesBroadSignificantOEs))),
                                                                 
                                                                 
                                                                 CorrectionOnSampleSizeOnlyInVitro = cbind(OE = c(No = dim(unique(LogRatioDifferenceTTestsForLipidSpeciesbwls[, c("NewLTPs", "ConvertedOELipidName")]))[1] - dim(BroadSubsetPairsUponOverexpression2)[1],
                                                                                                                  Yes = dim(BroadSubsetPairsUponOverexpression2)[1]),
                                                                                                           
                                                                                                           
                                                                                                           Overlap = c(No = sum(!(is.na(CombinedViewOfScreens2$InVitro)) &
                                                                                                                                  
                                                                                                                                  CombinedViewOfScreens2$LipidInOEData &
                                                                                                                                  CombinedViewOfScreens2$LTPProteinInOEData &
                                                                                                                                  
                                                                                                                                  !CombinedViewOfScreens2$MatchesBroadSignificantOEs), 
                                                                                                                       Yes = sum(!(is.na(CombinedViewOfScreens2$InVitro)) &
                                                                                                                                   
                                                                                                                                   CombinedViewOfScreens2$LipidInOEData &
                                                                                                                                   CombinedViewOfScreens2$LTPProteinInOEData &
                                                                                                                                   
                                                                                                                                   CombinedViewOfScreens2$MatchesBroadSignificantOEs))),
                                                                 
                                                                 
                                                                 CorrectionOnSampleSizeFocusOnScreenLTPSOnlyInVitro = cbind(OE = c(No = dim(unique(LogRatioDifferenceTTestsForLipidSpeciesbwls[LogRatioDifferenceTTestsForLipidSpeciesbwls[,"NewLTPs"] %in% unique(as.character(CombinedViewOfScreens2$LTPProtein)), c("NewLTPs", "ConvertedOELipidName")]))[1] - dim(BroadSubsetPairsUponOverexpression2[BroadSubsetPairsUponOverexpression2[,"NewLTPs"] %in% unique(as.character(CombinedViewOfScreens2$LTPProtein)),])[1],
                                                                                                                                   Yes = dim(BroadSubsetPairsUponOverexpression2[BroadSubsetPairsUponOverexpression2[,"NewLTPs"] %in% unique(as.character(CombinedViewOfScreens2$LTPProtein)),])[1]),
                                                                                                                            
                                                                                                                            
                                                                                                                            Overlap = c(No = sum(!(is.na(CombinedViewOfScreens2$InVitro)) &
                                                                                                                                                   
                                                                                                                                                   CombinedViewOfScreens2$LipidInOEData &
                                                                                                                                                   CombinedViewOfScreens2$LTPProteinInOEData &
                                                                                                                                                   
                                                                                                                                                   !CombinedViewOfScreens2$MatchesBroadSignificantOEs), 
                                                                                                                                        Yes = sum(!(is.na(CombinedViewOfScreens2$InVitro)) &
                                                                                                                                                    
                                                                                                                                                    CombinedViewOfScreens2$LipidInOEData &
                                                                                                                                                    CombinedViewOfScreens2$LTPProteinInOEData &
                                                                                                                                                    
                                                                                                                                                    CombinedViewOfScreens2$MatchesBroadSignificantOEs))),
                                                                 
                                                                 
                                                                 CorrectionOnSampleSizeBothInCelluloAndInVitro = cbind(OE = c(No = dim(unique(LogRatioDifferenceTTestsForLipidSpeciesbwls[, c("NewLTPs", "ConvertedOELipidName")]))[1] - dim(BroadSubsetPairsUponOverexpression2)[1],
                                                                                                                              Yes = dim(BroadSubsetPairsUponOverexpression2)[1]),
                                                                                                                       
                                                                                                                       
                                                                                                                       Overlap = c(No = sum(!(is.na(CombinedViewOfScreens2$InCellulo)) & !(is.na(CombinedViewOfScreens2$InVitro)) &
                                                                                                                                              
                                                                                                                                              CombinedViewOfScreens2$LipidInOEData &
                                                                                                                                              CombinedViewOfScreens2$LTPProteinInOEData &
                                                                                                                                              
                                                                                                                                              !CombinedViewOfScreens2$MatchesBroadSignificantOEs), 
                                                                                                                                   Yes = sum(!(is.na(CombinedViewOfScreens2$InCellulo)) & !(is.na(CombinedViewOfScreens2$InVitro)) &
                                                                                                                                               
                                                                                                                                               CombinedViewOfScreens2$LipidInOEData &
                                                                                                                                               CombinedViewOfScreens2$LTPProteinInOEData &
                                                                                                                                               
                                                                                                                                               CombinedViewOfScreens2$MatchesBroadSignificantOEs))),
                                                                 
                                                                 
                                                                 CorrectionOnSampleSizeFocusOnScreenLTPSBothInCelluloAndInVitro = cbind(OE = c(No = dim(unique(LogRatioDifferenceTTestsForLipidSpeciesbwls[LogRatioDifferenceTTestsForLipidSpeciesbwls[,"NewLTPs"] %in% unique(as.character(CombinedViewOfScreens2$LTPProtein)), c("NewLTPs", "ConvertedOELipidName")]))[1] - dim(BroadSubsetPairsUponOverexpression2[BroadSubsetPairsUponOverexpression2[,"NewLTPs"] %in% unique(as.character(CombinedViewOfScreens2$LTPProtein)),])[1],
                                                                                                                                               Yes = dim(BroadSubsetPairsUponOverexpression2[BroadSubsetPairsUponOverexpression2[,"NewLTPs"] %in% unique(as.character(CombinedViewOfScreens2$LTPProtein)),])[1]),
                                                                                                                                        
                                                                                                                                        
                                                                                                                                        Overlap = c(No = sum(!(is.na(CombinedViewOfScreens2$InCellulo)) & !(is.na(CombinedViewOfScreens2$InVitro)) &
                                                                                                                                                               
                                                                                                                                                               CombinedViewOfScreens2$LipidInOEData &
                                                                                                                                                               CombinedViewOfScreens2$LTPProteinInOEData &
                                                                                                                                                               
                                                                                                                                                               !CombinedViewOfScreens2$MatchesBroadSignificantOEs), 
                                                                                                                                                    Yes = sum(!(is.na(CombinedViewOfScreens2$InCellulo)) & !(is.na(CombinedViewOfScreens2$InVitro)) &
                                                                                                                                                                
                                                                                                                                                                CombinedViewOfScreens2$LipidInOEData &
                                                                                                                                                                CombinedViewOfScreens2$LTPProteinInOEData &
                                                                                                                                                                
                                                                                                                                                                CombinedViewOfScreens2$MatchesBroadSignificantOEs))))
# A part of the output below and above forms the basis for visualizations and analyses in Fig.1d                                             

ResultsFisherExactTestsOverlapVsOverexpressionSubclasses <- lapply(ListOfInputMatricesOverlapVsOverexpressionForSubClasses, fisher.test)


ResultsOverviewOfStatisticsMappingOfSubclasses <- cbind(do.call("rbind", lapply(ListOfInputMatricesOverlapVsOverexpressionForSubClasses, function(x){x[4:1]})),
                                                        do.call("rbind", lapply(ResultsFisherExactTestsOverlapVsOverexpressionSubclasses, function(y){c(y[["conf.int"]][1:2], y[["estimate"]], y[["p.value"]])}))) 

colnames(ResultsOverviewOfStatisticsMappingOfSubclasses) <- c("SignificantInTheOverlap", "NonsignificantInTheOverlap", "SignificantInOEData", "NonSignificantInOEData", "LowerBound95%ConfidenceInterval", "UpperBound95%ConfidenceInterval", "OddsRatio", "p-value")
write.table(ResultsOverviewOfStatisticsMappingOfSubclasses, file = "./Output/ResultsOverviewOfStatisticsMappingOfSubclasses10042025.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE) # Before: "D:/Rest/Rest5/ResultsOverviewOfStatisticsMappingOfSubclasses12032025.tsv"


LaterSubsetOfOverlap <- CombinedViewOfScreens2[!(is.na(CombinedViewOfScreens2$InVitro)) &
                                                 
                                                 CombinedViewOfScreens2$LipidInOEData &
                                                 CombinedViewOfScreens2$LTPProteinInOEData,]

sum(LaterSubsetOfOverlap$MatchesBroadSignificantOEs) # 21
sum(LaterSubsetOfOverlap$MatchesNarrowSignificantOEs) # 20


####

LipidClassLiteratureDataSetslc4hdrwan <- LipidClassLiteratureDataSetslc4hdr
rownames(LipidClassLiteratureDataSetslc4hdrwan) <- rownames(InVivoDataSetslc4hdr)

library(reshape2) #! Package dependency that is removable
OverviewOfLTPLipidSubclassesAndNovelties <- do.call("cbind", lapply(list(InVivoDataSetslc4hdr, InVitroDataSetslc4hdr, LipidClassLiteratureDataSetslc4hdrwan), melt))

OverviewOfLTPLipidSubclassesAndNovelties2 <- OverviewOfLTPLipidSubclassesAndNovelties[!(is.na(OverviewOfLTPLipidSubclassesAndNovelties[,3]) & is.na(OverviewOfLTPLipidSubclassesAndNovelties[,6])),c(1:3,6,9)]
colnames(OverviewOfLTPLipidSubclassesAndNovelties2) <- c("LipidSubclass", "LTPProtein", "InCellulo", "InVitro", "NovelPair")

OverviewOfLTPLipidSubclassesAndNovelties2OnlyLogicalsPresented <- cbind(OverviewOfLTPLipidSubclassesAndNovelties2[,c(1,2)], !is.na(OverviewOfLTPLipidSubclassesAndNovelties2[,3]), !is.na(OverviewOfLTPLipidSubclassesAndNovelties2[,4]), OverviewOfLTPLipidSubclassesAndNovelties2[,5])
colnames(OverviewOfLTPLipidSubclassesAndNovelties2OnlyLogicalsPresented) <- c("LipidSubclass", "LTPProtein", "InCellulo", "InVitro", "NovelPair")

write.table(OverviewOfLTPLipidSubclassesAndNovelties2OnlyLogicalsPresented, file = "./Output/469777_SupplementaryTable5b.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(OverviewOfLTPLipidSubclassesAndNovelties2, file = "./Output/AlternativeToSupplementaryTable5bWithNormalizedIntensities08052025.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# Aggregate the insights of the previous files to only information on presence or absence of a lipid detected by each of the LTPs
# This file was later updated in its headers, and STARD11 was converted to CERT naming, and some ORPs were removed as with the other data

write.table(cbind.data.frame(LTP = as.character(unique(OverviewOfLTPLipidSubclassesAndNovelties2OnlyLogicalsPresented$LTPProtein)), 
                             do.call("rbind", lapply(as.character(unique(OverviewOfLTPLipidSubclassesAndNovelties2OnlyLogicalsPresented$LTPProtein)), 
                                                     
                                                     function(x){colSums(OverviewOfLTPLipidSubclassesAndNovelties2OnlyLogicalsPresented[OverviewOfLTPLipidSubclassesAndNovelties2OnlyLogicalsPresented$LTPProtein == x, c("InCellulo", "InVitro")])>0}))),
            file = "./Output/469777_SupplementaryTable5b_v4.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# The above is the basis for part of Extended Data Table 6


#### (Re)make of figure panels for PI of PITPNA,B,C1 and SEC14L2, and PC of STARD2,10; also maybe cellular and liposomes one in new lines and LTPs not vertically combined
# PI-part: previous scripts used

LegendName <- "Legend"
LegendColor <- col_fung

WidthAdaptor <- 160
HeightAdaptor <- 160

# PI-entries and specific LTPs subset of all parts of the list
SubsetPITransLayersCondensedRowEntries <- lapply(NewPITransLayersCondensedRowEntries, function(x){x[(sapply(strsplit(rownames(x), "_"), "[[", 2) == "PI") & (sapply(strsplit(rownames(x), "_"), "[[", 1) %in% c("PITPNA", "PITPNB", "PITPNC1", "SEC14L2")),]})

CutOffPresenceProcent2 <- 0
SelectedColumnsForPIs <- colnames(SubsetPITransLayersCondensedRowEntries$InCellulo)[(colSums(SubsetPITransLayersCondensedRowEntries$InCellulo, na.rm = TRUE) > CutOffPresenceProcent2)|(colSums(SubsetPITransLayersCondensedRowEntries$InVitro, na.rm = TRUE) > CutOffPresenceProcent2)]

SubsetPITransLayersCondensedRowEntries2 <- lapply(SubsetPITransLayersCondensedRowEntries, function(x){x[,SelectedColumnsForPIs]})



# Bring in the part for the liposomes, and do other selection of species as the one above

SubselectedListOfPIPCPAInMobilizedRangexx <- setNames(lapply(match(c("PI", "PC", "PA"), names(PercentagesVesicleLipidomesSubclassesMatricesOfCombinedTissuesx)), function(y){unlist(lapply(list(NAsToZerosConverter(setNames(PercentagesVesicleLipidomesSubclassesMatricesOfCombinedTissuesx[[y]][SelectedColumnsForPIs], SelectedColumnsForPIs))),
                                                                                                                                                                                           function(x){x*100/max(x)}))}), c("PI", "PC", "PA"))

LiposomesSubsetForPI <- do.call("rbind", lapply(rownames(SubsetPITransLayersCondensedRowEntries2$InCellulo), function(x){SubselectedListOfPIPCPAInMobilizedRangexx$PI}))
rownames(LiposomesSubsetForPI) <- rownames(SubsetPITransLayersCondensedRowEntries2$InCellulo)

SubsetPITransLayersCondensedRowEntries2$Liposomes <- LiposomesSubsetForPI


TextDataframeWithLTPAnnotation2 <- data.frame(LTP = factor(sapply(strsplit(rownames(SubsetPITransLayersCondensedRowEntries2$InCellulo), "_"), "[[", 1), levels = unique(sapply(strsplit(rownames(SubsetPITransLayersCondensedRowEntries2$InCellulo), "_"), "[[", 1))))
AnnotationDataframeWithLTPAnnotation2 <- rowAnnotation(df = TextDataframeWithLTPAnnotation2, col = list(LTP = c("SEC14L2" = "#7E549F", "PITPNA" = "#FFCB3E", "PITPNB" = "#E0B01C", "PITPNC1" = "#A47C00")))

library(ComplexHeatmap)
library(RColorBrewer)

pdf("./Output/RemakeOfPISpecificSubsetOfLTPsForPanelsInFigure5NowTo6_12052025.pdf",
    width = unit(20, "mm"), height = unit(20, "mm")) # Additional part for figure 6c

Heatmap(SubsetPITransLayersCondensedRowEntries2$InCellulo, name = LegendName, col = LegendColor, rect_gp = gpar(type = "none"), column_title = "", 
        width = unit(WidthAdaptor*dim(SubsetPITransLayersCondensedRowEntries2$InCellulo)[2]/min(dim(SubsetPITransLayersCondensedRowEntries2$InCellulo)), "mm"), height = unit(HeightAdaptor*dim(SubsetPITransLayersCondensedRowEntries2$InCellulo)[1]/min(dim(SubsetPITransLayersCondensedRowEntries2$InCellulo)), "mm"),
        
        
        cell_fun = function(j, i, x, y, width, height, fill) {
          
          grid.rect(x = x, y = y-0.5*height, width = width, height = SubsetPITransLayersCondensedRowEntries2$InCellulo[i,j]/100*height, gp = gpar(col = NA, fill = brewer.pal(9,"Blues")[5], lwd = 1), just = c("center","bottom"))
          grid.rect(x = x, y = y-0.5*height, width = 0.5*width, height = SubsetPITransLayersCondensedRowEntries2$InVitro[i,j]/100*height, gp = gpar(col = NA, fill = brewer.pal(9,"Oranges")[5], lwd = 1), just = c("center","bottom"))
          
          grid.rect(x = x, y = y-0.5*height, width = width, height = SubsetPITransLayersCondensedRowEntries2$Cellular[i,j]/100*height, gp = gpar(col = "LightBlue", fill = NA, lwd = 3.2), just = c("center","bottom"))
          grid.rect(x = x, y = y-0.5*height, width = 0.8*width, height = SubsetPITransLayersCondensedRowEntries2$Liposomes[i,j]/100*height, gp = gpar(col = brewer.pal(9,"Oranges")[6], fill = NA, lwd = 2), just = c("center","bottom"))
          
        }, cluster_rows = FALSE, cluster_columns = FALSE, 
        cluster_row_slices = FALSE, cluster_column_slices = FALSE,
        
        row_title = " ", row_title_gp = gpar(fontsize = 10),
        left_annotation = AnnotationDataframeWithLTPAnnotation2,
        
        row_labels = sapply(strsplit(rownames(SubsetPITransLayersCondensedRowEntries2$InCellulo), "_"), "[[", 2),
        row_split = LETTERS[1:dim(SubsetPITransLayersCondensedRowEntries2$InCellulo)[1]]
        
)
dev.off() # Used as additional part for Fig.6c (combined in Adobe Illustrator)

# Get data from the final file for PC with STARD2 and STARD10
PCSubsetForSTARD2STARD10 <- lapply(c("in vivo", "in vitro"), function(x){LTPLipidConnectionsDataAggregated[(LTPLipidConnectionsDataAggregated$LTPProtein %in% c("STARD2","STARD10")) & (LTPLipidConnectionsDataAggregated$LikelySubclass == "PC") & (LTPLipidConnectionsDataAggregated$Screen == x),]})

PCSubsetForSTARD2STARD10b <- setNames(lapply(PCSubsetForSTARD2STARD10, function(x){as.matrix(as.data.frame(x[,-(1:3)])*100/do.call(pmax,as.data.frame(x[,-(1:3)])))}), nm = c("InCellulo", "InVitro"))
PCSubsetForSTARD2STARD10c <- lapply(PCSubsetForSTARD2STARD10b, function(x){rownames(x) <- c("STARD10", "STARD2"); return(x)})

CutOffPresenceProcent2 <- 5 
SelectedColumnsForPCs <- colnames(PCSubsetForSTARD2STARD10c$InCellulo)[(colSums(PCSubsetForSTARD2STARD10c$InCellulo, na.rm = TRUE) > CutOffPresenceProcent2)|(colSums(PCSubsetForSTARD2STARD10c$InVitro, na.rm = TRUE) > CutOffPresenceProcent2)]

PCSubsetForSTARD2STARD10d <- lapply(PCSubsetForSTARD2STARD10c, function(x){x[,SelectedColumnsForPCs]})



CellularBackgroundForPC <- LipidSubclassesAddedToBackground170620214[["PC"]][,"Cellular"]

CellularBackgroundForPC2 <- setNames(CellularBackgroundForPC, nm = sapply(strsplit(gsub("\\)", "", rownames(LipidSubclassesAddedToBackground170620214[["PC"]])), split = "\\("), "[[", 2))
CellularBackgroundForPC4 <- CellularBackgroundForPC2[SelectedColumnsForPCs]

CellularBackgroundForPC4[is.na(CellularBackgroundForPC4)] <- 0


CellularSubsetForPC <- do.call("rbind", lapply(rownames(PCSubsetForSTARD2STARD10d$InCellulo), function(x){CellularBackgroundForPC4}))
rownames(LiposomesSubsetForPC) <- rownames(PCSubsetForSTARD2STARD10d$InCellulo)

PCSubsetForSTARD2STARD10d$Cellular <- CellularSubsetForPC


LiposomesSubsetForPC <- do.call("rbind", lapply(rownames(PCSubsetForSTARD2STARD10d$InCellulo), function(x){PercentagesVesicleLipidomesSubclassesMatricesOfCombinedTissuesx[["PC"]][SelectedColumnsForPCs]}))
rownames(LiposomesSubsetForPC) <- rownames(PCSubsetForSTARD2STARD10d$InCellulo)

PCSubsetForSTARD2STARD10d$Liposomes <- LiposomesSubsetForPC


TextDataframeWithLTPAnnotation2ForPC <- data.frame(LTP = factor(c("STARD10", "STARD2"), levels = c("STARD10", "STARD2")))
AnnotationDataframeWithLTPAnnotation2ForPC <- rowAnnotation(df = TextDataframeWithLTPAnnotation2ForPC, col = list(LTP = c("STARD10" = "#7E549F", "STARD2" = "#FFCB3E")))

pdf("./Output/WorkdocForFigure6PanelPCToSTARD10And2ExtendedWithLiposomesAndCellularDataAt5PercentCutoff_13052025.pdf",
    width = unit(100, "mm"), height = unit(20, "mm")) # Additional part for figure 6c

Heatmap(PCSubsetForSTARD2STARD10d$InCellulo, name = LegendName, col = LegendColor, rect_gp = gpar(type = "none"), column_title = "", 
        width = unit(WidthAdaptor*dim(PCSubsetForSTARD2STARD10d$InCellulo)[2]/min(dim(PCSubsetForSTARD2STARD10d$InCellulo)), "mm"), height = unit(HeightAdaptor*dim(PCSubsetForSTARD2STARD10d$InCellulo)[1]/min(dim(PCSubsetForSTARD2STARD10d$InCellulo)), "mm"),
        
        
        cell_fun = function(j, i, x, y, width, height, fill) {
          
          grid.rect(x = x, y = y-0.5*height, width = width, height = PCSubsetForSTARD2STARD10d$InCellulo[i,j]/100*height, gp = gpar(col = NA, fill = brewer.pal(9,"Blues")[5], lwd = 1), just = c("center","bottom"))
          grid.rect(x = x, y = y-0.5*height, width = 0.5*width, height = PCSubsetForSTARD2STARD10d$InVitro[i,j]/100*height, gp = gpar(col = NA, fill = brewer.pal(9,"Oranges")[5], lwd = 1), just = c("center","bottom"))
          
          grid.rect(x = x, y = y-0.5*height, width = width, height = PCSubsetForSTARD2STARD10d$Cellular[i,j]/100*height, gp = gpar(col = "LightBlue", fill = NA, lwd = 3.2), just = c("center","bottom"))
          grid.rect(x = x, y = y-0.5*height, width = 0.8*width, height = PCSubsetForSTARD2STARD10d$Liposomes[i,j]/100*height, gp = gpar(col = brewer.pal(9,"Oranges")[6], fill = NA, lwd = 2), just = c("center","bottom"))
          
        }, cluster_rows = FALSE, cluster_columns = FALSE, 
        cluster_row_slices = FALSE, cluster_column_slices = FALSE,
        
        row_title = " ", row_title_gp = gpar(fontsize = 10),
        left_annotation = AnnotationDataframeWithLTPAnnotation2ForPC,
        
        row_labels = c("PC", "PC"),
        row_split = LETTERS[1:dim(PCSubsetForSTARD2STARD10d$InCellulo)[1]]
        
)
dev.off() # Also contributes to Fig.6c

# Export Extended Data Table 5 for data behind figure panel 2a (new article version): in cellulo, in vitro, and novelties
InCelluloInVitroNoveltySupplTable7Basis <- list(InVivoDataSetslc4hdr[,ReorderedLTPsByManualSeriation[c(1:19,21:24,28:43)]],
                                                
                                                InVitroDataSetslc4hdr[,ReorderedLTPsByManualSeriation[c(1:19,21:24,28:43)]],
                                                LipidClassLiteratureDataSetslc4hdr[,ReorderedLTPsByManualSeriation[c(1:19,21:24,28:43)]])

InCelluloInVitroNoveltySupplTable7Basis2 <- lapply(InCelluloInVitroNoveltySupplTable7Basis, function(x){
  colnames(x)[colnames(x) == "STARD11"] <- "CERT";
  
  rownames(x)[rownames(x) == "CH"] <- "Sterol";
  return(x)
  
})
InCelluloInVitroNoveltySupplTable7Basis2[[3]][is.na(InCelluloInVitroNoveltySupplTable7Basis2[[1]]) & is.na(InCelluloInVitroNoveltySupplTable7Basis2[[2]])] <- NA

# Write to Extended Data Table 5 files (in vitro was later updated to remove OSBPL9-PS connection because of discussion on the strength of the associated HPTLC datapoint)
write.table(InCelluloInVitroNoveltySupplTable7Basis2[[1]], file = "./Output/69777_SupplentaryTable7A_InCellulo.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

write.table(InCelluloInVitroNoveltySupplTable7Basis2[[2]], file = "./Output/69777_SupplentaryTable7B_InVitro.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
write.table(InCelluloInVitroNoveltySupplTable7Basis2[[3]], file = "./Output/69777_SupplentaryTable7C_PairNovelty.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

# THe above files form the basis for Extended Data Tables 5A, 5B and 5C.


write.table(OverexpressionHELA[c(1:2, 441:470),c(1:5,8:11,14:17)], file = "./Output/69777_SupplementaryTable11HeLaSubset.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
# The above table forms the basis for the Extended Data Table 11.


DetailsOfOverexpressionHEK2Split2 <- cbind(DetailsOfOverexpressionHEK2Split[,1:2], ControlMeans = rowMeans(DetailsOfOverexpressionHEK2Split[,c(4,10,16)])) 

# Make Supplementary Table for the cellular lipidome # The result is the basis for Extended Data Table 10B
SupplementaryTableForCellularLipidome <- OverexpressionHEK5[,c(1,3,9,15)]

colnames(SupplementaryTableForCellularLipidome) <- c("Lipid", "SampleAOfHEK293Lipidome", "SampleBOfHEK293Lipidome", "SampleCOfHEK293Lipidome")
write.table(SupplementaryTableForCellularLipidome, file = "./Output/469777_SupplementaryTableForCellularLipidome.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


# Make the Extended Data Table 7A & related steps

# Load in volume data from last version
StructuralVolumesLipidSpeciesLTPsx <- read.csv(file = "./InputData/469777_ExtendedDataTable7ALTPLipidVolumesColored.txt", header = TRUE, sep = "\t", as.is = TRUE)

StructuralVolumesLipidSpeciesLTPsxx <- cbind(StructuralVolumesLipidSpeciesLTPsx[,c("LTP", "lipid", "LTP.volume")], round(StructuralVolumesLipidSpeciesLTPsx[, "lipid.volume"], 2))
colnames(StructuralVolumesLipidSpeciesLTPsxx) <- c("Protein.Name", "Lipid", "Pocket.volume", "Lipid.volume")

StructuralVolumesLipidSpeciesLTPs2xx <- cbind(StructuralVolumesLipidSpeciesLTPsxx, LipidToPocketPercentage = StructuralVolumesLipidSpeciesLTPsxx$Lipid.volume*100/StructuralVolumesLipidSpeciesLTPsxx$Pocket.volume)
StructuralVolumesLipidSpeciesLTPs4xx <- unique(StructuralVolumesLipidSpeciesLTPs2xx)

any(duplicated(paste(StructuralVolumesLipidSpeciesLTPs4xx[,"Protein.Name"], StructuralVolumesLipidSpeciesLTPs4xx[,"Lipid"], sep = "_"))) # FALSE ==> OK


# Write out this file as cleaner Extended Data Table 7A
write.table(StructuralVolumesLipidSpeciesLTPs4xx, file = "./Output/469777_SupplementaryTable7ARoundedAndReducedAlternativeVersion06062025.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

CombinedDataOfScreens <- CleanedExtendedDataTable1Basis2 # Cleaner version with novelties
CombinedDataOfScreens[CombinedDataOfScreens == "STARD11"] <- "CERT"

CombinedDataOfScreens2 <- cbind(CombinedDataOfScreens, do.call("rbind", lapply(1:dim(CombinedDataOfScreens)[1], function(x){
  StructuralVolumesLipidSpeciesLTPs4[(StructuralVolumesLipidSpeciesLTPs4[,"Protein.Name"] %in% CombinedDataOfScreens[x,"LTPProtein"]) & (StructuralVolumesLipidSpeciesLTPs4[,"Lipid"] %in% CombinedDataOfScreens[x,"Lipid"]),3:5]
  
})))


CombinedDataOfScreens2xx <- cbind(CombinedDataOfScreens, do.call("rbind", lapply(1:dim(CombinedDataOfScreens)[1], function(x){
  StructuralVolumesLipidSpeciesLTPs4xx[(StructuralVolumesLipidSpeciesLTPs4xx[,"Protein.Name"] %in% CombinedDataOfScreens[x,"LTPProtein"]) & (StructuralVolumesLipidSpeciesLTPs4xx[,"Lipid"] %in% CombinedDataOfScreens[x,"Lipid"]),3:5]
  
}))) 
all(CombinedDataOfScreens2xx == CombinedDataOfScreens2) # Check: TRUE

CombinedDataOfScreens2 <- CombinedDataOfScreens2xx # Can stream through like this because equivalent
MaxKnownLipidToPocket <- max(CombinedDataOfScreens2[CombinedDataOfScreens2$NovelLinkLTPLipidSubclass == 0,"LipidToPocketPercentage"]) # 42.48058

CombinedDataOfScreens4 <- cbind(CombinedDataOfScreens2, AboveMaxOfKnown = CombinedDataOfScreens2[, "LipidToPocketPercentage"] > MaxKnownLipidToPocket)
CombinedDataOfScreens4[,c("LTPProtein", "Screen", "LikelySubclass","NovelLinkLTPLipidSubclass", "LipidToPocketPercentage")]

OverviewLTPLipidSubclassNoveltiesOverMaxSpace <- aggregate(CombinedDataOfScreens4$AboveMaxOfKnown, by = CombinedDataOfScreens4[,c("LTPProtein", "LikelySubclass", "Screen", "NovelLinkLTPLipidSubclass")], FUN = function(x){mean(x)*100})
colnames(OverviewLTPLipidSubclassNoveltiesOverMaxSpace) <- c("LTPProtein", "LipidSubclass", "Screen", "NovelLinkLTPLipidSubclass", "PercentageOfLipidSpeciesAboveKnownStructuralLimit")


OverviewLTPLipidSubclassNoveltiesOverMaxSpace2 <- cbind(OverviewLTPLipidSubclassNoveltiesOverMaxSpace, HigherLevelLipidSubclass = LipidNameConverterHigherSubclasses[match(OverviewLTPLipidSubclassNoveltiesOverMaxSpace$LipidSubclass, LipidNameConverterHigherSubclasses[,"StandardLevel"]), "HigherLevel"])

SignificantLTPLipidSubclassChangesInOE <- as.matrix(BroadSubsetPairsUponOverexpression2)
SignificantLTPLipidSubclassChangesInOE[SignificantLTPLipidSubclassChangesInOE == "STARD11"] <- "CERT" 

AllLTPLipidSubclassCombosInOE <- as.matrix(unique(LogRatioDifferenceTTestsForLipidSpeciesbwls[, c("NewLTPs", "ConvertedOELipidName")]))
AllLTPLipidSubclassCombosInOE[AllLTPLipidSubclassCombosInOE == "STARD11"] <- "CERT"


OverviewLTPLipidSubclassNoveltiesOverMaxSpace4 <- cbind(OverviewLTPLipidSubclassNoveltiesOverMaxSpace2, do.call("rbind", lapply(1:dim(OverviewLTPLipidSubclassNoveltiesOverMaxSpace2)[1], function(x){
  
  c(SignificantChangeOEData = any((SignificantLTPLipidSubclassChangesInOE[,"NewLTPs"] == as.character(OverviewLTPLipidSubclassNoveltiesOverMaxSpace2[x, "LTPProtein"])) &
                                    (SignificantLTPLipidSubclassChangesInOE[,"ConvertedOELipidName"] == as.character(OverviewLTPLipidSubclassNoveltiesOverMaxSpace2[x, "HigherLevelLipidSubclass"]))),
    
    CoveredInOEData = any((AllLTPLipidSubclassCombosInOE[,"NewLTPs"] == as.character(OverviewLTPLipidSubclassNoveltiesOverMaxSpace2[x, "LTPProtein"])) &
                            (AllLTPLipidSubclassCombosInOE[,"ConvertedOELipidName"] == as.character(OverviewLTPLipidSubclassNoveltiesOverMaxSpace2[x, "HigherLevelLipidSubclass"]))))
  
})))
write.table(OverviewLTPLipidSubclassNoveltiesOverMaxSpace4, file = "./Output/469777_SupplementaryTable8DataConnections.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# Previous table output is basis for the Extended Data Table 8.



# Analysis of CERT overexpression in HeLa cells

OverexpressionHELA <- read.csv(file = "./InputData/ExtractedOverexpressionHELAFromKenji15082019.txt", header = TRUE, sep = "\t", as.is = TRUE)
OverexpressionHELA2 <- OverexpressionHELA[-(1:2),]

colnames(OverexpressionHELA2) <- paste(colnames(OverexpressionHELA), OverexpressionHELA[1,], OverexpressionHELA[2,], sep = "_")
rownames(OverexpressionHELA2) <- 1:nrow(OverexpressionHELA2)

OverexpressionHELA4 <- apply(OverexpressionHELA2, 2, function(x){gsub(",",".",x)})
OverexpressionHELA5 <- cbind.data.frame(OverexpressionHELA4[,1], apply(OverexpressionHELA4[,-1], 2, function(x){as.numeric(x)}))

OverexpressionHELA5[,1] <- as.character(OverexpressionHELA5[,1])


OverexpressionHELARatiosMatrixAll <- do.call("cbind", lapply(1:9, function(x){(OverexpressionHELA5[,(2*x+1)] + 0.00001)/(OverexpressionHELA5[,(2*x)] + 0.00001)}))
colnames(OverexpressionHELARatiosMatrixAll) <- paste0(colnames(OverexpressionHELA5)[sapply(1:9,function(x){(2*x+1)})], "_ratio")

rownames(OverexpressionHELARatiosMatrixAll) <- OverexpressionHELA5[,1]


OverexpressionHELALogRatiosMatrix <- log10(OverexpressionHELARatiosMatrixAll)
OverexpressionHELALogRatiosGeneral2 <- t(OverexpressionHELALogRatiosMatrix[c("Cer", "CerP", "SM", "HexCer", "diHexCer", "triHexCer", "GM3", "DAG", "PA", "PA O-", "PC", "PC O-", "PE", "PE O-", "PS", "PI", "PG", "CL", "LPA", "LPC", "LPC O-", "LPE", "LPE O-", "LPS", "LPI", "LPG", "CE", "Chol :"),c(1,4,7,2,5,8,3,6,9)])

OverexpressionHELALogRatiosMelted <- melt(OverexpressionHELALogRatiosGeneral2)
colnames(OverexpressionHELALogRatiosMelted) <- c("Sample", "Lipid", "Ratio")


OverexpressionHELALogRatiosMelted2 <- cbind(OverexpressionHELALogRatiosMelted, LTPType = sapply(strsplit(as.character(OverexpressionHELALogRatiosMelted[,"Sample"]), split = "_"), "[[", 3))

StatHELALogRatiosMelted2 <- do.call("rbind",lapply(levels(OverexpressionHELALogRatiosMelted2[,"Lipid"]), function(y){c(y,sapply(c("CERT","Sec14L1"), function(x){unlist(t.test(OverexpressionHELALogRatiosMelted2[(OverexpressionHELALogRatiosMelted2[,"LTPType"] == "control") & (OverexpressionHELALogRatiosMelted2[,"Lipid"] == y),"Ratio"], OverexpressionHELALogRatiosMelted2[(OverexpressionHELALogRatiosMelted2[,"LTPType"] == x) & (OverexpressionHELALogRatiosMelted2[,"Lipid"] == y),"Ratio"])[c(3,5,4)])}))}))
colnames(StatHELALogRatiosMelted2) <- c("Lipid", "CERTp.value", "CERTControlMeanEstimate", "CERTProteinMeanEstimate", "CERTDiffConfidenceLow", "CERTDiffConfidenceHigh", "SEC14L1p.value", "SEC14L1ControlMeanEstimate", "SEC14L1ProteinMeanEstimate", "SEC14L1DiffConfidenceLow", "SEC14L1DiffConfidenceHigh")

StatHELALogRatiosMelted4 <- do.call("cbind", list(StatHELALogRatiosMelted2, CERTPercentDifference = 100*10^(as.numeric(StatHELALogRatiosMelted2[,4]) - as.numeric(StatHELALogRatiosMelted2[,3]))-100, SEC14L1PercentDifference = 100*10^(as.numeric(StatHELALogRatiosMelted2[,9]) - as.numeric(StatHELALogRatiosMelted2[,8]))-100))  


pdf("./Output/CERTAndSEC14L1LipidChangesHELAWelchsTtest150820192.pdf")

plot(as.numeric(StatHELALogRatiosMelted4[,"CERTPercentDifference"]), -log10(as.numeric(StatHELALogRatiosMelted4[,"CERTp.value"])), xlab = "Percent difference lipid abundance: protein vs. control", ylab = "-log10(p-value) [higher is better] (Welch's 2 sample t-test) ", pch = 16, col = "#FC4E07", main = "CERT (HELA)", xlim = c(-100,550))
abline(h = -log10(0.05), col = "lightgrey")

points(as.numeric(StatHELALogRatiosMelted4[as.numeric(StatHELALogRatiosMelted4[,"CERTp.value"]) <= 0.075,"CERTPercentDifference"]), -log10(as.numeric(StatHELALogRatiosMelted4[as.numeric(StatHELALogRatiosMelted4[,"CERTp.value"]) <= 0.075,"CERTp.value"])), pch = 16, col = "#E7B800")
points(as.numeric(StatHELALogRatiosMelted4[as.numeric(StatHELALogRatiosMelted4[,"CERTp.value"]) <= 0.05,"CERTPercentDifference"]), -log10(as.numeric(StatHELALogRatiosMelted4[as.numeric(StatHELALogRatiosMelted4[,"CERTp.value"]) <= 0.05,"CERTp.value"])), pch = 16, col = "#00AFBB")

text(as.numeric(StatHELALogRatiosMelted4[as.numeric(StatHELALogRatiosMelted4[,"CERTp.value"]) <= 0.075,"CERTPercentDifference"]), -log10(as.numeric(StatHELALogRatiosMelted4[as.numeric(StatHELALogRatiosMelted4[,"CERTp.value"]) <= 0.075,"CERTp.value"])), labels = paste0(StatHELALogRatiosMelted4[as.numeric(StatHELALogRatiosMelted4[,"CERTp.value"]) <= 0.075,"Lipid"], " (", round(as.numeric(StatHELALogRatiosMelted4[as.numeric(StatHELALogRatiosMelted4[,"CERTp.value"]) <= 0.075,"CERTPercentDifference"])), "%)"), pos = 4, cex = 0.8)
text(500, -log10(0.05), labels = "p-value 0.05", pos = 3, cex = 0.8, col = "darkgrey")

plot(as.numeric(StatHELALogRatiosMelted4[,"SEC14L1PercentDifference"]), -log10(as.numeric(StatHELALogRatiosMelted4[,"SEC14L1p.value"])), xlim = c(-100,550), xlab = "Percent difference lipid abundance: protein vs. control", ylab = "-log10(p-value) [higher is better] (Welch's 2 sample t-test) ", pch = 16, col = "#FC4E07", main = "SEC14L1 (HELA)")
abline(h = -log10(0.05), col = "lightgrey")

points(as.numeric(StatHELALogRatiosMelted4[as.numeric(StatHELALogRatiosMelted4[,"SEC14L1p.value"]) <= 0.075,"SEC14L1PercentDifference"]), -log10(as.numeric(StatHELALogRatiosMelted4[as.numeric(StatHELALogRatiosMelted4[,"SEC14L1p.value"]) <= 0.075,"SEC14L1p.value"])), pch = 16, col = "#E7B800")
points(as.numeric(StatHELALogRatiosMelted4[as.numeric(StatHELALogRatiosMelted4[,"SEC14L1p.value"]) <= 0.05,"SEC14L1PercentDifference"]), -log10(as.numeric(StatHELALogRatiosMelted4[as.numeric(StatHELALogRatiosMelted4[,"SEC14L1p.value"]) <= 0.05,"SEC14L1p.value"])), pch = 16, col = "#00AFBB")

text(as.numeric(StatHELALogRatiosMelted4[as.numeric(StatHELALogRatiosMelted4[,"SEC14L1p.value"]) <= 0.075,"SEC14L1PercentDifference"]), -log10(as.numeric(StatHELALogRatiosMelted4[as.numeric(StatHELALogRatiosMelted4[,"SEC14L1p.value"]) <= 0.075,"SEC14L1p.value"])), labels = paste0(StatHELALogRatiosMelted4[as.numeric(StatHELALogRatiosMelted4[,"SEC14L1p.value"]) <= 0.075,"Lipid"], " (", round(as.numeric(StatHELALogRatiosMelted4[as.numeric(StatHELALogRatiosMelted4[,"SEC14L1p.value"]) <= 0.075,"SEC14L1PercentDifference"])), "%)"), pos = 4, cex = 0.8)
text(500, -log10(0.05), labels = "p-value 0.05", pos = 3, cex = 0.8, col = "darkgrey")

dev.off()
# Basis for HeLa figure for CERT from Extended Data Fig.5a.

# Reset options to orginal ones from user
options(OriginalOptions)