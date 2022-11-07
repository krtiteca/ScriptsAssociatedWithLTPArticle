################ Document with the clean-up and calculations for figures of LTP article
## Legend:
# # Done and file was found 
# (#) Done in other context, or consequence of other files
#  #X Extra input
#  #. Relevant files for this part
  
################ Structure of this document

######## General functions 

######## MS-data input & clean-up

######## Supplementary material: quality tests
#### Supplementary material: ECN test
#### Supplementary material: adducts distribution

######## Figure 2
#### Fig.2a barchart
#### Fig.2a donut-diagram
#### Fig.2a circle-legend
#### (Fig.2b side-histogram)
#### Fig.2b circle heatmap visualization
#### Fig.2c unsaturation-distribution

######## Figure 3
#### Fig.3a
#### Supplementary material: Protein domain based ordering of LTPs
#### Fig.3b
#### Fig.3c

######## Figure 4
#### Fig.4a: Species version
#### (Fig. 4a: Subclasses version)
#### Fig.4b

######## Figure 5
#### Panel 5A
#### Panel 5B
#### Panel 5C
#### (Alternatives for panels of figure 5)
#### Supplementary material: results of Fisher exact tests

######## Figure 6
#### Fig.6b
#### Fig.6c


######## General functions
#### Function to substitute rownames by first column and then remove the first column

Col1ToRowNames <-function(x){
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


#### Function to convert the Zeros to NAs
ZerosToNAsConverter <- function(x){x[x == 0] <- NA; return(x)}

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

######## MS-data input 
#### MS results: in cellulo and in vitro: clean-up of input data 

# Note: although SEC14L1 data might be present at different stages, all entries associated with these were removed from the final results as a conservative precaution, 
# because of uncertainties about a potential experimental mix-up (in some cases) of CERT into the SEC14L1 entry (but not the inverse mix-up).

AntonellaData18062019 <- read.csv(file = "C:/Users/Kevin/Documents/RData/Rest2/ACGData/FiguresByKT/invivo_postproc_input_final_18062019.tsv", header = TRUE, sep = "\t", as.is = TRUE)
EnricData18062019 <- read.csv(file = "C:/Users/Kevin/Documents/RData/Rest2/ACGData/FiguresByKT/invitro_postproc_input_final_18062019.tsv", header = TRUE, sep = "\t", as.is = TRUE)

PureVersionAntonellaData18062019 <- AntonellaData18062019[AntonellaData18062019$cls == "I" | AntonellaData18062019$cls == "IV.7",]
PureVersionEnricData18062019 <- EnricData18062019[EnricData18062019$cls == "I" | EnricData18062019$cls == "IV.7",]

CombinedDataWithPureClasses18062019 <- rbind(cbind(PureVersionAntonellaData18062019, ScreenType = "in vivo"), cbind(PureVersionEnricData18062019 , ScreenType = "in vitro"))
write.table(CombinedDataWithPureClasses18062019, file="./RData/Rest2/ACGData/FiguresByKT/CombinedDataWithPureClasses18062019.csv", sep="\t", row.names = TRUE, quote = FALSE)

library(tidyr)
CombinedDataWithPureClasses180620192 <- unite(CombinedDataWithPureClasses18062019, "CompletelyUnique", c("protein","ionm","uid","mz","intensity"), sep = "_", remove = FALSE)

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

write.table(CombinedDataWithPureClasses1806201913, file="./RData/Rest2/ACGData/FiguresByKT/CleanConservativeDataWithoutFilters17072019.csv", sep="\t", row.names = FALSE, quote = FALSE)
#! Where supplementary table X is derived from


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
#### Supplementary material: ECN test


library(RColorBrewer)

BlueColorRangeUnsaturations <- setNames(rev(paste0(colorRampPalette(brewer.pal(9,"Blues")[-(1:2)])(11), "52")), as.character(0:10))
OrangeColorRangeUnsaturations <- setNames(rev(paste0(colorRampPalette(brewer.pal(9,"Oranges")[-1])(11), "52")), as.character(0:10))

BlueColorRangeUnsaturationsOpaque <- setNames(rev(colorRampPalette(brewer.pal(9,"Blues")[-(1:2)])(11)), as.character(0:10))
OrangeColorRangeUnsaturationsOpaque <- setNames(rev(colorRampPalette(brewer.pal(9,"Oranges")[-1])(11)), as.character(0:10))

#. LipidomicsQualityControlECNTestsUnsaturationInfoColored21012022.pdf #
pdf("C:/Users/Kevin/Documents/RData/Rest2/ACGData/FiguresByKT/LipidomicsQualityControlECNTestsUnsaturationInfoColored21012022.pdf")

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


#### Supplementary material: adducts distribution
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


pdf("C:/Users/Kevin/Documents/RData/Rest2/ACGData/FiguresByKT/LipidomicsQualityControlLTPLipidSpeciesCombinationsDistributionOfAdductsForTheScreensIndependentlyInOneGraph11032022.pdf")

barplot(MergedTablesForAdducts, beside = TRUE, las = 2, col = ColorMatrixTryOut["odd",], main = "Number of LTP - lipid species as different adducts", las = 1, xlab = "Number of adducts observed", ylab = "Number of LTPs - lipid species")
abline(h = seq(from = 50, to = 200, by = 50), col = "#FFFFFF33", lwd = 3.2)

dev.off()
# Previous: selected plot for inclusion in the figures

# What percentages covered with multiple adducts per screen #! #X
100-((MergedTablesForAdducts[,1]*100)/rowSums(MergedTablesForAdducts, na.rm = TRUE))

# in cellulo   in vitro 
# 26.07004   44.62617 

# What percentage covered with multiple adducts in both screens combined #! #X
100-(AdductTabelScreensCombined[1]*100)/sum(AdductTabelScreensCombined, na.rm = TRUE) # 40.94488%


# Independent visualization of previous figures

par()$mar # 5.1 4.1 4.1 2.1
pdf("C:/Users/Kevin/Documents/RData/Rest2/ACGData/FiguresByKT/ObservedNumberOfAdductsPerEachScreen15032022.pdf")


par(mar = c(12.8, 4.1, 4.1, 2))

barplot(sort(table(PureCombined32[PureCombined32[,"Screen"] == "in vivo","Adduct"]), decreasing = TRUE), col = ColorMatrixTryOut["odd",1], las = 2, ylab = "Observed number of LTP - lipid species pairs", main = "Distribution of adducts in cellulo")
barplot(sort(table(PureCombined32[PureCombined32[,"Screen"] == "in vitro","Adduct"]), decreasing = TRUE), col = ColorMatrixTryOut["odd",2], las = 2, ylab = "Observed number of LTP - lipid species pairs", main = "Distribution of adducts in vitro")

dev.off()
par(mar = c(5.3, 4.1, 4.1, 2))

# Combined visualization of lipid species as different adduct or in different screen
pdf("C:/Users/Kevin/Documents/RData/Rest2/ACGData/FiguresByKT/LipidomicsQualityControlLTPLipidSpeciesCombinationsDistribution20012022.pdf")

barplot(AdductTabelScreensCombined, col = "Grey", main = "Number of LTP - lipid species observations \n as different adduct or in different screen", las = 1, xlab = "Number of observations as adducts or in screens", ylab = "Number of LTPs - lipid species")
abline(h = seq(from = 50, to = 350, by = 50), col = "#FFFFFF33", lwd = 3.2)

dev.off()


######## Figure 2
#### Fig.2a barchart

#. TemporaryNewBarplotsForLipidomeOverviewGraph26042020.pdf #
#. Linear addition for things not in heatmaps some converted to percentages #!


# Carb distribution

CarbDistribution2403220 <- Col1ToRowNames(merge(aggregate(PureAntonella32b$Intensity, by = list(PureAntonella32b$TotalCarbonChainLength), FUN = sum), 
                                                aggregate(PureEnric32$Intensity, by = list(PureEnric32$TotalCarbonChainLength), FUN = sum), by = "Group.1", all = TRUE))

CarbDistribution2403220[setdiff(as.character(16:70), rownames(CarbDistribution2403220)),] <- 0
CarbDistribution2403220[is.na(CarbDistribution2403220)] <- 0

CarbDistribution2403220 <- CarbDistribution2403220[as.character(16:70),]
colnames(CarbDistribution2403220) <- c("in vivo", "in vitro")

ColorVectora <- c("white", brewer.pal(9,"Blues"))
ColorVectore <- c("white", brewer.pal(9,"Oranges"))


pdf("C:/Users/Kevin/Documents/RData/Rest2/ACGData/FiguresByKT/TemporaryNewBarplotsForLipidomeOverviewGraph26042020b.pdf")

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
# Previous figures: Different variants used to add upper part to figure 2a further in Adobe Illustrator

#### Fig.2a donut-diagram
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
pdf("C:/Users/Kevin/Documents/RData/Rest2/ACGData/FiguresByKT/Panel2CColorVsColorSidesFlippedOver09122021.pdf")

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

#### Fig.2a circle-legend
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


pdf("C:/Users/Kevin/Documents/RData/Rest2/ACGData/FiguresByKT/CirclesForTheLegendOfTheLipidomeGraph09102020.pdf",
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

#### (Fig.2b histograms at side) #X
#. LipidsIntensityFrationsOfEven26012022WithLysolipidsCollapsed27012022CleanedUpVersionWithTwoColorSchemesOddAndNoReducedToEmpty1403202222binsc.pdf (#)

#. LipidsIntensityFrationsOfEven26012022WithLysolipidsCollapsed27012022 #
library(reshape2)

CovertToHeatmapMatrixn <- function(InputDataStartMatrix, CompactionMatrixForTotal, HeadgroupOrder, CarbOrder, IonMode){
  reshapedstartofmatrix <- dcast(InputDataStartMatrix[InputDataStartMatrix$IonMode == IonMode,], LikelySubclass ~ TotalCarbonChainLength, value.var = "Intensity", fun.aggregate = sum)
  
  reshapedstartofmatrix2 <- reshapedstartofmatrix[(reshapedstartofmatrix$LikelySubclass != "VE" & reshapedstartofmatrix$LikelySubclass != "VA" & reshapedstartofmatrix$LikelySubclass != "P40"), colnames(reshapedstartofmatrix) != "NaN"]
  GenHead <- CompactionMatrixForTotal[match(reshapedstartofmatrix2$LikelySubclass, CompactionMatrixForTotal[,1]),2] # Changed CompactionMatrixForTotaln2 To CompactionMatrixForTotal #
  
  reshapedstartofmatrix4 <- rowsum(reshapedstartofmatrix2[,-1], GenHead)
  reshapedstartofmatrix4[,setdiff(CarbOrder, colnames(reshapedstartofmatrix4))] <- 0
  
  reshapedstartofmatrix4[setdiff(HeadgroupOrder, rownames(reshapedstartofmatrix4)),] <- 0
  return(reshapedstartofmatrix4[HeadgroupOrder, CarbOrder])
  
}



library("dplyr")

y <- PureAntonella32b # Changed from PureAntonella32 to PureAntonella32b #
Normalization1OfPureAntonella32 <- mutate(y, MinMaxRangeNorm = sapply(1:dim(y)[1], function(x){MinMaxNormFuncn(y$Intensity[x],y$IonMode[x],y)}))

y <- PureEnric32
Normalization1OfPureEnric32 <- mutate(y, MinMaxRangeNorm = sapply(1:dim(y)[1], function(x){MinMaxNormFuncn(y$Intensity[x],y$IonMode[x],y)}))

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
pdf("C:/Users/Kevin/Documents/RData/Rest2/ACGData/FiguresByKT/LipidsIntensityFrationsOfEven26012022WithLysolipidsCollapsed27012022.pdf")

barplot(t(EvenPercentagesForLipidsEmptycwl[rev(rownames(EvenPercentagesForLipidsEmptycwl)),2:1]), beside = TRUE, horiz = TRUE, las = 2, col = "white", border = "grey", xlim = c(0,1000), xaxt = "n") 
barplot(t(EvenPercentagesForLipidsAbsencescwl[rev(rownames(EvenPercentagesForLipidsAbsencescwl)),2:1]), add = TRUE, beside = TRUE, horiz = TRUE, las = 2, col = "grey", border = "grey", xlim = c(0,1000), xaxt = "n", yaxt = "n") 

barplot(t(EvenPercentagesForLipidscwl[rev(rownames(EvenPercentagesForLipidscwl)),2:1]), add = TRUE, beside = TRUE, horiz = TRUE, las = 2, col = ColorMatrixTryOut["odd",2:1], xlim = c(0,1000), xaxt = "n", yaxt = "n") 
axis(side = 1, at = c(0,25,50,75,100))

dev.off()
# Previous figure: basis for further clean-up in Adobe Illustrator: especially rastering and color schemes


#### Fig.2b circle heatmap visualization

#. CirclesToVisualizeLipidomeAndTheBackgroundEmptyCirclesForInVivoData150620205c2 #
#. CirclesToVisualizeLipidomeAndTheBackgroundEmptyCirclesForInVivoData150620205c2DifferentGreyColorsAndCircleOrderGridToBack2SameLineTicknessCombinedWorkDocument4PGPAbovePGGappedAndCollapsed4add.pdf (#)

library(RColorBrewer)
library(ComplexHeatmap)


# Import and cleaning of shotgun lipidomics data for HEK293 cells used

OverexpressionHEK <- read.csv(file = "C:/Users/Kevin/Documents/RData/Rest2/ACGData/FiguresByKT/ExtractedDataHEK293FromKenji08082019.txt", header = TRUE, sep = "\t", as.is = TRUE)
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
#! Maybe provide documents from this point on, instead of originals?

LipidConversionsBL <- cbind(unique(DetailsOfOverexpressionHEK2Split2$LipidClass), c("LPA", "LPC", "LPC", "LPE", "LPE", "LPG", "LPG", "LPI", "LPI", "LPS", "LPS", "GM3", "Chol", "SM", "PC", "PC", "HexCer", "CE", "DAG", "Cer", "Hex2Cer", "PE", "PE", "PG", "PI", "PI", "PS", "SHexCer", "CerP", "CL", "PA", "PA"))
DetailsOfOverexpressionHEK2Split2bl <- cbind(DetailsOfOverexpressionHEK2Split2, LipidConversionsBL[match(DetailsOfOverexpressionHEK2Split2$LipidClass, LipidConversionsBL[,1]),2])

DetailsOfOverexpressionHEK2Split2bl2 <- cbind(DetailsOfOverexpressionHEK2Split2bl, do.call("rbind", strsplit(x = DetailsOfOverexpressionHEK2Split2bl$CarbChainLengthAndUnsat, split = ":")))
colnames(DetailsOfOverexpressionHEK2Split2bl2)[4:6] <- c("Headgroup", "ChainLength", "Unsaturation")

library("reshape2")
DetailsOfOverexpressionHEK2Split2bl2WideVersion <- Col1ToRowNames(dcast(DetailsOfOverexpressionHEK2Split2bl2, Headgroup ~ ChainLength, value.var = "ControlMeans", fun.aggregate = sum))

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


pdf("C:/Users/Kevin/Documents/RData/Rest2/ACGData/FiguresByKT/CirclesToVisualizeLipidomeAndTheBackgroundEmptyCirclesForInVivoData150620205c2.pdf",
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
# Previous figure: Basis for putting rest of figure 2 together, with shades of grey changed in Adobe Illustrator afterwards, and PGP moved in location & gaps between groups

#### Fig.2c unsaturation-distribution
#. ComparisonOfUnsaturationPreferencesInCelluloVsInVitroVsCells140420222.pdf #

# Import of the HEK293 shotgun lipidomics data
LTPBackground1_2 <- read.csv(file = "C:/Users/Kevin/Documents/RData/Rest2/ACGData/FiguresByKT/Rest/Rest/DataKenjiShotgunLipidomicsOfLTPExpressingCells20032019/tableauOutput_20181023B_version2.txt", header = TRUE, sep = "\t", as.is = TRUE)

# Translate lipids manually outside R
write.table(LTPBackground1_2, file="./RData/Rest2/ACGData/FiguresByKT/LTPBackground1_217062021.tsv", sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Import again #! Maybe subset for the info and provide this?
LipidSubclassesAddedToBackground17062021 <- read.csv(file = "C:/Users/Kevin/Documents/RData/Rest2/ACGData/FiguresByKT/LipidSubclassesAddedToBackground17062021.txt", header = TRUE, sep = "\t", as.is = TRUE)

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

SpeciesOverviewLipidBodyConnections$Amount <- as.numeric(as.character(SpeciesOverviewLipidBodyConnections$Amount)) #! Added later as correction


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
pdf("C:/Users/Kevin/Documents/RData/Rest2/ACGData/FiguresByKT/ComparisonOfUnsaturationPreferencesInCelluloVsInVitroVsCells140420222.pdf")

for(x in names(ListLipidUnsatComp2)){
  y <- ListLipidUnsatComp2[[x]]
  
  barplot(t(y), beside = TRUE, col = c(brewer.pal(9,"Blues")[5], brewer.pal(9,"Oranges")[5], "grey"), las = 1, xlab = "", ylab = "") #, xaxt = "n", yaxt="n"
  abline(h = seq(5, ceiling(max(y)/10)*10, 5), col = "#FFFFFF52", lwd = 3.2)
  
  title(xlab = "Total unsaturation of lipids", line = 2.5)
  title(ylab = x, line = 2.5)
  
  legend("topright", legend = c("In Cellulo", "In Vitro", "Cellular"), col = c(brewer.pal(9,"Blues")[5], brewer.pal(9,"Oranges")[5], "grey"), pch = 15, bty = "n", cex = 0.95)
}

dev.off()
# Second page of previous figure-set: Used for panel in figure after some adaptations in Adobe Illustrator such as orientation and color scheme to fit with the other panels

######## Figure 3
#### Fig.3a

#. Panel3A09122021.pdf (#)
#. HeatmapOfTheLTPLocationsWithoutLipidInformationWithoutDomainsAdded251120214hdrBlackTrianglesAddedlwd1ZOrderEnhancedAndLegendAdded2.pdf (#)

#. (HeatmapOfTheLTPLocationsWithLipidInformationWithDomainsAdded11072020c.pdf # Domain and subcellular locations added: more extensive version that is not relevant for article, so not included) #!
#. HeatmapOfTheLTPLocationsWithoutLipidInformationWithoutDomainsAdded251120214hdrBlackTrianglesAddedlwd2 #


#### HPTLC-data import as dataframe with factors in columns #!

HPTLCDataInVivoAndInVitro <- read.csv(file = "C:/Users/Kevin/Documents/RData/Rest2/ACGData/FiguresByKT/KnownHPTLCResultsScraped31032020.tsv", header = FALSE, sep = "\t", as.is = FALSE)
colnames(HPTLCDataInVivoAndInVitro) <- c("LTPProtein", "Lipid", "Screen")


#### HPTLC-data: first clean-up and formatting #!

HPTLCDataInVivoAndInVitrob <- HPTLCDataInVivoAndInVitro[HPTLCDataInVivoAndInVitro[,1] != "SEC14L1",] # Removed because of a potential experimental issue
HPTLCDataInVivoAndInVitrob$LTPProtein <- factor(HPTLCDataInVivoAndInVitrob$LTPProtein, levels = levels(HPTLCDataInVivoAndInVitrob$LTPProtein)[levels(HPTLCDataInVivoAndInVitrob$LTPProtein) != "SEC14L1"])


HPTLCSpecificitiesPerScreen <- lapply(c("A","E"), function(x){Col1ToRowNames(dcast(HPTLCDataInVivoAndInVitrob[HPTLCDataInVivoAndInVitrob[,"Screen"] == x,], LTPProtein ~ Lipid, drop = FALSE, fun.aggregate = length)) })

library(reshape2)
HPTLCSpecificitiesPerScreen2 <- HPTLCSpecificitiesPerScreen

colnames(HPTLCSpecificitiesPerScreen2[[1]]) <- paste0(colnames(HPTLCSpecificitiesPerScreen2[[1]]), "*")
colnames(HPTLCSpecificitiesPerScreen2[[2]]) <- paste0(colnames(HPTLCSpecificitiesPerScreen2[[2]]), "*")

#### MS-data: clean-up and formatting #!
LTPProteins <- unique(c(unique(PureAntonella32b$LTPProtein), unique(PureEnric32$LTPProtein)))

HeadgroupOrderlnl <- c("d*Cer", "dCer", "DHCer", "DHOH*Cer", "tCer", "d*CerP", "d*HexCer", "t*HexCer", "t*Hex2Cer", "d*SHexCer", "d*SM", "DHSM", "t*SM", "FA", "FAL", "LPC",
                       "LPE", "LPE-O", "LPG", "DAG", "PA", "PC", "PC-O", "PE", "PE-O", "PI", "PS", "PG", "PG/BMP", "BMP", "TAG", "PGP", "CL", "VA")

CovertToHeatmapMatrixlnl <- function(InputDataStartMatrixlnl, HeadgroupOrderlnl, LTPOrder, IonMode){ # CompactionMatrixForTotal, 
  reshapedstartofmatrixlnl <- dcast(InputDataStartMatrixlnl[InputDataStartMatrixlnl$IonMode == IonMode,], LikelySubclass ~ LTPProtein, value.var = "Intensity", fun.aggregate = sum)
  
  reshapedstartofmatrixlnl2 <- reshapedstartofmatrixlnl[reshapedstartofmatrixlnl$LikelySubclass != "P40", colnames(reshapedstartofmatrixlnl) != "NaN"]
  
  
  reshapedstartofmatrixlnl4 <- Col1ToRowNames(reshapedstartofmatrixlnl2)
  reshapedstartofmatrixlnl4[,setdiff(LTPOrder, colnames(reshapedstartofmatrixlnl4))] <- 0
  
  reshapedstartofmatrixlnl4[setdiff(HeadgroupOrderlnl, rownames(reshapedstartofmatrixlnl4)),] <- 0
  return(reshapedstartofmatrixlnl4[HeadgroupOrderlnl, LTPOrder])
  
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

MeltedInVivoCombined <- rbind(melt(LTPMatrixTopInVivommn), melt(t(as.matrix(5*HPTLCSpecificitiesPerScreen2[[1]]))))
MeltedInVitroCombined <- rbind(melt(LTPMatrixTopInVitrommn), melt(t(as.matrix(5*HPTLCSpecificitiesPerScreen2[[2]]))))

CastInVivoCombined <- Col1ToRowNames(dcast(MeltedInVivoCombined, Var1 ~ Var2, value.var = "value"))
CastInVitroCombined <- Col1ToRowNames(dcast(MeltedInVitroCombined, Var1 ~ Var2, value.var = "value"))

CastInVivoCombined[is.na(CastInVivoCombined)] <- 0
CastInVitroCombined[is.na(CastInVitroCombined)] <- 0

# Make overview of LTP classes
MainDomainsOfTheLTPs <- unique(rbind(PureAntonella32b[,1:2], PureEnric32[,1:2])) 

MainDomainsOfTheLTPs2 <- cbind(LTPProtein = colnames(CastInVivoCombined), MainDomain = MainDomainsOfTheLTPs[match(colnames(CastInVivoCombined), MainDomainsOfTheLTPs[,1]),2])
MainDomainsOfTheLTPs2[36:43,2] <- "OSBP"

MainDomainsOfTheLTPs4 <- MainDomainsOfTheLTPs2[order(MainDomainsOfTheLTPs2[,2], MainDomainsOfTheLTPs2[,1]),][c(1:29,32:37,30:31,38:40,43,41,42),]


# Presence of LTP-lipid pairs for different technologies
MainDomainsOfTheLTPs5 <- as.data.frame(list(MainDomainsOfTheLTPs4,
                                            
                                            HPTLC = ifelse(MainDomainsOfTheLTPs4[,"LTPProtein"] %in% rownames(HPTLCSpecificitiesPerScreen2[[1]]),"Present","Absent"),
                                            LCMS = ifelse(MainDomainsOfTheLTPs4[,"LTPProtein"] %in% MainDomainsOfTheLTPs[,"LTPProtein"],"Present","Absent")))


# Extraction of information from Uniprot

BiocManager::install("UniprotR")
library(UniprotR) # Version 1.2.4

MainDomainsOfTheLTPs7 <- cbind(MainDomainsOfTheLTPs5, 
                               UniprotNames = c("Q86WG3", "Q7Z465", "P12271", "O76054", "Q9UDX3", "O43304", "B5MCN3", "P49638", "Q9BTX7", # "Q92503" removed
                                                
                                                "Q9NZD2", "Q5TA50", "Q00169", "P48739", "Q9UKF7", "P17213", "Q8N4F0", "P29373", "P07148", "Q01469",
                                                "O15540", "P31025", "Q6UWW0", "P02689", "P09455", "P02753", "P82980", "P17900", "P22059", "Q969R2",
                                                
                                                "Q9BXW6", "Q9H1P3", "Q9H0X9", "Q9BZF2", "Q9BZF1", "Q96SU4", "Q9BXB5", "Q9BXB4", "Q6YN16", "P22307",
                                                "Q9UJQ7", "Q9UKL6", "Q9Y365", "Q9Y5P4"))

MainDomainsOfTheLTPs8 <- cbind(MainDomainsOfTheLTPs7, GetFamily_Domains(MainDomainsOfTheLTPs7$UniprotNames))
# If issues with the GetFamily_Domains function: look at following commented out entries to work with same data.

# Make sure that we have an exact saved RDS-image of this file in time to avoid issues with Uniprot
# saveRDS(object = MainDomainsOfTheLTPs8, file = "C:/Users/Kevin/Documents/RData/Rest2/ACGData/FiguresByKT/BackUpOfMainDomainsOfTheLTPs8.rds")

# Load saved RDS-image of this file (if needed, otherwise this can be skipped)
# MainDomainsOfTheLTPs8FromTheStorage <- readRDS("C:/Users/Kevin/Documents/RData/Rest2/ACGData/FiguresByKT/BackUpOfMainDomainsOfTheLTPs8.rds")

# Check that nothing went wrong during the conversion (if needed, otherwise this can be skipped)
# identical(MainDomainsOfTheLTPs8FromTheStorage, MainDomainsOfTheLTPs8)

# Reassign reloaded variable to the old one (if needed, otherwise this can be skipped) (watch out to certainly do previous step first)
# MainDomainsOfTheLTPs8 <- MainDomainsOfTheLTPs8FromTheStorage

MainDomainsOfTheLTPs8LongVersionDomains <- do.call("rbind", lapply(1:dim(MainDomainsOfTheLTPs8)[1], function(i){do.call("cbind", list(as.character(MainDomainsOfTheLTPs8$LTPProtein[i]), "Domain", matrix(unlist(strsplit(as.character(MainDomainsOfTheLTPs8[i,"Domain..FT."]), split = "\\; ")), ncol = 3, byrow = TRUE)))}))
MainDomainsOfTheLTPs8LongVersionDomains2 <- cbind(MainDomainsOfTheLTPs8LongVersionDomains, do.call("rbind",strsplit(gsub("DOMAIN ", "", MainDomainsOfTheLTPs8LongVersionDomains[,3]), "\\.\\.")))

MainDomainsOfTheLTPs8LongVersionDomains4 <- cbind(MainDomainsOfTheLTPs8LongVersionDomains2, RegionName = gsub(" /note=", "", MainDomainsOfTheLTPs8LongVersionDomains2[,4]))  


MainDomainsOfTheLTPs8LongVersionCoils <- do.call("rbind", lapply(1:dim(MainDomainsOfTheLTPs8)[1], function(i){do.call("cbind", list(as.character(MainDomainsOfTheLTPs8$LTPProtein[i]), "CC", matrix(unlist(strsplit(as.character(MainDomainsOfTheLTPs8[i,"Coiled.coil"]), split = "\\; ")),ncol = 2, byrow = TRUE)))}))
AmountOfCoiledCoilsKnown <- sapply(strsplit(as.character(MainDomainsOfTheLTPs8[,"Coiled.coil"]), split = "\\; "), length)/2

MainDomainsOfTheLTPs8LongVersionCoils2 <- cbind(MainDomainsOfTheLTPs8LongVersionCoils, RegionName = unlist(sapply(AmountOfCoiledCoilsKnown, function(x){if(x != 0.5){paste0("CC",1:x)}else{NA}})))
MainDomainsOfTheLTPs8LongVersionCoils4 <- cbind(MainDomainsOfTheLTPs8LongVersionCoils2, do.call("rbind", strsplit(gsub("COILED ", "", MainDomainsOfTheLTPs8LongVersionCoils2[,3]), split = "\\.\\.")))

ListCompositionalBiasElements <- strsplit(MainDomainsOfTheLTPs8$Compositional.bias, "\\;")
MainDomainsOfTheLTPs8LongVersionBias <- do.call("cbind", list(as.character(MainDomainsOfTheLTPs8[,1]), "CompBias", do.call("rbind", lapply(ListCompositionalBiasElements,function(x){c(strsplit(gsub("COMPBIAS ", "", x[1]), split = "\\.\\.")[[1]], gsub("  /note=", "", x[2], fixed = TRUE))}))))

ListOfTheMotifParts <- strsplit(MainDomainsOfTheLTPs8$Motif, "\\;")
MainDomainsOfTheLTPs8LongVersionMotifs <- do.call("cbind", list(as.character(MainDomainsOfTheLTPs8[,1]), "Motif", do.call("rbind", lapply(ListOfTheMotifParts, function(x){if(is.na(x)){NA}else{c(strsplit(gsub("MOTIF ", "", x[1]), split = "\\.\\.")[[1]], gsub("  /note=", "", x[2], fixed =TRUE))}}))))

SequenceAndDomainHighlightsLTPs <- as.data.frame(do.call("rbind", list(MainDomainsOfTheLTPs8LongVersionDomains4[,c(1:2,6:8)], MainDomainsOfTheLTPs8LongVersionCoils4[,c(1,2,6,7,5)], MainDomainsOfTheLTPs8LongVersionBias, MainDomainsOfTheLTPs8LongVersionMotifs)))
colnames(SequenceAndDomainHighlightsLTPs) <- c("LTPProtein", "TypeRegion", "StartRegion", "StopRegion", "RegionName")

SequenceAndDomainHighlightsLTPs[,3] <- as.numeric(as.character(SequenceAndDomainHighlightsLTPs[,3]))
SequenceAndDomainHighlightsLTPs[,4] <- as.numeric(as.character(SequenceAndDomainHighlightsLTPs[,4]))

SequenceAndDomainHighlightsLTPs2 <- SequenceAndDomainHighlightsLTPs[!is.na(SequenceAndDomainHighlightsLTPs$StartRegion),]
library(reshape2)

CastProteinDomainsOfLTPs <- Col1ToRowNames(dcast(SequenceAndDomainHighlightsLTPs2, LTPProtein ~ RegionName, value.var = "StartRegion"))
CastProteinDomainsOfLTPs2 <- CastProteinDomainsOfLTPs

CastProteinDomainsOfLTPs2[setdiff(MainDomainsOfTheLTPs4[,"LTPProtein"], rownames(CastProteinDomainsOfLTPs2)),] <- NA
CastProteinDomainsOfLTPs4 <- CastProteinDomainsOfLTPs2[MainDomainsOfTheLTPs4[,"LTPProtein"],]


CastMainProteinDomainsLTPs <- Col1ToRowNames(dcast(as.data.frame(MainDomainsOfTheLTPs4), LTPProtein ~ MainDomain, fun.aggregate = length))[MainDomainsOfTheLTPs4[,"LTPProtein"],]

CastMainProteinDomainsLTPswnas <- CastMainProteinDomainsLTPs
CastMainProteinDomainsLTPswnas[CastMainProteinDomainsLTPswnas == 0] <- NA

CastManyProteinDomainsLTPs <- cbind(CastProteinDomainsOfLTPs4,CastMainProteinDomainsLTPswnas*round(mean(as.matrix(CastProteinDomainsOfLTPs4), na.rm = TRUE)))


OrderDecreasingPerMainDomain <- unlist(lapply(unique(MainDomainsOfTheLTPs4[,2]), function(x){names(sort(rowSums(CastManyProteinDomainsLTPs[MainDomainsOfTheLTPs4[MainDomainsOfTheLTPs4[,2] == x,1],], na.rm = TRUE), decreasing = TRUE))}))
MainDomainsOfTheLTPs4x <- as.data.frame(cbind(OrderDecreasingPerMainDomain, MainDomainsOfTheLTPs4[,2]))

# Sort by type of domains and split in the types
DomainTypes <- rbind(unique(SequenceAndDomainHighlightsLTPs2[,c(2,5)])[c(2,1,3,6,4,5,7:8,11:9,13,14,12),], cbind(TypeRegion = "MainDomain", RegionName = levels(MainDomainsOfTheLTPs4x[,2])[c(1:3, 9, 7, 4:6,8)]))

DomainTypes2 <- as.matrix(DomainTypes)
DomainTypes2[c(15,18),2] <- paste0(DomainTypes2[c(15,18),2],".1")


ListDomainsOfTheLTPs <- lapply(levels(MainDomainsOfTheLTPs4x[,2])[c(1:3, 9, 7, 4:6,8)], function(x){MainDomainsOfTheLTPs4x[MainDomainsOfTheLTPs4x[,2] == x,]})

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

#### HPTLC-data: further clean-up and formatting #!
HPTLCSpecificitiesPerScreen2hdr <- HPTLCSpecificitiesPerScreen2

colnames(HPTLCSpecificitiesPerScreen2hdr[[1]]) <- c("Cer*", "Sterol", "DAG", "PC", "PE", "PG", "PIPs", "PS")
colnames(HPTLCSpecificitiesPerScreen2hdr[[2]]) <- c("Cer*", "Sterol", "DAG", "PC", "PE", "PG", "PIPs", "PS")


library(reshape2)

HPTLCSpecificitiesPerScreen2hdrm <- lapply(HPTLCSpecificitiesPerScreen2hdr, function(x){melt(t(as.matrix(5*x)))})
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
}),] # Only leaves PS_OSBPL9: makes sense.

MeltedInVivoCombinedslchdr <- rbind(melt(LTPMatrixTopInVivommnslc), HPTLCSpecificitiesPerScreen2hdrm4[[1]]) #!
MeltedInVitroCombinedslchdr <- rbind(melt(LTPMatrixTopInVitrommnslc), HPTLCSpecificitiesPerScreen2hdrm4[[2]]) #!

CastInVivoCombinedslchdr <- as.matrix(Col1ToRowNames(dcast(MeltedInVivoCombinedslchdr, Var1 ~ Var2, value.var = "value", fun.aggregate = function(x){max(x,na.rm = TRUE)})))
CastInVitroCombinedslchdr <- as.matrix(Col1ToRowNames(dcast(MeltedInVitroCombinedslchdr, Var1 ~ Var2, value.var = "value", fun.aggregate = function(x){max(x,na.rm = TRUE)})))

CastInVivoCombinedslchdr[is.infinite(CastInVivoCombinedslchdr)] <- 0
CastInVitroCombinedslchdr[is.infinite(CastInVitroCombinedslchdr)] <- 0

all(colnames(CastInVitroCombinedslchdr) %in% rownames(DomainsByRowColumnAssociationsReordered)) # TRUE
all(rownames(DomainsByRowColumnAssociationsReordered) %in% colnames(CastInVitroCombinedslchdr)) # FALSE

rownames(DomainsByRowColumnAssociationsReordered)[!(rownames(DomainsByRowColumnAssociationsReordered) %in% colnames(CastInVitroCombinedslchdr))]
# Missing parts: "OSBP"    "OSBPL1A" "OSBPL10" "OSBP2"   "OSBPL8"  "OSBPL11" "OSBPL7"  "OSBPL2"


CastInVitroCombinedslchdr2 <- cbind(CastInVitroCombinedslchdr, matrix(data = 0, 
                                                                      
                                                                      nrow = nrow(CastInVitroCombinedslchdr),
                                                                      ncol = length(rownames(DomainsByRowColumnAssociationsReordered)[!(rownames(DomainsByRowColumnAssociationsReordered) %in% colnames(CastInVitroCombinedslchdr))]),
                                                                      
                                                                      dimnames = list(rownames(CastInVitroCombinedslchdr),
                                                                                      rownames(DomainsByRowColumnAssociationsReordered)[!(rownames(DomainsByRowColumnAssociationsReordered) %in% colnames(CastInVitroCombinedslchdr))])))

InVivoDataSetTotalslchdr <- as.matrix(CastInVivoCombinedslchdr[, rownames(DomainsByRowColumnAssociationsReordered)])
InVitroDataSetTotalslchdr <- as.matrix(CastInVitroCombinedslchdr2[, rownames(DomainsByRowColumnAssociationsReordered)]) 

c(rownames(InVivoDataSetTotalslchdr), rownames(InVitroDataSetTotalslchdr))[!unique(c(rownames(InVivoDataSetTotalslchdr), rownames(InVitroDataSetTotalslchdr))) %in% rownames(InVivoDataSetTotalslchdr)]
# character(0)

MissingLipidEntriesInVitro <- c(rownames(InVivoDataSetTotalslchdr), rownames(InVitroDataSetTotalslchdr))[!unique(c(rownames(InVivoDataSetTotalslchdr), rownames(InVitroDataSetTotalslchdr))) %in% rownames(InVitroDataSetTotalslchdr)]
# "Sterol" "PIPs"

MeltedInVivoCombinedslc <- rbind(melt(LTPMatrixTopInVivommnslc), melt(t(as.matrix(5*HPTLCSpecificitiesPerScreen2[[1]]))))
MeltedInVitroCombinedslc <- rbind(melt(LTPMatrixTopInVitrommnslc), melt(t(as.matrix(5*HPTLCSpecificitiesPerScreen2[[2]]))))

CastInVivoCombinedslc <- as.matrix(Col1ToRowNames(dcast(MeltedInVivoCombinedslc, Var1 ~ Var2, value.var = "value", fun.aggregate = function(x){max(x,na.rm = TRUE)})))
CastInVitroCombinedslc <- as.matrix(Col1ToRowNames(dcast(MeltedInVitroCombinedslc, Var1 ~ Var2, value.var = "value", fun.aggregate = function(x){max(x,na.rm = TRUE)})))

CastInVivoCombinedslc[is.infinite(CastInVivoCombinedslc)] <- 0
CastInVitroCombinedslc[is.infinite(CastInVitroCombinedslc)] <- 0

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

# Literature coverage #!
LiteratureDataLTPsFurtherIntegrated26052020 <- read.csv(file = "C:/Users/Kevin/Documents/RData/Rest2/ACGData/FiguresByKT/LiteratureDataLTPsFurtherIntegrated26052020.txt", header = TRUE, sep = "\t", as.is = TRUE)

LiteratureDataLTPsFurtherIntegrated26052020$InVivo <- as.numeric(gsub(",",".", LiteratureDataLTPsFurtherIntegrated26052020$InVivo))
LiteratureDataLTPsFurtherIntegrated26052020$InVitro <- as.numeric(gsub(",",".", LiteratureDataLTPsFurtherIntegrated26052020$InVitro))

# 29.6% Literature consensus
sum(LiteratureDataLTPsFurtherIntegrated26052020$LiteratureConsensus)/length(LiteratureDataLTPsFurtherIntegrated26052020$LiteratureConsensus) # 0.296

# Import of literature hits overview
ExpandedLiteratureHitsFromClassesOfLipidsWithoutRedundantStarsStored <- read.table(file="C:/Users/Kevin/Documents/RData/Rest2/ACGData/FiguresByKT/ExpandedLiteratureHitsFromClassesOfLipidsWithoutRedundantStarsStored05112022.csv", sep = "\t")

# Previous: LipidClassLiteratureDataSet2 <- !((InVivoDataSet == 0) & (InVitroDataSet == 0)) & (LipidClassLiteratureDataSet == 0) # Changed to unique entries avoid looping and avoid possible confusion
LipidClassLiteratureDataSet2 <- !((InVivoDataSetWithoutRedundantStars == 0) & (InVitroDataSetWithoutRedundantStars == 0)) & (as.matrix(ExpandedLiteratureHitsFromClassesOfLipidsWithoutRedundantStarsStored) == 0)

library(reshape2)
LiteratureConsensusLinksWideVersion <- Col1ToRowNames(dcast(data = LiteratureDataLTPsFurtherIntegrated26052020, Lipid ~ Protein, value.var = "LiteratureConsensus",fun.aggregate = sum))

LiteratureConsensusLinksWideVersion2 <- LiteratureConsensusLinksWideVersion[rownames(LipidClassLiteratureDataSet2), colnames(LipidClassLiteratureDataSet2)]
LiteratureConsensusLinksWideVersion4 <- !((InVivoDataSetWithoutRedundantStars == 0) & (InVitroDataSetWithoutRedundantStars == 0)) & (LiteratureConsensusLinksWideVersion2 == 0) # InVivoDataSetWithoutRedundantStars instead of InVivoDataSet to avoid loop and same for in vitro data


LipidClassLiteratureDataSet <- as.matrix(LiteratureConsensusLinksWideVersion4)[c(1:19,21:28,32,29:31,33,20,34,35,36),]

LipidClassLiteratureDataSetslc <- LipidClassLiteratureDataSet[c(1,6,7,9:11,14:36),]
rownames(LipidClassLiteratureDataSetslc) <- rownames(InVivoDataSetslc)

# Correction on original heatmap: literature data for OSBPL2-PIPs known now & PC-O and PE-O should not be seen as known only because PC and PE are known #!
# Also immediately correct the rownames

LipidClassLiteratureDataSetslc4hdr <- LipidClassLiteratureDataSetslc
rownames(LipidClassLiteratureDataSetslc) <- rownames(InVivoDataSetslc4hdr)

LipidClassLiteratureDataSetslc4hdr["PIPs","OSBPL2"] <- FALSE
LipidClassLiteratureDataSetslc4hdr["PC-O","STARD2"] <- TRUE

LipidClassLiteratureDataSetslc4hdr["PC-O","STARD10"] <- TRUE
LipidClassLiteratureDataSetslc4hdr["PC-O","GM2A"] <- TRUE

LipidClassLiteratureDataSetslc4hdr["PC-O","LCN1"] <- TRUE
LipidClassLiteratureDataSetslc4hdr["PE-O","STARD10"] <- TRUE


ORPDomains <- read.csv(file = "C:/Users/Kevin/Documents/RData/Rest2/ACGData/FiguresByKT/ORPDomains270520202.txt", header = TRUE, sep = "\t", as.is = TRUE)

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
SplitByColsVector <- factor(MainDomainsOfTheLTPs4xx[,2], levels = levels(MainDomainsOfTheLTPs4x[,2])[c(1:3, 9, 7, 4:6,8)])

GraphNameslc <- "White circles with black borders for novelty (vs. consensus literature knowledge)(sphingolipids aggregated)"
library(ComplexHeatmap)

pdf("C:/Users/Kevin/Documents/RData/Rest2/ACGData/FiguresByKT/HeatmapOfTheLTPLocationsWithoutLipidInformationWithoutDomainsAdded251120214hdrBlackTrianglesAddedlwd2.pdf",
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

# Previous figure: Basis for Figure 3a, with several stylistic updates in Adobe Illustrator and the inclusion of a legend.



#### Protein-domain-based ordering of LTPs

#. SupplementaryInformationOverviewSeriationOfProteinDomains061120222.pdf (#) #!
#. HeatmapListDomainsReorderedTypesDomainsEnhanced06112022.pdf #

library(ComplexHeatmap)
LegendName <- "Legend"

WidthAdaptor <- 160
HeightAdaptor <- 160

SplitByColsVector <- factor(MainDomainsOfTheLTPs4xx[,2], levels = levels(MainDomainsOfTheLTPs4x[,2])[c(1:3, 9, 7, 4:6,8)])
DomainsLTPsReorderedAfterSeriation <- colnames(DomainsByRowColumnAssociationsReordered5)[c(1:13,15,14,17,16,18:31)]

ReorderedDomainsInputMatrix <- t(as.matrix(DomainsByRowColumnAssociationsReordered5[ReorderedLTPsByManualSeriation, DomainsLTPsReorderedAfterSeriation]))
rownames(ReorderedDomainsInputMatrix)[rownames(ReorderedDomainsInputMatrix) %in% c("CRAL-TRIO", "Ankyrin", "lipocalin", "adh_short", "TMpd", "TMpd2")] <-  c("CRAL_TRIO", "Ankyrin repeats", "Lipocalin", "Adh_short", "TM1", "TM2")

ColFunctionForDomains <- colorRamp2(1:max(CastManyProteinDomainsLTPs, na.rm = TRUE), colorRampPalette(c("#fafa6e", "#2a4858"))(max(CastManyProteinDomainsLTPs, na.rm = TRUE)), space = "RGB")


DomainTypes5 <- cbind(TypeRegion = c(rep("Domain", 18), rep("Region", 9), rep("Motif", 4)),
                      RegionName = c("PRELI/MSF1", "CRAL_TRIO_N", "CRAL-TRIO", "GOLD", "BNIP2", "CRAL_TRIO_2", "GLTP", "IP_trans", "START", "PH", "ORD", "Ankyrin", "LBP_BPI_CETP", "lipocalin", "ML", "Thiolase", "adh_short","SCP2", "TMpd", "TMpd2", "CC1", "CC2", "CC3", "CC4", "Ser-rich", "Ala/Gly-rich", "Poly-Leu", "FFAT", "Signal peptide", "Nuclear localization signal", "PTS1"))

pdf("C:/Users/Kevin/Documents/RData/Rest2/ACGData/FiguresByKT/HeatmapListDomainsReorderedTypesDomainsEnhanced06112022.pdf",
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

#### Fig.3b
#. LinearBarplotsForNoveltiesLTPLipidClassPairsUpdatedVersionsAfterGroupMeetingCommentsPercentagesMovedAxisExtendedTextUpdated140320223b2.pdf 

# Largely manual. Adobe Illustator was used for most of visualization. #!


#### Fig.3c
#. PITPFocussedPIPCPARemakeGraph25022022LipidsClusterdTogetherWithThinAndThickLinesOverlappingAndWithGuideLinesAddedToDifferentLevelsSeveralFigures28022022.pdf (#)

#. PITPFocussedPIPCPARemakeStripedVerticalLinesGraph140320223d.pdf (#)
#. BarplotHeatmapPITransportersCondensedRowEntriesOnlyScreenColumnsAt5ProcentCutoffAndAnnotationRows12102021.pdf #


library(reshape2)

LTPLipidConnectionsDataSet <- rbind(PureAntonella32b, PureEnric32)
LipidSpeciesLTPSpecificList <- setNames(lapply(unique(LTPLipidConnectionsDataSet[,"LikelySubclass"]),
                                               
                                               function(y){dcast(aggregate(LTPLipidConnectionsDataSet[LTPLipidConnectionsDataSet[,"LikelySubclass"] == y, "Intensity"], 
                                                                           by=list(paste0(LTPLipidConnectionsDataSet[LTPLipidConnectionsDataSet[,"LikelySubclass"] == y, "LTPProtein"], "(", LTPLipidConnectionsDataSet[LTPLipidConnectionsDataSet[,"LikelySubclass"] == y, "Screen"], ")"),
                                                                                   
                                                                                   LTPLipidConnectionsDataSet[LTPLipidConnectionsDataSet[,"LikelySubclass"] == y, "Lipid"]), 
                                                                           FUN=function(x){sum(x, na.rm = TRUE)}),
                                                                 
                                                                 formula = Group.2 ~ Group.1, value.var = "x")}),
                                        nm = unique(LTPLipidConnectionsDataSet[,"LikelySubclass"]))

rowSumsDFS <- function(x){if(class(x) == "data.frame"){rowSums(x, na.rm = TRUE)}else{x}}
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

library(reshape2)
LTPLipidConnectionsDataAggregated <- dcast(LTPLipidConnectionsDataSubset, LTPProtein + LikelySubclass + Screen ~ CarbonChain, value.var="Intensity", fun.aggregate = function(x){sum(x, na.rm = TRUE)})

PITransportersComobilizedPatterns <- LTPLipidConnectionsDataAggregated[LTPLipidConnectionsDataAggregated[, "LTPProtein"] %in% unique(LTPLipidConnectionsDataAggregated[LTPLipidConnectionsDataAggregated$LikelySubclass == "PI", "LTPProtein"])[c(1,2,6,3,4,5)],]
PITransportersComobilizedPatterns2 <- cbind(PITransportersComobilizedPatterns[,1:3], PITransportersComobilizedPatterns[,4:100][,colSums(PITransportersComobilizedPatterns[,4:100]) != 0])

PITransportersComobilizedPatterns4 <- PITransportersComobilizedPatterns2
PITransportersComobilizedPatterns4[,-(1:3)] <- PITransportersComobilizedPatterns4[,-(1:3)]*100/do.call(pmax, PITransportersComobilizedPatterns4[,-(1:3)])

PITransportersComobilizedPatterns4 <- cbind(LTP_Lipid = paste(PITransportersComobilizedPatterns4[,1], PITransportersComobilizedPatterns4[,2], sep = "_"), PITransportersComobilizedPatterns4)
NewRownamesPITransporters <- unique(PITransportersComobilizedPatterns4$LTP_Lipid)

GroupingListPITransporters <- list(list(22, 21, 4 , c(1:3,5), 10, c(6:9), 13, 12, 11, 18, 15, 14, c(16:17), 20, 19),
                                   list(NA, NA, NA, c("BPI", NA, "in vivo"), NA, c("GM2A", NA, "in vivo"), NA, NA, NA, NA, NA, NA, c("PITPNB", "PG/BMP", "in vitro"), NA, NA))

PITransportersComobilizedPatterns2CondensedRowEntries <- as.matrix(do.call("rbind", lapply(1:length(GroupingListPITransporters[[1]]), function(x){
  if(length(GroupingListPITransporters[[1]][[x]]) == 1){
    
    PITransportersComobilizedPatterns2[GroupingListPITransporters[[1]][[x]],]
  }else{
    
    c(GroupingListPITransporters[[2]][[x]], colSums(PITransportersComobilizedPatterns2[GroupingListPITransporters[[1]][[x]],-(1:3)]))
    
    
  }
})))

PITransportersComobilizedPatterns2CondensedRowEntries[4,2] <- "Other"
PITransportersComobilizedPatterns2CondensedRowEntries[6,2] <- "Other"

PITransportersComobilizedPatterns2CondensedRowEntriesNumericSubsetData <- PITransportersComobilizedPatterns2CondensedRowEntries[,-(1:3)]
mode(PITransportersComobilizedPatterns2CondensedRowEntriesNumericSubsetData) <- "numeric" 

PITransportersComobilizedPatterns4CondensedRowEntries <- cbind(PITransportersComobilizedPatterns2CondensedRowEntries[,1:3], as.data.frame(PITransportersComobilizedPatterns2CondensedRowEntriesNumericSubsetData))
PITransportersComobilizedPatterns4CondensedRowEntries[,-(1:3)] <- PITransportersComobilizedPatterns4CondensedRowEntries[,-(1:3)]*100/do.call(pmax, PITransportersComobilizedPatterns4CondensedRowEntries[,-(1:3)])

PITransportersComobilizedPatterns4CondensedRowEntries <- cbind(LTP_Lipid = paste(PITransportersComobilizedPatterns4CondensedRowEntries[,1], PITransportersComobilizedPatterns4CondensedRowEntries[,2], sep = "_"), PITransportersComobilizedPatterns4CondensedRowEntries)
NewRownamesPITransportersCondensedRowEntries <- unique(PITransportersComobilizedPatterns4CondensedRowEntries$LTP_Lipid)


library("plyr")

CellularListPITransporters <- lapply(sapply(strsplit(sapply(strsplit(as.character(NewRownamesPITransporters), "_"), "[[", 2), "/"), "[[", 1), function(x){rbind(rownames(LipidSubclassesAddedToBackground170620214[[x]]), 
                                                                                                                                                            ZerosToNAsConverter(LipidSubclassesAddedToBackground170620214[[x]][,"Cellular"]))})

CellularListPITransportersVectorsWithoutTheNAs <- lapply(1:length(CellularListPITransporters), function(x){as.data.frame(t(setNames(as.numeric(CellularListPITransporters[[x]][2,!is.na(CellularListPITransporters[[x]][2,])]),
                                                                                                                                    gsub("O-","",sapply(strsplit(gsub(")","",CellularListPITransporters[[x]][1,!is.na(CellularListPITransporters[[x]][2,])]), "\\("), "[[", 2)))))})

CarbonChainColNames <- sort(unique(c(colnames(PITransportersComobilizedPatterns4)[-(1:4)], unlist(lapply(CellularListPITransportersVectorsWithoutTheNAs, names)))))[c(1:54,56:60,55)]
PITransportersComobilizedPatterns4WithCellularAdded <- ZerosToNAsConverter(rbind.fill(PITransportersComobilizedPatterns4, do.call("cbind", list(LTP_Lipid = "cellular",
                                                                                                                                                
                                                                                                                                                LTPProtein = "cellular",
                                                                                                                                                LikelySubclass = sapply(strsplit(sapply(strsplit(as.character(NewRownamesPITransporters), "_"), "[[", 2), "/"), "[[", 1),
                                                                                                                                                
                                                                                                                                                Screen = "cellular",
                                                                                                                                                do.call("rbind.fill", CellularListPITransportersVectorsWithoutTheNAs)))))[c("LTP_Lipid", "LTPProtein", "LikelySubclass", "Screen", CarbonChainColNames)]


PITransportersComobilizedPatterns4WithCellularAddedCondensedRowEntries <- ZerosToNAsConverter(rbind.fill(PITransportersComobilizedPatterns4CondensedRowEntries, do.call("cbind", list(LTP_Lipid = "cellular",
                                                                                                                                                                                      
                                                                                                                                                                                      LTPProtein = "cellular",
                                                                                                                                                                                      LikelySubclass = sapply(strsplit(sapply(strsplit(as.character(NewRownamesPITransporters), "_"), "[[", 2), "/"), "[[", 1),
                                                                                                                                                                                      
                                                                                                                                                                                      Screen = "cellular",
                                                                                                                                                                                      do.call("rbind.fill", CellularListPITransportersVectorsWithoutTheNAs)))))[c("LTP_Lipid", "LTPProtein", "LikelySubclass", "Screen", CarbonChainColNames)]

PITransportersComobilizedPatterns4WithCellularAdded2 <- unique(PITransportersComobilizedPatterns4WithCellularAdded)
PITransportersComobilizedPatterns4WithCellularAddedCondensedRowEntries2 <- unique(PITransportersComobilizedPatterns4WithCellularAddedCondensedRowEntries)


EmptyLayerPITransporters <- matrix(NA, 
                                   
                                   nrow = length(NewRownamesPITransporters), 
                                   ncol = length(CarbonChainColNames), 
                                   
                                   dimnames = list(NewRownamesPITransporters,
                                                   CarbonChainColNames))

NewPITransLayers <- list(InCellulo = EmptyLayerPITransporters,
                         InVitro = EmptyLayerPITransporters,
                         
                         Cellular = EmptyLayerPITransporters)



EmptyLayerPITransportersCondensedRowEntries <- matrix(NA, 
                                                      
                                                      nrow = length(NewRownamesPITransportersCondensedRowEntries), 
                                                      ncol = length(CarbonChainColNames), 
                                                      
                                                      dimnames = list(NewRownamesPITransportersCondensedRowEntries,
                                                                      CarbonChainColNames))

NewPITransLayersCondensedRowEntries <- list(InCellulo = EmptyLayerPITransportersCondensedRowEntries,
                                            InVitro = EmptyLayerPITransportersCondensedRowEntries,
                                            
                                            Cellular = EmptyLayerPITransportersCondensedRowEntries)


for(i in rownames(NewPITransLayers$InCellulo)){
  if(any((PITransportersComobilizedPatterns4WithCellularAdded2$Screen == "in vivo") & (PITransportersComobilizedPatterns4WithCellularAdded2$LTP_Lipid == i))){
    
    NewPITransLayers$InCellulo[i,] <- as.numeric(ZerosToNAsConverter(PITransportersComobilizedPatterns4WithCellularAdded2[(PITransportersComobilizedPatterns4WithCellularAdded2$Screen == "in vivo") & (PITransportersComobilizedPatterns4WithCellularAdded2$LTP_Lipid == i),-(1:4)]))
  }}

for(i in rownames(NewPITransLayers$InVitro)){
  if(any((PITransportersComobilizedPatterns4WithCellularAdded2$Screen == "in vitro") & (PITransportersComobilizedPatterns4WithCellularAdded2$LTP_Lipid == i))){
    
    NewPITransLayers$InVitro[i,] <- as.numeric(ZerosToNAsConverter(PITransportersComobilizedPatterns4WithCellularAdded2[(PITransportersComobilizedPatterns4WithCellularAdded2$Screen == "in vitro") & (PITransportersComobilizedPatterns4WithCellularAdded2$LTP_Lipid == i),-(1:4)]))
  }}

for(i in rownames(NewPITransLayers$Cellular)){
  NewPITransLayers$Cellular[i,] <- as.numeric(ZerosToNAsConverter(PITransportersComobilizedPatterns4WithCellularAdded2[(PITransportersComobilizedPatterns4WithCellularAdded2$Screen == "cellular") & (PITransportersComobilizedPatterns4WithCellularAdded2$LikelySubclass == strsplit(strsplit(i,"_")[[1]][2],"/")[[1]][1]),-(1:4)]))
  
}


for(i in rownames(NewPITransLayersCondensedRowEntries$InCellulo)){
  if(any((PITransportersComobilizedPatterns4WithCellularAddedCondensedRowEntries2$Screen == "in vivo") & (PITransportersComobilizedPatterns4WithCellularAddedCondensedRowEntries2$LTP_Lipid == i))){
    
    NewPITransLayersCondensedRowEntries$InCellulo[i,] <- as.numeric(ZerosToNAsConverter(PITransportersComobilizedPatterns4WithCellularAddedCondensedRowEntries2[(PITransportersComobilizedPatterns4WithCellularAddedCondensedRowEntries2$Screen == "in vivo") & (PITransportersComobilizedPatterns4WithCellularAddedCondensedRowEntries2$LTP_Lipid == i),-(1:4)]))
  }}

for(i in rownames(NewPITransLayersCondensedRowEntries$InVitro)){
  if(any((PITransportersComobilizedPatterns4WithCellularAddedCondensedRowEntries2$Screen == "in vitro") & (PITransportersComobilizedPatterns4WithCellularAddedCondensedRowEntries2$LTP_Lipid == i))){
    
    NewPITransLayersCondensedRowEntries$InVitro[i,] <- as.numeric(ZerosToNAsConverter(PITransportersComobilizedPatterns4WithCellularAddedCondensedRowEntries2[(PITransportersComobilizedPatterns4WithCellularAddedCondensedRowEntries2$Screen == "in vitro") & (PITransportersComobilizedPatterns4WithCellularAddedCondensedRowEntries2$LTP_Lipid == i),-(1:4)]))
  }}

GroupingListPITransportersAtLaterStage <- list(20, 19, 4 , c(1:3,5), 10, c(6:9), 12, 11, 16, 13, c(14:15), 18, 17)
NewPITransLayersCondensedRowEntries$Cellular <- do.call("rbind", lapply(GroupingListPITransportersAtLaterStage, function(x){
  
  if(length(unlist(x)) == 1){
    NewPITransLayers$Cellular[unlist(x),]
    
  }else{
    colSums(NewPITransLayers[["Cellular"]][unlist(x),], na.rm = TRUE)*100/max(colSums(NewPITransLayers[["Cellular"]][unlist(x),], na.rm = TRUE))
    
  }
}))

rownames(NewPITransLayersCondensedRowEntries$Cellular) <- rownames(NewPITransLayersCondensedRowEntries$InCellulo)
NewPITransLayersCondensedRowEntries$Cellular <- ZerosToNAsConverter(NewPITransLayersCondensedRowEntries$Cellular)
  
LegendName <- "Legend"
LegendColor <- col_fung

WidthAdaptor <- 160
HeightAdaptor <- 160

TextDataframeWithLTPAnnotation <- data.frame(LTP = factor(sapply(strsplit(rownames(NewPITransLayersCondensedRowEntries$InCellulo), "_"), "[[", 1), levels = unique(sapply(strsplit(rownames(NewPITransLayersCondensedRowEntries$InCellulo), "_"), "[[", 1))))
AnnotationDataframeWithLTPAnnotation <- rowAnnotation(df = TextDataframeWithLTPAnnotation, col = list(LTP = c("SEC14L2" = "#7E549F", "BPI" = "#FB836F", "GM2A" = "#C1549C", "PITPNA" = "#FFCB3E", "PITPNB" = "#E0B01C", "PITPNC1" = "#A47C00")))

CutOffPresenceProcent <- 5
CSPIDCRE <- (colSums(NewPITransLayersCondensedRowEntries$InCellulo, na.rm = TRUE) > CutOffPresenceProcent)|(colSums(NewPITransLayersCondensedRowEntries$InVitro, na.rm = TRUE) > CutOffPresenceProcent)

library(ComplexHeatmap)
library(RColorBrewer)

pdf(paste0("C:/Users/Kevin/Documents/RData/Rest2/ACGData/FiguresByKT/BarplotHeatmapPITransportersCondensedRowEntriesOnlyScreenColumnsAt", CutOffPresenceProcent,"ProcentCutoffAndAnnotationRows12102021.pdf"),
    width = unit(20, "mm"), height = unit(20, "mm"))

Heatmap(NewPITransLayersCondensedRowEntries$InCellulo[,CSPIDCRE], name = LegendName, col = LegendColor, rect_gp = gpar(type = "none"), column_title = "", 
        width = unit(WidthAdaptor*dim(NewPITransLayersCondensedRowEntries$InCellulo[,CSPIDCRE])[2]/min(dim(NewPITransLayersCondensedRowEntries$InCellulo[,CSPIDCRE])), "mm"), height = unit(HeightAdaptor*dim(NewPITransLayersCondensedRowEntries$InCellulo[,CSPIDCRE])[1]/min(dim(NewPITransLayersCondensedRowEntries$InCellulo[,CSPIDCRE])), "mm"),
        
        
        cell_fun = function(j, i, x, y, width, height, fill) {
          
          grid.rect(x = x, y = y-0.5*height, width = width, height = NewPITransLayersCondensedRowEntries$InCellulo[,CSPIDCRE][i,j]/100*height, gp = gpar(col = NA, fill = brewer.pal(9,"Blues")[5], lwd = 1), just = c("center","bottom"))
          grid.rect(x = x, y = y-0.5*height, width = 0.5*width, height = NewPITransLayersCondensedRowEntries$InVitro[,CSPIDCRE][i,j]/100*height, gp = gpar(col = NA, fill = brewer.pal(9,"Oranges")[5], lwd = 1), just = c("center","bottom"))
          
          grid.rect(x = x, y = y-0.5*height, width = width, height = NewPITransLayersCondensedRowEntries$Cellular[,CSPIDCRE][i,j]/100*height, gp = gpar(col = "LightGrey", fill = NA, lwd = 2), just = c("center","bottom"))
        }, cluster_rows = FALSE, cluster_columns = FALSE, 
        

        cluster_row_slices = FALSE, cluster_column_slices = FALSE,
        
        row_title = " ", row_title_gp = gpar(fontsize = 10),
        left_annotation = AnnotationDataframeWithLTPAnnotation,
        
        row_labels = sapply(strsplit(rownames(NewPITransLayersCondensedRowEntries$InCellulo), "_"), "[[", 2),
        row_split = c(rep("A", 2), rep("B",2), rep("C",2), rep("D",2), rep("E",3), rep("F",2))
        
)
dev.off()

# Previous graph further in Adobe Illustrator used to subset for panel figure 3c and also for figure of SEC14L2
# Note: only anything above 5% is visualized here to reduce the complexity of the figures

######## Figure 4
#### Fig.4a

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



library(reshape2)

HPTLCaggregated <- rbind(cbind(melt(as.matrix(5*HPTLCSpecificitiesPerScreen2[[1]])), ProteinDomain = NA, Screen = "in vivo", TotalCarbonChainLength = NA, TotalCarbonChainUnsaturations = NA),
                         cbind(melt(as.matrix(5*HPTLCSpecificitiesPerScreen2[[2]])), ProteinDomain = NA, Screen = "in vitro", TotalCarbonChainLength = NA, TotalCarbonChainUnsaturations = NA))

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


AggregatedAllScreenDataNormalized4 <- dcast(AggregatedAllScreenDataNormalized2, LTPProtein + ProteinDomain + TotalCarbonChainLength + TotalCarbonChainUnsaturations + LikelySubclass ~ Screen, value.var = "x")
AggregatedAllScreenDataNormalized4[,1] <- factor(AggregatedAllScreenDataNormalized4[,1], levels = ReorderedLTPsByManualSeriation)

AggregatedAllScreenDataNormalized4$LikelySubclass <- factor(AggregatedAllScreenDataNormalized4$LikelySubclass, levels = c("Cer*", "d*Cer", "dCer", "DHCer", "DHOH*Cer", "tCer", "d*CerP", "d*HexCer", "t*HexCer", "t*Hex2Cer", "d*SHexCer", "d*SM", "DHSM", "t*SM", "FA",       
                                                                                                                          "FAL", "LPC", "LPE", "LPE-O", "LPG", "PA", "PC", "PC-O", "PE", "PE-O", "PI", "PIPs", "PS", "PGP", "PG", "PG/BMP", "BMP", "CL", "DAG", "TAG", "Sterol", "VA") )

AggregatedAllScreenDataNormalized4$BandSize <- 1 


HPTLCRowsInData <- which(AggregatedAllScreenDataNormalized4$TotalCarbonChainLength == 0)
AggregatedAllScreenDataNormalized4ifh <- AggregatedAllScreenDataNormalized4

HPTLCRowsInDataNAsToBeIntroduced <- do.call("cbind", list(AggregatedAllScreenDataNormalized4ifh[HPTLCRowsInData,],
                                                          RowToBeEvaluated = HPTLCRowsInData,
                                                          
                                                          InVitroRemovalData = sapply(HPTLCRowsInData, function(x){paste(AggregatedAllScreenDataNormalized4ifh[x,1], AggregatedAllScreenDataNormalized4ifh[x,5], sep = "_") %in% paste(AggregatedAllScreenDataNormalized4ifh[-x,][!is.na(AggregatedAllScreenDataNormalized4ifh[-x,][,"in vitro"]), ][,1], AggregatedAllScreenDataNormalized4ifh[-x,][!is.na(AggregatedAllScreenDataNormalized4ifh[-x,][,"in vitro"]), ][,5], sep = "_")}),
                                                          InVivoRemovalData = sapply(HPTLCRowsInData, function(x){paste(AggregatedAllScreenDataNormalized4ifh[x,1], AggregatedAllScreenDataNormalized4ifh[x,5], sep = "_") %in% paste(AggregatedAllScreenDataNormalized4ifh[-x,][!is.na(AggregatedAllScreenDataNormalized4ifh[-x,][,"in vivo"]), ][,1], AggregatedAllScreenDataNormalized4ifh[-x,][!is.na(AggregatedAllScreenDataNormalized4ifh[-x,][,"in vivo"]), ][,5], sep = "_")})))

AggregatedAllScreenDataNormalized4ifh[HPTLCRowsInDataNAsToBeIntroduced[HPTLCRowsInDataNAsToBeIntroduced[,"InVitroRemovalData"],"RowToBeEvaluated"], "in vitro"] <- NA
AggregatedAllScreenDataNormalized4ifh[HPTLCRowsInDataNAsToBeIntroduced[HPTLCRowsInDataNAsToBeIntroduced[,"InVivoRemovalData"],"RowToBeEvaluated"], "in vivo"] <- NA

AggregatedAllScreenDataNormalized4hdr <- AggregatedAllScreenDataNormalized4ifh[!(is.na(AggregatedAllScreenDataNormalized4ifh[,"in vivo"]) & is.na(AggregatedAllScreenDataNormalized4ifh[,"in vitro"])),]


InCelluloInVitroPresence <- (InVivoDataSetslc != 0)|(InVitroDataSetslc != 0)
mode(InCelluloInVitroPresence) <- "numeric"


# Input of lipid colors

LipidColorsForCircos <- read.csv(file = "C:/Users/Kevin/Documents/RData/Rest2/ACGData/FiguresByKT/LipidColorsUpdatesForCircos011120212.txt", header = TRUE, sep = "\t", as.is = TRUE)[,1:2]
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

library(circlize)
pdf("C:/Users/Kevin/Documents/RData/Rest2/ACGData/FiguresByKT/CircosExtendedWithSpecies051120217WithTracksSwitchedDecreasedTextSize0412WhiteRingRemovedSmallerScreenTracksWithoutInnerHalfringByDiffHeight0ColorSchemeUpdatedWithoutSomeTLCEntriesWithBlackOutsides7b.pdf")

circos.clear()
circos.par(start.degree = -70, clock.wise = FALSE, cell.padding = c(0,0,0,0))

testchord2 <- chordDiagram(AggregatedAllScreenDataNormalized4hdr[,c(5,1,8)], annotationTrack = NULL, grid.col = LipidColorsForCircos4b, directional = -1, diffHeight = mm_h(0),
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

cdmresaddedvalues <- cbind(testchord2, AggregatedAllScreenDataNormalized4hdr)
colnames(cdmresaddedvalues)[colnames(cdmresaddedvalues) == "in vivo"] <- "in cellulo"


i <- 10

for(i in seq_len(nrow(cdmresaddedvalues))) {
  if(cdmresaddedvalues$value1[i] > 0) {
    
    circos.rect(cdmresaddedvalues[i, "x1"], y1, cdmresaddedvalues[i, "x1"] - abs(cdmresaddedvalues[i, "value1"]), y1 + (y2-y1)*0.45, 
                col = "#ededed", 
                
                border = "#ededed",
                sector.index = cdmresaddedvalues$rn[i], track.index = 3)
    
    circos.rect(cdmresaddedvalues[i, "x1"], y1 + (y2-y1)*0.55, cdmresaddedvalues[i, "x1"] - abs(cdmresaddedvalues[i, "value1"]), y2, 
                col = "#ededed", 
                
                border = "#ededed",
                sector.index = cdmresaddedvalues$rn[i], track.index = 3)
    
    circos.rect(cdmresaddedvalues[i, "x2"], y1, cdmresaddedvalues[i, "x2"] - abs(cdmresaddedvalues[i, "value1"]), y1 + (y2-y1)*0.45, 
                col = "#ededed", 
                
                border = "#ededed",
                sector.index = cdmresaddedvalues$cn[i], track.index = 3)
    
    circos.rect(cdmresaddedvalues[i, "x2"], y1 + (y2-y1)*0.55, cdmresaddedvalues[i, "x2"] - abs(cdmresaddedvalues[i, "value1"]), y2, 
                col = "#ededed", 
                
                border = "#ededed",
                sector.index = cdmresaddedvalues$cn[i], track.index = 3)
    
    circos.rect(cdmresaddedvalues[i, "x1"], y1, cdmresaddedvalues[i, "x1"] - abs(cdmresaddedvalues[i, "value1"]), y1 + (y2-y1)*0.45, 
                col = col_funo(cdmresaddedvalues[i, "in vitro"]), 
                
                border = col_funo(cdmresaddedvalues[i, "in vitro"]),
                sector.index = cdmresaddedvalues$rn[i], track.index = 3)
    
    circos.rect(cdmresaddedvalues[i, "x1"], y1 + (y2-y1)*0.55, cdmresaddedvalues[i, "x1"] - abs(cdmresaddedvalues[i, "value1"]), y2, 
                col = col_funb(cdmresaddedvalues[i, "in cellulo"]), 
                
                border = col_funb(cdmresaddedvalues[i, "in cellulo"]),
                sector.index = cdmresaddedvalues$rn[i], track.index = 3)
    
    circos.rect(cdmresaddedvalues[i, "x2"], y1, cdmresaddedvalues[i, "x2"] - abs(cdmresaddedvalues[i, "value1"]), y1 + (y2-y1)*0.45, 
                col = col_funo(cdmresaddedvalues[i, "in vitro"]), 
                
                border = col_funo(cdmresaddedvalues[i, "in vitro"]),
                sector.index = cdmresaddedvalues$cn[i], track.index = 3)
    
    circos.rect(cdmresaddedvalues[i, "x2"], y1 + (y2-y1)*0.55, cdmresaddedvalues[i, "x2"] - abs(cdmresaddedvalues[i, "value1"]), y2, 
                col = col_funb(cdmresaddedvalues[i, "in cellulo"]), 
                
                border = col_funb(cdmresaddedvalues[i, "in cellulo"]),
                sector.index = cdmresaddedvalues$cn[i], track.index = 3)
    
    
    if((cdmresaddedvalues[i, "TotalCarbonChainLength"] == 0) & !is.na(cdmresaddedvalues[i, "in vitro"]) & (cdmresaddedvalues[i, "rn"] != "VA")){
      
      circos.rect(cdmresaddedvalues[i, "x1"], y1-(y2-y1)*0.2, cdmresaddedvalues[i, "x1"] - abs(cdmresaddedvalues[i, "value1"]), y1, 
                  col = "black", 
                  
                  border = "black",
                  sector.index = cdmresaddedvalues$rn[i], track.index = 3)
      
      circos.rect(cdmresaddedvalues[i, "x2"], y1-(y2-y1)*0.2, cdmresaddedvalues[i, "x2"] - abs(cdmresaddedvalues[i, "value1"]), y1, 
                  col = "black", 
                  
                  border = "black",
                  sector.index = cdmresaddedvalues$cn[i], track.index = 3)
      
    }
    
    
    
    if((cdmresaddedvalues[i, "TotalCarbonChainLength"] == 0) & !is.na(cdmresaddedvalues[i, "in cellulo"]) & (cdmresaddedvalues[i, "rn"] != "VA")){
      
      circos.rect(cdmresaddedvalues[i, "x1"], (y2+0.2), cdmresaddedvalues[i, "x1"] - abs(cdmresaddedvalues[i, "value1"]), y2, 
                  col = "black", 
                  
                  border = "black",
                  sector.index = cdmresaddedvalues$rn[i], track.index = 3)
      
      circos.rect(cdmresaddedvalues[i, "x2"], (y2+0.2), cdmresaddedvalues[i, "x2"] - abs(cdmresaddedvalues[i, "value1"]), y2, 
                  col = "black", 
                  
                  border = "black",
                  sector.index = cdmresaddedvalues$cn[i], track.index = 3)
      
    }
  }  
  
}


dev.off()
# Previous figure: Basis for figure panel 4a, after some stylistic enhancements in Adobe Illustrator, such as changes (of orientation of) some labels and optimization of some of the colors

#### (Subclasses version: not used in final figure.)
#. CircosExtended031120215SterolChangedAndTracksSwappedDifferentTrackColorsLinkColorsUpdatedRemovalWhiteSectionsWithoutSomeTLCEntries24112021BlackOuterBlocksWithLegendAdded.pdf (#)

#. CircosExtended031120215SterolChangedAndTracksSwappedDifferentTrackColorsLinkColorsUpdatedRemovalWhiteSectionsWithoutSomeTLCEntries24112021BlackOuterBlocks.pdf #


InCelluloInVitroPresence2 <- InCelluloInVitroPresence
rownames(InCelluloInVitroPresence2)[rownames(InCelluloInVitroPresence2) == "CH"] <- "Sterol"

InCelluloInVitroPresence2hdr <- InCelluloInVitroPresence2


# Updated version of sub-class - LTP connections circos with changed colors and removed white circles for TLC-reduced dataset, based on data that only has TLC if no other data there for the specific combo
pdf("C:/Users/Kevin/Documents/RData/Rest2/ACGData/FiguresByKT/CircosExtended031120215SterolChangedAndTracksSwappedDifferentTrackColorsLinkColorsUpdatedRemovalWhiteSectionsWithoutSomeTLCEntries24112021BlackOuterBlocks.pdf")

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

for(i in seq_len(nrow(cdm_res))) {
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
}

dev.off()
# The previous figure is the alternative version of the circos focussed on only subclasses, instead of species, which was not used in the final figures of the article

#### Fig.4b
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
pdf("C:/Users/Kevin/Documents/RData/Rest2/ACGData/FiguresByKT/LTPSubclassDistributionMSAndTLCBothScreensInBargraphPlot170120222021AmountToNumber18012022.pdf")

barplot(LTPSubclassDistributionMSAndTLCBothScreens2, beside = TRUE, col = ColorMatrixTryOut["odd",], main = "For MS with HPTLC data", las = 1, xlab = "Number of lipid subclasses bound", ylab = "Number of LTPs") 
abline(h = 1:max(LTPSubclassDistributionMSAndTLCBothScreens2), col = "#FFFFFF33", lwd = 3.2)

legend("topright", legend = c("in cellulo", "in vitro"), col = ColorMatrixTryOut["odd",], pch = 15, bty = "n", cex = 1)
dev.off()

# Previous figure was used for panel 4b: only y-axis further extended


######## Figure 5
# Comparison of lipid co-mobilization with co-regulation and co-localization in tissues and subcellularly

# Note: first we describe the figures used for making the panels of figure 5, 
# afterwards we also describe an alternative version that is not part of the final figure 5 with Nrd0 in bandwidth selection for the Gaussian kernels for the densities.

# The final sets of figures that were used for the panels of figure 5 are: Panel5AWithTailoredScriptAndMoreDetailedVisualizationNowIncludesAllTooWithBinSize00117022022NoSegmentStacking18032022.pdf;
# Panel5BWithTailoredScriptAndMoreDetailedVisualizationNowIncludesAllTooWithBinSize001170220222NoSegmentStacking18032022.pdf; Panel5CWithTailoredScriptAndMoreDetailedVisualizationBinSize001NoSegmentStacking18032022.pdf


#### Panel 5A

#. Panel5AWithTheNrd0WdAdaptedWithScalesImplemented.pdf (#) (Non-used version of final file)
#. Panel5AWithoutTheNrd0WdAdaptedWithScalesImplemented.pdf (#) (Intermediate cleaned-up version of final file in Illustrator)

#. Panel5AWithTailoredScriptAndMoreDetailedVisualizationNowIncludesAllTooWithBinSize00117022022NoSegmentStacking18032022.pdf # (Basis for work to final version of article: see after part for panel 5C.)
#. (OverlayOfTheViolinPlotsWithBeanlineMediansAttemptToScaleInInverseOrderNrd0CorrectionDone11092020RemakeWithAllIncludedAndOverplottingWithReferenceUnity21032022.pdf # Not used version with Nrd0 and unity line.)


KoeberlinCorrelations <- read.csv(file = "C:/Users/Kevin/Documents/RData/Rest2/ACGData/Rest2/Rest2/Koeberlin_Snijder_Cell_2015_lipid_lipid_correlations.txt", header = TRUE, sep = "\t", as.is = TRUE)

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


library(reshape2)

KoeberlinCorrelationsConsensusNames2 <- KoeberlinCorrelationsConsensusNames
KoeberlinCorrelationsConsensusNames2[upper.tri(KoeberlinCorrelationsConsensusNames2)] <- NA

KoeberlinCorrelationsConsensusNames2Long <- melt(cbind(rownames(KoeberlinCorrelationsConsensusNames2), KoeberlinCorrelationsConsensusNames2))
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
pdf("C:/Users/Kevin/Documents/RData/Rest2/ACGData/FiguresByKT/Panel5AWithTailoredScriptAndMoreDetailedVisualizationNowIncludesAllTooWithBinSize00117022022NoSegmentStacking18032022.pdf")

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

# All of the figures in the second set of figure pages in the previous document were combined into Panel 5a and further esthetically enhanced and integrated with the other parts of Figure 5 in Adobe Illustrator.
# For the statistical significances the results of the Fisher exact tests on the difference in the medians was used (see further), and also added in Adobe Illustrator.


#### Panel 5B

#. Panel5BWithoutTheNrd0WdAdaptedWithScalesImplemented22032022.pdf (#)
#. Panel5BWithTheNrd0WdAdaptedWithScalesImplemented22032022.pdf (#)

#. Panel5BWithTailoredScriptAndMoreDetailedVisualizationNowIncludesAllTooWithBinSize001170220222NoSegmentStacking18032022.pdf # (Basis for used version (without the Nrd0))
#. Panel5BWithGreenUnityDistributionWithNrd0AndWd02322032022b # See after the 5C part. Not finally used version.

Coexp_Mander_NoTreshold_Rework21022019 <- read.csv(file = "C:/Users/Kevin/Documents/RData/Rest2/ACGData/FiguresByKT/Rest/Rest/Coexp_Mander_NoTreshold_Rework21022019.csv", header = TRUE, sep = ";", as.is = TRUE)
Lipid_Classes_Sergio_Subsets_210220194 <- read.csv(file = "C:/Users/Kevin/Documents/RData/Rest2/ACGData/FiguresByKT/Rest/Rest/Lipid_Classes_Sergio_Subsets_210220194.txt", header = TRUE, sep = "\t", as.is = TRUE)

# table(Lipid_Classes_Sergio_Subsets_210220194$SuperClass)
# 1032 lipid-like molecules; 179 other molecules; 1 (224): wrong annotation: -hydroxyprednisolone"

# -hydroxyprednisolone" Lipids and lipid-like molecules                           Other 
#                               1                            1032                             179 

# Eliminate non-lipids from METASPACE subsetting
SubsetOfLipidClassesSergio <- Lipid_Classes_Sergio_Subsets_210220194[Lipid_Classes_Sergio_Subsets_210220194$SuperClass != "Other",]

Coexp_Mander_NoTreshold_Rework210220192 <- Coexp_Mander_NoTreshold_Rework21022019[,3:47]
rownames(Coexp_Mander_NoTreshold_Rework210220192) <- paste(Coexp_Mander_NoTreshold_Rework21022019[,1], Coexp_Mander_NoTreshold_Rework21022019[,2], sep = "_")

MList <- lapply(1:dim(Coexp_Mander_NoTreshold_Rework210220192)[2],function(y){sapply(Coexp_Mander_NoTreshold_Rework210220192[,y],function(x){if(is.na(x)){NA}else{strsplit(GetStuffBetweenBrackets(x)[[1]], ",")[[1]][1]}})})
MMatrix <- do.call("cbind", MList)

mode(MMatrix) <- "numeric"
MMeans <- cbind(Coexp_Mander_NoTreshold_Rework21022019[,1:2], rowMeans(MMatrix, na.rm = TRUE))

MMeans2 <- cbind(MMeans,cbind(sapply(MMeans[,1], function(x){x %in% SubsetOfLipidClassesSergio[,1]}),sapply(MMeans[,2], function(x){x %in% SubsetOfLipidClassesSergio[,1]})))
LipidSubsetMMeans <- MMeans2[MMeans2[,4] & MMeans2[,5],]

LipidSubsetMMeans_2 <- do.call("cbind", list(LipidSubsetMMeans,
                                             "From" = gsub(pattern = "([A-Z](?![0-9]))", replacement = "\\11", x = LipidSubsetMMeans[,1], perl = TRUE),
                                             
                                             "To" = gsub(pattern = "([A-Z](?![0-9]))", replacement = "\\11", x = LipidSubsetMMeans[,2], perl = TRUE)))


LipidSubsetMMeans_4 <- data.frame(LipidSubsetMMeans_2[,1:5], 
                                  FromTo = paste(LipidSubsetMMeans_2[,"From"], LipidSubsetMMeans_2[,"To"], sep = "_"), 
                                  
                                  LipidSubsetMMeans_2[,6:7], 
                                  stringsAsFactors = FALSE)

HeadGroupConversionMatrix <- do.call("cbind", list("HeadGroup" = c("PC", "LPC", "PE", "LPE", "PG", "LPG", "PG/BMP", "BMP", "PS", "PI", "FA", "PA", "TAG", "DAG", "VE", "VA"),
                                                   "C" = c(8, 8, 5, 5, 6, 6, 6, 6, 6, 9, 0, 3, 3, 3, 26, 20),
                                                   
                                                   "H" = c(18, 18, 12, 12, 13, 13, 13, 13, 12, 16, 1, 7, 5, 6, 50, 30),
                                                   "N" = c(1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0),
                                                   
                                                   "O" = c(6, 6, 6, 6, 8, 8, 8, 8, 8, 11, 1, 6, 3, 3, 2, 1),
                                                   "P" = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0)))


LipidSubsetMMeans_4Subset <- data.frame(Cooccurrence = LipidSubsetMMeans_4[,3], LipidPairMatchingType = "All")

AggregatedLTPLipidPairsCombined_7wvi <- unique(rbind(cbind(PureAntonella32b[,c("LTPProtein", "Lipid", "LikelySubclass")], Screen = "A"), cbind(PureEnric32[,c("LTPProtein", "Lipid", "LikelySubclass")], Screen = "E")))
colnames(AggregatedLTPLipidPairsCombined_7wvi) <- c("LTPProtein", "LipidSpecies", "LipidSubclass", "Screen")

# Convertion: a: acyl; e: ether; d: d sphingosine; h: DH sphingosine; t: t sphingosine)
# Alternative nomenclature: similar to Koeberlin but not exactly the same: lysolipids indicated by headgroup, so PC_a_C and not LPC_a_C

ConversionMatrixForLipidsForFullDataAlternativeKoeberlinNomenclaturewvi <- cbind(as.character(levels(AggregatedLTPLipidPairsCombined_7wvi[,"LipidSubclass"])), c("BMP_aa_C", "CL_aaaa_C", "Cer_da_C", "CerP_da_C", "HexCer_da_C", "SHexCer_da_C", "SM_da_C", "DAG_aa_C", "Cer_da_C", "Cer_ha_C", "Cer_ta_C", "SM_ha_C",
                                                                                                                                                                 "FA_a_C", "FA_e_C", "PC_a_C", "PE_a_C", "PE_e_C", "PG_a_C", "PA_aa_C", "PC_aa_C", "PC_ae_C", "PE_aa_C", "PE_ae_C", "PG_aa_C", "PG/BMP_aa_C", "PGP_aa_C", 
                                                                                                                                                                 
                                                                                                                                                                 "PI_aa_C",  "PS_aa_C", "Hex2Cer_ta_C", "HexCer_ta_C", "SM_ta_C", "TAG_aaa_C", "Cer_ta_C", "VA_nn_C"))


AggregatedLTPLipidPairsCombined_8wvi <- cbind(AggregatedLTPLipidPairsCombined_7wvi, 
                                              LipidSpeciesAlternativeNomenclature = paste0(ConversionMatrixForLipidsForFullDataAlternativeKoeberlinNomenclaturewvi[match(AggregatedLTPLipidPairsCombined_7wvi$LipidSubclass, ConversionMatrixForLipidsForFullDataAlternativeKoeberlinNomenclaturewvi[,1]),2], 
                                                                                           
                                                                                           unlist(lapply(GetStuffBetweenBrackets(AggregatedLTPLipidPairsCombined_7wvi$LipidSpecies), function(x){if(length(x) != 0){mgsub::mgsub(x[length(x)], c("d", "DH", "t", "-OH", "\\*", "O-"), rep("",6))}else{"nn:nn"}}))))
write.table(AggregatedLTPLipidPairsCombined_8wvi, file="./RData/Rest2/ACGData/FiguresByKT/LTPLipidConnectionsWithSphingolipidsIncludedAndNewestDataKevinTiteca24092020.csv", sep="\t", row.names = FALSE, quote = FALSE)


# Corrected version of PI Headgroup #!

HeadGroupConversionMatrixc <- HeadGroupConversionMatrix
HeadGroupConversionMatrixc[HeadGroupConversionMatrixc[,"HeadGroup"] == "PI","H"] <- "17"


HeadGroupConversionMatrix2 <- rbind(cbind(HeadGroupConversionMatrixc, "S" = rep(0,16)),
                                    
                                    do.call("cbind", list("HeadGroup" = c("PGP", "CL", "Cer", "CerP", "SM", "HexCer", "Hex2Cer", "SHexCer"),
                                                          "C" = c(6,9,3,3,8,9,15,9), # Double check result here #!
                                                          
                                                          "H" = c(14,18,7,8,19,17,27,17), # Double check result here #!
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
                                "HAmount" = c(2,2,2,2,2,2,2,2,2,2,1,2,3,2,0,0,2,4,2,2,2,2,2,2))) # Similar to just counting the amount of linkages, so maybe replace #!

AggregatedLTPLipidPairsCombined_8wvi2 <- cbind(AggregatedLTPLipidPairsCombined_8wvi, do.call("rbind",strsplit(as.character(AggregatedLTPLipidPairsCombined_8wvi$LipidSpeciesAlternativeNomenclature), "\\_C|\\:")))
AggregatedLTPLipidPairsCombined_8wvi2 <- cbind(AggregatedLTPLipidPairsCombined_8wvi2, do.call("rbind", strsplit(as.character(AggregatedLTPLipidPairsCombined_8wvi2[,6]),"\\_")))

colnames(AggregatedLTPLipidPairsCombined_8wvi2)[6:10] <- c("SubclassKoeberlinlike", "TotalCarb", "TotalUnsat", "Headgroup", "Linkages")


library(stringr)
AggregatedLTPLipidPairsCombined_8wvi4 <- cbind(AggregatedLTPLipidPairsCombined_8wvi2, do.call("cbind", structure(lapply(ConnectorPieceConversionMatrix2[,1], function(x){str_count(AggregatedLTPLipidPairsCombined_8wvi2$Linkages, x)}), names = paste0(ConnectorPieceConversionMatrix2[,1],"Counts"))))

Carb <- as.numeric(as.character(AggregatedLTPLipidPairsCombined_8wvi4[, "TotalCarb"]))
MissingCarbByLinkage <- colSums(t(AggregatedLTPLipidPairsCombined_8wvi4[,paste0(ConnectorPieceConversionMatrix2[,1],"Counts")])*c(1,1,5,5,5,0))

# C residuals fatty acids
AggregatedLTPLipidPairsCombined_8wvi4c <- cbind(AggregatedLTPLipidPairsCombined_8wvi4, CResidualFattyAcids = (ifelse(!is.na(Carb), Carb, 0) - MissingCarbByLinkage))

# H residuals fatty acids
AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms <- cbind(AggregatedLTPLipidPairsCombined_8wvi4c, HResidualFattyAcids = 2*AggregatedLTPLipidPairsCombined_8wvi4c$CResidualFattyAcids + as.numeric(HTips2[match(as.character(AggregatedLTPLipidPairsCombined_8wvi4c$Headgroup), HTips2[,"HeadGroup"]), "HAmount"]) - 2*as.numeric(as.character(AggregatedLTPLipidPairsCombined_8wvi4c$TotalUnsat)))

AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms[is.na(AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms$HResidualFattyAcids), "HResidualFattyAcids"] <- 0


# Correction for the da sphingolipids, because otherwise 2H to much subtracted #! Leave it or correct it?
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

ObservedPossibleChemicalCombinationsForLipidsOfLTPs <- do.call("rbind", lapply(1:dim(LTPsScreens)[1], function(y){
  if(length(unique(as.character(DataSet[(DataSet$LTPProtein == LTPsScreens[y,]$LTPProtein) & (DataSet$Screen == LTPsScreens[y,]$Screen),"ChemicalFormulaWithOnes"]))) > 1){
    
    cbind(LTPsScreens[y,], t(combn(unique(as.character(DataSet[(DataSet$LTPProtein == LTPsScreens[y,]$LTPProtein) & (DataSet$Screen == LTPsScreens[y,]$Screen),"ChemicalFormulaWithOnes"])), 2)))}}))


AllPossibleChemicalCombinationsForLipidsOfLTPs <- t(combn(as.character(unique(DataSet$ChemicalFormulaWithOnes)),2)) 
#All possible pairs: 30381 <=> observed pairs: 10256

dim(ObservedPossibleChemicalCombinationsForLipidsOfLTPs[ObservedPossibleChemicalCombinationsForLipidsOfLTPs$Screen == "A",])[1] #2777
dim(ObservedPossibleChemicalCombinationsForLipidsOfLTPs[ObservedPossibleChemicalCombinationsForLipidsOfLTPs$Screen == "E",])[1] #7479

ObservedPossibleChemicalCombinationsForLipidsOfLTPs2 <- do.call("cbind", list(ObservedPossibleChemicalCombinationsForLipidsOfLTPs, 
                                                                              ChemicalPairs1 = apply(ObservedPossibleChemicalCombinationsForLipidsOfLTPs[,3:4], 1, function(x){paste0(x,collapse = "_")}),
                                                                              
                                                                              ChemicalPairs2 = apply(ObservedPossibleChemicalCombinationsForLipidsOfLTPs[,4:3], 1, function(x){paste0(x,collapse = "_")})))


AllPossibleChemicalCombinationsForLipidsOfLTPs2 <- do.call("cbind", list(AllPossibleChemicalCombinationsForLipidsOfLTPs, 
                                                                         ChemicalPairs1 = apply(AllPossibleChemicalCombinationsForLipidsOfLTPs[,1:2], 1, function(x){paste0(x,collapse = "_")}),
                                                                         
                                                                         ChemicalPairs2 = apply(AllPossibleChemicalCombinationsForLipidsOfLTPs[,2:1], 1, function(x){paste0(x,collapse = "_")})))
# Maybe separate the background calculations for the screens

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


CooccurrenceOfChemicalPairsList <- lapply(list(ObservedPossibleChemicalCombinationsForLipidsOfLTPs2, 
                                               as.data.frame(AllPossibleChemicalCombinationsForLipidsOfLTPs2), 
                                               
                                               as.data.frame(AllPossibleChemicalCombinationsForLipidsOfLTPs2a), 
                                               as.data.frame(AllPossibleChemicalCombinationsForLipidsOfLTPs2e)), 
                                          
                                          function(y){cbind(y, Cooccurrence = LipidSubsetMMeans_4[match(as.character(y[,"ChemicalPairs1"]), as.character(LipidSubsetMMeans_4[,"FromTo"])),3])})


# Write to do in Excel #! Entered earlier to avoid repetition and looping
write.table(t(t(levels(AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms4$LipidSubclass))), file="./RData/Rest2/ACGData/FiguresByKT/SubclassesForObservedLipids02102020.csv", sep="\t", row.names = FALSE, quote = FALSE)

# Input again with the consensus subclasses and classes, BMP seen as a type of PG, as well as PG/BMP #! Entered earlier to avoid repetition and looping
SubclassesMatchingToConsensusSubclassesAndClasses <- read.csv(file = "./RData/Rest2/ACGData/FiguresByKT/SubclassesForObservedLipidsConsensusSubclassesAndClasses02102020.txt", header = TRUE, sep = "\t", as.is = TRUE, quote = "")


CooccurrenceOfChemicalPairsList2 <- list()

y <- CooccurrenceOfChemicalPairsList[[1]]
CooccurrenceOfChemicalPairsList2[[1]] <- do.call("cbind", list(y,
                                                               
                                                               Subclass1 = AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms4[match(apply(y[,1:3], 1, function(x){paste0(x,collapse = "_")}),
                                                                                                                                    apply(AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms4[,c("LTPProtein", "Screen","ChemicalFormulaWithOnes")], 1, function(x){paste0(x,collapse = "_")})), "LipidSubclass"],
                                                               
                                                               Subclass2 = AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms4[match(apply(y[,c(1,2,4)], 1, function(x){paste0(x,collapse = "_")}),
                                                                                                                                    apply(AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms4[,c("LTPProtein", "Screen","ChemicalFormulaWithOnes")], 1, function(x){paste0(x,collapse = "_")})), "LipidSubclass"]
                                                               
))



y <- CooccurrenceOfChemicalPairsList[[2]]

CooccurrenceOfChemicalPairsList2[[2]] <- do.call("cbind", list(y,
                                                               Subclass1 = AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms4[match(as.character(y$V1), as.character(AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms4$ChemicalFormulaWithOnes)), "LipidSubclass"],
                                                               
                                                               Subclass2 = AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms4[match(as.character(y$V2), as.character(AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms4$ChemicalFormulaWithOnes)), "LipidSubclass"],
                                                               SubclassesMatchingToConsensusSubclassesAndClasses[match(AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms4[match(as.character(y$V1), as.character(AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms4$ChemicalFormulaWithOnes)), "LipidSubclass"], SubclassesMatchingToConsensusSubclassesAndClasses[,1]),2:3],
                                                               
                                                               SubclassesMatchingToConsensusSubclassesAndClasses[match(AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms4[match(as.character(y$V2), as.character(AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms4$ChemicalFormulaWithOnes)), "LipidSubclass"], SubclassesMatchingToConsensusSubclassesAndClasses[,1]),2:3]))



y <- CooccurrenceOfChemicalPairsList[[3]]

CooccurrenceOfChemicalPairsList2[[3]] <- do.call("cbind", list(y,
                                                               Subclass1 = AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms4[match(as.character(y$V1), as.character(AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms4$ChemicalFormulaWithOnes)), "LipidSubclass"],
                                                               
                                                               Subclass2 = AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms4[match(as.character(y$V2), as.character(AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms4$ChemicalFormulaWithOnes)), "LipidSubclass"],
                                                               SubclassesMatchingToConsensusSubclassesAndClasses[match(AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms4[match(as.character(y$V1), as.character(AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms4$ChemicalFormulaWithOnes)), "LipidSubclass"], SubclassesMatchingToConsensusSubclassesAndClasses[,1]),2:3],
                                                               
                                                               SubclassesMatchingToConsensusSubclassesAndClasses[match(AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms4[match(as.character(y$V2), as.character(AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms4$ChemicalFormulaWithOnes)), "LipidSubclass"], SubclassesMatchingToConsensusSubclassesAndClasses[,1]),2:3]))



y <- CooccurrenceOfChemicalPairsList[[4]]

CooccurrenceOfChemicalPairsList2[[4]] <- do.call("cbind", list(y,
                                                               Subclass1 = AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms4[match(as.character(y$V1), as.character(AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms4$ChemicalFormulaWithOnes)), "LipidSubclass"],
                                                               
                                                               Subclass2 = AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms4[match(as.character(y$V2), as.character(AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms4$ChemicalFormulaWithOnes)), "LipidSubclass"],
                                                               SubclassesMatchingToConsensusSubclassesAndClasses[match(AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms4[match(as.character(y$V1), as.character(AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms4$ChemicalFormulaWithOnes)), "LipidSubclass"], SubclassesMatchingToConsensusSubclassesAndClasses[,1]),2:3],
                                                               
                                                               SubclassesMatchingToConsensusSubclassesAndClasses[match(AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms4[match(as.character(y$V2), as.character(AggregatedLTPLipidPairsCombined_8wvi4WithAllAtoms4$ChemicalFormulaWithOnes)), "LipidSubclass"], SubclassesMatchingToConsensusSubclassesAndClasses[,1]),2:3]))


for(i in 2:4){
  CooccurrenceOfChemicalPairsList2[[i]] <- cbind(CooccurrenceOfChemicalPairsList2[[i]], LipidPairMatchingType = (CooccurrenceOfChemicalPairsList2[[i]][,8] == CooccurrenceOfChemicalPairsList2[[i]][,10]) + (CooccurrenceOfChemicalPairsList2[[i]][,9] == CooccurrenceOfChemicalPairsList2[[i]][,11]))
  
}


CooccurrenceOfChemicalPairsListExpandedWithAll <- CooccurrenceOfChemicalPairsList2[[1]]
CooccurrenceOfChemicalPairsListExpandedWithAll$LipidPairMatchingType <- "All"

CooccurrenceOfChemicalPairsListExpandedWithAll2 <- rbind(CooccurrenceOfChemicalPairsList2[[1]], CooccurrenceOfChemicalPairsListExpandedWithAll)



CooccurrenceOfChemicalPairsList2[[1]] <- do.call("cbind", list(CooccurrenceOfChemicalPairsList2[[1]],
                                                               
                                                               SubclassesMatchingToConsensusSubclassesAndClasses[match(CooccurrenceOfChemicalPairsList2[[1]]$Subclass1, SubclassesMatchingToConsensusSubclassesAndClasses[,1]),2:3],
                                                               SubclassesMatchingToConsensusSubclassesAndClasses[match(CooccurrenceOfChemicalPairsList2[[1]]$Subclass2, SubclassesMatchingToConsensusSubclassesAndClasses[,1]),2:3]))

CooccurrenceOfChemicalPairsList2[[1]] <- cbind(CooccurrenceOfChemicalPairsList2[[1]], LipidPairMatchingType = (CooccurrenceOfChemicalPairsList2[[1]][,10] == CooccurrenceOfChemicalPairsList2[[1]][,12]) + (CooccurrenceOfChemicalPairsList2[[1]][,11] == CooccurrenceOfChemicalPairsList2[[1]][,13]))



PossibleCooccurrenceOfChemicalPairs <- rbind(cbind(CooccurrenceOfChemicalPairsList2[[3]], Screen = "A"), cbind(CooccurrenceOfChemicalPairsList2[[4]], Screen = "E"))

PossibleCooccurrenceOfChemicalPairs2 <- PossibleCooccurrenceOfChemicalPairs
PossibleCooccurrenceOfChemicalPairs2$LipidPairMatchingType <- "All"

PossibleCooccurrenceOfChemicalPairs2 <- rbind(PossibleCooccurrenceOfChemicalPairs, PossibleCooccurrenceOfChemicalPairs2)



pdf("C:/Users/Kevin/Documents/RData/Rest2/ACGData/FiguresByKT/Panel5BWithTailoredScriptAndMoreDetailedVisualizationNowIncludesAllTooWithBinSize001170220222NoSegmentStacking18032022.pdf")

BinSize <- 0.01
MaxList2TissuesColoc <- list()

StatListTissuesColoc <- list()
StatList2TissuesColoc <- list()

StatListTissuesColocX2 <- list()
StatList2TissuesColocX2 <- list()

ScreenColorsForFilling <- cbind(c("A","E"),c(brewer.pal(9,"Blues")[4], brewer.pal(9,"Oranges")[4]))
ScreenColorsForFillingPot <- cbind(c("A","E"),c(brewer.pal(9,"Blues")[7], brewer.pal(9,"Oranges")[7]))


for(SimilarityIndex in c("All", as.character(0:2))){
  
  x1a <- na.omit(as.numeric(CooccurrenceOfChemicalPairsListExpandedWithAll2[(CooccurrenceOfChemicalPairsListExpandedWithAll2[,"Screen"] == "A") & (CooccurrenceOfChemicalPairsListExpandedWithAll2[,"LipidPairMatchingType"] == SimilarityIndex),"Cooccurrence"]))
  x1e <- na.omit(as.numeric(CooccurrenceOfChemicalPairsListExpandedWithAll2[(CooccurrenceOfChemicalPairsListExpandedWithAll2[,"Screen"] == "E") & (CooccurrenceOfChemicalPairsListExpandedWithAll2[,"LipidPairMatchingType"] == SimilarityIndex),"Cooccurrence"]))
  
  x1apot <- na.omit(as.numeric(PossibleCooccurrenceOfChemicalPairs2[(PossibleCooccurrenceOfChemicalPairs2[,"Screen"] == "A") & (PossibleCooccurrenceOfChemicalPairs2[,"LipidPairMatchingType"] == SimilarityIndex),"Cooccurrence"]))
  x1epot <- na.omit(as.numeric(PossibleCooccurrenceOfChemicalPairs2[(PossibleCooccurrenceOfChemicalPairs2[,"Screen"] == "E") & (PossibleCooccurrenceOfChemicalPairs2[,"LipidPairMatchingType"] == SimilarityIndex),"Cooccurrence"]))
  
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
  SimilarityIndex <- "2"
  
  x1a <- na.omit(as.numeric(CooccurrenceOfChemicalPairsListExpandedWithAll2[(CooccurrenceOfChemicalPairsListExpandedWithAll2[,"Screen"] == "A") & (CooccurrenceOfChemicalPairsListExpandedWithAll2[,"LipidPairMatchingType"] == SimilarityIndex),"Cooccurrence"]))
  x1e <- na.omit(as.numeric(CooccurrenceOfChemicalPairsListExpandedWithAll2[(CooccurrenceOfChemicalPairsListExpandedWithAll2[,"Screen"] == "E") & (CooccurrenceOfChemicalPairsListExpandedWithAll2[,"LipidPairMatchingType"] == SimilarityIndex),"Cooccurrence"]))
  
  x1apot <- na.omit(as.numeric(PossibleCooccurrenceOfChemicalPairs2[(PossibleCooccurrenceOfChemicalPairs2[,"Screen"] == "A") & (PossibleCooccurrenceOfChemicalPairs2[,"LipidPairMatchingType"] == SimilarityIndex),"Cooccurrence"]))
  x1epot <- na.omit(as.numeric(PossibleCooccurrenceOfChemicalPairs2[(PossibleCooccurrenceOfChemicalPairs2[,"Screen"] == "E") & (PossibleCooccurrenceOfChemicalPairs2[,"LipidPairMatchingType"] == SimilarityIndex),"Cooccurrence"]))
  
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


# All of the figures in the second set of figure pages in the previous document were combined into Panel 5b and further esthetically enhanced and integrated with the other parts of Figure 5 in Adobe Illustrator.
# For the statistical significances the results of the Fisher exact tests on the difference in the medians was used (see further), and also added in Adobe Illustrator.


#### Panel 5C

#. Panel5CWithoutTheNrd0WdAdaptedWithScalesImplemented22032022.pdf (#) (Name of the finally used version for inclusion in figure 5.)
#. Panel5CWithTheNrd0WdAdaptedWithScalesImplemented22032022.pdf (#) (Not the basis of the final version.)

#. Panel5CWithTailoredScriptAndMoreDetailedVisualizationBinSize001NoSegmentStacking18032022.pdf # (Finally used version for inclusion in figure 5 based on this.)
#. Panel5CWithGreenUnityDistributionAndBwWithNrd0Wd01422032022.pdf # (Not the basis of the final version.)

MolPctOfMads <- Col1ToRowNames(read.csv(file = "C:/Users/Kevin/Documents/RData/Rest2/ACGData/FiguresByKT/processed.molpct_generatedbyMads.csv", header = TRUE, sep = ",", as.is = TRUE))[,1:38] # Get rid of empty space at end columns
ExperimentTypesOfMads <- sapply(strsplit(colnames(MolPctOfMads), "\\."), "[",1)


AggregateMPOfMads <- sapply(unique(ExperimentTypesOfMads), function(x){rowMeans(MolPctOfMads[,ExperimentTypesOfMads == x], na.rm = TRUE)})

AggregateMPOfMads0 <- AggregateMPOfMads
AggregateMPOfMads0[is.na(AggregateMPOfMads0)] <- 0


AggregateMPOfMads0RTL <- setNames(split(AggregateMPOfMads0, row(AggregateMPOfMads0)), rownames(AggregateMPOfMads0))

AggregateMPOfMads0RTLMOC <- do.call("cbind", setNames(lapply(AggregateMPOfMads0RTL, function(y){unlist(lapply(AggregateMPOfMads0RTL, function(x){MOCFunction(y, x)}))}), rownames(AggregateMPOfMads0)))
write.table(AggregateMPOfMads0RTLMOC, file="./RData/Rest2/ACGData/FiguresByKT/AggregateMPOfMads0RTLMOC.tsv", sep="\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

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

gc()


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
    cbind(unique(InCelluloReducedLinks[InCelluloReducedLinks[,1] == x, 1:2]), t(combn(as.character(InCelluloReducedLinks[InCelluloReducedLinks[,1] == x, 3]), m = 2)))
    
  }}))



InVitroReducedLinksCouples <- do.call("rbind", lapply(unique(InVitroReducedLinks[,1]), function(x){
  
  if(sum(InVitroReducedLinks[,1] == x) > 1){
    cbind(unique(InVitroReducedLinks[InVitroReducedLinks[,1] == x, 1:2]), t(combn(as.character(InVitroReducedLinks[InVitroReducedLinks[,1] == x, 3]), m = 2)))
    
  }}))


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


pdf("C:/Users/Kevin/Documents/RData/Rest2/ACGData/FiguresByKT/Panel5CWithTailoredScriptAndMoreDetailedVisualizationBinSize001NoSegmentStacking18032022.pdf")

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

# All of the figures in the second set of figure pages in the previous document were combined into Panel 5b and further esthetically enhanced and integrated with the other parts of Figure 5 in Adobe Illustrator.
# For the statistical significances the results of the Fisher exact tests on the difference in the medians was used (see further), and also added in Adobe Illustrator.

#### Similar figures for panels from figure 5, but with Nrd0 and with unity lines to enable internal comparisons, y-axes, and correct fitting with eachother during integration of figures in Adobe Illustrator.
#### Not the finally used version of these panels (see parts before, for the finally used versions).


#### Panel 5A Alternative With Nrd0

# After this "all" gets integrated in 5a in a similar set-up as for 5b and 5c too
CorrelationsMatrixae8_8g <- CorrelationsMatrixae8_8

CorrelationsMatrixae8_8g[,"MatchLevel"] <- "All"
CorrelationsMatrixae8_8combo <- rbind(CorrelationsMatrixae8_8g, CorrelationsMatrixae8_8)

KoeberlinCorrelationsConsensusNames2LongReducedVersion7g <- KoeberlinCorrelationsConsensusNames2LongReducedVersion7
KoeberlinCorrelationsConsensusNames2LongReducedVersion7g[,"MatchingNumber"] <- "All"

KoeberlinCorrelationsConsensusNames2LongReducedVersion7combo <- rbind(KoeberlinCorrelationsConsensusNames2LongReducedVersion7g, KoeberlinCorrelationsConsensusNames2LongReducedVersion7) 


pdf("C:/Users/Kevin/Documents/RData/Rest2/ACGData/FiguresByKT/OverlayOfTheViolinPlotsWithBeanlineMediansAttemptToScaleInInverseOrderNrd0CorrectionDone11092020RemakeWithAllIncludedAndOverplotting18032022.pdf")
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
# Previous figure is Nrd0-version for panel 5A that was not finally used in the article, but without an unit line for the beanplots as reference for further steps

# Steps to include an unit line for all the beanplots
UnitBean <- do.call("rbind.data.frame", lapply(c(as.character(2:0),"All"), function(x){cbind.data.frame(as.numeric(seq(-1,1, by = 0.01)), x)}))

colnames(UnitBean) <- c("correlation", "MatchingNumber")


pdf("C:/Users/Kevin/Documents/RData/Rest2/ACGData/FiguresByKT/OverlayOfTheViolinPlotsWithBeanlineMediansAttemptToScaleInInverseOrderNrd0CorrectionDone11092020RemakeWithAllIncludedAndOverplottingWithReferenceUnity21032022.pdf")
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
         col = c(NA,NA, NA,NA),
         
         axes=F, log = "",
         horizontal = TRUE,  boxwex = 1,
         
         cutmin = -1,
         cutmax = 1,
         
         border = "green", wd = TestBeanForSpecial218032022$wd,
         overalline = "median", beanlines = "median",
         
         add = TRUE)
dev.off()

# Previous figure is the Nrd0 corrected version for panel 5A that was not finally used in the article, with an unit line for the beanplots as reference for further steps
rm(TestBeanForGeneral2standard18032022)


#### Panel 5C Alternative With Nrd0 (and unit line)

UnitBeanSubcellLoc <- do.call("rbind.data.frame", lapply(c("Species", "Sub-Class", "Class", "All"), function(x){cbind.data.frame(as.numeric(seq(0,1, by = 0.01)), x)}))
colnames(UnitBeanSubcellLoc) <- c("SubcellularColocalization", "PairType")

pdf("C:/Users/Kevin/Documents/RData/Rest2/ACGData/FiguresByKT/Panel5CWithGreenUnityDistributionAndBwWithNrd0Wd01422032022.pdf")
beanplot(as.numeric(as.character(SubcellularLocalizationOverlapOurDatax[,"SubcellularColocalization"])) ~ factor(SubcellularLocalizationOverlapOurDatax[,"Screen"], levels = c("In Vitro", "In Cellulo")) * factor(SubcellularLocalizationOverlapOurDatax[,"PairType"], levels = c("Species", "Sub-Class", "Class", "All")),
         
         main = "Co-transport: lipid hierarchy preference by co-occurrence in cells", side = "both", xlab="Sub-cellular co-occurrence", ll = 0.04, bw = "nrd0", wd = 0.14, # wd = 0.23,
         col = list(brewer.pal(9,"Oranges")[4], c(brewer.pal(9,"Blues")[4], "black")), method = "overplot", 
         
         axes=F,
         horizontal = TRUE,
         
         cutmin = 0,
         cutmax = 1,
         
         border = FALSE,
         overalline = "median")

AllPairsInSubcellularData <- rbind(AggregateMPOfMads0RTLMOC11x[,c("PairType", "SubcellularColocalization", "Pairs")], cbind(PairType = "All", AggregateMPOfMads0RTLMOC11x[,c("SubcellularColocalization", "Pairs")]))
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
         col = c(NA,NA, NA,NA),
         
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
# Previous figure is the Nrd0 corrected version for panel 5C that was not finally used in the article, with an unit line for the beanplots as reference for further steps


#### Panel 5B Alternative With Nrd0 (and unit line)

UnitBeanTissueLoc <- do.call("rbind.data.frame", lapply(c("Species", "Sub-Class", "Class", "All"), function(x){cbind.data.frame(as.numeric(seq(0,1, by = 0.01)), x)}))
colnames(UnitBeanTissueLoc) <- c("Cooccurrence", "LipidPairMatchingType")

pdf("C:/Users/Kevin/Documents/RData/Rest2/ACGData/FiguresByKT/Panel5BWithGreenUnityDistributionWithNrd0AndWd02322032022b.pdf")
beanplot(as.numeric(as.character(CooccurrenceOfChemicalPairsListExpandedWithAll2[,"Cooccurrence"])) ~ factor(CooccurrenceOfChemicalPairsListExpandedWithAll2[,"Screen"], levels = c("E", "A")) * factor(CooccurrenceOfChemicalPairsListExpandedWithAll2[,"LipidPairMatchingType"], levels = c(as.character(2:0), "All")), 
         
         main = "Co-transport: lipid hierarchy preference by co-occurrence in tissues", side = "both", xlab="Co-occurrence", ll = 0.04, wd = 0.23,
         col = list(brewer.pal(9,"Oranges")[4], c(brewer.pal(9,"Blues")[4], "black")), method = "overplot",
         
         axes=F, bw = "nrd0", 
         horizontal = TRUE,
         
         cutmin = 0,
         cutmax = 1,
         
         border = FALSE,
         overalline = "median")


beanplot(as.numeric(as.character(PossibleCooccurrenceOfChemicalPairs2[,"Cooccurrence"])) ~ factor(PossibleCooccurrenceOfChemicalPairs2[,"Screen"], levels = c("E", "A")) * factor(PossibleCooccurrenceOfChemicalPairs2[,"LipidPairMatchingType"], levels = c(as.character(2:0), "All")), 
         
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


beanplot(as.numeric(UnitBeanSubcellLoc[,"SubcellularColocalization"]) ~ factor(UnitBeanSubcellLoc[,"PairType"], levels = c("Species", "Sub-Class", "Class", "All")), ll = 0, #!
         
         # bw = "nrd0", #alpha = 0.5, #!
         col = c(NA,NA, NA,NA),
         
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
# Previous figure is the Nrd0 corrected version for panel 5A that was not finally used in the article, with an unit line for the beanplots as reference for further steps


#### Results of Fisher exact tests

#. FisherExactTestResultsStandardCorrectedWithAll18082022.csv # Not finally used version because not all-round applicable & least strict
#. FisherExactTestResultsWithEqualizedBackgroundCorrectedWithAll18082022.csv # Finally used version of comparisons: applicable to all & even more strict

write.table(structure(do.call("rbind", StatList), dimnames = list(paste0(rep(c("A","E"),3), unlist(lapply(c("All",as.character(0:2)), function(x){rep(x, 2)}))), c("NoNeg", "NoPos", "YesNeg", "YesPos", "LowConf", "HighConf", "OddsRatio", "pValue"))),
            file="./RData/Rest2/ACGData/FiguresByKT/FisherExactTestResultsStandardCorrectedWithAll18082022.csv", sep="\t", row.names = TRUE, quote = FALSE)

write.table(structure(do.call("rbind", StatList2), dimnames = list(paste0(rep(c("A","E"),3), unlist(lapply(c("All",as.character(0:2)), function(x){rep(x, 2)}))), c("NoNeg", "NoPos", "YesNeg", "YesPos", "LowConf", "HighConf", "OddsRatio", "pValue"))),
            file="./RData/Rest2/ACGData/FiguresByKT/FisherExactTestResultsWithEqualizedBackgroundCorrectedWithAll18082022.csv", sep="\t", row.names = TRUE, quote = FALSE)

#. FisherExactTestResultsStandardForTheCooccurrences09102020.csv # Not finally used version because not all-round applicable & least strict
#. FisherExactTestResultsWithEqualizedBackgroundForTheCooccurrencesCorrected18082022.csv # Finally used version of comparisons: applicable to all & even more strict


SimilarityIndex <- "All"

x1a <- RemoveNAsFromVector(as.numeric(CooccurrenceOfChemicalPairsListExpandedWithAll2[(CooccurrenceOfChemicalPairsListExpandedWithAll2[,"Screen"] == "A") & (CooccurrenceOfChemicalPairsListExpandedWithAll2[,"LipidPairMatchingType"] == SimilarityIndex),"Cooccurrence"]))
x1e <- RemoveNAsFromVector(as.numeric(CooccurrenceOfChemicalPairsListExpandedWithAll2[(CooccurrenceOfChemicalPairsListExpandedWithAll2[,"Screen"] == "E") & (CooccurrenceOfChemicalPairsListExpandedWithAll2[,"LipidPairMatchingType"] == SimilarityIndex),"Cooccurrence"]))

pdForx1a <- RemoveNAsFromVector(as.numeric(as.character(PossibleCooccurrenceOfChemicalPairs2[(PossibleCooccurrenceOfChemicalPairs2[,"Screen"] == "A") & (PossibleCooccurrenceOfChemicalPairs2[,"LipidPairMatchingType"] == SimilarityIndex),"Cooccurrence"])))
pdForx1e <- RemoveNAsFromVector(as.numeric(as.character(PossibleCooccurrenceOfChemicalPairs2[(PossibleCooccurrenceOfChemicalPairs2[,"Screen"] == "E") & (PossibleCooccurrenceOfChemicalPairs2[,"LipidPairMatchingType"] == SimilarityIndex),"Cooccurrence"])))

StatListForCooccurrences[[SimilarityIndex]] <- structure(rbind(FishersExactTestBelowAbove05(SubsetX = x1a, SupersetY = c(pdForx1a,pdForx1a)), FishersExactTestBelowAbove05(SubsetX = x1e, SupersetY = c(pdForx1e,pdForx1e))), 
                                                         dimnames = list(c("A","E"), c("NoNeg", "NoPos", "YesNeg", "YesPos", "LowConf", "HighConf", "OddsRatio", "pValue")))

StatList2ForCooccurrences[[SimilarityIndex]] <- structure(rbind(FishersExactTestNegPosEqualBackground2(SubsetX = x1a, SupersetY = c(pdForx1a,pdForx1a)), FishersExactTestNegPosEqualBackground2(SubsetX = x1e, SupersetY = c(pdForx1e,pdForx1e))), 
                                                          dimnames = list(c("A","E"), c("NoNeg", "NoPos", "YesNeg", "YesPos", "LowConf", "HighConf", "OddsRatio", "pValue")))

Stat1FullAllForCooccurrences <- structure(rbind(FishersExactTestBelowAbove05(SubsetX = x1a, SupersetY = c(as.numeric(as.character(LipidSubsetMMeans_4Subset[,1])),as.numeric(as.character(LipidSubsetMMeans_4Subset[,1])))), FishersExactTestBelowAbove05(SubsetX = x1e, SupersetY = c(as.numeric(as.character(LipidSubsetMMeans_4Subset[,1])),as.numeric(as.character(LipidSubsetMMeans_4Subset[,1]))))), 
                                          dimnames = list(c("A","E"), c("NoNeg", "NoPos", "YesNeg", "YesPos", "LowConf", "HighConf", "OddsRatio", "pValue")))

Stat2FullAllForCooccurrencesCorrectedVersion <- structure(rbind(FishersExactTestNegPosEqualBackground2(SubsetX = x1a, SupersetY = c(as.numeric(as.character(LipidSubsetMMeans_4Subset[,1])),as.numeric(as.character(LipidSubsetMMeans_4Subset[,1])))), FishersExactTestNegPosEqualBackground2(SubsetX = x1e, SupersetY = c(as.numeric(as.character(LipidSubsetMMeans_4Subset[,1])),as.numeric(as.character(LipidSubsetMMeans_4Subset[,1]))))), 
                                                          dimnames = list(c("A","E"), c("NoNeg", "NoPos", "YesNeg", "YesPos", "LowConf", "HighConf", "OddsRatio", "pValue")))


for(SimilarityIndex in as.character(0:2)){
  
  x1a <- RemoveNAsFromVector(as.numeric(CooccurrenceOfChemicalPairsListExpandedWithAll2[(CooccurrenceOfChemicalPairsListExpandedWithAll2[,"Screen"] == "A") & (CooccurrenceOfChemicalPairsListExpandedWithAll2[,"LipidPairMatchingType"] == SimilarityIndex),"Cooccurrence"]))
  x1e <- RemoveNAsFromVector(as.numeric(CooccurrenceOfChemicalPairsListExpandedWithAll2[(CooccurrenceOfChemicalPairsListExpandedWithAll2[,"Screen"] == "E") & (CooccurrenceOfChemicalPairsListExpandedWithAll2[,"LipidPairMatchingType"] == SimilarityIndex),"Cooccurrence"]))
  
  pdForx1a <- RemoveNAsFromVector(as.numeric(as.character(PossibleCooccurrenceOfChemicalPairs2[(PossibleCooccurrenceOfChemicalPairs2[,"Screen"] == "A") & (PossibleCooccurrenceOfChemicalPairs2[,"LipidPairMatchingType"] == SimilarityIndex),"Cooccurrence"])))
  pdForx1e <- RemoveNAsFromVector(as.numeric(as.character(PossibleCooccurrenceOfChemicalPairs2[(PossibleCooccurrenceOfChemicalPairs2[,"Screen"] == "E") & (PossibleCooccurrenceOfChemicalPairs2[,"LipidPairMatchingType"] == SimilarityIndex),"Cooccurrence"])))
  
  StatListForCooccurrences[[SimilarityIndex]] <- structure(rbind(FishersExactTestBelowAbove05(SubsetX = x1a, SupersetY = c(pdForx1a,pdForx1a)), FishersExactTestBelowAbove05(SubsetX = x1e, SupersetY = c(pdForx1e,pdForx1e))), 
                                                           dimnames = list(c("A","E"), c("NoNeg", "NoPos", "YesNeg", "YesPos", "LowConf", "HighConf", "OddsRatio", "pValue")))
  
  StatList2ForCooccurrences[[SimilarityIndex]] <- structure(rbind(FishersExactTestNegPosEqualBackground2(SubsetX = x1a, SupersetY = c(pdForx1a,pdForx1a)), FishersExactTestNegPosEqualBackground2(SubsetX = x1e, SupersetY = c(pdForx1e,pdForx1e))), 
                                                            dimnames = list(c("A","E"), c("NoNeg", "NoPos", "YesNeg", "YesPos", "LowConf", "HighConf", "OddsRatio", "pValue")))
  
}


StatListForCooccurrences[["AllAll"]] <- Stat1FullAllForCooccurrences


StatList2ForCooccurrencesCorrected <- StatList2ForCooccurrences
StatList2ForCooccurrencesCorrected[["AllAll"]] <- Stat2FullAllForCooccurrencesCorrectedVersion

write.table(structure(do.call("rbind", StatListForCooccurrences[c("AllAll", "All", as.character(0:2))]), dimnames = list(paste0(rep(c("A","E"),5), unlist(lapply(c("AllAll", "All", as.character(0:2)), function(x){rep(x, 2)}))), c("NoNeg", "NoPos", "YesNeg", "YesPos", "LowConf", "HighConf", "OddsRatio", "pValue"))),
            file="./RData/Rest2/ACGData/FiguresByKT/FisherExactTestResultsStandardForTheCooccurrences09102020.csv", sep="\t", row.names = TRUE, quote = FALSE)

write.table(structure(do.call("rbind", StatList2ForCooccurrencesCorrected[c("AllAll", "All", as.character(0:2))]), dimnames = list(paste0(rep(c("A","E"),5), unlist(lapply(c("AllAll", "All", as.character(0:2)), function(x){rep(x, 2)}))), c("NoNeg", "NoPos", "YesNeg", "YesPos", "LowConf", "HighConf", "OddsRatio", "pValue"))),
            file="./RData/Rest2/ACGData/FiguresByKT/FisherExactTestResultsWithEqualizedBackgroundForTheCooccurrencesCorrected18082022.csv", sep="\t", row.names = TRUE, quote = FALSE)

# Subcellular Fisher Exact Tests # Finally used version of comparisons: only the median versions, for full reference (grey lines in graphs) and for what possible (colored lines in graphs) and not the others because applicable to all & even more strict, and double references used for comparible size in comparisons.
# They are all named ...21022022.csv . #! Double-check the correct versions!

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
              file= paste0("./RData/Rest2/ACGData/FiguresByKT/",names(FisherTestsForSubcellularColocalizations)[i],"21022022.csv"), quote = FALSE)
  
}



######## Figure 6

#### Fig.6b
#. BarplotHeatmapTestSphingolipidTransporters06092021WithSidesAdded2CorrecteddSMZeroToNA.pdf # (R-basis of final figure.)

#. BarplotHeatmapTestSphingolipidTransporters06092021WithSidesAdded2CorrecteddSMZeroToNA10CorrectMetabolism5MoreAdvancedVersions4WithHighlighLinesWithStrongerCrosses.pdf (#) 
# (Last page is further intermediate figure after integration in Adobe Illustrator: Final figure has LCN1 information removed, legends simplified, CERS information added to the top, and highlights of most important part focusses.)


LTPLipidConnectionsDataSubset <- cbind(LTPLipidConnectionsDataSet[,c("LTPProtein", "LikelySubclass", "Screen", "Intensity")], CarbonChain = paste(LTPLipidConnectionsDataSet[,"TotalCarbonChainLength"],LTPLipidConnectionsDataSet[,"TotalCarbonChainUnsaturations"], sep = ":"))

library(reshape2)
LTPLipidConnectionsDataAggregated <- dcast(LTPLipidConnectionsDataSubset, LTPProtein + LikelySubclass + Screen ~ CarbonChain, value.var="Intensity", fun.aggregate = function(x){sum(x, na.rm = TRUE)})

SphingolipidSubclasses <- c("d*Cer", "dCer", "DHCer", "DHOH*Cer", "tCer", "d*CerP", "d*SM", "DHSM", "t*SM", "d*HexCer", "t*HexCer", "d*SHexCer", "t*Hex2Cer")
LTPLipidConnectionsDataSphingolipids <- LTPLipidConnectionsDataAggregated[LTPLipidConnectionsDataAggregated$LikelySubclass %in% SphingolipidSubclasses,]

LTPLipidConnectionsDataSphingolipids2 <- LTPLipidConnectionsDataSphingolipids[,c(rep(TRUE,3), colSums(LTPLipidConnectionsDataSphingolipids[,4:101]) != 0)]



CellularSphingolipidDistributionInfo <- Col1ToRowNames(dcast(do.call("rbind", lapply(SphingolipidSubclasses[SphingolipidSubclasses %in% names(LipidSubclassesAddedToBackground170620214)], function(x){
  
  data.frame(CarbonChain = mgsub::mgsub(sapply(strsplit(rownames(LipidSubclassesAddedToBackground170620214[[x]]), "\\("), "[[", 2), c(")", "d", "\\*"), rep("",3)),
             Cellular = LipidSubclassesAddedToBackground170620214[[x]][,"Cellular"],
             
             Lipid = x)})),
  Lipid ~ CarbonChain, value.var = "Cellular"))

CellularSphingolipidDistributionInfo2 <- CellularSphingolipidDistributionInfo[rowSums(CellularSphingolipidDistributionInfo, na.rm = TRUE) > 0,
                                                                              colSums(CellularSphingolipidDistributionInfo, na.rm = TRUE) > 0]

CellularSphingolipidDistributionInfo4 <- do.call("cbind", list(LTPProtein = "Cellular", 
                                                               LikelySubclass = rownames(CellularSphingolipidDistributionInfo2),
                                                               
                                                               Screen = "Cellular",
                                                               CellularSphingolipidDistributionInfo2))

LTPLipidConnectionsDataSphingolipidsWithCellularAdded <- as.data.frame(t(Col1ToRowNames(merge(t(LTPLipidConnectionsDataSphingolipids2), t(CellularSphingolipidDistributionInfo4), by = 0, all = TRUE))))
LTPLipidConnectionsDataSphingolipidsWithCellularMatrix <- as.matrix(LTPLipidConnectionsDataSphingolipidsWithCellularAdded[,1:22])

mode(LTPLipidConnectionsDataSphingolipidsWithCellularMatrix) <- "numeric"
LTPLipidConnectionsDataSphingolipidsWithCellularMatrix[is.na(LTPLipidConnectionsDataSphingolipidsWithCellularMatrix)] <- 0

LTPLipidConnectionsDataSphingolipidsWithCellularMatrixToRowMax <- do.call("rbind", lapply(1:dim(LTPLipidConnectionsDataSphingolipidsWithCellularMatrix)[1], 
                                                                                          function(x){LTPLipidConnectionsDataSphingolipidsWithCellularMatrix[x,]*100/max(LTPLipidConnectionsDataSphingolipidsWithCellularMatrix[x,], na.rm = TRUE)}))

rownames(LTPLipidConnectionsDataSphingolipidsWithCellularMatrixToRowMax) <- apply(LTPLipidConnectionsDataSphingolipidsWithCellularAdded[,23:25], 1 , paste , collapse = "_" )
LTPLipidConnectionsDataSphingolipidsWithCellularMatrixToRowMax2 <- LTPLipidConnectionsDataSphingolipidsWithCellularMatrixToRowMax[c(19,14:18,20,7,21,9,8,11,10,12,13,22,2,1,6,5,3,4),]

OriginalRowsSphingoTrans <- do.call("rbind", strsplit(rownames(LTPLipidConnectionsDataSphingolipidsWithCellularMatrixToRowMax2), split = "_"))
NewRowsSphingoTrans <- unique(OriginalRowsSphingoTrans[OriginalRowsSphingoTrans[,3] != "Cellular", 1:2])


EmptyLayer <- matrix(NA, 
                     
                     nrow = dim(NewRowsSphingoTrans)[1], 
                     ncol = dim(LTPLipidConnectionsDataSphingolipidsWithCellularMatrixToRowMax2)[2], 
                     
                     dimnames = list(paste(NewRowsSphingoTrans[,1],NewRowsSphingoTrans[,2],sep = "_"),
                                     colnames(LTPLipidConnectionsDataSphingolipidsWithCellularMatrixToRowMax2)))

NewSphingoTransLayers <- list(InCellulo = EmptyLayer,
                              InVitro = EmptyLayer,
                              
                              Cellular = EmptyLayer)


#! Maybe make more compact in for loop over a dataframe?

for(j in which(OriginalRowsSphingoTrans[,3] == "in vivo")){
  for(i in which(paste0(NewRowsSphingoTrans[,1],NewRowsSphingoTrans[,2]) == paste0(OriginalRowsSphingoTrans[j,1],OriginalRowsSphingoTrans[j,2]))){
    
    NewSphingoTransLayers$InCellulo[i,] <- ZerosToNAsConverter(LTPLipidConnectionsDataSphingolipidsWithCellularMatrixToRowMax2[j,])
  }}

for(j in which(OriginalRowsSphingoTrans[,3] == "in vitro")){
  for(i in which(paste0(NewRowsSphingoTrans[,1],NewRowsSphingoTrans[,2]) == paste0(OriginalRowsSphingoTrans[j,1],OriginalRowsSphingoTrans[j,2]))){
    
    NewSphingoTransLayers$InVitro[i,] <- ZerosToNAsConverter(LTPLipidConnectionsDataSphingolipidsWithCellularMatrixToRowMax2[j,])
  }}

for(j in which(OriginalRowsSphingoTrans[,3] == "Cellular")){
  for(i in which(NewRowsSphingoTrans[,1] == OriginalRowsSphingoTrans[j,1])){
    
    NewSphingoTransLayers$Cellular[i,] <- ZerosToNAsConverter(LTPLipidConnectionsDataSphingolipidsWithCellularMatrixToRowMax2[j,])
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

pdf("C:/Users/Kevin/Documents/RData/Rest2/ACGData/FiguresByKT/BarplotHeatmapTestSphingolipidTransporters06092021WithSidesAdded2CorrecteddSMZeroToNA.pdf",
    width = unit(20, "mm"), height = unit(25, "mm"))

Heatmap(NewSphingoTransLayers$InCellulo, name = LegendName, col = LegendColor, rect_gp = gpar(type = "none"), column_title = "", 
        width = unit(WidthAdaptor*dim(NewSphingoTransLayers$InCellulo)[2]/min(dim(NewSphingoTransLayers$InCellulo)), "mm"), height = unit(HeightAdaptor*dim(NewSphingoTransLayers$InCellulo)[1]/min(dim(NewSphingoTransLayers$InCellulo)), "mm"),
        
        
        cell_fun = function(j, i, x, y, width, height, fill) {

          
          grid.rect(x = x, y = y, width = width, height = NewSphingoTransLayers$InCellulo[i,j]/100*height, gp = gpar(col = NA, fill = brewer.pal(9,"Blues")[5], lwd = 1), just = c("center","bottom"))
          grid.rect(x = x, y = y, width = 0.5*width, height = NewSphingoTransLayers$InVitro[i,j]/100*height, gp = gpar(col = NA, fill = brewer.pal(9,"Oranges")[5], lwd = 1), just = c("center","bottom"))
          
          grid.rect(x = x, y = y, width = width, height = NewSphingoTransLayers$Cellular[i,j]/100*height, gp = gpar(col = "LightGrey", fill = NA, lwd = 2), just = c("center","bottom"))
        }, cluster_rows = FALSE, cluster_columns = FALSE, 
        

        cluster_row_slices = FALSE, cluster_column_slices = FALSE,
        
        row_title = " ", row_title_gp = gpar(fontsize = 10),
        left_annotation = AnnotationDataframeLeftLipidsLTPs2,
        
        row_labels = sapply(strsplit(rownames(NewSphingoTransLayers$InCellulo), "_"), "[[", 1),
        row_split = c(rep("A", 5), rep("B",1), rep("C",3), "D", rep("E",4))
        
)
dev.off() 

# Previous figure is R-basis of final figure. 
# Final figure has LCN1 information removed, legends added, CERS information added to the top, and highlights of most important part focusses by putting black rectagles in Adobe Illustrator.

#### Fig.6c
#. CERTAndSEC14L1LipidChangesHEK293WelchsTTest120820192.pdf #

min(OverexpressionHEK5[,2:ncol(OverexpressionHEK5)][OverexpressionHEK5[,2:ncol(OverexpressionHEK5)] != 0])
# 0.000100165 --> 0.00001 should be fine as background (to avoid 0 issues in the calculations); corrected the x+1 in lines below 

OverexpressionHEKRatiosMatrixAll <- do.call("cbind", lapply(1:9, function(x){(OverexpressionHEK5[,(2*x+1)] + 0.00001)/(OverexpressionHEK5[,(2*x)] + 0.00001)}))
colnames(OverexpressionHEKRatiosMatrixAll) <- paste0(colnames(OverexpressionHEK5)[sapply(1:9,function(x){(2*x+1)})], "_ratio")

rownames(OverexpressionHEKRatiosMatrixAll) <- OverexpressionHEK5[,1]
OverexpressionHEKLogRatiosMatrix <- log10(OverexpressionHEKRatiosMatrixAll)

OverexpressionHEKLogRatiosGeneral2 <- t(OverexpressionHEKLogRatiosMatrix[c("Cer", "CerP", "SM", "HexCer", "SHexCer", "diHexCer", "GM3", "DAG", "PA", "PA O-", "PC", "PC O-", "PE", "PE O-", "PS", "PI", "PI O-", "PG", "CL", "LPA", "LPC", "LPC O-", "LPE", "LPE O-", "LPS", "LPS O-", "LPI", "LPI O-", "LPG", "LPG O-", "CE", "Chol :"),c(1,4,7,2,5,8,3,6,9)])
library("reshape2")

OverexpressionHEKLogRatiosMelted <- melt(OverexpressionHEKLogRatiosGeneral2)
colnames(OverexpressionHEKLogRatiosMelted) <- c("Sample", "Lipid", "Ratio")


OverexpressionHEKLogRatiosMelted2 <- cbind(OverexpressionHEKLogRatiosMelted, LTPType = sapply(strsplit(as.character(OverexpressionHEKLogRatiosMelted[,"Sample"]), split = "_"), "[[", 3))

StatHEKLogRatiosMelted2 <- do.call("rbind",lapply(levels(OverexpressionHEKLogRatiosMelted2[,"Lipid"]), function(y){c(y,sapply(c("CERT","Sec14L1"), function(x){unlist(t.test(OverexpressionHEKLogRatiosMelted2[(OverexpressionHEKLogRatiosMelted2[,"LTPType"] == "control") & (OverexpressionHEKLogRatiosMelted2[,"Lipid"] == y),"Ratio"], OverexpressionHEKLogRatiosMelted2[(OverexpressionHEKLogRatiosMelted2[,"LTPType"] == x) & (OverexpressionHEKLogRatiosMelted2[,"Lipid"] == y),"Ratio"])[c(3,5,4)])}))}))
colnames(StatHEKLogRatiosMelted2) <- c("Lipid", "CERTp.value", "CERTControlMeanEstimate", "CERTProteinMeanEstimate", "CERTDiffConfidenceLow", "CERTDiffConfidenceHigh", "SEC14L1p.value", "SEC14L1ControlMeanEstimate", "SEC14L1ProteinMeanEstimate", "SEC14L1DiffConfidenceLow", "SEC14L1DiffConfidenceHigh")

StatHEKLogRatiosMelted4 <- do.call("cbind", list(StatHEKLogRatiosMelted2, CERTPercentDifference = 100*10^(as.numeric(StatHEKLogRatiosMelted2[,4]) - as.numeric(StatHEKLogRatiosMelted2[,3]))-100, SEC14L1PercentDifference = 100*10^(as.numeric(StatHEKLogRatiosMelted2[,9]) - as.numeric(StatHEKLogRatiosMelted2[,8]))-100))  


pdf("C:/Users/Kevin/Documents/RData/Rest2/ACGData/FiguresByKT/CERTAndSEC14L1LipidChangesHEK293WelchsTTest120820192.pdf")

plot(as.numeric(StatHEKLogRatiosMelted4[,"CERTPercentDifference"]), -log10(as.numeric(StatHEKLogRatiosMelted4[,"CERTp.value"])), xlab = "Percent difference lipid abundance: protein vs. control", ylab = "-log10(p-value) [higher is better] (Welch's 2 sample t-test) ", pch = 16, col = "#FC4E07", main = "CERT (HEK293)", xlim = c(-100,550))
abline(h = -log10(0.05), col = "lightgrey")

points(as.numeric(StatHEKLogRatiosMelted4[as.numeric(StatHEKLogRatiosMelted4[,"CERTp.value"]) <= 0.075,"CERTPercentDifference"]), -log10(as.numeric(StatHEKLogRatiosMelted4[as.numeric(StatHEKLogRatiosMelted4[,"CERTp.value"]) <= 0.075,"CERTp.value"])), pch = 16, col = "#E7B800")
points(as.numeric(StatHEKLogRatiosMelted4[as.numeric(StatHEKLogRatiosMelted4[,"CERTp.value"]) <= 0.05,"CERTPercentDifference"]), -log10(as.numeric(StatHEKLogRatiosMelted4[as.numeric(StatHEKLogRatiosMelted4[,"CERTp.value"]) <= 0.05,"CERTp.value"])), pch = 16, col = "#00AFBB")

text(as.numeric(StatHEKLogRatiosMelted4[as.numeric(StatHEKLogRatiosMelted4[,"CERTp.value"]) <= 0.075,"CERTPercentDifference"]), -log10(as.numeric(StatHEKLogRatiosMelted4[as.numeric(StatHEKLogRatiosMelted4[,"CERTp.value"]) <= 0.075,"CERTp.value"])), labels = paste0(StatHEKLogRatiosMelted4[as.numeric(StatHEKLogRatiosMelted4[,"CERTp.value"]) <= 0.075,"Lipid"], " (", round(as.numeric(StatHEKLogRatiosMelted4[as.numeric(StatHEKLogRatiosMelted4[,"CERTp.value"]) <= 0.075,"CERTPercentDifference"])), "%)"), pos = 4, cex = 0.8)
text(500, -log10(0.05), labels = "p-value 0.05", pos = 3, cex = 0.8, col = "darkgrey")

dev.off()
# The above figure is the basis for Figure panel 6c, and was optimized further and integrated by Adobe Illustrator use.

