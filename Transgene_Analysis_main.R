require(data.table)
require(dplyr)
require(purrr)
require(reshape2)
require(openxlsx)
require(ggplot2)


##### Modify per run ######

path<-'/camhpc/ngs/TechDev/Transgene_Seq/TST11287_Transgene_Seq_method_optimization_detection_levels/snv/Nextera' #input file path
outpath<-'/camhpc/ngs/TechDev/Transgene_Seq/TST11287_Transgene_Seq_method_optimization_detection_levels/R_analysis/Nextera/Eric_script_test' #output file path
excelFilename<-'mAB_Seq_method_optimization_detection_Nextera.xlsx' #will deposit in output filepath
path_source_file <-'/camhpc/ngs/TechDev/Transgene_Seq/TST11287_Transgene_Seq_method_optimization_detection_levels/R_analysis/Nextera/Eric_script_test/transgene_analysis_fxns.R'

#which mutation% thresholds to use for summarizing samples and mutations in excel output
thresholdMuts<-c(1,2,5)
# define Coding Sequence Starts and Stops for Regions (note the end coordinate is last nt of CDS stop codon!)
CDS1_start<-1263
CDS1_stop_end<-2660
# CDS2_start <- NULL
# CDS2_stop_end <- NULL
if ( all( exists("CDS1_start"), exists("CDS2_start"), exists("CDS1_stop_end"), exists("CDS2_stop_end") ) ) {  # if 2 coding regions
  starts<-c(CDS1_start, CDS2_start)      # list of cds starts
  ends<-c(CDS1_stop_end, CDS2_stop_end)  # list of cds ends 
} else {   # if 1 coding region
  starts<-c(CDS1_start)  # list of cds starts
  ends<-c(CDS1_stop_end) # list of cds ends 
}

# define names for CDS regions
CDS1_region_name <- "aSyn_Hc"   # e.g. mAb_HC
# CDS2_region_name <- "aSyn_Lc" # e.g. mAb_LC
if ( all( exists("CDS1_region_name"), exists("CDS2_region_name")) ) {  # if 2 coding regions
  regions<-c(CDS1_region_name, CDS2_region_name)
} else {   # if only 1 coding region
  regions<-c(CDS1_region_name)
}

# provide coverage thresholds for reporting out data in excel file
thresholdCoverage<-500  # coverage minimum across any vector position below which samples will be listed as low coverage
thresholdDeletions<-50  # coverage minimum for deletions across any vector position above wich samples and their deletions
thresholdInsertions<-45 # coverage minimum for insertions across any vector position above wich samples and their insertions

########


source(path_source_file) # source file with functions transgene_analysis_fxns.R

args<-list.files(path, pattern = 'callstats.txt', full.names = T)

if (length(args)==0){
  stop('No callstats files found in this directory')
}

dir.create(outpath, showWarnings = FALSE)
setwd(outpath)

Analysis<-lapply(args, FUN=transgene_analysis, starts=starts, ends=ends, thresholdCoverage=thresholdCoverage, thresholdIns=thresholdInsertions, thresholdDel=thresholdDeletions, thresholdMuts=thresholdMuts)

mutNames<-lapply(Analysis, function(x){return(x[[1]])})
samplesBelowCoverage<-lapply(Analysis, function(x){if(x[[2]]==T){return(x[[1]])}})
samplesWinsertions<-lapply(Analysis, function(x){if(dim(x[[3]])[1]>0){return(x[[1]])}})
samplesWdeletions<-lapply(Analysis, function(x){if(dim(x[[4]])[1]>0){return(x[[1]])}})

insertions<-lapply(Analysis, function(x){return(x[[3]])})
insertions<-insertions %>% reduce(full_join, by = c("CHROM", "pos", "ref"))%>%arrange(pos)

deletions<-lapply(Analysis, function(x){return(x[[4]])})
deletions<-deletions %>% reduce(full_join, by = c("CHROM", "pos", "ref"))%>%arrange(pos)

mutations<-lapply(Analysis, function(x){return(x[[5]])})
mutations<-lapply(as.character(thresholdMuts), function(x){mutations %>% map(x)})
mutations<-lapply(mutations, function(x){x%>%reduce(full_join, by=c("CHROM", "pos", "ref"))%>%arrange(pos)})
names(mutations)<-as.character(thresholdMuts)

samplesWmutations<-lapply(Analysis, function(x){return(x[[6]])})
samplesWmutations<-lapply(as.character(thresholdMuts), function(x){samplesWmutations %>% map(x)})
names(samplesWmutations)<-as.character(thresholdMuts)

CDStables<-lapply(Analysis, function(x){return(x[[7]])})
names(CDStables)<-mutNames

createOpenXlsx(excelFilename, samplesBelowCoverage, samplesWdeletions, samplesWinsertions, insertions, deletions, mutations, samplesWmutations)

maxMut<-max(unlist(lapply(CDStables, function(x){max(x[,Mutant])})))
mat<-cbind(starts, ends, regions)
apply(cbind(names(CDStables), CDStables), 1,FUN=regionPlots, mat=mat, Max=maxMut, interactive=F)

apply(mat, 1, FUN=depthPlotToFile, y=cbind(names(CDStables), CDStables))

list.files()
