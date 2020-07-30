# Databricks notebook source
library(data.table)
library(dplyr)
library(purrr)
library(reshape2)
library(openxlsx)
library(ggplot2)


# subfunction of CDS_subset - used to subset mutation data to only covering CDS regions
start_end_CDS<-function(mutTable, row){
  start=row[[1]]
  end=row[[2]]
  mutTable<-data.table(mutTable)
  mutTable_CDS<-mutTable[(pos>=start)&(pos<=end)]
  return(mutTable_CDS)
}

# subfunction subsetting mutation data from callstats.txt over whole vector to data only from CDS regions (per sample)
CDS_subset<-function(mutTable, starts, ends){
  mat<-cbind(starts, ends)
  list_of_out<-apply(mat, MARGIN = 1, FUN = start_end_CDS, mutTable=mutTable)
  combined<-do.call(rbind, list_of_out)
  return(combined)
}

# subfunction summarizing point mutations above threshold from callstats.txt file (per sample)
mutCount<-function(mutTable, threshold, Name){
  mutTableThresh<-mutTable[Mutant>threshold]
  mutTableThresh<-mutTableThresh[,.(CHROM, pos, ref, G, A, T, C)]
  names(mutTableThresh)[4:7]<-unlist(lapply(names(mutTableThresh)[4:7], function(x){paste(x, Name, sep='_')}))
  return(mutTableThresh)
}

# subfunction summarizing either deletions or insertions from callstats.txt file (per sample)
inDelCount<-function(mutTable, threshold, inDel, Name){
  if(inDel=='inserted'){count<-'ins_count'}
  else{count<-'del_count'}
  bools<-mutTable[,..count]>threshold
  inDelTable<-mutTable[bools[,],.(CHROM, pos, ref)]
  inDelTable<-cbind(inDelTable, mutTable[bools[,],..inDel],mutTable[bools[,],..count])
  names(inDelTable)[4]<-paste(names(inDelTable)[4], Name, sep='_')
  names(inDelTable)[5]<-paste(names(inDelTable)[5], Name, sep='_')
  return(inDelTable)
}

# MAIN ANALYSIS FUNCTION -  wrapper function of subfunctions per sample and summarizes outputs into a list format
transgene_analysis<-function(mutFile, starts, ends, regions, thresholdCoverage, thresholdIns, thresholdDel, thresholdMuts){
  mutTable<-fread(mutFile)
  mutName<-gsub('.callstats.txt', "", basename(mutFile))
  
  mutTable_CDS<-CDS_subset(mutTable, starts, ends)
  
  if(dim(mutTable_CDS)[1]>0){
    MinimumCoverage<-min(mutTable_CDS$depth)
    if(MinimumCoverage<thresholdCoverage){belowCoverage=T}
    else{belowCoverage=F} #this is a possible place to split this function if manual sample removal is necessary
    #could also see how many reads are below thresh and filter by that, would need to ask Joost
    insertions<-inDelCount(mutTable_CDS, thresholdIns, 'inserted',mutName)
    deletions<-inDelCount(mutTable_CDS, thresholdDel, 'deleted', mutName)
    thresholdMutTables<-lapply(thresholdMuts, FUN=mutCount, mutTable=mutTable_CDS, Name=mutName)
    aboveMutCoverage<-lapply(thresholdMutTables, FUN = function(x){if(dim(x)[1]>0){return(mutName)}})
    names(thresholdMutTables)<-as.character(thresholdMuts)
    names(aboveMutCoverage)<-as.character(thresholdMuts)
    return(list(mutName,belowCoverage, insertions, deletions, thresholdMutTables, aboveMutCoverage,mutTable_CDS))
  }
}

# previous scripts used xlsx read package to create and output mutation data to excel which depends on Java and caused problems (versioning) - here this package has been replaced with java-free package OpenXlsx 
# createXlsx<-function(xlsxName, samplesBelowCoverage, samplesWdeletions, samplesWinsertions, insertions, deletions, mutations, samplesWmutations){
#   addMutationSheets<-function(mutation, sampleWmutation){
#     cov=names(sampleWmutation)
#     addDataFrame(data.frame(unlist(sampleWmutation)), sheet=createSheet(wb=WBdata, sheetName = paste0('Smpls', cov, '%')), col.names=F, row.names=F)
#     addDataFrame(mutation, sheet=createSheet(WBdata, sheetName = paste0('Muts>=', cov, '%')), row.names=F)
#   }  
#   WBdata <- createWorkbook();
#   addDataFrame(data.frame(unlist(samplesBelowCoverage)), sheet=createSheet(wb=WBdata, sheetName = 'SmplsTooLowCov'), col.names=F, row.names=F)
#   for(i in length(mutations):1){
#     addMutationSheets(mutation = mutations[i][[1]], sampleWmutation = samplesWmutations[i])  }
#   addDataFrame(data.frame(unlist(samplesWdeletions)), sheet=createSheet(wb=WBdata, sheetName = 'SmplsDEL'), col.names=F, row.names=F)
#   addDataFrame(deletions, sheet=createSheet(wb=WBdata, sheetName = 'MutsDEL'), row.names=F)
#   addDataFrame(data.frame(unlist(samplesWinsertions)), sheet=createSheet(wb=WBdata, sheetName = 'SmplsINS'), col.names=F, row.names=F)
#   addDataFrame(insertions, sheet=createSheet(wb=WBdata, sheetName = 'MutsINS'), row.names=F)
#   saveWorkbook(WBdata, xlsxName)
# }

createOpenXlsx<-function(xlsxName, samplesBelowCoverage, samplesWdeletions, samplesWinsertions, insertions, deletions, mutations, samplesWmutations){
  addMutationSheets<-function(mutation, sampleWmutation){
    cov=names(sampleWmutation)
    myAddWorkSheet(wb=WBdata, sheetName = paste0('Smpls', cov, '%'))
    writeData(wb = WBdata,data.frame(unlist(sampleWmutation)), sheet=paste0('Smpls', cov, '%'), colNames=F, rowNames=F)
    myAddWorkSheet(wb=WBdata, sheetName = paste0('Muts>=', cov, '%'))
    mutation=Filter(function(x) !(all(x=="")), mutation)
    writeData(wb = WBdata,mutation, sheet=paste0('Muts>=', cov, '%'), rowNames=F)
  }  
  WBdata <- createWorkbook()

  myAddWorkSheet(wb=WBdata, sheetName = 'SmplsTooLowCov')
  writeData(wb = WBdata, data.frame(unlist(samplesBelowCoverage)), sheet= 'SmplsTooLowCov', colNames=F, rowNames=F)
  for(i in length(mutations):1){
    addMutationSheets(mutation = mutations[i][[1]], sampleWmutation = samplesWmutations[i])  }
  myAddWorkSheet(wb=WBdata, sheetName = 'SmplsDEL')
  writeData(wb = WBdata,data.frame(unlist(samplesWdeletions)), sheet='SmplsDEL', colNames=F, rowNames=F)
  myAddWorkSheet(wb=WBdata, sheetName = 'MutsDEL')
  deletions=Filter(function(x) !(all(x=="")), deletions)
  writeData(wb = WBdata,deletions, sheet='MutsDEL', rowNames=F)
  myAddWorkSheet(wb = WBdata, sheetName = 'SmplsINS')
  writeData(wb = WBdata,data.frame(unlist(samplesWinsertions)), sheet= 'SmplsINS', colNames=F, rowNames=F)
  myAddWorkSheet(wb=WBdata, sheetName = 'MutsINS')
  insertions=Filter(function(x) !(all(x=="")), insertions)
  writeData(wb = WBdata,insertions, sheet='MutsINS', rowNames=F)
  saveWorkbook(WBdata, xlsxName, overwrite = T)
}

myAddWorkSheet<-copy(addWorksheet) #addWorksheet isn't properly exported or something

regionPlots<-function(x, mat, Max, interactive){
  if(interactive){    apply(mat, 1, FUN = regionPlot, y=x, Max=Max)  }
  else{    apply(mat, 1, FUN = regionPlotToFile, y=x, Max=Max)  }
}

# for creating interactive plots in Databricks - commented out here
# regionPlot<-function(y, row, Max){
  # sample<-y[[1]]
  # start<-as.numeric(row[[1]])
  # end<-as.numeric(row[[2]])
  # regionName<-row[[3]]
  # toPlot<-y[[2]]
  # plot(toPlot[pos>=start&pos<=end, pos],toPlot[pos>=start&pos<=end, Mutant], main=paste(regionName, sample, sep = '_'), xlab='NTposition', ylab='%MutatedBases', ylim=c(0,Max))
  # lines(toPlot[pos>=start&pos<=end, pos],toPlot[pos>=start&pos<=end, Mutant], type='h')
# }

regionPlotToFile<-function(y, row, Max){
  sample<-y[[1]]
  start<-as.numeric(row[[1]])
  end<-as.numeric(row[[2]])
  regionName<-row[[3]]
  toPlot<-y[[2]]
  print(paste('Plotting region', regionName, 'for sample', sample, 'to file', sep=' '))
  png(paste0(regionName, '_',sample, '.png'))
  plot(toPlot[pos>=start&pos<=end, pos],toPlot[pos>=start&pos<=end, Mutant], main=paste(regionName, sample, sep = '_'), xlab='NTposition', ylab='%MutatedBases', ylim=c(0,Max))
  lines(toPlot[pos>=start&pos<=end, pos],toPlot[pos>=start&pos<=end, Mutant], type='h')
  dev.off()
}

# for creating interactive plots in Databricks - commented out here
# depthPlots<-function(z, mat, interactive){
  # if(interactive){ apply(mat, 1, FUN=depthPlot, y=z)  }
  # else{ apply(mat, 1, FUN=depthPlotToFile, y=z)  }
# }
# depthPlot<-function(y, row){
  # start<-as.numeric(row[1])
  # end<-as.numeric(row[2])
  # regionName<-row[3]
  # tabNames<-y[,1]
  # tabs<-y[,2]
  # tabs<-lapply(tabs, function(x){x[pos>=start&pos<=end,.(CHROM,pos,ref,depth)]})
  # maxDepth<-max(unlist(lapply(tabs, function(x){max(x[,depth])})))
  
  # merged<-tabs %>% reduce(full_join, by = c("CHROM", "pos", "ref"))%>%arrange(pos)
  # names(merged)[4:length(names(merged))]<-unlist(tabNames)
  # merged<-melt(merged, id=c('CHROM', 'pos', 'ref'))
  # ggplot(merged)+
    # aes(x = pos, y = value, color = variable)+
    # ylab('depth')+
    # geom_point()
# }

depthPlotToFile<-function(y, row){
  start<-as.numeric(row[1])
  end<-as.numeric(row[2])
  regionName<-row[3]
  tabNames<-y[,1]
  tabs<-y[,2]
  tabs<-lapply(tabs, function(x){x[pos>=start&pos<=end,.(CHROM,pos,ref,depth)]})
  maxDepth<-max(unlist(lapply(tabs, function(x){max(x[,depth])})))
  
  merged<-tabs %>% reduce(full_join, by = c("CHROM", "pos", "ref"))%>%arrange(pos)
  names(merged)[4:length(names(merged))]<-unlist(tabNames)
  merged<-melt(merged, id=c('CHROM', 'pos', 'ref'))
  png(paste0(regionName, '_Coverage_AllSamples.png'))
  print({
    #png(paste0(regionName, '_Coverage_AllSamples.png'), width = '240', height = '240')
    ggplot(merged)+
      aes(x = pos, y = value, color = variable)+
      ylab('depth')+
      geom_point()+
      theme(axis.text=element_text(size=10))
    #ggsave(paste0(regionName, '_Coverage_AllSamples.png'), width = 4, height=2.5)
  })
  dev.off()
}