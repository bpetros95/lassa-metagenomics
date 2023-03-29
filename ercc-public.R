## version history -------------------------------------------------------------
# v1 by B.A.P. 29 Mar 2023

## setup -----------------------------------------------------------------------
rm(list=ls(all.names=TRUE))
setwd("/insert/wd/here")
suppressPackageStartupMessages(library(tidyverse))
options(stringsAsFactors=FALSE)
theme_set(theme_classic())

## inputs ----------------------------------------------------------------------
ercc <- read.csv("ercc_data.csv")

## assess for sample contamination ---------------------------------------------
assignments = data.frame(sample = colnames(ercc)[2:length(ercc)], ercc = NA, 
                         purity = NA, reads = NA, keep = FALSE)

# for each sample, calculate purity as max % of ERCC reads assigned to a specific ERCC
for (s in assignments$sample){
  reads = sum(ercc[,s])
  assignments$reads[which(assignments$sample == s)] = reads
  assignments$purity[which(assignments$sample == s)] = round(max(ercc[,s])/reads, 3)
  if (reads > 0){
    assignments$ercc[which(assignments$sample == s)] = 
      ercc$seq_mapped_to[which(ercc[,s] == max(ercc[,s]))]}
  else {
    assignments$ercc[which(assignments$sample == s)] = "no reads"}}
rm(reads, s)

## identify samples w/ low ERCC purity------------------------------------------
spiked = data.frame(ercc = unique(assignments$ercc), mean = NA, sd = NA)

# parameterize the distribution of purity scores on per-ERCC basis
for (e in spiked$ercc){
  dt = assignments[assignments$ercc == e & assignments$reads >= 100,]
  spiked$mean[which(spiked$ercc == e)] = mean(dt$purity, na.rm = TRUE)
  spiked$sd[which(spiked$ercc == e)] = sd(dt$purity, na.rm = TRUE)}
rm(dt, e)

# keep samples with high purity scores
for (s in assignments$sample){
  ercc = assignments$ercc[which(assignments$sample == s)]
  purity = assignments$purity[which(assignments$sample == s)]
  # keep samples with < 100 reads regardless of purity (too little data)
  if (assignments$reads[which(assignments$sample == s)] < 100){
    assignments$keep[which(assignments$sample == s)] = TRUE}
  # keep samples with >= 99% purity
  else if (assignments$purity[which(assignments$sample == s)] >= 0.99){
    assignments$keep[which(assignments$sample == s)] = TRUE}
  # keep samples with purity > 3 sd below ERCC's mean purity
  else if (purity > (spiked$mean[which(spiked$ercc == ercc)] - 
                     3*spiked$sd[which(spiked$ercc == ercc)])){
    assignments$keep[which(assignments$sample == s)] = TRUE}}
rm (ercc, s, purity)

## plot per-ERCC purity scores -------------------------------------------------
purity_plot = ggplot(assignments[assignments$reads > 100,], aes(ercc, purity)) + 
  geom_boxplot() +
  scale_y_continuous(limits = c(0, 1.05)) +
  xlab('ERCC') + ylab('Purity') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
purity_plot

