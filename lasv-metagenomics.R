## version history -------------------------------------------------------------
# v1 by B.A.P. 29 Mar 2023

## setup -----------------------------------------------------------------------
rm(list=ls(all.names=TRUE))
setwd("/insert/wd/here")
suppressPackageStartupMessages(library(tidyverse))
options(stringsAsFactors=FALSE)
theme_set(theme_classic())

## inputs ----------------------------------------------------------------------
samples = read.csv("data/in/samples.csv", na.strings="", row.names = 1) %>%
  mutate_at(c('LF_sx', 'siddle_2018', 'cluster', 'unknown', 'tree'), as.numeric)
samples = samples[samples$ERCC_flag != 1,] #rm samples w/ low ERCC purity scores

# reads per virus & per sample
reads <- read.csv("data/in/reads.csv", na.strings="", row.names = 1)
reads = data.frame(rows=rownames(reads)[row(reads)],
                   cols=colnames(reads)[col(reads)],
                   value=as.numeric(unlist(reads)))
colnames(reads) <- c('sample', 'virus', 'reads')

# reads per virus & per control
ctrl <-read.csv("data/in/neg_control_reads.csv")

## filtering--------------------------------------------------------------------
# determine max percent of reads assigned to a given virus across controls per 
# seq batch, among viruses that surpass the 5 read cut-off
ctrl <- ctrl %>% 
  filter(reads >= 5) %>%
  mutate(prct = reads/total_reads) %>%
  group_by(virus, batch) %>%
  summarise(threshold = max(prct)) %>%
  as.data.frame()
rownames(ctrl) <- paste(ctrl$virus, ctrl$batch, sep = "_")

# flag samples with at least 5 reads and a greater % of reads assigned to virus 
# than was assigned to that virus in batch-specific controls
reads = reads %>% 
  mutate(prct = reads/samples[sample, "reads"]) %>%
  mutate(batch = samples[sample, "batch"]) %>%
  mutate(threshold = ifelse(paste(virus, batch, sep = "_") %in% row.names(ctrl), 
                            ctrl[paste(virus, batch, sep = "_"), 'threshold'], 0)) %>%
  mutate(keep = ((reads >= 5) & (prct > threshold)))
rm(ctrl)

# filter reads for those meeting thresholds
reads_filt = reads %>% filter(keep)

# assembly performed for all rows in reads_filt in DNA Nexus via the viral-ngs
# software, filtered to rows where assembly length > 10% of reference genome 
# length, and manually grouped by viral family
reads_asbl = read.csv("data/in/viruses_family.csv")

## characterize samples w/ symptoms of Lassa Fever -----------------------------
asbl_lasv = reads_asbl %>% filter(virus == 'LASV')

# isolate samples with LF symptoms
lassa <- reads %>% 
  mutate(case_ctrl = samples[sample, "case_ctrl"]) %>% 
  mutate(LF_sx = samples[sample, "LF_sx"]) %>% 
  filter(!is.na(LF_sx)) %>%
  filter(virus == 'Lassa_mammarenavirus') %>%
  mutate(PCR = ifelse(case_ctrl == "control", "neg", "pos")) 
lassa$asbl = lassa$sample %in% asbl_lasv$sample
rownames(lassa) = lassa$sample
rm(asbl_lasv)

# compare RT-qPCR positivity to metagenomic positivity
results <- data.frame(PCR = c('positive', 'positive', 'negative', 'negative'),
                      Metagenomics = c('positive', 'negative', 'positive', 'negative'),
                      Counts = c(sum(lassa$PCR == 'pos' & lassa$asbl),
                                 sum(lassa$PCR == 'pos' & !(lassa$asbl)),
                                 sum(lassa$PCR == 'neg' & lassa$asbl),
                                 sum(lassa$PCR == 'neg' & !(lassa$asbl))))
p <- fisher.test(matrix(results$Counts, nrow = 2, ncol = 2))$p.value

# heatmap of PCR vs metagenomic positivity
pcrplot <- ggplot(results, aes(PCR, Metagenomics)) + 
  geom_tile(aes(fill = Counts)) +
  geom_text(aes(label = Counts)) +
  scale_fill_gradient(low = "grey", high = "blue") +
  labs(x="PCR", y="Metagenomics", fill = 'Number of Cases')
pcrplot
rm(p, results)

# count number of LASV+ and LASV- samples containing each virus
counts = reads_asbl %>% mutate(LASV = lassa[sample, "PCR"]) %>%
  mutate(case_ctrl = lassa[sample, "case_ctrl"]) %>% 
  filter(!is.na(LASV)) %>%
  group_by(LASV, virus) %>%
  summarise(cases=n()) %>%
  filter(virus != "LASV") %>% ungroup %>%
  complete(LASV, virus, fill = list(cases = 0)) %>%
  group_by(LASV, virus, cases) %>%
  summarise(tot = if (LASV == 'pos') sum(lassa$PCR == 'pos') 
            else sum(lassa$PCR == 'neg')) %>%
  mutate(prct = cases/tot) %>%
  mutate(lab = ifelse(LASV == "pos", "positive", "negative"))

# heatmap of viruses by LASV RT-qPCR status
hm <- ggplot(counts, aes(virus, lab)) + 
  geom_tile(aes(fill = counts$prct)) +
  geom_text(aes(label = cases)) +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.75, hjust = 0.5)) +
  scale_fill_gradient(low = "light gray", high = "blue") +
  labs(y="Lassa Virus PCR", x="Assembled Viral Genome", fill = 'Percent of Cases')
hm   
