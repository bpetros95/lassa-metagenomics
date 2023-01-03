
## setup -----------------------------------------------------------------------
rm(list=ls(all.names=TRUE))
setwd("/insert/path/here")
lapply(list('reshape2', 'tidyverse'), library, character.only = TRUE)
options(stringsAsFactors=FALSE)

## wrangle inputs --------------------------------------------------------------

# sample metadata
samples = read.csv("data/in/samples.csv", na.strings="", row.names = 1) %>%
  mutate_at(c('LF_sx', 'siddle_2018', 'cluster', 'unknown', 'tree'), as.numeric)

# sequencing metadata
meta = read.csv("data/in/cleaned-meta.csv", na.strings = "", row.names = 1)

# reads per virus, per sample
reads <- read.csv("data/in/reads.csv", na.strings="", row.names = 1)

# isolate reads from samples in samples df (NB: removes technical or biological 
# replicates, retaining sample with greatest number of reads)
reads = reads[rownames(reads) %in% rownames(samples), ]
reads = data.frame(rows=rownames(reads)[row(reads)], 
                   cols=colnames(reads)[col(reads)],
                   value=as.numeric(unlist(reads)))
colnames(reads) <- c('sample', 'virus', 'reads')

# maximum percent of reads assigned to a given virus across water (NTC) controls 
# in seq batch, among viruses assigned minimally 5 reads
water <-read.csv("data/in/ntc-counts.csv") %>% 
  filter(count >= 5) %>%
  mutate(prct = count/reads) %>%
  group_by(virus, batch) %>%
  summarise(threshold = max(prct))
water = as.data.frame(water)
rownames(water) <- paste(water$virus, water$batch, sep = "_")

# flag taxa that have at least 5 reads and a greater % of reads assigned to them 
# than maximal % assigned among batch NTCs 
reads = reads %>% 
  mutate(prct = reads/samples[sample, "reads"]) %>%
  mutate(batch = samples[sample, "batch"]) %>%
  mutate(threshold = ifelse(paste(virus, batch, sep = "_") %in% row.names(water), 
                    water[paste(virus, batch, sep = "_"), 'threshold'], 0)) %>%
  mutate(keep = ((reads >= 5) & (prct > threshold)))

# filter reads based on flag and write to file
reads_filt = reads %>% filter(keep)
write.csv(reads_filt, "data/out/filt-reads.csv", row.names=TRUE)

# import assemblies performed on all sample-taxon pairs in reads_filt using the
# viral-ngs software (https://viral-ngs.readthedocs.io/en/latest/) in DNA Nexus 
reads_asbl = read.csv("data/in/filt-assembled.csv", na.strings="") %>%
  mutate(prct = assembly_length/ref_genome_length) %>%
  filter(prct >= 0.1) %>%
  mutate(isth = meta[sample, "isth"])

# write file with sample-taxon pairs that surpassed 10% genome threshold
write.csv(reads_asbl, "data/out/viruses_present.csv", row.names=FALSE)

# import sample-taxon pairs that were manually parsed to remove known lab 
# contaminants (e.g., MLV) and collapse distinct viral species at genus/family levels
reads_asbl = read.csv("data/in/viruses_family.csv")

# isolate samples from patients suspected of Lassa Fever
lassa <- reads %>% 
  mutate(case_ctrl = samples[sample, "case_ctrl"]) %>% 
  mutate(LF_sx = samples[sample, "LF_sx"]) %>% 
  mutate(isth = meta[sample, "isth"]) %>% 
  filter(!is.na(LF_sx)) %>%
  filter(virus == 'Lassa_mammarenavirus') %>%
  mutate(PCR = ifelse(case_ctrl == "control", "neg", "pos"))
rownames(lassa) = lassa$sample

## RT-qPCR vs metagenomics ------------------------------------------------------
results <- data.frame(PCR = c('pos', 'pos', 'neg', 'neg'),
                      Metagenomics = c('pos', 'neg', 'pos', 'neg'),
                      Counts = c(sum(lassa$PCR == 'pos' & lassa$reads >=5),
                                 sum(lassa$PCR == 'pos' & lassa$reads <5),
                                 sum(lassa$PCR == 'neg' & lassa$reads >=5),
                                 sum(lassa$PCR == 'neg' & lassa$reads <5)))

p <- fisher.test(matrix(results$Counts, nrow = 2, ncol = 2))$p.value

# contingency table of PCR vs. metagenomic positivity
pcr <- ggplot(results, aes(PCR, Metagenomics)) + 
  geom_tile(aes(fill = Counts)) +
  geom_text(aes(label = Counts)) +
  scale_fill_gradient(low = "grey", high = "blue") +
  labs(x="PCR", y="Metagenomics (5+ Reads)", fill = 'Number of Cases')
pcr

# compare number of reads mapped to LASV for PCR-positive vs. -negative samples
p <- wilcox.test(lassa$reads[lassa$PCR == 'pos'], 
                 lassa$reads[lassa$PCR == 'neg'], alternative = "greater")

## plot metagenomic reads as function of PCR positivity
rd <- ggplot(lassa, aes(PCR, reads)) + 
  geom_boxplot(fill = 'grey', color = 'black') +
  scale_y_log10(limits=c(1,10000000)) +
  labs(x="PCR", y="LASV Reads") +
  geom_jitter(shape=16, position=position_jitter(0.2), color = 'blue') +
  annotate("text",x=1.5, y=10000000,label= "p < 0.001", size = 4) +
  geom_segment(aes(x = 1, y = 6000000, xend = 2, yend = 6000000))
rd

## viruses in LASV+ and LASV- samples ------------------------------------------

# count number of LASV+ and LASV- samples containing each virus
counts = reads_asbl %>% mutate(LASV = lassa[sample, "PCR"]) %>%
  mutate(case_ctrl = lassa[sample, "case_ctrl"]) %>% 
  filter(!is.na(LASV)) %>%
  group_by(LASV, virus) %>%
  summarise(cases=n()) %>%
  filter(virus != "LASV") %>%
  ungroup %>%
  complete(LASV, virus, fill = list(cases = 0)) %>%
  group_by(LASV, virus, cases) %>%
  summarise(tot = if (LASV == 'pos') sum(lassa$PCR == 'pos') 
            else sum(lassa$PCR == 'neg')) %>%
  mutate(prct = cases/tot)

# heatmap of identified viruses
hm <- ggplot(counts, aes(virus, LASV)) + 
  geom_tile(aes(fill = counts$prct)) +
  geom_text(aes(label = cases)) +
  theme(axis.text.x = element_text(angle = 20, vjust = 0.75, hjust = 0.5)) +
  scale_fill_gradient(low = "light gray", high = "blue") + 
  labs(y="Lassa Virus PCR", x="Assembled Viral Genome", fill = 'Percent of Cases')
hm   

# number of LASV+ samples with 1 or more viral co-infections
sum(lassa$sample[lassa$PCR == 'pos'] %in% filter(reads_asbl, virus!='LASV')$sample)

# number of LASV- samples with 1 or more viral infections
sum(lassa$sample[lassa$PCR == 'neg'] %in% reads_asbl$sample)

## taxa enrichment in LASV+ vs. LASV- samples collected from ISTH in 2018-------

# count number of LASV+ and LASV- samples containing each virus
counts_enrich = reads_asbl %>%
  mutate(LASV = lassa[sample, "PCR"]) %>%
  mutate(case_ctrl = lassa[sample, "case_ctrl"]) %>%
  filter(case_ctrl != 'NA')  %>%
  group_by(LASV, virus) %>%
  summarise(cases=n()) %>%
  filter(virus != "LASV") %>%
  ungroup %>%
  complete(LASV, virus, fill = list(cases = 0)) %>%
  group_by(LASV, virus, cases) %>%
  summarise(tot = ifelse((LASV == 'pos'), sum(lassa$case_ctrl == 'case'),
                         sum(lassa$case_ctrl == 'control'))) %>%
  mutate(prct = cases/tot)%>%
  filter(virus != 'LASV')

# calculate p-vals via Fisher's exact test
pvals = data.frame(virus = unique(counts_enrich$virus), p = 1)

for (vir in pvals$virus){
  pos_cases = counts_enrich$cases[counts_enrich$LASV == 'pos' & counts_enrich$virus == vir]
  neg_cases = counts_enrich$cases[counts_enrich$LASV == 'neg' & counts_enrich$virus == vir]
  mat = matrix(c(pos_cases, sum((lassa$PCR == 'pos') & lassa$isth) - pos_cases, neg_cases,
                 sum((lassa$PCR == 'neg') & lassa$isth) - neg_cases), nrow = 2, ncol = 2)
  pvals$p[pvals$virus == vir] = fisher.test(mat)$p.value
}
rm(mat, neg_cases, pos_cases, vir)

# adjust p-vals via Benjamini-Yekutieli method
pvals$fdr = p.adjust(pvals$p, method = "BY")

# bar chart of identified viruses in 2018 LASV+ vs LASV- ISTH samples
bc <- ggplot(counts_enrich, aes(fill=LASV, y=prct, x=virus)) +
  geom_bar(position="dodge", stat="identity") +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  scale_fill_manual(values=c("neg"="grey","pos"="purple")) +
  labs(x="Virus", y="Percent of Cases", fill = 'Lassa Virus PCR')
bc

## figure 2 --------------------------------------------------------------------
fig2 = cowplot::plot_grid(cowplot::plot_grid(pcr, rd, ncol = 2, 
                labels=c("A", "B")), hm, bc, nrow = 3,labels=c(" ", "C", "D"))
fig2
