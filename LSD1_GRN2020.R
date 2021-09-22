################
## GRN LSD1 4 studies
## Files aligned and prepared by Andrew
## The GRN generated with Christian's script
## 04.05.2020
## Daria
################

if (file.exists("/g/scb/zaugg/bunina/RNAseq/bodyPcrispr.R")) {
  dir = "/g/scb"
} else {
  dir = "/Volumes"
}
file.exists(file = file.path(dir, "zaugg/bunina/RNAseq", "bodyPcrispr.R"))

if (file.exists("/g/scb2/zaugg/bunina/WT129/crisprRNA/data/bams/Elavl18cl1.sorted.bam")) {
  dir1 = "/g/scb2/zaugg"
} else {
  dir1 = "/Volumes/zaugg-1"
}
file.exists(file = file.path(dir1, "bunina/WT129/crisprRNA/data/bams", "Elavl18cl1.sorted.bam"))

showMemoryUse <- function(sort="size", decreasing=TRUE, limit) {
  
  objectList <- ls(parent.frame())
  
  oneKB <- 1024
  oneMB <- 1048576
  oneGB <- 1073741824
  
  memoryUse <- sapply(objectList, function(x) as.numeric(object.size(eval(parse(text=x)))))
  
  memListing <- sapply(memoryUse, function(size) {
    if (size >= oneGB) return(paste(round(size/oneGB,2), "GB"))
    else if (size >= oneMB) return(paste(round(size/oneMB,2), "MB"))
    else if (size >= oneKB) return(paste(round(size/oneKB,2), "kB"))
    else return(paste(size, "bytes"))
  })
  
  memListing <- data.frame(objectName=names(memListing),memorySize=memListing,row.names=NULL)
  
  if (sort=="alphabetical") memListing <- memListing[order(memListing$objectName,decreasing=decreasing),] 
  else memListing <- memListing[order(memoryUse,decreasing=decreasing),] #will run if sort not specified or "size"
  
  if(!missing(limit)) memListing <- memListing[1:limit,]
  
  print(memListing, row.names=FALSE)
  return(invisible(memListing))
}

library("DESeq2")
library("gplots")
library("RColorBrewer")
library("pheatmap")
#library("plyr")
library("ggplot2")
# library("ComplexHeatmap")#
library("dplyr")
library(tidyverse)
library(reshape2)
library(IRanges)
library(GenomicRanges)
library(GenomicAlignments)
library(biomaRt)
library(rtracklayer)
library(ggrepel)
library(viridis)
library(ChIPQC)
library(ChIPseeker)
# R version 3.5.1


# Initial filter by R>0.4 -------------------------------------------------


#check out GRN output files:
GRNtable = read.delim(file = file.path(dir1, "bunina/LSD1/GRN/output/new_run_2020", "GRN_withGenes_TF-PeakFDR0.3_peakGene_rawP0.3.tsv.gz"))
colnames(GRNtable)

##select sigCorrelated peaks-genes and work with those:
#only take positive correlations, as these are peak-gene correlations and we assume peak accessibility has to be positively correlated with target express
GRNtable_sig = GRNtable %>%
  subset(r_peak_gene > 0.4)
summary(GRNtable_sig$r_peak_gene)
summary(GRNtable_sig$p.adj_peak_gene)
# saveRDS(GRNtable_sig, file = file.path(dir1, "bunina/LSD1/GRN/output/new_run_2020", "GRNtable_sig.rds")) #r > 0.4
GRNtable_sig = readRDS(file.path(dir1, "bunina/LSD1/GRN/output/new_run_2020", "GRNtable_sig.rds"))

###
## also filter by abs.corr values!!!
###

GRNtable_sig = GRNtable_sig %>%
  subset(abs(r_TF_peak) > 0.4)

#more stringent corr threshold:
# GRNtable_sig = GRNtable_sig %>%
#   subset(r_peak_gene > 0.6)

#explore the table a bit:
GRNtable_stats = as.data.frame(table(droplevels(GRNtable_sig[,c(1,6)])))
# length(GRNtable_stats$Var1[GRNtable_stats$Freq > 100]) #204
colnames(GRNtable_stats)

ggplot(data = GRNtable_stats, aes(x = TF, y = Freq)) + geom_bar(stat = "identity", aes(fill = fdr_direction_TF_peak)) + theme_classic()

#wide format:

GRNtable_stats_wide = GRNtable_stats %>%
  spread(fdr_direction_TF_peak, Freq)

sum(GRNtable_stats_wide$neg) + sum(GRNtable_stats_wide$pos)

library(ggrepel)
ggplot(data = GRNtable_stats_wide, aes(x = neg, y = pos)) + geom_point() + theme_classic() + scale_x_log10() + scale_y_log10() + 
  geom_text_repel(data = subset(GRNtable_stats_wide, pos > 0 & neg > 0), aes(label = TF), 
                  box.padding = 0.35, size = 2, segment.size = 0.2)
  #geom_label(data = subset(GRNtable_stats_wide, pos > 0 & neg > 0), size = 0.6, aes(label = TF))

GRNtable_stats_wide$type = ifelse(GRNtable_stats_wide$neg > 0 & GRNtable_stats_wide$pos, "both", ifelse(GRNtable_stats_wide$neg > 0, "neg", "pos"))
table(GRNtable_stats_wide$type)



# Filter GRN by protein coding & lncRNAs ----------------------------------

library(org.Hs.eg.db)
library(biomaRt)
mart_hg19 = useMart(host='grch37.ensembl.org', biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")

filt_hg19 = listFilters(mart_hg19)
filt_hg19[filt_hg19$name == "biotype",]
listFilterValues(mart = mart_hg19, filter = "biotype")

#check which attribute values are biotype (to find lncRNAs):
attr_hg19 = listAttributes(mart_hg19)
egs_hg19_all = getBM(attributes = c('ensembl_gene_id', 'chromosome_name','start_position', 'external_gene_name', 
                                    'end_position','strand', 'gene_biotype', 'description'), 
                     mart=mart_hg19)
unique(egs_hg19_all$gene_biotype) #of course, its "lincRNA" instead of lncRNA like on their website! what a mess...

# ds_hg19 = useDataset('hsapiens_gene_ensembl', mart=mart_hg19)
egs_hg19 = getBM(attributes = c('ensembl_gene_id', 'chromosome_name','start_position', 'external_gene_name', 'end_position','strand', 'entrezgene_id'), 
                 mart=mart_hg19, filters = "biotype", values = c("protein_coding", "lincRNA")) #30k instead of 22k protein-coding
egs_hg19$TSS = ifelse(egs_hg19$strand == 1, egs_hg19$start_position, egs_hg19$end_position)

GRNtable_sig = GRNtable_sig %>%
  subset(ENSEMBL %in% egs_hg19$ensembl_gene_id)

saveRDS(GRNtable_sig[,1:12], file = file.path(dir1, "bunina/LSD1/GRN/output/new_run_2020", "GRNtable_peakGene_ProtCod.rds"))
GRNtable_sig = readRDS(file.path(dir1, "bunina/LSD1/GRN/output/new_run_2020", "GRNtable_peakGene_ProtCod.rds"))

###
### liftover to hg38 for Annique's GWAS:
###
library(liftOver)
hg19to38_chain = import.chain(file.path(dir1, "bunina/LSD1/GRN/data", "hg19ToHg38.over.chain.gz"))
#gives an error of the chain file. save as bed and liftover via the website!
GRNtable_sig_df = as.data.frame(GRNtable_sig)

GRNtable_sig_df$chr = sapply(strsplit(as.character(GRNtable_sig_df$peak), ":"), "[", 1)
GRNtable_sig_df$range = sapply(strsplit(as.character(GRNtable_sig_df$peak), ":"), "[", 2)
GRNtable_sig_df$start = sapply(strsplit(as.character(GRNtable_sig_df$range), "-"), "[", 1)
GRNtable_sig_df$end = sapply(strsplit(as.character(GRNtable_sig_df$range), "-"), "[", 2)

GRNtable_sig_df_cc = GRNtable_sig_df[!duplicated(GRNtable_sig_df$peak),]
GRNtable_sig_df_cc$row = 1:nrow(GRNtable_sig_df_cc)
write.table(GRNtable_sig_df_cc[,c("chr", "start", "end", "strand", "row")], file = file.path(dir1, "bunina/LSD1/GRN/output/new_run_2020", "GRNtable_peakGene_ProtCod.bed"), 
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE) #has to be row name, not the peak id - otherwise does not output read position!

GRNtable_sig_df_hg38 = read.table(file = file.path(dir1, "bunina/LSD1/GRN/output/new_run_2020", "GRNtable_peakGene_ProtCod_hg38.bed"))
GRNtable_sig_df_hg38$newPeak = paste0(GRNtable_sig_df_hg38$V1, ":", GRNtable_sig_df_hg38$V2, "-", GRNtable_sig_df_hg38$V3)

GRNtable_sig_df_cc$newPeak = GRNtable_sig_df_hg38$newPeak[match(GRNtable_sig_df_cc$row, GRNtable_sig_df_hg38$V5)]
GRNtable_sig_df$newPeak = GRNtable_sig_df_cc$newPeak[match(GRNtable_sig_df$peak, GRNtable_sig_df_cc$peak)]
GRNtable_sig_df_hg38m = GRNtable_sig_df[complete.cases(GRNtable_sig_df$newPeak),]
GRNtable_sig_df_hg38m$peak = GRNtable_sig_df_hg38m$newPeak
write.csv(GRNtable_sig_df_hg38m[,1:8], file = file.path(dir1, "bunina/LSD1/GRN/output/new_run_2020", "GRNtable_peakGene_ProtCod_hg38forGWAS.csv"), 
          row.names = FALSE, quote = FALSE)

#add single TSS position to calculate one distance to peak:
GRNtable_sig$geneTSS = egs_hg19$TSS[match(GRNtable_sig$ENSEMBL, egs_hg19$ensembl_gene_id)]

#add minimum distance from either of peak ends to the gene's TSS:
GRNtable_sig$peak_left = sapply(strsplit(sapply(strsplit(as.character(GRNtable_sig$peak), ":"), "[", 2), "-"), "[",1)
GRNtable_sig$peak_right = sapply(strsplit(sapply(strsplit(as.character(GRNtable_sig$peak), ":"), "[", 2), "-"), "[",2)
GRNtable_sig$peak_gene_dist = rowMins(cbind(abs(as.numeric(GRNtable_sig$peak_left) - GRNtable_sig$geneTSS), abs(as.numeric(GRNtable_sig$peak_right) - GRNtable_sig$geneTSS)))

#clean up the table:
GRNtable_sig_old = GRNtable_sig
# GRNtable_sig = GRNtable_sig[,-c(9:11,18:19)]
GRNtable_sig = GRNtable_sig[,1:12]


# GRN stats for Figure2B - NEW 2021 ---------------------------------------

GRNtable_sig = readRDS(file.path(dir1, "bunina/LSD1/GRN/output/new_run_2020", "GRNtable_peakGene_ProtCod.rds"))
# GRNtable_sig_df = as.data.frame(GRNtable_sig)
#see above to get the object below:
GRNtable_sig_df_cc$geneTSS = egs_hg19$TSS[match(GRNtable_sig_df_cc$ENSEMBL, egs_hg19$ensembl_gene_id)]
GRNtable_sig_df_cc$start = as.numeric(GRNtable_sig_df_cc$start)
GRNtable_sig_df_cc$end = as.numeric(GRNtable_sig_df_cc$end)
GRNtable_sig_df_cc$distToPeak = min(GRNtable_sig_df_cc$start - GRNtable_sig_df_cc$geneTSS, GRNtable_sig_df_cc$end - GRNtable_sig_df_cc$geneTSS)

GRNtable_sig_df_cc = GRNtable_sig_df_cc %>%
  mutate(dist1 = abs(start - geneTSS), dist2 = abs(end - geneTSS) ) %>%
  rowwise() %>%
  mutate(distToPeak = min(dist1, dist2))

GRNtable_sig_df_cc_gr = GRanges(seqnames = GRNtable_sig_df_cc$chr, ranges = IRanges(start = GRNtable_sig_df_cc$start, end = GRNtable_sig_df_cc$end), 
                                peak = GRNtable_sig_df_cc$peak, gene = GRNtable_sig_df_cc$ENSEMBL)

egs_hg19_gr = GRanges(seqnames = paste0("chr", egs_hg19$chromosome_name), ranges = IRanges(start = egs_hg19$TSS, 
                                                                                           end = egs_hg19$TSS), 
                      gene = egs_hg19$ensembl_gene_id)

nearestGene = distanceToNearest(GRNtable_sig_df_cc_gr, egs_hg19_gr)
mcols(nearestGene)$egs_gene = mcols(egs_hg19_gr)$gene[subjectHits(nearestGene)]

mcols(GRNtable_sig_df_cc_gr)$nearestGene = mcols(nearestGene)$egs_gene[queryHits(nearestGene)]

GRNtable_sig_df_cc$nearestGene = mcols(GRNtable_sig_df_cc_gr)$nearestGene[match(GRNtable_sig_df_cc$peak, mcols(GRNtable_sig_df_cc_gr)$peak)]
GRNtable_sig_df_cc$sameGene = ifelse(GRNtable_sig_df_cc$ENSEMBL == GRNtable_sig_df_cc$nearestGene, TRUE, FALSE)
table(GRNtable_sig_df_cc$sameGene) # 13114 / 96050

pdf(file = file.path(dir1, "bunina/LSD1/GRN/output/plots/newGRN2020", "BarplotPeakToNearbyGene.pdf"), 
    useDingbats = FALSE, width = 2.2, height = 2)
as.data.frame(table(GRNtable_sig_df_cc$sameGene)) %>%
  ggplot(aes(x = "GRN", y = Freq, fill = Var1)) + geom_bar(stat = "identity") +
  theme_classic() + scale_fill_brewer(palette = "Dark2") + ylab("Number of GRN peaks") + xlab("Peaks-gene links")
dev.off()

#check:
GRNtable_sig_df_cc %>%
  subset(distToPeak < 1000000) %>%
  ggplot(aes(x = distToPeak )) + geom_histogram(bins = 100) + theme_classic()

#now filter out > 500kb ones:
pdf(file = file.path(dir1, "bunina/LSD1/GRN/output/plots/newGRN2020", "BarplotDistToGene.pdf"), 
    useDingbats = FALSE, width = 4, height = 2)
GRNtable_sig_df_cc %>%
  subset(distToPeak < 500000) %>%
  ggplot(aes(x = distToPeak, fill = sameGene )) + geom_histogram(bins = 100) + theme_classic() + scale_fill_brewer(palette = "Dark2")
dev.off()

nrow(subset(GRNtable_sig_df_cc, distToPeak < 1000000)) #95512, so 538 regions with higher dist
96050 - nrow(subset(GRNtable_sig_df_cc, distToPeak < 500000)) #7020 regions! which genes are these:

GRNtable_sig_df_cc %>%
  subset(distToPeak > 500000) %>%
  group_by(sameGene) %>%
  distinct(ENSEMBL, .keep_all = TRUE) %>%
  tally()

2662+63 # N of unique genes

###
### check in the new GRN:
###
library(GRNdev)
GRNnew = readRDS(file = "/g/scb/zaugg/zaugg_shared/data/GRN/data/Daria_iPSCs_mesoderm_cardiomyocytes/outputNew/GRN.rds")
GRNnew = GRN::add_TF_gene_correlation(GRNnew, corMethod = "pearson", nCores = 1, forceRerun = TRUE)
GRNnew.all = GRN::getGRNConnections(GRNnew, type = "TF_peaks", include_TF_gene_correlations = TRUE)
GRNnew.all

#now peak-gene:
GRNnew.all_genes = GRN::getGRNConnections(GRNnew, type = "peak_genes", include_TF_gene_correlations = TRUE)
GRNnew.all_genes

GRNnew.expr = read.table("/g/scb/zaugg/zaugg_shared/data/GRN/data/Daria_iPSCs_mesoderm_cardiomyocytes/outputNew/GRN_TF-PeakFDR0.3_expression_peakGeneFDR0.3_allowMissingGenes_perm0.tsv.gz")
colnames(GRNnew.expr)
GRNnew.expr[1,]
colnames(GRNnew.expr) <- GRNnew.expr[1,]
GRNnew.expr = GRNnew.expr[-1,] #24021238
# GRNnew.expr = subset(GRNnew.expr, !is.na(gene.ENSEMBL) )
unique(GRNnew.expr$gene.type)

GRNnew.expr_cc = GRNnew.expr[!duplicated(GRNnew.expr$peak.ID),]
GRNnew.expr_cc$peak_gene.distance = as.numeric(GRNnew.expr_cc$peak_gene.distance)
summary(GRNnew.expr_cc$peak_gene.distance)

nrow(subset(GRNnew.expr_cc, peak_gene.distance < 500000))

###
## Add significant TFs ====
###

#new run of diffTF is here: /g/scb2/zaugg/bunina/LSD1/ATACseq/output/extension50
TFresults = read.table(file = file.path(dir1, "bunina/LSD1/ATACseq/output/extension50/LSD1_2019.summary.tsv"), 
                       header = T, sep = "\t")
hocoIDs = read.table(file = file.path(dir1, "bunina/LSD1/ATACseq/data", "HOCOTFID2ENSEMBL.txt"), header = T)
TFresults$ensembl = hocoIDs$ENSEMBL[match(TFresults$TF, hocoIDs$HOCOID)]
TFresults$type = ifelse(TFresults$weighted_meanDifference < 0, "UP_kids", "UP_dads")

# Filter out TFs not expressed in the fibroblasts:
TFresults$exprFib = ifelse(TFresults$ensembl %in% res_fibOrdered_df$ensembl, TRUE, FALSE) #see below for the object
table(TFresults$exprFib)
write.table(TFresults, file = file.path(dir1, "bunina/LSD1/ATACseq/output/extension50/LSD1_2019.summary_ensembl.tsv"), sep = "\t", quote = F)
TFresults = read.table(file.path(dir1, "bunina/LSD1/ATACseq/output/extension50/LSD1_2019.summary_ensembl.tsv"))

TFresults_sig = TFresults[TFresults$Cohend_factor != 1 & TFresults$pvalueAdj < 0.05 & TFresults$exprFib == TRUE,]

#now add this info to the GRN table:
GRNtable_sig$sigDiffTF = ifelse(GRNtable_sig$TF %in% TFresults_sig$TF, TRUE, FALSE)
#for stats check here:
GRNtable_stats_wide$sigDiffTF = ifelse(GRNtable_stats_wide$TF %in% TFresults_sig$TF, TRUE, FALSE)
table(GRNtable_stats_wide$sigDiffTF)

#is GRN enriched for diffTFs?
length(GRNtable_stats_wide$TF[GRNtable_stats_wide$TF %in% TFresults$TF])
fisherDF = data.frame(sigTF = c(59,26), nonSig = c(339,215)) #for r>0.4 sig GRN
# fisherDF = data.frame(sigTF = c(52,33), nonSig = c(293,260)) #for r>0.6 sig GRN
fisher.test(fisherDF) #nope!


# sigGRN stats ------------------------------------------------------------

#first drop empty factor levels (TFs):
GRNtable_sig$TF = factor(GRNtable_sig$TF)

#how many unique TFs in the network:
length(unique(GRNtable_sig$TF)) #393

#unique peaks:
length(unique(GRNtable_sig$peak)) #96050

#unique genes:
length(unique(GRNtable_sig$ENSEMBL)) #15128

#total TF-peak connections:
nrow(GRNtable_sig[!duplicated(GRNtable_sig[,c(1,5)]),]) #506564, total rows 1.546.958

#median number of peaks per TF:
grn_sigStats = as.data.frame(table(GRNtable_sig$TF[!duplicated(GRNtable_sig[,c("TF", "peak")])]))
median(grn_sigStats$Freq) #116
#kept old factor levels - remove these!
grn_sigStats = grn_sigStats[grn_sigStats$Freq > 0,]
ggplot(data = grn_sigStats, aes(x = log10(Freq) )) + geom_histogram(bins = 100) + theme_classic()

#median number of genes per TF:(genes can be multiple times! not unique!)
grn_sigStats$NgenesAll = as.data.frame(table(GRNtable_sig$TF))$Freq
median(grn_sigStats$NgenesAll) #274

grn_sigStats$NgenesUnique = as.data.frame(table(GRNtable_sig$TF[!duplicated(GRNtable_sig[,c(1,7)])]))$Freq
median(grn_sigStats$NgenesUnique) #246



###
## add fibr gene expression values to the table ====
###

res_fib10ordered = readRDS(file.path(dir1, "bunina/LSD1/RNAseq/myRNAseqFibr/output/Robjects", "res_fib10ordered.Rds"))
res_fibOrdered_df = as.data.frame(res_fib10ordered)
res_fibOrdered_df$ensembl = row.names(res_fibOrdered_df)
res_fibOrdered_df_red = cbind(res_fibOrdered_df[,c(7,2)], TF = NA, type = "all")
colnames(res_fibOrdered_df_red)[[1]] <- "ENSEMBL"

GRNtable_sig_fib = merge(GRNtable_sig, res_fibOrdered_df, by.x = "ENSEMBL", by.y = "ensembl", all.x = TRUE) #take all peaks per gene (R > 0.4)

###
### plot expression of all GRN genes in fibroblasts (see if not biased):
###

colnames(GRNtable_sig_LSDp_3tissues_sel_k27ac)
GRNtable_sig_LSDp_3tissues_sel_k27ac$log2FC_fibro = res_fibOrdered_df_red$log2FoldChange[match(GRNtable_sig_LSDp_3tissues_sel_k27ac$ENSEMBL, 
                                                                                               res_fibOrdered_df_red$ENSEMBL)]

pdf(file = file.path(dir1, "bunina/LSD1/GRN/output/plots/newGRN2020", "DensityBoxAllGRNgenesFibroExpr.pdf"), 
    useDingbats = FALSE, width = 4, height = 2)
GRNtable_sig_LSDp_3tissues_sel_k27ac %>%
  distinct(ENSEMBL, .keep_all = TRUE) %>%
  ggplot(aes(x = log2FC_fibro)) + geom_density() + theme_classic() + xlab("Log2FC in LSD1mut/WT, all GRN genes") +
  geom_boxplot(aes(x = log2FC_fibro, y = 0.5)) + geom_vline(xintercept = 0, linetype = 2)
dev.off()

###
# calculate log2FC of TF targets vs. all expressed genes, plot per-TF volcano plot:
## TF regulons in fibroblasts: ====
###
df_red_fib = GRNtable_sig_fib[,c("ENSEMBL", "TF", "sigDiffTF", "log2FoldChange", "padj")]
# df_red_fib = GRNtable_sig_LSDp_3tissues_sel_k27ac[,c("ENSEMBL", "TF", "sigDiffTF", "log2FoldChange", "padj")]
df_red_fib = df_red_fib[!duplicated(df_red_fib[,c(1,2)]) & !is.na(df_red_fib$log2FoldChange),]

#how many are diffTFs:
length(unique(df_red_fib$TF[df_red_fib$sigDiffTF == TRUE])) #58 - only sigDiffTFs expressed in fibroblasts

df_red_fib = GRNtable_sig_LSDp_3tissues_sel_k27ac %>%
  dplyr::select(c("ENSEMBL", "TF", "sigDiffTF", "log2FoldChange")) %>%
  distinct(ENSEMBL, TF, .keep_all = TRUE) %>%
  drop_na(log2FoldChange)


df_red_stats_fib = as.data.frame(table(df_red_fib$TF))
df_red_stats_fib$test.pval = NA
df_red_stats_fib$meanDiff = NA
df_red_stats_fib = df_red_stats_fib[df_red_stats_fib$Freq > 0,]
which(duplicated(df_red_stats_fib$Var1)) #none, good

for (i in as.character(df_red_stats_fib$Var1)) {
  sel_df = df_red_fib[df_red_fib$TF == i,c(1,2,4)]
  sel_df = sel_df[!is.na(sel_df[,3]),]
  if (nrow(sel_df) > 2) {
    sel_df$type = "targets"
    df_merge_fib = rbind(sel_df[sel_df$TF == i,], res_fibOrdered_df_red)
    # df_merge_fib = df_merge_fib[!is.na(df_merge_fib$log2FC_kidsVsDads),]
    df_red_stats_fib[df_red_stats_fib$Var1 == i,]$test.pval = t.test(df_merge_fib$log2FoldChange[df_merge_fib$type == "targets"], 
                                                                            df_merge_fib$log2FoldChange[df_merge_fib$type == "all"], 
                                                                       var.equal = FALSE)$p.value
    df_red_stats_fib[df_red_stats_fib$Var1 == i,]$meanDiff = mean(df_merge_fib$log2FoldChange[df_merge_fib$type == "targets"], na.rm = TRUE) - 
      mean(df_merge_fib$log2FoldChange[df_merge_fib$type == "all"], na.rm = TRUE)
  }
  
}

df_red_stats_fib$p.adjust = p.adjust(df_red_stats_fib$test.pval, method = "BH")
df_red_stats_fib$sig = ifelse(df_red_stats_fib$p.adjust < 0.05, TRUE, FALSE)
table(df_red_stats_fib$sig)
df_red_stats_fib[df_red_stats_fib$Var1 == "ZEB1",]

#add ensembl TF names for that first:
hocoIDs = read.table(file = file.path(dir1, "bunina/LSD1/ATACseq/data", "HOCOTFID2ENSEMBL.txt"), header = T)
df_red_stats_fib$TFname = hocoIDs$ENSEMBL[match(df_red_stats_fib$Var1, hocoIDs$HOCOID)]

#the line below keeps only those cases where TF targets are expressed in fibroblasts and TFs have > 2 targets, which makes sense to do:
df_red_stats_fib = df_red_stats_fib[complete.cases(df_red_stats_fib$test.pval),] #it might have removed some of the significantly diffTFs, as there are now only 43 of them left out of 58

###
#add TF expression itself in addition to targets:
###
df_red_stats_fib$TFexpr = res_fibOrdered_df$log2FoldChange[match(df_red_stats_fib$TFname, res_fibOrdered_df$ensembl)]
df_red_stats_fib$TFexpr_padj = res_fibOrdered_df$padj[match(df_red_stats_fib$TFname, res_fibOrdered_df$ensembl)]
df_red_stats_fib$TFexpr_sig = ifelse(df_red_stats_fib$TFexpr_padj <= 0.05, TRUE, FALSE)

#some TFs are not expressed in fibroblasts, for the stats below makes sense to remove them:
df_red_stats_fib = df_red_stats_fib[complete.cases(df_red_stats_fib$TFexpr),]

#save the old fib object and rename the new one to fit the code in the stats chunk:
# df_red_stats_fib_old = df_red_stats_fib
# df_red_stats_fib = df_red_stats_fib1

###
### plot overall violin plots for an example tfs targets vs. all genes changes: (take zeb1)
###

df_merge_fib_zeb1 = rbind( cbind(df_red_fib[df_red_fib$TF == "ZEB1",c(1,2,4)], type = "target"), res_fibOrdered_df_red)
colnames(df_merge_fib_zeb1)
df_merge_fib_zeb1$type = factor(df_merge_fib_zeb1$type, levels = c("all", "target"))

give.n <- function(x){
  return(c(y = -2, label = length(x))) 
  # experiment with the multiplier to find the perfect position: median(x)*1.2
}

pdf(file = file.path(dir1, "bunina/LSD1/GRN/output/plots/newGRN2020", "ViolinZEB1_log2FC_fibrobl_minDistGene.pdf"), 
    width = 3, height = 4, useDingbats = FALSE)
ggplot(data = df_merge_fib_zeb1, aes(x = type, y = log2FoldChange, fill = type)) + 
  geom_violin(draw_quantiles = 0.5) + 
  theme_classic(base_size = 9) + ggtitle(label = "Zeb1 targets log2FC in fibroblasts") + 
  stat_summary(fun.data = give.n, geom = "text", fun = median, 
               position = position_dodge(width = 0.75))
dev.off()

pdf(file = file.path(dir1, "bunina/LSD1/GRN/output/plots/newGRN2020", "DensityZEB1_log2FC_fibrobl_sigGRN.pdf"), 
    width = 4, height = 1.5, useDingbats = FALSE)
ggplot(data = df_merge_fib_zeb1, aes(x = log2FoldChange, fill = type)) + 
  geom_density(alpha = 0.4) + 
  geom_vline(xintercept = mean(df_merge_fib_zeb1$log2FoldChange[df_merge_fib_zeb1$type == "all"]), linetype = 2, size = 0.4, color = "orange") +
  geom_vline(xintercept = mean(df_merge_fib_zeb1$log2FoldChange[df_merge_fib_zeb1$type == "target"]), linetype = 2, size = 0.4, color = "blue") +
  theme_classic(base_size = 9) + ggtitle(label = "Zeb1 targets log2FC in fibroblasts")
dev.off()


# Fibroblast GRN stats ----------------------------------------------------

###
#the same as above but with removed ageing genes:
# (originals are in the script "LSD1fibrobl_myRNAseq.Rmd)
###

ageGenes_sig = readRDS(file.path(dir1, "bunina/LSD1/RNAseq/myRNAseqFibr/output", "ageGenes_sig_pval.Rds"))
ageGenes300m = readRDS(file.path(dir1, "bunina/LSD1/RNAseq/myRNAseqFibr/output", "ageGenes300m.Rds"))
ageGenes_sig$TFhoco = hocoIDs$HOCOID[match(ageGenes_sig$ensembl_gene_id, hocoIDs$ENSEMBL)]
ageGenes300m$TFhoco = hocoIDs$HOCOID[match(ageGenes300m$ensembl_gene_id, hocoIDs$ENSEMBL)]

df_red_stats_fib$ageingG = ifelse(df_red_stats_fib$TFname %in% ageGenes_sig$ensembl_gene_id | 
                                    df_red_stats_fib$TFname %in% ageGenes300m$ensembl_gene_id, TRUE, FALSE)
df_red_stats_fib[,12] = factor(df_red_stats_fib[,12],levels = c(TRUE, FALSE) ) #lapply(df_red_stats_fib[,c(5:6,8)], factor, levels = c(TRUE, FALSE))
table(df_red_stats_fib[,c(11,12)])
table(df_red_stats_fib[,c(6,12)])

df_red_stats_fib_noAge = df_red_stats_fib[df_red_stats_fib$ageingG == FALSE,]

fisher.test(table(df_red_stats_fib_noAge[,c(6,11)])) #nope! oddsR 1.34, pval 0.53

library(ggrepel)
pdf(file = file.path(dir1, "bunina/LSD1/GRN/output/plots/newGRN2020", "PerTFtest_meanDiff_fibrobl_GRNtableSig.pdf"), 
    width = 5, height = 4, useDingbats = FALSE)
ggplot(data = df_red_stats_fib_noAge, aes(x = Freq, y = -log10(p.adjust))) + geom_point(size = 0.4, aes(colour = sigDiffTF)) + 
  theme_classic() + geom_hline(yintercept = -log10(0.05), color = "red", linetype = 2) +
  geom_text_repel(data = subset(df_red_stats_fib_noAge, p.adjust < 0.05 & sigDiffTF == TRUE), aes(label = Var1), 
                  box.padding = 0.35, size = 2, segment.size = 0.2)
ggplot(data = df_red_stats_fib_noAge, aes(x = meanDiff, y = -log10(p.adjust))) + geom_point(size = 0.4, aes(colour = sigDiffTF)) + 
  theme_classic() + geom_hline(yintercept = -log10(0.05), color = "red", linetype = 2) + geom_vline(xintercept = 0, linetype = 2, size = 0.3) +
  geom_text_repel(data = subset(df_red_stats_fib_noAge, p.adjust < 0.05 & sigDiffTF == TRUE), aes(label = Var1), 
                  box.padding = 0.35, size = 3, segment.size = 0.2, segment.alpha = 0.3)
#now plot TF expression on Y-axis, and targets LFC on the x as before:
ggplot(data = df_red_stats_fib_noAge, aes(x = meanDiff, y = TFexpr)) + geom_point(size = 0.4, aes(colour = sigDiffTF)) + 
  theme_classic() + geom_hline(yintercept = 0, size = 0.3, linetype = 2) + geom_vline(xintercept = 0, linetype = 2, size = 0.3) +
  geom_text_repel(data = subset(df_red_stats_fib_noAge, p.adjust < 0.05 & sigDiffTF == TRUE), aes(label = Var1), 
                  box.padding = 0.35, size = 2, segment.size = 0.2, segment.alpha = 0.3)
#color by p.adj:
ggplot(data = df_red_stats_fib_noAge, aes(x = meanDiff, y = TFexpr)) + geom_point(size = 0.4, aes(colour = p.adjust)) + 
  theme_classic() + geom_hline(yintercept = 0, size = 0.3, linetype = 2) + geom_vline(xintercept = 0, linetype = 2, size = 0.3) +
  geom_text_repel(data = subset(df_red_stats_fib, p.adjust < 0.05 & sigDiffTF == TRUE), aes(label = Var1), 
                  box.padding = 0.35, size = 3, segment.size = 0.2, segment.alpha = 0.3)
dev.off()

#same but without labels:

pdf(file = file.path(dir1, "bunina/LSD1/GRN/output/plots/newGRN2020", "PerTFtest_meanDiff_fibrobl_GRNtableSig_NoLabels.pdf"), 
    width = 4, height = 3.5, useDingbats = FALSE)
ggplot(data = df_red_stats_fib_noAge, aes(x = meanDiff, y = -log10(p.adjust))) + geom_point(size = 0.2, aes(colour = sig)) + 
  theme_classic() + geom_hline(yintercept = -log10(0.05), linetype = 2) + geom_vline(xintercept = 0, linetype = 2, size = 0.3)
dev.off()

#with labels, color by padj:
pdf(file = file.path(dir1, "bunina/LSD1/GRN/output/plots/newGRN2020", "PerTFtest_meanDiff_fibrobl_GRNtableSigColorPadj.pdf"), 
    width = 4, height = 3.5, useDingbats = FALSE)
ggplot(data = df_red_stats_fib_noAge, aes(x = meanDiff, y = -log10(p.adjust))) + geom_point(size = 0.6, aes(colour = sig, shape = sigDiffTF)) + 
  theme_classic() + geom_hline(yintercept = -log10(0.05), linetype = 2) + geom_vline(xintercept = 0, linetype = 2, size = 0.5) + 
  geom_text_repel(data = subset(df_red_stats_fib_noAge, p.adjust < 0.05 & sigDiffTF == TRUE & ageingG == FALSE), aes(label = Var1), 
                  box.padding = 0.4, size = 2.5, segment.size = 0.2, segment.alpha = 0.3) + scale_color_manual(values = c("#377EB8", "darkgrey"))
ggplot(data = df_red_stats_fib_noAge, aes(x = meanDiff, y = -log10(p.adjust))) + geom_point(size = 0.6, aes(colour = sig)) + 
  theme_classic() + geom_hline(yintercept = -log10(0.05), linetype = 2) + geom_vline(xintercept = 0, linetype = 2, size = 0.5) + 
  geom_text_repel(data = subset(df_red_stats_fib_noAge, p.adjust < 0.05 & sigDiffTF == TRUE & ageingG == FALSE), aes(label = Var1), 
                  box.padding = 0.4, size = 2.5, segment.size = 0.2, segment.alpha = 0.3) + scale_color_manual(values = c("#377EB8", "darkgrey"))
ggplot(data = df_red_stats_fib_noAge, aes(x = meanDiff, y = -log10(p.adjust))) + geom_point(size = 0.6, aes(colour = sig, shape = sigDiffTF)) + 
  theme_classic() + geom_hline(yintercept = -log10(0.05), linetype = 2) + geom_vline(xintercept = 0, linetype = 2, size = 0.5) + 
  geom_text_repel(data = subset(df_red_stats_fib_noAge, p.adjust < 0.05 & sigDiffTF == TRUE & ageingG == FALSE), aes(label = Var1), 
                  box.padding = 0.4, size = 2.5, segment.size = 0.2, segment.alpha = 0.3) 
dev.off()

### Final new figure 2 scatter fibro:
pdf(file = file.path(dir1, "bunina/LSD1/GRN/output/plots/newGRN2020", "FibroMeanDiff_finalFigureScatter.pdf"), 
    width = 4, height = 3.5, useDingbats = FALSE)
ggplot(data = df_red_stats_fib_noAge, aes(x = meanDiff, y = -log10(p.adjust))) + geom_point(size = 0.6, aes(colour = sig)) + 
  theme_classic() + geom_hline(yintercept = -log10(0.05), linetype = 2) + geom_vline(xintercept = 0, linetype = 2, size = 0.5) + 
  geom_text_repel(data = subset(df_red_stats_fib_noAge, p.adjust < 0.05 & sigDiffTF == TRUE), aes(label = Var1), 
                  box.padding = 0.4, size = 2.5, segment.size = 0.2, segment.alpha = 0.3) + scale_color_manual(values = c("#377EB8", "darkgrey"))
dev.off()

###
### as boxplots:
###
pdf(file = file.path(dir1, "bunina/LSD1/GRN/output/plots/newGRN2020", "BoxplMeanDiff_fibroblColorSig.pdf"), 
    width = 2.5, height = 4, useDingbats = FALSE)
ggplot(data = df_red_stats_fib, aes(x = sig, y = meanDiff)) + geom_boxplot(aes(fill = sig), alpha = 0.4) + 
  theme_classic() + geom_hline(yintercept = 0, linetype = 2)
dev.off()


###
### is diffTF and GRN mean diff correlated?
### diffTF score is opposite sign of GRN! reverse here
df_red_stats_fib$diffTF_score = -TFresults$weighted_meanDifference[match(df_red_stats_fib$Var1, TFresults$TF)]
df_red_stats_fib$sig_both = ifelse(df_red_stats_fib$sig == TRUE & df_red_stats_fib$sigDiffTF == TRUE, TRUE, FALSE)
df_red_stats_fib$sig_both = factor(df_red_stats_fib$sig_both, levels = c(TRUE, FALSE))

pdf(file = file.path(dir1, "bunina/LSD1/GRN/output/plots/newGRN2020", "DiffTFvsGRNtargetsExprScatter.pdf"), 
    width = 6, height = 5, useDingbats = FALSE)
df_red_stats_fib %>%
  subset(ageingG == FALSE) %>%
  ggplot(aes(x = meanDiff, y = diffTF_score, color = sig_both)) + geom_point(size = 0.5) + 
  theme_classic() + geom_vline(xintercept = 0, linetype = 2) + geom_hline(yintercept = 0, linetype = 2) + 
  geom_text_repel(data = subset(df_red_stats_fib, sig_both == TRUE & ageingG == FALSE), aes(label = Var1), 
                  segment.size = 0.2, size = 3, color = "red", segment.alpha = 0.2) + geom_smooth(method = "lm", linetype = 2, size = 0.5, alpha = 0.5) +
  ylab("TF activity, LSD1mut/WT") + xlab("TF target mean expression, LSD1mut/WT")
# all sig targets in fib:
df_red_stats_fib %>%
  subset(ageingG == FALSE) %>%
  ggplot(aes(x = meanDiff, y = diffTF_score, color = sig)) + geom_point(size = 0.5) + 
  theme_classic() + geom_vline(xintercept = 0, linetype = 2) + geom_hline(yintercept = 0, linetype = 2) + 
  geom_text_repel(data = subset(df_red_stats_fib, sig_both == TRUE & ageingG == FALSE), aes(label = Var1), 
                  segment.size = 0.2, size = 3, color = "red", segment.alpha = 0.2) + geom_smooth(method = "lm", linetype = 2, size = 0.5, alpha = 0.5) +
  ylab("TF activity, LSD1mut/WT") + xlab("TF target mean expression, LSD1mut/WT")
dev.off()

cor(df_red_stats_fib$meanDiff[df_red_stats_fib$ageingG == FALSE], df_red_stats_fib$diffTF_score[df_red_stats_fib$ageingG == FALSE]) #-0.04
cor(df_red_stats_fib$meanDiff[df_red_stats_fib$sig_both == TRUE & df_red_stats_fib$ageingG == FALSE], 
    df_red_stats_fib$diffTF_score[df_red_stats_fib$sig_both == TRUE & df_red_stats_fib$ageingG == FALSE]) #0.53
cor(df_red_stats_fib$meanDiff[df_red_stats_fib$sig == TRUE & df_red_stats_fib$ageingG == FALSE], 
    df_red_stats_fib$diffTF_score[df_red_stats_fib$sig == TRUE & df_red_stats_fib$ageingG == FALSE]) #0.04

table(df_red_stats_fib[df_red_stats_fib$ageingG == FALSE,c("sig_both", "TFexpr_sig")])

df_red_stats_fib %>%
  subset(ageingG == FALSE & sig_both == TRUE) %>%
  select(Var1, TFexpr_sig)




# GRN overlap with iPSC inhibitor/KO hits -------------------------------------

#add target gene info: 
all_sig_genes_df_full = read.csv(file.path(dir1, "bunina/LSD1/RNAseq/iPSCs/data/geneClusters", "lsd_ipscAllSigGenes.csv"))
all_sig_genes_df_full_05 = all_sig_genes_df_full[all_sig_genes_df_full$geneFDR == 5,]

#add to the table: (taken from below, next time re-run to have this table first maybe)
GRNtable_sig_LSDp$giuseppe_hits = ifelse(GRNtable_sig_LSDp$ENSEMBL %in% all_sig_genes_df_full_05$ensembl, TRUE, FALSE)
GRNtable_sig_LSDp[,c("giuseppe_hits")] <- lapply(GRNtable_sig_LSDp[,c("giuseppe_hits")], factor, levels = c(TRUE, FALSE))


# Add iPSCs gene expression data ------------------------------------------

dds_ipsc_mod = readRDS(file.path(dir1, "bunina/LSD1/RNAseq/iPSCs/output/Robjects",
                                 "dds_iPSCsDadsKidsRNAseq_noPro3.Rds"))
dds_ipsc_mod <- dds_ipsc_mod[ rowSums(counts(dds_ipsc_mod)) > 10, ]
res_ipsc_mod = results(dds_ipsc_mod)
summary(res_ipsc_mod)
res_ipsc_df = as.data.frame(res_ipsc_mod)
res_ipsc_df$ensembl = sapply(strsplit(row.names(res_ipsc_df), "[.]"), "[", 1)

#modify a bit for later:
res_ipsc_df_red = res_ipsc_df[,c(7,2)]
res_ipsc_df_red$TF = NA
res_ipsc_df_red$type = "all"
colnames(res_ipsc_df_red)[c(1,2)] <- c("ENSEMBL", "log2FC_ipscs")

#now merge with the network table:
GRNtable_sig_LSDp_3tissues_ipsc = GRNtable_sig_LSDp_3tissues
GRNtable_sig_LSDp_3tissues_ipsc$log2FC_ipscs = res_ipsc_df$log2FoldChange[match(GRNtable_sig_LSDp_3tissues_ipsc$ENSEMBL, 
                                                                                res_ipsc_df$ensembl)]
df_red_ipsc = GRNtable_sig_LSDp_3tissues_ipsc[,c("ENSEMBL", "TF", "sigDiffTF", "log2FC_ipscs")]
df_red_ipsc = df_red_ipsc[!duplicated(df_red_ipsc[,c(1,2)]) & !is.na(df_red_ipsc$log2FC_ipscs),]

df_red_ipsc = GRNtable_sig_LSDp_3tissues_sel_k27ac %>%
  dplyr::select(c("ENSEMBL", "TF", "sigDiffTF", "log2FC_ipscs")) %>%
  distinct(ENSEMBL, TF, .keep_all = TRUE) %>%
  drop_na(log2FC_ipscs)

### density plots TFtargets meanDiff fibroblasts+iPSCs ====


#how many are diffTFs:
length(unique(df_red_ipsc$TF[df_red_ipsc$sigDiffTF == TRUE])) #58 - only sigDiffTFs expressed in fibroblasts

df_red_stats_ipsc = as.data.frame(table(df_red_ipsc$TF))
df_red_stats_ipsc$test.pval = NA
df_red_stats_ipsc$meanDiff = NA
df_red_stats_ipsc = df_red_stats_ipsc[df_red_stats_ipsc$Freq > 0,]
which(duplicated(df_red_stats_ipsc$Var1)) #none, good

for (i in as.character(df_red_stats_ipsc$Var1)) {
  sel_df = df_red_ipsc[df_red_ipsc$TF == i,c(1,2,4)]
  sel_df = sel_df[!is.na(sel_df[,3]),]
  if (nrow(sel_df) > 2) {
    sel_df$type = "targets"
    df_merge_ipsc = rbind(sel_df[sel_df$TF == i,], res_ipsc_df_red)
    df_red_stats_ipsc[df_red_stats_ipsc$Var1 == i,]$test.pval = t.test(df_merge_ipsc$log2FC_ipscs[df_merge_ipsc$type == "targets"], 
                                                                     df_merge_ipsc$log2FC_ipscs[df_merge_ipsc$type == "all"], 
                                                                     var.equal = FALSE)$p.value
    df_red_stats_ipsc[df_red_stats_ipsc$Var1 == i,]$meanDiff = mean(df_merge_ipsc$log2FC_ipscs[df_merge_ipsc$type == "targets"], na.rm = TRUE) - 
      mean(df_merge_ipsc$log2FC_ipscs[df_merge_ipsc$type == "all"], na.rm = TRUE)
  }
  
}

df_red_stats_ipsc$p.adjust = p.adjust(df_red_stats_ipsc$test.pval, method = "BH")
df_red_stats_ipsc$sig = ifelse(df_red_stats_ipsc$p.adjust < 0.05, TRUE, FALSE)
df_red_stats_ipsc$sig[is.na(df_red_stats_ipsc$sig)] <- FALSE
table(df_red_stats_ipsc$sig)
df_red_stats_ipsc$sigDiffTF = ifelse(df_red_stats_ipsc$Var1 %in% TFresults_sig$TF, TRUE, FALSE)
df_red_stats_ipsc$sigDiffTF = factor(df_red_stats_ipsc$sigDiffTF, levels = c(TRUE, FALSE))



library(ggrepel)
pdf(file = file.path(dir1, "bunina/LSD1/GRN/output/plots/newGRN2020", "PerTFtest_meanDiff_iPSCdata.pdf"), 
    width = 4, height = 3.5, useDingbats = FALSE)
ggplot(data = df_red_stats_ipsc, aes(x = meanDiff, y = -log10(p.adjust))) + geom_point(size = 0.4, aes(colour = sigDiffTF)) + 
  theme_classic() + geom_hline(yintercept = -log10(0.05), color = "red", linetype = 2) + geom_vline(xintercept = 0, linetype = 2, size = 0.3) +
  geom_text_repel(data = subset(df_red_stats_ipsc, p.adjust < 0.05 & sigDiffTF == TRUE), aes(label = Var1), 
                  box.padding = 0.35, size = 2, segment.size = 0.2, segment.alpha = 0.3)
ggplot(data = df_red_stats_ipsc, aes(x = Freq, y = -log10(p.adjust))) + geom_point(size = 0.4, aes(colour = sigDiffTF)) + 
  theme_classic() + geom_hline(yintercept = -log10(0.05), color = "red", linetype = 2) +
  geom_text_repel(data = subset(df_red_stats_ipsc, p.adjust < 0.05 & sigDiffTF == TRUE), aes(label = Var1), 
                  box.padding = 0.35, size = 2, segment.size = 0.2)

dev.off()

#plot as density:

pdf(file = file.path(dir1, "bunina/LSD1/GRN/output/plots/newGRN2020", "DensityMeanDiff_fibro_iPSCdata.pdf"),
    width = 4, height = 2, useDingbats = FALSE) #using non-LSD1filtered data
# pdf(file = file.path(dir1, "bunina/LSD1/GRN/output/plots/newGRN2020", "DensityMeanDiff_fibro_iPSCdata_LSD1pFilter.pdf"), 
#     width = 4, height = 2, useDingbats = FALSE) #with re-done LSD1-filtered data:
df_red_stats_ipsc$sig = factor(df_red_stats_ipsc$sig, levels = c(TRUE, FALSE))
ggplot(data = df_red_stats_ipsc, aes(x = meanDiff)) + geom_density(aes(fill = sig), alpha = 0.4) + 
  theme_classic() + geom_vline(xintercept = 0, linetype = 2, size = 0.3)
#same for fibroblasts data from above:
df_red_stats_fib_noAge$sig = factor(df_red_stats_fib_noAge$sig, levels = c(TRUE, FALSE))
ggplot(data = df_red_stats_fib_noAge, aes(x = meanDiff)) + geom_density(aes(fill = sig), alpha = 0.4) + 
  theme_classic() + geom_vline(xintercept = 0, linetype = 2, size = 0.3)
dev.off()





# 28.12.20 VPA endo data --------------------------------------------------

### add VPA mut/ctrl expression changes to the network: if the negative shift reduced?

### first add it with unique name to the main GRN table and save for later:
GRNtable_sig_LSDp_3tissues_sel_k27ac$endoCTRL = res_CTRLmut_wt$log2FoldChange[match(GRNtable_sig_LSDp_3tissues_sel_k27ac$ENSEMBL,
                                                                                    res_CTRLmut_wt$ensembl)]
GRNtable_sig_LSDp_3tissues_sel_k27ac$endoVPA = res_VPAmut_wt$log2FoldChange[match(GRNtable_sig_LSDp_3tissues_sel_k27ac$ENSEMBL,
                                                                                  res_VPAmut_wt$ensembl)]
saveRDS(GRNtable_sig_LSDp_3tissues_sel_k27ac, file = file.path(dir1, "bunina/LSD1/GRN/output/Robjects", "GRNtable_sig_LSDp_3tissues_wVPAdata.rds"))
GRNtable_sig_LSDp_3tissues_sel_k27ac = readRDS(file.path(dir1, "bunina/LSD1/GRN/output/Robjects", "GRNtable_sig_LSDp_3tissues_wVPAdata.rds"))


###
### now add it with mut/wt column name as with smarca2 etc. mut data to make densty plots:
###

res_VPAmut_wt = readRDS(file.path(dir1, "bunina/LSD1/VPAtreatment/output", "res_VPAmut_wt.rds"))
res_VPAmut_wt = as.data.frame(res_VPAmut_wt)
res_VPAmut_wt$ensembl = sapply(strsplit(row.names(res_VPAmut_wt), "[.]"), "[", 1)

res_VPA_df_red = res_VPAmut_wt[,c("ensembl", "log2FoldChange")]
res_VPA_df_red$TF = NA
res_VPA_df_red$type = "all"
colnames(res_VPA_df_red)[c(1,2)] <- c("ENSEMBL", "mut_vs_wt_log" )

#now merge with the network table:
GRNtable_sig_LSDp_3tissues_sel_k27ac$mut_vs_wt_log = res_VPAmut_wt$log2FoldChange[match(GRNtable_sig_LSDp_3tissues_sel_k27ac$ENSEMBL,
                                                                                       res_VPAmut_wt$ensembl)]

df_red_VPA = GRNtable_sig_LSDp_3tissues_sel_k27ac[,c("ENSEMBL", "TF", "sigDiffTF", "mut_vs_wt_log")]
df_red_VPA = df_red_VPA[!duplicated(df_red_VPA[,c(1,2)]) & !is.na(df_red_VPA$mut_vs_wt_log),]


#how many are diffTFs:
length(unique(df_red_VPA$TF[df_red_VPA$sigDiffTF == TRUE])) #58 - only sigDiffTFs expressed in fibroblasts

df_red_stats_VPA = as.data.frame(table(df_red_VPA$TF))
df_red_stats_VPA$test.pval = NA
df_red_stats_VPA$meanDiff = NA
df_red_stats_VPA = df_red_stats_VPA[df_red_stats_VPA$Freq > 0,]
which(duplicated(df_red_stats_VPA$Var1)) #none, good

for (i in as.character(df_red_stats_VPA$Var1)) {
  sel_df = df_red_VPA[df_red_VPA$TF == i,c(1,2,4)]
  sel_df = sel_df[!is.na(sel_df[,3]),]
  if (nrow(sel_df) > 2) {
    sel_df$type = "targets"
    df_merge_VPA = rbind(sel_df[sel_df$TF == i,], res_VPA_df_red)
    df_red_stats_VPA[df_red_stats_VPA$Var1 == i,]$test.pval = wilcox.test(df_merge_VPA$mut_vs_wt_log[df_merge_VPA$type == "targets"],
                                                                              df_merge_VPA$mut_vs_wt_log[df_merge_VPA$type == "all"])$p.value
    df_red_stats_VPA[df_red_stats_VPA$Var1 == i,]$meanDiff = mean(df_merge_VPA$mut_vs_wt_log[df_merge_VPA$type == "targets"], na.rm = TRUE) - 
      mean(df_merge_VPA$mut_vs_wt_log[df_merge_VPA$type == "all"], na.rm = TRUE, trim = 0.1)
  }
  
}

df_red_stats_VPA = df_red_stats_VPA[!is.na(df_red_stats_VPA$meanDiff),]
df_red_stats_VPA$p.adjust = p.adjust(df_red_stats_VPA$test.pval, method = "BH")
df_red_stats_VPA$sig = ifelse(df_red_stats_VPA$p.adjust < 0.05, TRUE, FALSE)
df_red_stats_VPA$sig[is.na(df_red_stats_VPA$sig)] <- FALSE
table(df_red_stats_VPA$sig)
df_red_stats_VPA$type = "VPA"

###
### same for endo CTRL differentiation for consistency:
###

res_CTRLmut_wt = readRDS(file = file.path(dir1, "bunina/LSD1/VPAtreatment/output", "res_CTRLmut_wt.rds"))
res_CTRLmut_wt = as.data.frame(res_CTRLmut_wt)
res_CTRLmut_wt$ensembl = sapply(strsplit(row.names(res_CTRLmut_wt), "[.]"), "[", 1)

res_CTRL_df_red = res_CTRLmut_wt[,c("ensembl", "log2FoldChange")]
res_CTRL_df_red$TF = NA
res_CTRL_df_red$type = "all"
colnames(res_CTRL_df_red)[c(1,2)] <- c("ENSEMBL", "mut_vs_wt_log" )

#now merge with the network table:
GRNtable_sig_LSDp_3tissues_sel_k27ac$mut_vs_wt_log = res_CTRLmut_wt$log2FoldChange[match(GRNtable_sig_LSDp_3tissues_sel_k27ac$ENSEMBL,
                                                                                        res_CTRLmut_wt$ensembl)]

df_red_CTRL = GRNtable_sig_LSDp_3tissues_sel_k27ac[,c("ENSEMBL", "TF", "sigDiffTF", "mut_vs_wt_log")]
df_red_CTRL = df_red_CTRL[!duplicated(df_red_CTRL[,c(1,2)]) & !is.na(df_red_CTRL$mut_vs_wt_log),]


#how many are diffTFs:
length(unique(df_red_CTRL$TF[df_red_CTRL$sigDiffTF == TRUE])) #58 - only sigDiffTFs expressed in fibroblasts

df_red_stats_CTRL = as.data.frame(table(df_red_CTRL$TF))
df_red_stats_CTRL$test.pval = NA
df_red_stats_CTRL$meanDiff = NA
df_red_stats_CTRL = df_red_stats_CTRL[df_red_stats_CTRL$Freq > 0,]
which(duplicated(df_red_stats_CTRL$Var1)) #none, good

for (i in as.character(df_red_stats_CTRL$Var1)) {
  sel_df = df_red_CTRL[df_red_CTRL$TF == i,c(1,2,4)]
  sel_df = sel_df[!is.na(sel_df[,3]),]
  if (nrow(sel_df) > 2) {
    sel_df$type = "targets"
    df_merge_CTRL = rbind(sel_df[sel_df$TF == i,], res_CTRL_df_red)
    df_red_stats_CTRL[df_red_stats_CTRL$Var1 == i,]$test.pval = wilcox.test(df_merge_CTRL$mut_vs_wt_log[df_merge_CTRL$type == "targets"],
                                                                          df_merge_CTRL$mut_vs_wt_log[df_merge_CTRL$type == "all"])$p.value
    df_red_stats_CTRL[df_red_stats_CTRL$Var1 == i,]$meanDiff = mean(df_merge_CTRL$mut_vs_wt_log[df_merge_CTRL$type == "targets"], na.rm = TRUE) - 
      mean(df_merge_CTRL$mut_vs_wt_log[df_merge_CTRL$type == "all"], na.rm = TRUE, trim = 0.1)
  }
  
}

df_red_stats_CTRL = df_red_stats_CTRL[!is.na(df_red_stats_CTRL$meanDiff),]
df_red_stats_CTRL$p.adjust = p.adjust(df_red_stats_CTRL$test.pval, method = "BH")
df_red_stats_CTRL$sig = ifelse(df_red_stats_CTRL$p.adjust < 0.05, TRUE, FALSE)
df_red_stats_CTRL$sig[is.na(df_red_stats_CTRL$sig)] <- FALSE
df_red_stats_CTRL$type = "CTRL"

###
### merge VPA and ctrl dfs:
###

df_red_stats_VPA_Ctrl = rbind(df_red_stats_VPA, df_red_stats_CTRL)
df_red_stats_VPA_Ctrl$sigDiffTF = ifelse(df_red_stats_VPA_Ctrl$Var1 %in% TFresults_sig$TF, TRUE, FALSE)
df_red_stats_VPA_Ctrl[,c(6,8)] = lapply(df_red_stats_VPA_Ctrl[,c(6,8)], factor, levels = c(TRUE, FALSE))
saveRDS(df_red_stats_VPA_Ctrl, file = file.path(dir1, "bunina/LSD1/VPAtreatment/output", "df_red_stats_VPA_Ctrl.rds"))

library(ggrepel)
pdf(file = file.path(dir1, "bunina/LSD1/GRN/output/plots/newGRN2020", "PerTFtest_meanDiff_VPAendo.pdf"), 
    width = 4, height = 3.5, useDingbats = FALSE)
ggplot(data = df_red_stats_VPA_Ctrl, aes(x = meanDiff, y = -log10(p.adjust))) + geom_point(size = 0.4, aes(colour = type)) + 
  theme_classic() + geom_hline(yintercept = -log10(0.05), color = "red", linetype = 2) + geom_vline(xintercept = 0, linetype = 2, size = 0.3) 
dev.off()

#label TFs:
ggplot(data = df_red_stats_VPA_Ctrl, aes(x = meanDiff, y = -log10(p.adjust), colour = type)) + geom_point(size = 0.4) + 
  theme_classic() + geom_hline(yintercept = -log10(0.05), color = "red", linetype = 2) + 
  geom_vline(xintercept = 0, linetype = 2, size = 0.3) + scale_color_manual(values = c("#999999", "#CC79A7")) +
  geom_text_repel(data = subset(df_red_stats_VPA_Ctrl, p.adjust < 0.01 & abs(meanDiff) > 0.1), aes(label = Var1), 
                box.padding = 0.35, size = 2.5, segment.alpha = 0.3)

#label only diffTFs:
ggplot(data = df_red_stats_VPA_Ctrl, aes(x = meanDiff, y = -log10(p.adjust), colour = type)) + geom_point(size = 0.4) + 
  theme_classic() + geom_hline(yintercept = -log10(0.05), color = "red", linetype = 2) + 
  geom_vline(xintercept = 0, linetype = 2, size = 0.3) + scale_color_manual(values = c("#999999", "#CC79A7")) +
  geom_text_repel(data = subset(df_red_stats_VPA_Ctrl, p.adjust < 0.01 & sigDiffTF == TRUE), aes(label = Var1), 
                  box.padding = 0.35, size = 3, segment.alpha = 0.3)

#plot as density:

#match sig TFs - sig in ctrl or vpa:
df_red_stats_VPA_Ctrl$sig_ctrl = ifelse(df_red_stats_VPA_Ctrl$Var1 %in% df_red_stats_VPA_Ctrl$Var1[df_red_stats_VPA_Ctrl$sig == TRUE & df_red_stats_VPA_Ctrl$type == "CTRL"], TRUE, FALSE)
table(df_red_stats_VPA_Ctrl[,c("sig_ctrl", "type")])

#add means by group:
mean_vpa_ctrl <- df_red_stats_VPA_Ctrl %>%
  group_by(sig_ctrl, type) %>%
  summarize(mean=mean(meanDiff))

pdf(file = file.path(dir1, "bunina/LSD1/GRN/output/plots/newGRN2020", "DensityMeanDiff_VPAendo.pdf"), 
    width = 5, height = 1.6, useDingbats = FALSE)
#sig in VPA/CTRL separately defined:
df_red_stats_VPA_Ctrl$sig = factor(df_red_stats_VPA_Ctrl$sig, levels = c(TRUE, FALSE))
ggplot(data = df_red_stats_VPA_Ctrl, aes(x = meanDiff)) + geom_density(aes(fill = type), alpha = 0.7) + 
  scale_fill_manual(values = c("#999999", "#CC79A7")) +
  theme_classic() + geom_vline(xintercept = 0, linetype = 2, size = 0.3) + facet_wrap(~sig, scales = "free")
#sig in CTRL defined, same labeled in VPA:
df_red_stats_VPA_Ctrl$sig_ctrl = factor(df_red_stats_VPA_Ctrl$sig_ctrl, levels = c(TRUE, FALSE))
ggplot(data = df_red_stats_VPA_Ctrl, aes(x = meanDiff)) + geom_density(aes(fill = type), alpha = 0.5) + 
  scale_fill_manual(values = c("#999999", "#CC79A7")) + 
  geom_vline(data = mean_vpa_ctrl, aes(xintercept = mean, color = type), size=0.6) +
  scale_color_manual(values = c("#999999", "#CC79A7")) +
  theme_classic() + geom_vline(xintercept = 0, linetype = 2, size = 0.3) + facet_wrap(~sig_ctrl, scales = "free")
#remove 2 TF outliers for visualization:
df_red_stats_VPA_Ctrl %>%
  subset(Var1 != "ATF7" & Var1 != "PO3F4") %>%
  ggplot(aes(x = meanDiff)) + geom_density(aes(fill = type), alpha = 0.7) + 
  scale_fill_manual(values = c("#999999", "#CC79A7")) + 
  theme_classic() + geom_vline(xintercept = 0, linetype = 2, size = 0.3) + facet_wrap(~sig_ctrl, scales = "free")
dev.off()

### same but vertical:
pdf(file = file.path(dir1, "bunina/LSD1/GRN/output/plots/newGRN2020", "DensityMeanDiff_VPAendoV.pdf"), 
    width = 3, height = 3, useDingbats = FALSE)
#sig in VPA/CTRL separately defined:
df_red_stats_VPA_Ctrl$sig = factor(df_red_stats_VPA_Ctrl$sig, levels = c(TRUE, FALSE))
ggplot(data = df_red_stats_VPA_Ctrl, aes(x = meanDiff)) + geom_density(aes(fill = type), alpha = 0.7) + 
  scale_fill_manual(values = c("#999999", "#CC79A7")) +
  theme_classic() + geom_vline(xintercept = 0, linetype = 2, size = 0.3) + facet_wrap(~sig, ncol = 1, scales = "free")
#sig in CTRL defined, same labeled in VPA:
df_red_stats_VPA_Ctrl$sig_ctrl = factor(df_red_stats_VPA_Ctrl$sig_ctrl, levels = c(TRUE, FALSE))
ggplot(data = df_red_stats_VPA_Ctrl, aes(x = meanDiff)) + geom_density(aes(fill = type), alpha = 0.5) + 
  scale_fill_manual(values = c("#999999", "#CC79A7")) + 
  geom_vline(data = mean_vpa_ctrl, aes(xintercept = mean, color = type), size=0.6) +
  scale_color_manual(values = c("#999999", "#CC79A7")) +
  theme_classic() + geom_vline(xintercept = 0, linetype = 2, size = 0.3) + facet_wrap(~sig_ctrl, ncol = 1, scales = "free")
#remove 2 TF outliers for visualization:
df_red_stats_VPA_Ctrl %>%
  subset(Var1 != "ATF7" & Var1 != "PO3F4") %>%
  ggplot(aes(x = meanDiff)) + geom_density(aes(fill = type), alpha = 0.7) + 
  scale_fill_manual(values = c("#999999", "#CC79A7")) + 
  theme_classic() + geom_vline(xintercept = 0, linetype = 2, size = 0.3) + facet_wrap(~sig_ctrl, ncol = 1, scales = "free")
dev.off()

### as dotplots:
df_red_stats_VPA_Ctrlmerge = merge(df_red_stats_VPA, df_red_stats_CTRL, by = "Var1", suffixes = c(".vpa", ".ctrl"))
df_red_stats_VPA_Ctrlmerge$sig.ctrl = factor(df_red_stats_VPA_Ctrlmerge$sig.ctrl, levels = c(TRUE, FALSE))

#exclude TFs that are differential themselve in LSD1mut:
df_red_stats_VPA_Ctrlmerge$TF_DEG = ifelse(df_red_stats_VPA_Ctrlmerge$Var1 %in% 
                                              diff_genes_3layers$name[diff_genes_3layers$significant.endo == TRUE], TRUE, FALSE) #if sig_any - 8 TFs more to filter out

df_red_stats_VPA_Ctrlmerge$maxDE_TF = diff_genes_3layers$maxDE[match(df_red_stats_VPA_Ctrlmerge$Var1, diff_genes_3layers$name)]
df_red_stats_VPA_Ctrlmerge = subset(df_red_stats_VPA_Ctrlmerge, TF_DEG == FALSE & abs(maxDE_TF) < 0.58)
# table(df_red_stats_VPA_Ctrlmerge$TF_DEG)
df_red_stats_VPA_Ctrlmerge$isTop30TF = ifelse(df_red_stats_VPA_Ctrlmerge$Var1 %in% subnetw_30TFs$TF, "top30", "other")
df_red_stats_VPA_Ctrlmerge$isTop10TF = ifelse(df_red_stats_VPA_Ctrlmerge$Var1 %in% subnetw_30TFs$TF[1:10], TRUE, FALSE)

ggplot(data = df_red_stats_VPA_Ctrlmerge, aes(x = meanDiff.ctrl, y = meanDiff.vpa)) + geom_point(size = 0.5) + 
  geom_vline(xintercept = 0, linetype = 2) + geom_hline(yintercept = 0, linetype = 2) + 
  geom_abline(slope = 1, intercept = 0) + 
  geom_text_repel(data = subset(df_red_stats_VPA_Ctrlmerge, meanDiff.ctrl < -0.3), aes(label = Var1), 
                  box.padding = 0.35, size = 2.5, segment.alpha = 0.3) +
  theme_classic() + facet_wrap(~sig.ctrl, ncol = 1, scales = "free")

#color by sigTF in CTRL:
ggplot(data = df_red_stats_VPA_Ctrlmerge, aes(x = meanDiff.ctrl, y = meanDiff.vpa, color = sig.ctrl)) + geom_point(size = 0.5) + 
  geom_vline(xintercept = 0, linetype = 2) + geom_hline(yintercept = 0, linetype = 2) + 
  geom_abline(slope = 1, intercept = 0) + 
  geom_text_repel(data = subset(df_red_stats_VPA_Ctrlmerge, meanDiff.ctrl < -0.3), aes(label = Var1), 
                  box.padding = 0.35, size = 2.5, segment.alpha = 0.3) +
  theme_classic()

ggplot(data = df_red_stats_VPA_Ctrlmerge, aes(x = meanDiff.ctrl, y = VPAvsCTRL, color = sig.ctrl)) + geom_point(size = 0.5) + 
  geom_vline(xintercept = 0, linetype = 2) + geom_hline(yintercept = 0, linetype = 2) + 
  geom_text_repel(data = subset(df_red_stats_VPA_Ctrlmerge, meanDiff.ctrl < -0.3), aes(label = Var1), 
                  box.padding = 0.35, size = 2.5, segment.alpha = 0.3) +
  theme_classic()

#the object is saved below as csv!
pdf(file = file.path(dir1, "bunina/LSD1/GRN/output/plots/newGRN2020", "ScatterVPAmostRescuedTFsNonDE.pdf"), 
    width = 4, height = 2, useDingbats = FALSE)
df_red_stats_VPA_Ctrlmerge %>%
  subset(sig.ctrl == TRUE) %>%
  ggplot(aes(x = meanDiff.ctrl, y = VPAvsCTRL, color = isTop30TF)) + geom_point(size = 0.5) + 
  scale_color_manual(values = c("#666666", "#E7298A", "#7570B3")) +
  geom_vline(xintercept = 0, linetype = 2) + geom_hline(yintercept = 0, linetype = 2) + 
  geom_text_repel(data = subset(df_red_stats_VPA_Ctrlmerge, isTop10TF == TRUE), aes(label = Var1), 
                  box.padding = 0.35, size = 2.5, segment.alpha = 0.3) +
  theme_classic()
dev.off()

pdf(file = file.path(dir1, "bunina/LSD1/GRN/output/plots/newGRN2020", "DensityVPAmostRescuedTFsNonDE.pdf"), 
    width = 4, height = 1.5, useDingbats = FALSE)
df_red_stats_VPA_Ctrlmerge %>%
  subset(sig.ctrl == TRUE) %>%
  ggplot(aes(x = meanDiff.ctrl, fill = isTop30TF)) + geom_density() + 
  scale_fill_manual(values = c("#666666", "#E7298A", "#7570B3")) +
  geom_vline(xintercept = 0, linetype = 2) + theme_classic()
dev.off()


pdf(file = file.path(dir1, "bunina/LSD1/GRN/output/plots/newGRN2020", "DensityVPAmostRescuedTFsNonDE_Yaxis.pdf"), 
    width = 3, height = 1, useDingbats = FALSE)
df_red_stats_VPA_Ctrlmerge %>%
  subset(sig.ctrl == TRUE) %>%
  ggplot(aes(x = VPAvsCTRL, fill = isTop30TF)) + geom_density() + 
  scale_fill_manual(values = c("#666666", "#E7298A", "#7570B3")) +
  geom_vline(xintercept = 0, linetype = 2) + theme_classic() + scale_x_reverse()
dev.off()

# as 2D density:
colorscale = scale_fill_gradientn(
  colors = rev(brewer.pal(9, "YlGnBu")),
  values = c(0, exp(seq(-5, 0, length.out = 100))))

ggplot(data = df_red_stats_VPA_Ctrlmerge, aes(x = meanDiff.ctrl, y = meanDiff.vpa)) + 
  stat_density2d(h = 0.2, bins = 60, aes( color = sig.ctrl), alpha = 0.2, geom = "polygon") +
  coord_fixed() +
  geom_vline(xintercept = 0, linetype = 2) + geom_hline(yintercept = 0, linetype = 2) + 
  geom_abline(slope = 1, intercept = 0) +
  theme_classic() #+ facet_wrap(~sig.ctrl, ncol = 1)

df_red_stats_VPA_Ctrlmerge$VPAvsCTRL = df_red_stats_VPA_Ctrlmerge$meanDiff.vpa - df_red_stats_VPA_Ctrlmerge$meanDiff.ctrl
ggplot(data = df_red_stats_VPA_Ctrlmerge, aes(x = VPAvsCTRL, fill = sig.ctrl)) + 
  geom_density(alpha = 0.5) + theme_classic() + geom_vline(xintercept = 0, linetype = 2)

df_red_stats_VPA_Ctrlmerge$TFname = hocoIDs$name[match(df_red_stats_VPA_Ctrlmerge$Var1, hocoIDs$HOCOID)]
df_red_stats_VPA_Ctrlmerge$Descr = egs_hg19_all$description[match(df_red_stats_VPA_Ctrlmerge$TFname, egs_hg19_all$external_gene_name)]
write.csv(df_red_stats_VPA_Ctrlmerge, file = file.path(dir1, "bunina/LSD1/GRN/output", "TFtargetsMeanDiff_endoVPA.csv"))


ggplot(data = df_red_stats_VPA_Ctrlmerge, aes(x = meanDiff.ctrl, y = VPAvsCTRL)) + 
  stat_density2d(h = 0.2, bins = 60, aes( color = sig.ctrl), alpha = 0.2, geom = "polygon") +
  coord_fixed() +
  geom_vline(xintercept = 0, linetype = 2) + geom_hline(yintercept = 0, linetype = 2) + 
  theme_classic()

###
### violin plots:
###

pdf(file = file.path(dir1, "bunina/LSD1/GRN/output/plots/newGRN2020", "ViolinMeanDiff_VPAendo_CTRLsig.pdf"), 
    width = 3, height = 3, useDingbats = FALSE)
ggplot(data = df_red_stats_VPA_Ctrl, aes(x = sig_ctrl, y = meanDiff)) + geom_violin(draw_quantiles = 0.5, aes(fill = type)) + 
  scale_fill_manual(values = c("#999999", "#CC79A7")) +
  theme_classic() + geom_hline(yintercept = 0, linetype = 2, size = 0.3)
#or boxes:
ggplot(data = df_red_stats_VPA_Ctrl, aes(x = sig_ctrl, y = meanDiff)) + geom_boxplot(aes(fill = type)) + 
  scale_fill_manual(values = c("#999999", "#CC79A7")) +
  theme_classic() + geom_hline(yintercept = 0, linetype = 2, size = 0.3)
#remove 2 outlier TFs:
df_red_stats_VPA_Ctrl %>%
  subset(Var1 != "ATF7" & Var1 != "PO3F4") %>%
  ggplot(aes(x = sig_ctrl, y = meanDiff)) + geom_violin(draw_quantiles = 0.5, aes(fill = type)) + 
  scale_fill_manual(values = c("#999999", "#CC79A7")) +
  theme_classic() + geom_hline(yintercept = 0, linetype = 2, size = 0.3)
df_red_stats_VPA_Ctrl %>%
  subset(Var1 != "ATF7" & Var1 != "PO3F4") %>%
  ggplot(aes(x = sig_ctrl, y = meanDiff)) + geom_boxplot(aes(fill = type)) + 
  scale_fill_manual(values = c("#999999", "#CC79A7")) +
  theme_classic() + geom_hline(yintercept = 0, linetype = 2, size = 0.3)
dev.off()

pdf(file = file.path(dir1, "bunina/LSD1/GRN/output/plots/newGRN2020", "ECDFMeanDiff_VPAendo_CTRLsig.pdf"), 
    width = 4, height = 2.5, useDingbats = FALSE)
ggplot(data = df_red_stats_VPA_Ctrl, aes(x = meanDiff, color = type)) + stat_ecdf(aes(linetype = sig_ctrl)) + 
  scale_color_manual(values = c("#999999", "#CC79A7")) +
  theme_classic()
dev.off()

### is the difference significant?
t.test(df_red_stats_VPA_Ctrl$meanDiff[df_red_stats_VPA_Ctrl$sig_ctrl == TRUE & df_red_stats_VPA_Ctrl$type == "CTRL"], 
            df_red_stats_VPA_Ctrl$meanDiff[df_red_stats_VPA_Ctrl$sig_ctrl == TRUE & df_red_stats_VPA_Ctrl$type == "VPA"]) #p-value = 0.0001904

t.test(df_red_stats_VPA_Ctrl$meanDiff[df_red_stats_VPA_Ctrl$sig_ctrl == FALSE & df_red_stats_VPA_Ctrl$type == "CTRL"], 
            df_red_stats_VPA_Ctrl$meanDiff[df_red_stats_VPA_Ctrl$sig_ctrl == FALSE & df_red_stats_VPA_Ctrl$type == "VPA"]) #p-value = 0.03153



### are affected TFs more connected to sig DE endo-Vpa genes? ====

GRNtable_sig_LSDp_3tissues_sel_k27ac$sigEndo_CTRL = ifelse(GRNtable_sig_LSDp_3tissues_sel_k27ac$ENSEMBL %in% resSigEndoCtrl$ensembl, TRUE, FALSE)
table(df_red_stats_VPA_Ctrl[,c("sig", "type")])
GRNtable_sig_LSDp_3tissues_sel_k27ac$sigTF_endoCTRL = ifelse(GRNtable_sig_LSDp_3tissues_sel_k27ac$TF %in% 
                                        df_red_stats_VPA_Ctrl$Var1[df_red_stats_VPA_Ctrl$sig == TRUE & 
                                                                     df_red_stats_VPA_Ctrl$type == "CTRL"], TRUE, FALSE)
GRN_TFgene_u = GRNtable_sig_LSDp_3tissues_sel_k27ac %>%
  distinct(TF, ENSEMBL, .keep_all = TRUE)

GRN_TFgene_u$sigEndo_CTRL = factor(GRN_TFgene_u$sigEndo_CTRL, levels = c(TRUE, FALSE))
GRN_TFgene_u$sigTF_endoCTRL = factor(GRN_TFgene_u$sigTF_endoCTRL, levels = c(TRUE, FALSE))

table(GRN_TFgene_u[,c("sigEndo_CTRL", "sigTF_endoCTRL")])
fisher.test(table(GRN_TFgene_u[,c("sigEndo_CTRL", "sigTF_endoCTRL")])) #oddsR 1.35, p-value < 0.00000000000000022

vpaDEG_TFs = as.data.frame(table(subset(GRN_TFgene_u, sigTF_endoCTRL == TRUE)[,c("TF", "sigEndo_CTRL")]))
vpaDEG_TFs = vpaDEG_TFs %>%
  pivot_wider(names_from = sigEndo_CTRL, values_from = Freq)
colnames(vpaDEG_TFs)[2:3] <- c("sigDEG", "other")
vpaDEG_TFs = vpaDEG_TFs %>% 
  subset(other > 0) %>%
  mutate(ratio = 100*sigDEG/(sigDEG + other)) %>%
  arrange(desc(sigDEG))

vpaDEG_TFs = vpaDEG_TFs %>% 
  subset(other > 100) %>%
  mutate(ratio = 100*sigDEG/(sigDEG + other)) %>%
  arrange(desc(ratio))

vpaDEG_TFs$GRNclassif = classif_list_GRN$classification_q0.1_final[match(vpaDEG_TFs$TF, classif_list_GRN$TF)]
table(vpaDEG_TFs$GRNclassif) #not more repressors!
vpaDEG_TFs$TFname = hocoIDs$name[match(vpaDEG_TFs$TF, hocoIDs$HOCOID)]
write.csv(vpaDEG_TFs, file = file.path(dir1, "bunina/LSD1/GRN/output", "vpaRescuedDEG_TFs.csv"))

# Add lineage expression data ---------------------------------------------

diff_genes_3layers = readRDS(file.path(dir1, "oheachte/LSD1/Tissue_comparison/MergedDiffGenesData", "MergedDiffGenesData.rds"))
length(unique(diff_genes_3layers$gene))
diff_genes_3layers_long = gather(diff_genes_3layers[,c(2,4,8,12)], key = "tissue", value = "log2FC_kidsVsDads", -gene)
diff_genes_3layers_long$tissue_m = sapply(strsplit(diff_genes_3layers_long$tissue, "[.]"), "[",2)
diff_genes_3layers_long$type = "all"
colnames(diff_genes_3layers_long)[[1]] <- "ENSEMBL"
diff_genes_3layers_long$TF = NA

length(diff_genes_3layers$gene[diff_genes_3layers$log2FoldChange.endo < 0 & diff_genes_3layers$significant.endo == TRUE])
diff_genes_3layers$sig_any = ifelse(diff_genes_3layers$significant.card == TRUE | diff_genes_3layers$significant.endo == TRUE | 
                                      diff_genes_3layers$significant.neur == TRUE, TRUE, FALSE)

diff_genes_3layers_long %>% 
  group_by(tissue_m) %>%
  summarise_at(.vars = "log2FC_kidsVsDads", .funs = mean, na.rm = TRUE)

#merge the two tables by targets:
GRNtable_sig_LSDp_3tissues = merge(GRNtable_sig_LSDp, diff_genes_3layers, by.x = "ENSEMBL", by.y = "gene", all.x = TRUE)
length(unique(GRNtable_sig_LSDp_3tissues$ENSEMBL)) #13k out of 28k

#before the loop run this:
#colnames(GRNtable_sig_LSDp_3tissues[,c(1,2,13,19,23,27)])

# df_red = GRNtable_sig_LSDp_3tissues[,c("ENSEMBL", "TF", "sigDiffTF", "log2FoldChange.endo", "log2FoldChange.card", "log2FoldChange.neur")]
df_red = GRNtable_sig_LSDp_3tissues_sel_k27ac[,c("ENSEMBL", "TF", "sigDiffTF", "log2FoldChange.endo", "log2FoldChange.card", "log2FoldChange.neur")]

df_red = df_red[!duplicated(df_red[,c(1,2)]),] #here: keep unique TF targets for all stats and plots!

#for stats:
df_red_stats = as.data.frame(table(df_red$TF))

#plot histogram of N of targets per TF:
ggplot(data = df_red_stats, aes(x = Freq)) + geom_histogram(binwidth = 100) + theme_classic()

df_long = gather(df_red, key = "tissue", value = "log2FC_kidsVsDads", -ENSEMBL, -TF, -sigDiffTF)
df_long$tissue_m = sapply(strsplit(df_long$tissue, "[.]"), "[",2)
df_long$type = "targets"

# pdf(file = file.path(dir1, "bunina/LSD1/GRN/output/plots", "BoxplotsLog2FCtargets_byTF.pdf"), width = 4, height = 4, useDingbats = FALSE) #1:365
pdf(file = file.path(dir1, "bunina/LSD1/GRN/output/plots/newGRN2020", "BoxplotsLog2FCtargets_byTF.pdf"), width = 4, height = 4, useDingbats = FALSE)
for (i in unique(GRNtable_sig_LSDp_3tissues$TF)) {
  sel_df = df_long[df_long$TF == i,c(1,2,5:7)]
  if(length(sel_df$log2FC_kidsVsDads[!is.na(sel_df$log2FC_kidsVsDads)]) > 0) {
    df_long_merge = rbind(df_long[df_long$TF == i,c(1,2,5:7)], diff_genes_3layers_long[,-2])
    df_long_merge = df_long_merge[!is.na(df_long_merge$log2FC_kidsVsDads),]
    # df_long_merge = df_long_merge %>%
    #   distinct(ENSEMBL,  tissue_m, type, .keep_all = TRUE) #already done see above for df_red!
    xlabs_boxes = paste(names(table(df_long_merge[,c(5,4)])[1,]), "\nn=", table(df_long_merge[,c(5,4)])[2,], 
                        sep = "")
    plot(ggplot(data = df_long_merge, aes(x = tissue_m, y = log2FC_kidsVsDads)) + 
           geom_violin(aes(fill = type), draw_quantiles = 0.5) + 
           theme_classic() + ylim(-3,3) + ggtitle(i) + scale_x_discrete(labels = xlabs_boxes))
  }
  
} #stat_summary(data = df_long_merge[df_long_merge$type == "target" & df_long_merge$tissue_m == "endo",], fun.data = give.n, geom = "text") + 
dev.off()

#same but all TFs together on one plot:

pdf(file = file.path(dir1, "bunina/LSD1/GRN/output/plots/newGRN2020", "BoxplotsLog2FCtargets_allTFonePlot.pdf"), 
    width = 4, height = 4, useDingbats = FALSE)

df_long_merge = rbind(df_long[,c(1,2,5:7)], diff_genes_3layers_long[,-2])
df_long_merge = df_long_merge[!is.na(df_long_merge$log2FC_kidsVsDads),]
xlabs_boxes = paste(names(table(df_long_merge[,c(5,4)])[1,]), "\nn=", table(df_long_merge[,c(5,4)])[2,], 
                    sep = "")
ggplot(data = df_long_merge, aes(x = tissue_m, y = log2FC_kidsVsDads)) + 
  geom_violin(aes(fill = type), draw_quantiles = 0.5) + 
  theme_classic() + ylim(-3,3) + ggtitle("All TFs in the network") + scale_x_discrete(labels = xlabs_boxes)
#take only unique genes:
df_long_merge = df_long_merge[!duplicated(df_long_merge[,c(1,4,5)]),]
xlabs_boxes = paste(names(table(df_long_merge[,c(5,4)])[1,]), "\nn=", table(df_long_merge[,c(5,4)])[2,], 
                    sep = "")
ggplot(data = df_long_merge, aes(x = tissue_m, y = log2FC_kidsVsDads)) + 
  geom_violin(aes(fill = type), draw_quantiles = 0.5) + 
  theme_classic() + ylim(-3,3) + ggtitle("All TFs in the network, unique targets") + scale_x_discrete(labels = xlabs_boxes)
dev.off()

#test for significance:
shapiro.test(df_long_merge$log2FC_kidsVsDads[df_long_merge$type == "targets" & df_long_merge$tissue_m == j]) #not normally distr.! use wilcox test!

#not normally distributed is not too bad, more important is the variance between groups, check it:
var(df_long_merge$log2FC_kidsVsDads[df_long_merge$type == "targets" & df_long_merge$tissue_m == "endo"]) #only one TF already in this df!
var(df_long_merge$log2FC_kidsVsDads[df_long_merge$type == "all" & df_long_merge$tissue_m == "endo"]) #not equal, so use welch t-test for this!

df_red_stats$endo = NA
df_red_stats$card = NA
df_red_stats$neur = NA
df_red_stats$endo_meanDiff = NA
df_red_stats$card_meanDiff = NA
df_red_stats$neur_meanDiff = NA

for (i in as.character(df_red_stats$Var1)) {
  for (j in unique(df_long$tissue_m)) {
    sel_df = df_long[df_long$TF == i,c(1,2,5:7)]
    if(length(sel_df$log2FC_kidsVsDads[!is.na(sel_df$log2FC_kidsVsDads)]) > 6) {
      df_long_merge = rbind(df_long[df_long$TF == i,c(1,2,5:7)], diff_genes_3layers_long[,-2])
      df_long_merge = df_long_merge[!is.na(df_long_merge$log2FC_kidsVsDads),]
      df_long_merge = df_long_merge[!duplicated(df_long_merge[,c(1,4,5)]),] #filter for unique genes in each tissue!
      df_red_stats[df_red_stats$Var1 == i,j] = t.test(df_long_merge$log2FC_kidsVsDads[df_long_merge$type == "targets" & df_long_merge$tissue_m == j],
                                                           df_long_merge$log2FC_kidsVsDads[df_long_merge$type == "all" & df_long_merge$tissue_m == j],
                                                      var.equal = FALSE)$p.value
      df_red_stats[df_red_stats$Var1 == i,paste0(j, "_meanDiff")] = mean(df_long_merge$log2FC_kidsVsDads[df_long_merge$type == "targets" & df_long_merge$tissue_m == j], na.rm = TRUE) - 
        mean(df_long_merge$log2FC_kidsVsDads[df_long_merge$type == "all" & df_long_merge$tissue_m == j], na.rm = TRUE)
    }
    
  }
}
df_red_stats$sigDiffTF = ifelse(df_red_stats$Var1 %in% TFresults_sig$TF, TRUE, FALSE)

#make a long format of it:
df_red_stats_l = gather(df_red_stats[,1:5], endo:neur, key = "tissue", value = "p.value")
df_red_stats_l = df_red_stats_l[!is.na(df_red_stats_l$p.value),]
df_red_stats_l$p.adjust = p.adjust(df_red_stats_l$p.value, method = "BH")
df_red_stats_l$sig = ifelse(df_red_stats_l$p.adjust < 0.01, TRUE, FALSE)
table(df_red_stats_l[,c(3,6)])
df_red_stats_l[df_red_stats_l$Var1 == "ZEB1",]

#now the same long format but for meanDiff values:
df_red_stats_l1 = gather(df_red_stats[,c(1:2,6:8)], endo_meanDiff:neur_meanDiff, key = "tissue_m", value = "meanDiff")
df_red_stats_l1$tissue = sapply(strsplit(df_red_stats_l1$tissue_m, "_"), "[", 1)
df_red_stats_l = merge(df_red_stats_l, df_red_stats_l1[,-3], by = c("Var1", "Freq", "tissue"))
df_red_stats_l_cc = df_red_stats_l[complete.cases(df_red_stats_l$p.adjust),]

#are TFs changing gene expression enriched for diffTFs?
df_red_stats_l_cc$sigDiffTF = ifelse(df_red_stats_l_cc$Var1 %in% TFresults_sig$TF, TRUE, FALSE)
df_red_stats_l_cc[,c(6,8)] = lapply(df_red_stats_l_cc[,c(6,8)], factor, levels = c(TRUE, FALSE))
table(df_red_stats_l_cc[df_red_stats_l_cc$tissue == "neur",c(6,8)])
fisher.test(table(df_red_stats_l_cc[df_red_stats_l_cc$tissue == "neur",c(6,8)])) #no, for GRNtable_sig: oddsR 1.6, p-value = 0.09652

unique(df_red_stats_l_cc$tissue)
fisher.test(table(df_red_stats_l_cc[df_red_stats_l_cc$tissue == "card",c(6,8)])) #a little for cardio: oddsR 2.3, p-value = 0.01982; not for GRNtable_sig p-value = 0.2418
fisher.test(table(df_red_stats_l_cc[df_red_stats_l_cc$tissue == "endo",c(6,8)])) #nope, neither for GRNtable_sig: p-value = 0.7629, oddsR 1.12

###
### add TFs expression data and if they are sig. differential in tissues:
###

df_red_stats_l_cc$TF_ensembl = TFresults$ensembl[match(df_red_stats_l_cc$Var1, TFresults$TF)]
df_red_stats_l_cc_merge = merge(df_red_stats_l_cc, diff_genes_3layers_long[,c(1,3,4)], by.x = c("tissue", "TF_ensembl"), 
                                  by.y = c("tissue_m", "ENSEMBL"), suffixes = c("_target", "_TF"), all.x = TRUE)
colnames(df_red_stats_l_cc_merge)

# is TF differential in any tissue?
diff_genes_3layers_long_pval = gather(diff_genes_3layers[,c(2,5,9,13)], key = "tissue", value = "padj_tissues", -gene)
diff_genes_3layers_long_pval$tissue_m = sapply(strsplit(diff_genes_3layers_long_pval$tissue, "[.]"), "[",2)
#merge:
df_red_stats_l_cc_merge2 = merge(df_red_stats_l_cc_merge, diff_genes_3layers_long_pval[,c(1,3,4)], by.x = c("tissue", "TF_ensembl"), 
                                  by.y = c("tissue_m", "gene"), suffixes = c("_target", "_TF"), all.x = TRUE)
df_red_stats_l_cc_merge2$sigTFexpr = ifelse(df_red_stats_l_cc_merge2$padj_tissues < 0.1, "bold", "plain")
#remove NAs, otherwise ggplot complains!!
df_red_stats_l_cc_merge2$sigTFexpr[is.na(df_red_stats_l_cc_merge2$sigTFexpr)] <- "plain"
colnames(df_red_stats_l_cc_merge2)

df_red_stats_l_merge_lsd$sigDiffTF = ifelse(df_red_stats_l_merge_lsd$Var1 %in% TFresults_sig$TF, TRUE, FALSE)
table(df_red_stats_l_merge_lsd$sigDiffTF)


#plot and color by TF expression:
library(wesanderson)

df_red_stats_l_cc_merge2_sel = df_red_stats_l_cc_merge2[df_red_stats_l_cc_merge2$Var1 %in% GRNtable_sig_LSDp_3tissues_sel_all$TF,]


#now plot - all TFs and only sigDiffTFs too:
pdf(file = file.path(dir1, "bunina/LSD1/GRN/output/plots/newGRN2020", "PerTFtest_log2FCdiffTargetsToAll.pdf"), width = 7, height = 4, useDingbats = FALSE)
ggplot(data = df_red_stats_l_cc, aes(x = meanDiff, y = -log10(p.adjust))) + geom_point(size = 0.4, alpha = 0.5, aes(color = sig)) + 
  theme_classic() + geom_hline(yintercept = -log10(0.01), color = "grey", linetype = 2) + facet_wrap(~tissue) + 
  geom_vline(xintercept = 0, color = "grey", linetype = 2) +
  geom_text_repel(data = subset(df_red_stats_l_cc, p.adjust < 0.01 & abs(meanDiff) > 0.1), aes(label = Var1), 
                  box.padding = 0.35, size = 2, segment.alpha = 0.3)
ggplot(data = df_red_stats_l_cc[df_red_stats_l_cc$sigDiffTF == TRUE,], aes(x = meanDiff, y = -log10(p.adjust))) + geom_point(size = 0.4, alpha = 0.5, aes(color = sig)) + 
  theme_classic() + geom_hline(yintercept = -log10(0.01), color = "grey", linetype = 2) + facet_wrap(~tissue) + 
  geom_vline(xintercept = 0, color = "grey", linetype = 2) +
  geom_text_repel(data = subset(df_red_stats_l_cc, sigDiffTF == TRUE & p.adjust < 0.01), aes(label = Var1), 
                  box.padding = 0.35, size = 2, segment.alpha = 0.3)
ggplot(data = df_red_stats_l_cc_merge2_sel, aes(x = meanDiff, y = -log10(p.adjust))) + 
  geom_point(size = 0.3, alpha = 0.99, aes(color = log2FC_kidsVsDads)) + 
  scale_color_gradientn(colors = wes_palette("Zissou1", 15, type = "continuous")) +
  geom_hline(yintercept = -log10(0.05), color = "red", linetype = 2, size = 0.2) +
  geom_text_repel(data = subset(df_red_stats_l_cc_merge2_sel, p.adjust < 0.05), 
                  aes(label = Var1, fontface = sigTFexpr), box.padding = 0.35, size = 2, segment.size = 0.2) + labs(color = "TF expr") +
  theme_classic() + facet_wrap(~tissue) + geom_vline(xintercept = 0, color = "blue", linetype = 2, size = 0.2)
# ggplot(data = df_red_stats_l, aes(x = Freq, y = -log10(p.adjust))) + geom_point(size = 0.4) + 
#   theme_classic() + geom_hline(yintercept = -log10(0.01), color = "red", linetype = 2) + facet_wrap(~tissue)
dev.off()

###
### plot corr of TF activ and targets expression:
###

df_red_stats_l_cc$TF_weighMeanDiff = TFresults$weighted_meanDifference[match(df_red_stats_l_cc$Var1, TFresults$TF)]
df_red_stats_l_cc$sig_both = ifelse(df_red_stats_l_cc$sig == TRUE & df_red_stats_l_cc$sigDiffTF == TRUE, TRUE, FALSE)
df_red_stats_l_cc$sig_both = factor(df_red_stats_l_cc$sig_both, levels = c(TRUE, FALSE))

pdf(file = file.path(dir1, "bunina/LSD1/GRN/output/plots/newGRN2020", "DiffTFvsGRNtargetsExprScatter_lineages.pdf"), 
    width = 15, height = 4, useDingbats = FALSE)
ggplot(data = df_red_stats_l_cc, aes(x = meanDiff, y = TF_weighMeanDiff, color = sig)) + geom_point(size = 0.5) + 
  theme_classic() + geom_vline(xintercept = 0, linetype = 2) + geom_hline(yintercept = 0, linetype = 2) + 
  geom_text_repel(data = subset(df_red_stats_l_cc, sig == TRUE & sigDiffTF == TRUE), aes(label = Var1), 
                  segment.size = 0.2, size = 3, color = "red", segment.alpha = 0.2) + geom_smooth(method = "lm", linetype = 2, size = 0.5, alpha = 0.5) +
  ylab("TF activity, LSD1mut/WT") + xlab("TF target mean expression, LSD1mut/WT") + facet_wrap(~tissue, scales = "free")
ggplot(data = df_red_stats_l_cc, aes(x = meanDiff, y = TF_weighMeanDiff, color = sig_both)) + geom_point(size = 0.5) + 
  theme_classic() + geom_vline(xintercept = 0, linetype = 2) + geom_hline(yintercept = 0, linetype = 2) + 
  geom_text_repel(data = subset(df_red_stats_l_cc, sig == TRUE & sigDiffTF == TRUE), aes(label = Var1), 
                  segment.size = 0.2, size = 3, color = "red", segment.alpha = 0.2) + geom_smooth(method = "lm", linetype = 2, size = 0.5, alpha = 0.5) +
  ylab("TF activity, LSD1mut/WT") + xlab("TF target mean expression, LSD1mut/WT") + facet_wrap(~tissue, scales = "free")
dev.off()

# 21.08.20 compare meanDiff in fibr, iPSCs, lineages ----------------------

df_red_compare = merge(df_red_stats_fib[,c(1,4,6)], df_red_stats_ipsc[,c(1,4,6)], by = "Var1", suffixes = c(".fib", ".ipsc"), all = TRUE)
# df_red_compare = merge(df_red_compare, df_red_lsd_stats[,c(1,6:8)], by = "Var1", all = TRUE) #non-lsd peak filtered df:df_red_stats and df_red_tissues_stats
df_red_compare = merge(df_red_compare, df_red_stats[,c(1,6:8)], by = "Var1", all = TRUE)

#tissue data are not multiple testing corrected! instead of testing again, merge with a df containing sig info (from a chunk below):
df_red_stats_l_merge_lsd_w # or this for non-LSD1filtered data: df_red_tissues_stats_l_w
# df_red_compare = merge(df_red_compare, df_red_stats_l_merge_lsd_w[,c("Var1","card", "endo", "neur")], by = "Var1", all.x = TRUE) #  df_red_tissues_stats_l_w
df_red_compare = merge(df_red_compare, df_red_tissues_stats_l_w[,c("Var1","card", "endo", "neur")], by = "Var1", all.x = TRUE)
df_red_compare[,9:11] <- apply(df_red_compare[,9:11], 2, as.logical)
colnames(df_red_compare)[9:11] <- paste0("sig.", colnames(df_red_compare)[9:11])

colnames(df_red_compare)[6:8] <- paste0("meanDiff.", c("endo", "card", "neur"))

#long format:
df_red_compare_l = gather(df_red_compare[,c(1,2,4,6:8)], key = "tissue", value = "meanDiff", -Var1)
df_red_compare_l$tissue_m = sapply(strsplit(df_red_compare_l$tissue, "[.]"), "[", 2) 

df_red_compare_l1 = gather(df_red_compare[,c(1,3,5,9:11)], key = "tissue", value = "sig", -Var1)
df_red_compare_l1$tissue_m = sapply(strsplit(df_red_compare_l1$tissue, "[.]"), "[", 2) 

df_red_compare_all = merge(df_red_compare_l[,c(1,3:4)], df_red_compare_l1[,c(1,3:4)], by = c("Var1", "tissue_m"))

df_red_compare_all_sig = df_red_compare_all[df_red_compare_all$sig == TRUE,]
df_red_compare_all_sig = df_red_compare_all_sig[!is.na(df_red_compare_all_sig$Var1),]
df_red_compare_all_sig$tissue_m = factor(df_red_compare_all_sig$tissue_m, levels = c("fib", "ipsc", "endo", "card", "neur"))


ggplot(data = df_red_compare_all_sig, aes(x = tissue_m, y = meanDiff, group = Var1)) + 
  geom_point(stat='summary', fun = sum) +
  stat_summary(fun=sum, geom="line")

### not too visual, do heatmap instead:
df_red_compare_all_sig_w = spread(df_red_compare_all_sig[,1:3], key = "tissue_m", value = "meanDiff")

breaks_comp = c(seq(min(df_red_compare_all_sig$meanDiff),-0.05, length.out = 12), 
                seq(-0.049,0.049, length.out = 6), seq(0.05,max(df_red_compare_all_sig$meanDiff), length.out = 12))

#cluster rows before plotting:
library(cluster)
df_red_compare_all_sig_w_c = df_red_compare_all_sig_w
df_red_compare_all_sig_w_c[is.na(df_red_compare_all_sig_w_c)] <- 0 #otherwise doesn't cluster
pam_comp = pam(df_red_compare_all_sig_w_c[,2:6], 10)

#now add to the original data frame:
df_red_compare_all_sig_w$pam = pam_comp$clustering
df_red_compare_all_sig_w = df_red_compare_all_sig_w[order(df_red_compare_all_sig_w$pam),]

gapRows = as.character(aggregate(row.names(df_red_compare_all_sig_w), list(df_red_compare_all_sig_w$pam), tail, 1)[,2])
gapRowsN = which(row.names(df_red_compare_all_sig_w) %in% gapRows)

pdf(file = file.path(dir1, "bunina/LSD1/GRN/output/plots/newGRN2020", "heatmap_meanDiff_compareAllcells.pdf"), 
    width = 3, height = 13, useDingbats = FALSE)
pheatmap(df_red_compare_all_sig_w[,2:6], cluster_cols = FALSE, cluster_rows = FALSE, cellwidth = 10, 
         color = colorRampPalette(rev(brewer.pal(n = 5, name = "RdBu")))(30), na_col = "grey",
         labels_row = df_red_compare_all_sig_w$Var1, fontsize = 4, breaks = breaks_comp, gaps_row = gapRowsN )
pheatmap(df_red_compare_all_sig_w[,2:6], cluster_cols = FALSE, cluster_rows = FALSE, cellwidth = 10, 
         color = colorRampPalette(rev(brewer.pal(n = 5, name = "RdBu")))(30), na_col = "grey", legend = FALSE,
         labels_row = df_red_compare_all_sig_w$Var1, fontsize = 4, breaks = breaks_comp, gaps_row = gapRowsN )

dev.off()

###
### same but only for diffTFs and non-ageing ones:

# df_red_compare_all_sig_w$sigDiffTF = ifelse(df_red_compare_all_sig_w$Var1 %in% TFresults_sig$TF, TRUE, FALSE)
df_red_compare_all_sig_w$finalTFs = ifelse(df_red_compare_all_sig_w$Var1 %in% unique(as.character(GRNtable_sig_LSDp_3tissues_selSig$TF)), TRUE, FALSE)

df_red_compare_all_sig_w_red = df_red_compare_all_sig_w[df_red_compare_all_sig_w$finalTFs == TRUE,]
df_red_compare_all_sig_w_red_c = df_red_compare_all_sig_w_red
df_red_compare_all_sig_w_red_c[is.na(df_red_compare_all_sig_w_red_c)] <- 0 #otherwise doesn't cluster
set.seed(22)
pam_comp = kmeans(df_red_compare_all_sig_w_red_c[,c(2:6)], 3) #pam

#now add to the original data frame:
df_red_compare_all_sig_w_red$pam = pam_comp$cluster # pam_comp$clustering
df_red_compare_all_sig_w_red = df_red_compare_all_sig_w_red[order(df_red_compare_all_sig_w_red$pam),]

gapRows1 = as.character(aggregate(row.names(df_red_compare_all_sig_w_red), list(df_red_compare_all_sig_w_red$pam), tail, 1)[,2])
gapRowsN1 = which(row.names(df_red_compare_all_sig_w_red) %in% gapRows1)

breaks_comp1 = c(seq(min(df_red_compare_all_sig_w_red[,2:6], na.rm = TRUE),-0.05, length.out = 12), 
                seq(-0.049,0.049, length.out = 6), seq(0.05,max(df_red_compare_all_sig_w_red[,2:6], na.rm = TRUE), length.out = 12))


pdf(file = file.path(dir1, "bunina/LSD1/GRN/output/plots/newGRN2020", "heatmap_meanDiff_compareAllcells_diffTFs.pdf"), 
    width = 3, height = 5, useDingbats = FALSE)
pheatmap(df_red_compare_all_sig_w_red[,2:6], cluster_cols = FALSE, cluster_rows = FALSE, cellwidth = 10, 
         color = colorRampPalette(rev(brewer.pal(n = 5, name = "RdBu")))(30), na_col = "grey",
         labels_row = df_red_compare_all_sig_w_red$Var1, fontsize = 7, breaks = breaks_comp1, gaps_row = gapRowsN1)
pheatmap(df_red_compare_all_sig_w_red[,2:6], cluster_cols = FALSE, cluster_rows = FALSE, cellwidth = 10, 
         color = colorRampPalette(rev(brewer.pal(n = 5, name = "RdBu")))(30), na_col = "grey", legend = FALSE,
         labels_row = df_red_compare_all_sig_w_red$Var1, fontsize = 7, breaks = breaks_comp1, gaps_row = gapRowsN1)
dev.off()

###
### same but add TFs own log2FC in all tissues:
###
df_red_compare_all_sig_w_red$ensembl = TFresults_sig$ensembl[match(df_red_compare_all_sig_w_red$Var1, TFresults_sig$TF)]

df_red_compare_all_sig_w_red$fibr_TF = res_fibOrdered_df$log2FoldChange[match(df_red_compare_all_sig_w_red$ensembl, res_fibOrdered_df$ensembl)]
df_red_compare_all_sig_w_red$ipsc_TF = res_ipsc_df$log2FoldChange[match(df_red_compare_all_sig_w_red$ensembl, res_ipsc_df$ensembl)]
df_red_compare_all_sig_w_red$endo_TF = diff_genes_3layers$log2FoldChange.endo[match(df_red_compare_all_sig_w_red$ensembl, diff_genes_3layers$gene)]
df_red_compare_all_sig_w_red$card_TF = diff_genes_3layers$log2FoldChange.card[match(df_red_compare_all_sig_w_red$ensembl, diff_genes_3layers$gene)]
df_red_compare_all_sig_w_red$neur_TF = diff_genes_3layers$log2FoldChange.neur[match(df_red_compare_all_sig_w_red$ensembl, diff_genes_3layers$gene)]

breaks_comp1m = c(seq(min(df_red_compare_all_sig_w_red[,c(2:6,9:13)], na.rm = TRUE),-0.05, length.out = 12), 
                 seq(-0.049,0.049, length.out = 6), seq(0.05,max(df_red_compare_all_sig_w_red[,c(2:6,9:13)], na.rm = TRUE), length.out = 12))

breaks_comp2 = c(seq(min(df_red_compare_all_sig_w_red[,c(2:6)], na.rm = TRUE),-0.01, length.out = 12), 
                  seq(-0.009,0.009, length.out = 6), seq(0.01,max(df_red_compare_all_sig_w_red[,c(2:6)], na.rm = TRUE), length.out = 12))

df_red_compare_all_sig_w_red$GRNclassif = classif_list_GRN$classification_q0.1_final[match(df_red_compare_all_sig_w_red$Var1, classif_list_GRN$TF)]
df_red_compare_all_sig_w_red$GRNclassif = factor(df_red_compare_all_sig_w_red$GRNclassif, levels = unique(df_red_compare_all_sig_w_red$GRNclassif))

df_red_compare_all_sig_w_red_reorder = df_red_compare_all_sig_w_red[order(df_red_compare_all_sig_w_red$GRNclassif),]

### Apr2021 heatmap TF targets all cell types Figure 2 ====

# pdf(file = file.path(dir1, "bunina/LSD1/GRN/output/plots/newGRN2020", "heatmap_meanDiff_compareAllcells_diffTFs_wTFexpr.pdf"), 
#     width = 4, height = 5, useDingbats = FALSE)
pdf(file = file.path(dir1, "bunina/LSD1/GRN/output/plots/newGRN2020", "heatmap_meanDiff_20TFsFinal.pdf"), 
    width = 4, height = 5, useDingbats = FALSE)
#order by activ/repress/nd:
pheatmap(df_red_compare_all_sig_w_red_reorder[,c(2:6,9:13)], cluster_cols = FALSE, cluster_rows = FALSE, cellwidth = 10, 
         color = colorRampPalette(rev(brewer.pal(n = 5, name = "PuOr")))(30), na_col = "grey", gaps_col = 5, 
         annotation_row = df_red_compare_all_sig_w_red_reorder[,"GRNclassif",drop=F],
         labels_row = df_red_compare_all_sig_w_red_reorder$Var1, fontsize = 7, breaks = breaks_comp1m, gaps_row = c(6,18))
pheatmap(df_red_compare_all_sig_w_red_reorder[,c(2:6,9:13)], cluster_cols = FALSE, cluster_rows = FALSE, cellwidth = 10, 
         color = colorRampPalette(rev(brewer.pal(n = 5, name = "PuOr")))(30), na_col = "grey", gaps_col = 5, 
         annotation_row = df_red_compare_all_sig_w_red_reorder[,"GRNclassif",drop=F],
         labels_row = df_red_compare_all_sig_w_red_reorder$Var1, fontsize = 7, breaks = breaks_comp2, gaps_row = c(6,18))
#old ordering:
pheatmap(df_red_compare_all_sig_w_red[,c(2:6,9:13)], cluster_cols = FALSE, cluster_rows = FALSE, cellwidth = 10, 
         color = colorRampPalette(rev(brewer.pal(n = 5, name = "PuOr")))(30), na_col = "grey", gaps_col = 5, 
         annotation_row = df_red_compare_all_sig_w_red[,"GRNclassif",drop=F],
         labels_row = df_red_compare_all_sig_w_red$Var1, fontsize = 7, breaks = breaks_comp1m, gaps_row = gapRowsN1)
pheatmap(df_red_compare_all_sig_w_red[,c(2:6,9:13)], cluster_cols = FALSE, cluster_rows = FALSE, cellwidth = 10, 
         color = colorRampPalette(rev(brewer.pal(n = 5, name = "PuOr")))(30), na_col = "grey", gaps_col = 5, 
         annotation_row = df_red_compare_all_sig_w_red[,"GRNclassif",drop=F],
         labels_row = df_red_compare_all_sig_w_red$Var1, fontsize = 7, breaks = breaks_comp1m, gaps_row = gapRowsN1)
pheatmap(df_red_compare_all_sig_w_red[,c(2:6,9:13)], cluster_cols = FALSE, cluster_rows = FALSE, cellwidth = 10, 
         color = colorRampPalette(rev(brewer.pal(n = 5, name = "PuOr")))(30), na_col = "grey", legend = FALSE, gaps_col = 5,
         labels_row = df_red_compare_all_sig_w_red$Var1, fontsize = 7, breaks = breaks_comp1m, gaps_row = gapRowsN1)
dev.off()

### try boxplots of targets expr and tf expr:
df_red_compare_all_sig_w_long1 = gather(df_red_compare_all_sig_w_red[,c(1:6,10:15)], key = "cellType", value = "meanExpr", fib:neur_TF)
unique(df_red_compare_all_sig_w_long1$cellType)
df_red_compare_all_sig_w_long1$type = ifelse(df_red_compare_all_sig_w_long1$cellType %in% c("fib", "ipsc", "endo", "card", "neur"), "target_gene", "TF")
df_red_compare_all_sig_w_long1$cellTypeM = sapply(strsplit(df_red_compare_all_sig_w_long1$cellType, "_"),"[", 1)
df_red_compare_all_sig_w_long1[df_red_compare_all_sig_w_long1=="fibr"] <- "fib"
df_red_compare_all_sig_w_long1$tissue = ifelse(df_red_compare_all_sig_w_long1$cellTypeM == "fib", "fib", 
                                               ifelse(df_red_compare_all_sig_w_long1$cellTypeM == "ipsc", "ipsc", "differentiation"))
df_red_compare_all_sig_w_long1$tissue = factor(df_red_compare_all_sig_w_long1$tissue, levels = c("fib", "ipsc", "differentiation"))

library(ggbeeswarm)
pdf(file = file.path(dir1, "bunina/LSD1/GRN/output/plots/newGRN2020", "BoxplmeanDiff_20TFsFinal.pdf"), 
    width = 5, height = 4, useDingbats = FALSE)
df_red_compare_all_sig_w_long1 %>%
  subset(tissue != "ipsc") %>%
  ggplot(aes(x = GRNclassif, y = meanExpr)) +
  geom_boxplot(aes(fill = tissue)) + facet_wrap(~type, scales = "free_y") + theme_classic() + 
  theme(axis.text.x = element_text(angle=45, hjust=1)) + geom_hline(yintercept = 0, linetype = 2)
df_red_compare_all_sig_w_long1 %>%
  subset(tissue != "ipsc") %>%
  ggplot(aes(x = GRNclassif, y = meanExpr, color = tissue)) + geom_boxplot() +
  geom_beeswarm(size = 1.5, dodge.width = 0.7) + facet_wrap(~type, scales = "free_y") + theme_classic() + 
  theme(axis.text.x = element_text(angle=45, hjust=1)) + geom_hline(yintercept = 0, linetype = 2)
dev.off()

###
### same but all values, not only sig!
###

df_red_compare_all_w = spread(df_red_compare_all[,1:3], key = "tissue_m", value = "meanDiff")
df_red_compare_all_w$finalTFs = ifelse(df_red_compare_all_w$Var1 %in% unique(as.character(GRNtable_sig_LSDp_3tissues_selSig$TF)), TRUE, FALSE)

df_red_compare_all_w_red = df_red_compare_all_w[df_red_compare_all_w$finalTFs == TRUE,]
df_red_compare_all_w_red$ensembl = TFresults_sig$ensembl[match(df_red_compare_all_w_red$Var1, TFresults_sig$TF)]

df_red_compare_all_w_red$fibr_TF = res_fibOrdered_df$log2FoldChange[match(df_red_compare_all_w_red$ensembl, res_fibOrdered_df$ensembl)]
df_red_compare_all_w_red$ipsc_TF = res_ipsc_df$log2FoldChange[match(df_red_compare_all_w_red$ensembl, res_ipsc_df$ensembl)]
df_red_compare_all_w_red$endo_TF = diff_genes_3layers$log2FoldChange.endo[match(df_red_compare_all_w_red$ensembl, diff_genes_3layers$gene)]
df_red_compare_all_w_red$card_TF = diff_genes_3layers$log2FoldChange.card[match(df_red_compare_all_w_red$ensembl, diff_genes_3layers$gene)]
df_red_compare_all_w_red$neur_TF = diff_genes_3layers$log2FoldChange.neur[match(df_red_compare_all_w_red$ensembl, diff_genes_3layers$gene)]

df_red_compare_all_w_red$GRNclassif = classif_list_GRN$classification_q0.1_final[match(df_red_compare_all_w_red$Var1, classif_list_GRN$TF)]
df_red_compare_all_w_red$GRNclassif = factor(df_red_compare_all_w_red$GRNclassif, levels = unique(df_red_compare_all_w_red$GRNclassif))

df_red_compare_all_w_red_reorder = df_red_compare_all_w_red[order(df_red_compare_all_w_red$GRNclassif),]

pdf(file = file.path(dir1, "bunina/LSD1/GRN/output/plots/newGRN2020", "heatmap_meanDiff_20TFsFinalfull.pdf"), 
    width = 4, height = 5, useDingbats = FALSE)
#order by activ/repress/nd:
pheatmap(df_red_compare_all_w_red_reorder[,c(4,5,3,2,6,9:13)], cluster_cols = FALSE, cluster_rows = FALSE, cellwidth = 10, 
         color = colorRampPalette(rev(brewer.pal(n = 5, name = "PuOr")))(30), na_col = "grey", gaps_col = 5, 
         annotation_row = df_red_compare_all_w_red_reorder[,"GRNclassif",drop=F],
         labels_row = df_red_compare_all_w_red_reorder$Var1, fontsize = 7, breaks = breaks_comp1m, gaps_row = c(6,18)) #breaks_comp2
pheatmap(df_red_compare_all_w_red_reorder[,c(4,5,3,2,6)], cluster_cols = FALSE, cluster_rows = FALSE, cellwidth = 8, 
         color = colorRampPalette(rev(brewer.pal(n = 5, name = "PuOr")))(30), na_col = "grey", gaps_col = 5, 
         annotation_row = df_red_compare_all_w_red_reorder[,"GRNclassif",drop=F],
         labels_row = df_red_compare_all_w_red_reorder$Var1, fontsize = 7, breaks = breaks_comp2, gaps_row = c(6,18)) #
dev.off()

###
### same but with absolute expression of TFs (dirty, better use TPM values next time!!!)
###

df_red_compare_all_sig_w_red$fibr_TF_count = res_fibOrdered_df$baseMean[match(df_red_compare_all_sig_w_red$ensembl, res_fibOrdered_df$ensembl)]
df_red_compare_all_sig_w_red$ipsc_TF_count = res_ipsc_df$baseMean[match(df_red_compare_all_sig_w_red$ensembl, res_ipsc_df$ensembl)]
df_red_compare_all_sig_w_red$endo_TF_count = diff_genes_3layers$baseMean.endo[match(df_red_compare_all_sig_w_red$ensembl, diff_genes_3layers$gene)]
df_red_compare_all_sig_w_red$card_TF_count = diff_genes_3layers$baseMean.card[match(df_red_compare_all_sig_w_red$ensembl, diff_genes_3layers$gene)]
df_red_compare_all_sig_w_red$neur_TF_count = diff_genes_3layers$baseMean.neur[match(df_red_compare_all_sig_w_red$ensembl, diff_genes_3layers$gene)]

pheatmap(df_red_compare_all_sig_w_red[,c(16:20)], cluster_cols = FALSE, cluster_rows = FALSE, cellwidth = 10, 
         color = colorRampPalette((brewer.pal(n = 5, name = "PuRd")))(30), na_col = "grey",
         labels_row = df_red_compare_all_sig_w_red$Var1, fontsize = 7, gaps_row = gapRowsN1)




### TF regulons expression densities in differentiated iPSCs ====

df_red_tissues_stats = as.data.frame(table(GRNtable_sig_LSDp_3tissues_sel_k27ac$TF))
df_red_tissues_stats$endo = NA
df_red_tissues_stats$card = NA
df_red_tissues_stats$neur = NA
df_red_tissues_stats$endo_meanDiff = NA
df_red_tissues_stats$card_meanDiff = NA
df_red_tissues_stats$neur_meanDiff = NA
df_red_tissues_stats = df_red_tissues_stats[df_red_tissues_stats$Freq > 0,]

#plot histogram of N of targets per TF:
ggplot(data = df_red_tissues_stats, aes(x = Freq)) + geom_histogram(binwidth = 100) + theme_classic()

for (i in as.character(df_red_tissues_stats$Var1)) {
  for (j in unique(df_long$tissue_m)) {
    sel_df = df_long[df_long$TF == i,c(1,2,5:7)]
    if(length(sel_df$log2FC_kidsVsDads[!is.na(sel_df$log2FC_kidsVsDads)]) > 6) {
      df_long_merge = rbind(df_long[df_long$TF == i,c(1,2,5:7)], diff_genes_3layers_long[,-2])
      df_long_merge = df_long_merge[!is.na(df_long_merge$log2FC_kidsVsDads),]
      # df_long_merge = df_long_merge %>%
      #   distinct(ENSEMBL,  tissue_m, type, .keep_all = TRUE) #unique genes are already picked, see df_red object few chunks above!
      df_red_tissues_stats[df_red_tissues_stats$Var1 == i,j] = t.test(df_long_merge$log2FC_kidsVsDads[df_long_merge$type == "targets" & df_long_merge$tissue_m == j], 
                                                              df_long_merge$log2FC_kidsVsDads[df_long_merge$type == "all" & df_long_merge$tissue_m == j], 
                                                              var.equal = FALSE)$p.value
      df_red_tissues_stats[df_red_tissues_stats$Var1 == i,paste0(j, "_meanDiff")] = mean(df_long_merge$log2FC_kidsVsDads[df_long_merge$type == "targets" & df_long_merge$tissue_m == j], na.rm = TRUE) - 
        mean(df_long_merge$log2FC_kidsVsDads[df_long_merge$type == "all" & df_long_merge$tissue_m == j], na.rm = TRUE)
      
    }
    
  }
}

#make a long format of it:
df_red_tissues_stats_l = gather(df_red_tissues_stats[,1:5], endo:neur, key = "tissue", value = "p.value")
df_red_tissues_stats_l = df_red_tissues_stats_l[!is.na(df_red_tissues_stats_l$p.value),]
df_red_tissues_stats_l$p.adjust = p.adjust(df_red_tissues_stats_l$p.value, method = "BH")
df_red_tissues_stats_l$sig = ifelse(df_red_tissues_stats_l$p.adjust < 0.05, TRUE, FALSE) #changed to 0.05 FDR! for sig. diffTFs!
table(df_red_tissues_stats_l[,c(3,6)])

df_red_tissues_stats_l_w = spread(df_red_tissues_stats_l[,c(1:3,6)], key = "tissue", value = "sig")

#now the same long format but for meanDiff values:
df_red_tissues_stats_l1 = gather(df_red_tissues_stats[,c(1:2,6:8)], endo_meanDiff:neur_meanDiff, key = "tissue_m", value = "meanDiff")
df_red_tissues_stats_l1$tissue = sapply(strsplit(df_red_tissues_stats_l1$tissue_m, "_"), "[", 1)
df_red_stats_l_merge_tissues = merge(df_red_tissues_stats_l, df_red_tissues_stats_l1[,-3], by = c("Var1", "Freq", "tissue"))

#plot:
pdf(file = file.path(dir1, "bunina/LSD1/GRN/output/plots/newGRN2020", "DensityMeanDiff_lineagesNoLSDfilter.pdf"), 
    width = 8, height = 2, useDingbats = FALSE)
df_red_stats_l_merge_tissues$sig = factor(df_red_stats_l_merge_tissues$sig, levels = c(TRUE, FALSE))
ggplot(data = df_red_stats_l_merge_tissues, aes(x = meanDiff)) + geom_density(aes(fill = sig), alpha = 0.4) + 
  theme_classic() + geom_vline(xintercept = 0, linetype = 2, size = 0.3) + facet_wrap(~tissue, scales = "free")
dev.off()

saveRDS(df_red_stats_l_merge_tissues, file = file.path(dir1, "bunina/LSD1/GRN/output/Robjects", "df_red_stats_l_merge_tissuesNOlsdFilt.rds"))


# Explore lsd-filtered diffTF sigLineages ---------------------------------

#filter the 3tissues network to keep only lsd-bound peaks and TFs that are diffTFs and sig. changing lineage expression (see heatmap above):
GRNtable_sig_LSDp_3tissues_sel = GRNtable_sig_LSDp_3tissues[GRNtable_sig_LSDp_3tissues$LSDpeak == TRUE & 
                                                              GRNtable_sig_LSDp_3tissues$TF %in% df_red_stats_l_merge_lsd_w[rowSums(df_red_stats_l_merge_lsd_w[,4:6]) > 0 & 
                                                                                                                              df_red_stats_l_merge_lsd_w$sigDiffTF == TRUE,]$Var1 ,]
length(unique(GRNtable_sig_LSDp_3tissues_sel$TF)) #29 (with FDR 0.05)
length(unique(GRNtable_sig_LSDp_3tissues_sel$ENSEMBL)) #5714

GRNtable_sig_LSDp_3tissues_sel[,c("significant.endo", "significant.card", "significant.neur")] <- lapply(GRNtable_sig_LSDp_3tissues_sel[,c("significant.endo", "significant.card", "significant.neur")], 
                                                                                                     factor, levels = c(TRUE, FALSE))
table(GRNtable_sig_LSDp_3tissues_sel$significant.endo)

###
### add ageing info to the table (of the TFs themselves, to filter out those that might be ageing-related):
###
GRNtable_sig_LSDp_3tissues_sel$TFname = hocoIDs$ENSEMBL[match(GRNtable_sig_LSDp_3tissues_sel$TF, hocoIDs$HOCOID)]
GRNtable_sig_LSDp_3tissues_sel$ageingTF = ifelse(GRNtable_sig_LSDp_3tissues_sel$TFname %in% ageGenes_sig$ensembl_gene_id | 
                                                   GRNtable_sig_LSDp_3tissues_sel$TFname %in% ageGenes300m$ensembl_gene_id, TRUE, FALSE)
GRNtable_sig_LSDp_3tissues_sel[,"ageingTF"] = factor(GRNtable_sig_LSDp_3tissues_sel[,"ageingTF"],levels = c(TRUE, FALSE) )

unique(GRNtable_sig_LSDp_3tissues_sel$TF[GRNtable_sig_LSDp_3tissues_sel$ageingTF == FALSE]) #20
unique(GRNtable_sig_LSDp_3tissues_sel$TF[GRNtable_sig_LSDp_3tissues_sel$ageingTF == TRUE])

###
#ADDED NEW 29thMay: SELECT ONLY NON-AGEING TFs!
###
GRNtable_sig_LSDp_3tissues_sel = GRNtable_sig_LSDp_3tissues_sel[GRNtable_sig_LSDp_3tissues_sel$ageingTF == FALSE,]

###
### subnetwork stats:
###

#unique peaks:
length(unique(GRNtable_sig_LSDp_3tissues_sel$peak)) #4160

#unique genes:
length(unique(GRNtable_sig_LSDp_3tissues_sel$ENSEMBL)) #5484

#total TF-peak connections:
nrow(GRNtable_sig_LSDp_3tissues_sel[!duplicated(GRNtable_sig_LSDp_3tissues_sel[,c(1,5)]),]) #25203, total rows 29165

#median number of peaks per TF:
grn_sigStatsSel = as.data.frame(table(GRNtable_sig_LSDp_3tissues_sel$TF[!duplicated(GRNtable_sig_LSDp_3tissues_sel[,c(1,5)])]))
median(grn_sigStatsSel$Freq) #777
ggplot(data = grn_sigStatsSel, aes(x = Freq)) + geom_histogram(bins = 100) + theme_classic()

#median number of genes per TF:(genes can be multiple times! not unique!)
grn_sigStatsSel$NgenesAll = as.data.frame(table(GRNtable_sig_LSDp_3tissues_sel$TF))$Freq
median(grn_sigStatsSel$NgenesAll) #850

###
#save table for cytoscape etc.:
###
write.csv(GRNtable_sig_LSDp_3tissues_sel[,-3], file = file.path(dir1, "bunina/LSD1/GRN/output", "GRNsubnetwork_lsdFilter_27052020.csv"), quote = F, row.names = F)

###
### same but only signif genes in any of the lineages:
###
GRNtable_sig_LSDp_3tissues_selSig = GRNtable_sig_LSDp_3tissues_sel[GRNtable_sig_LSDp_3tissues_sel$significant.card == TRUE | 
                                                                           GRNtable_sig_LSDp_3tissues_sel$significant.endo == TRUE | 
                                                                           GRNtable_sig_LSDp_3tissues_sel$significant.neur == TRUE,]
write.csv(GRNtable_sig_LSDp_3tissues_selSig, file = file.path(dir1, "bunina/LSD1/GRN/output", "GRNsubnetwork_lsdFilter_sigDiffTargets_27052020.csv"), quote = F, row.names = F)
saveRDS(GRNtable_sig_LSDp_3tissues_selSig, file = file.path(dir1, "bunina/LSD1/GRN/output/Robjects", "GRNtable_sig_LSDp_3tissues_selSig.rds"))
GRNtable_sig_LSDp_3tissues_selSig = readRDS(file.path(dir1, "bunina/LSD1/GRN/output/Robjects", "GRNtable_sig_LSDp_3tissues_selSig.rds"))

#compare TF IDs to non-LSD-filtered:
unique(GRNtable_sig_LSDp_3tissues_sel_all$TF[GRNtable_sig_LSDp_3tissues_sel_all$TF %in% 
                                               GRNtable_sig_LSDp_3tissues_sel$TF]) #all 29 are also in that set!
unique(GRNtable_sig_LSDp_3tissues_sel_all$TF[!(GRNtable_sig_LSDp_3tissues_sel_all$TF %in% 
                                               GRNtable_sig_LSDp_3tissues_sel$TF)])


#are TFs in the subnetwork enriched in sth?
library(msigdbr)
library(clusterProfiler)
# msigdbr_show_species()

#select set C2 - curated gene set, for the other options see this vignette: https://yulab-smu.github.io/clusterProfiler-book/chapter3.html#msigdb-analysis
m_t2g <- msigdbr(species = "Homo sapiens", category = "C2") %>% 
  dplyr::select(gs_name, gene_symbol)
head(m_t2g)

m_t2g_H <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, gene_symbol)

hocoIDs$name = egs_hg19_all$external_gene_name[match(hocoIDs$ENSEMBL, egs_hg19_all$ensembl_gene_id)]
length(hocoIDs$name[!is.na(hocoIDs$name)])

msig_tfs = enricher(gene = as.character(hocoIDs$name[match(unique(GRNtable_sig_LSDp_3tissues_sel$TF), hocoIDs$HOCOID)]), 
                    TERM2GENE = m_t2g, qvalueCutoff = 0.2, universe = as.character(hocoIDs$name))
head(msig_tfs)
msig_tfs = as.data.frame(msig_tfs)
#Zeb1 etc. are in this pathway: CHARAFE_BREAST_CANCER_LUMINAL_VS_MESENCHYMAL_DN

#save TFs in a file for string network etc.:
# write.table(as.character(hocoIDs$name[match(unique(GRNtable_sig_LSDp_3tissues_sel$TF), hocoIDs$HOCOID)]), 
#                                             file.path(dir1, "bunina/LSD1/GRN/output", "subnetworkTFs_diffTFs_sigLineageChanges.txt"), 
#                                        quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

#for GRNtable_sig unfiltered:
write.table(as.character(hocoIDs$name[match(unique(GRNtable_sig_LSDp_3tissues_sel$TF), hocoIDs$HOCOID)]), 
            file.path(dir1, "bunina/LSD1/GRN/output", "subnetworkTFs_diffTFs_sigLineageChanges_GRNunfiltered.txt"), 
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

###
#are all genes in the subnetwork enriched in sth?
###

msig_subnetw = enricher(gene = unique(GRNtable_sig_LSDp_3tissues_sel$SYMBOL), TERM2GENE = m_t2g, qvalueCutoff = 0.2, universe = unique(GRNtable_sig$SYMBOL))
head(msig_subnetw)
dotplot(msig_subnetw, showCategory = 20)
msig_subnetw_Alldf = as.data.frame(msig_subnetw)
msig_subnetw_Alldf = msig_subnetw_Alldf[order(msig_subnetw_Alldf$Count, decreasing = TRUE),]
msig_subnetw_df = msig_subnetw_df[msig_subnetw_df$qvalue < 0.1,] #all below 0.1

#take Hallmark genes dataset from MsigDB:
msig_subnetw_H = enricher(gene = unique(GRNtable_sig_LSDp_3tissues_sel$SYMBOL), TERM2GENE = m_t2g_H, qvalueCutoff = 0.2, 
                        universe = unique(GRNtable_sig$SYMBOL))
dotplot(msig_subnetw_H, showCategory = 20)

#GO terms instead:
go_subnetw = enrichGO(gene = unique(GRNtable_sig_LSDp_3tissues_sel$SYMBOL), OrgDb = org.Hs.eg.db, universe = unique(GRNtable_sig$SYMBOL), 
                      ont = "BP", keyType = "SYMBOL", pAdjustMethod = "BH", pvalueCutoff  = 0.05, qvalueCutoff  = 0.1)
dotplot(go_subnetw, showCategory = 12)

### for each TF separately:
GRNtable_sig_LSDp_3tissues_sel$TF = factor(GRNtable_sig_LSDp_3tissues_sel$TF) #this drops unused levels!
table(GRNtable_sig_LSDp_3tissues_sel$TF)

length(unique(GRNtable_sig_LSDp_3tissues_sel$peak)) #approx. half of peaks are unique, check the ovelap?
length(unique(GRNtable_sig_LSDp_3tissues_sel$ENSEMBL[GRNtable_sig_LSDp_3tissues_sel$TF == "ZEB1"]))

# do it with comparecluster function:
msig_subnetw_eachTF = compareCluster(SYMBOL ~ TF, data = GRNtable_sig_LSDp_3tissues_sel, fun = "enricher", TERM2GENE = m_t2g, universe = unique(GRNtable_sig$SYMBOL), 
                                     qvalueCutoff  = 0.1 )

#same for hallmark genes dataset: 
msig_subnetw_eachTF_H = compareCluster(SYMBOL ~ TF, data = GRNtable_sig_LSDp_3tissues_sel, fun = "enricher", TERM2GENE = m_t2g_H, universe = unique(GRNtable_sig$SYMBOL), 
                                     qvalueCutoff  = 0.1 )#

library(org.Hs.eg.db)
go_subnetw_eachTF = compareCluster(SYMBOL ~ TF, data = GRNtable_sig_LSDp_3tissues_sel, fun = "enrichGO", OrgDb = org.Hs.eg.db, universe = unique(GRNtable_sig$SYMBOL), 
                                   ont = "BP", keyType = "SYMBOL", pAdjustMethod = "BH", pvalueCutoff  = 0.05, qvalueCutoff  = 0.1)

pdf(file = file.path(dir1, "bunina/LSD1/GRN/output/plots/newGRN2020", "MsigDB_subnetworkTargetsEachTFdotplot.pdf"), 
    width = 12, height = 8, useDingbats = FALSE)
dotplot(go_subnetw_eachTF, font.size = 6, showCategory = 5)
dotplot(msig_subnetw_eachTF, font.size = 6, showCategory = 5)
dev.off()

pdf(file = file.path(dir1, "bunina/LSD1/GRN/output/plots/newGRN2020", "MsigDB_subnetworkTargetsEachTFdotplot_HallmarkGenes.pdf"), 
    width = 8, height = 5, useDingbats = FALSE)
dotplot(msig_subnetw_H, showCategory = 20, title = "Enrichments of all targets of all TFs")
dotplot(msig_subnetw_eachTF_H, font.size = 6, showCategory = 10, title = "Enrichments of targets per TF")
dev.off()

###
### now GO terms:
###
go_subnetw_eachTF = compareCluster(SYMBOL ~ TF, data = GRNtable_sig_LSDp_3tissues_sel, fun = "enrichGO", OrgDb = org.Hs.eg.db, universe = unique(GRNtable_sig$SYMBOL), 
                                   ont = "BP", keyType = "SYMBOL", pAdjustMethod = "BH", pvalueCutoff  = 0.05, qvalueCutoff  = 0.1 )
dotplot(go_subnetw_eachTF) #doesn't look very interesting, MsigDB had nicer terms!


###
### Reactome pathways:
###

#first add entrez gene IDs to the table:

GRNtable_sig_LSDp_3tissues_sel$entrez_id = egs_hg19$entrezgene_id[match(GRNtable_sig_LSDp_3tissues_sel$ENSEMBL, egs_hg19$ensembl_gene_id)]
length(GRNtable_sig_LSDp_3tissues_sel$entrez_id[is.na(GRNtable_sig_LSDp_3tissues_sel$entrez_id)]) #143
GRNtable_sig$entrez_id = egs_hg19$entrezgene_id[match(GRNtable_sig$ENSEMBL, egs_hg19$ensembl_gene_id)]

library(ReactomePA)
library(clusterProfiler)
x <- enrichPathway(gene = unique(GRNtable_sig_LSDp_3tissues_sel$entrez_id), pvalueCutoff=0.05, readable=T)
head(as.data.frame(x))
dotplot(x, showCategory = 20)
emapplot(x)
cnetplot(x, showCategory = 3)

go_subnetw_eachTF_reactome = compareCluster(entrez_id ~ TF, data = GRNtable_sig_LSDp_3tissues_sel, fun = "enrichPathway",  
                                   pvalueCutoff  = 0.05, qvalueCutoff  = 0.1 )
dotplot(go_subnetw_eachTF_reactome, showCategory = 5)

#save some of these plots:
pdf(file = file.path(dir1, "bunina/LSD1/GRN/output/plots/newGRN2020", "SubnetworkTargetsEachTF_ReactomePA.pdf"), 
    width = 15, height = 10, useDingbats = FALSE, pointsize = 8)
dotplot(x, showCategory = 20)
emapplot(x)
cnetplot(x, showCategory = 3)
dotplot(go_subnetw_eachTF_reactome, showCategory = 5, font.size = 5)
dev.off()

pdf(file = file.path(dir1, "bunina/LSD1/GRN/output/plots/newGRN2020", "SubnetworkTargetsEachTF_ReactomePA_cc.pdf"), 
    width = 8, height = 6, useDingbats = FALSE, pointsize = 8)
dotplot(go_subnetw_eachTF_reactome, showCategory = 5, font.size = 5)
dev.off()



# Heatmaps of selected TF expression across samples -----------------------

dds <- readRDS(file = file.path(dir1, "oheachte/GRN_new/prep_RNA", "dds.rds"))
sample_ids_grn = as.data.frame(colData(dds))
table(sample_ids_grn$cond)

vsd_grn <- varianceStabilizingTransformation(dds, blind=FALSE)
vsd_grn_df = as.data.frame(assay(vsd_grn))

#now plot the selected 20 TFs (non-ageing, from LSD-bound network):
sel_TFs_noAge = unique(as.character(GRNtable_sig_LSDp_3tissues_sel$TFname))

vsd_grn_df_tfs = vsd_grn_df[row.names(vsd_grn_df) %in% sel_TFs_noAge,]
vsd_grn_df_tfs$TF = GRNtable_sig_LSDp_3tissues_sel$TF[match(row.names(vsd_grn_df_tfs), GRNtable_sig_LSDp_3tissues_sel$TFname)]
row.names(vsd_grn_df_tfs) <- vsd_grn_df_tfs$TF

# gapRows = as.character(aggregate(row.names(protRNAdf_sel), list(protRNAdf_sel$Cluster.number), tail, 1)[,2])
# gapRowsN = which(row.names(protRNAdf_sel) %in% gapRows)

#add TFs significance in lineages+GRN data from above:
annoDF_tfs = df_red_stats_l_merge_lsd_w[rowSums(df_red_stats_l_merge_lsd_w[,4:6]) > 0 & 
                                          df_red_stats_l_merge_lsd_w$sigDiffTF == TRUE & 
                                          df_red_stats_l_merge_lsd_w$ageing == FALSE,4:6]

pdf(file = file.path(dir1, "bunina/LSD1/GRN/output/plots/newGRN2020", "PheatmapSelTFsExprGRNsamples.pdf"),
    width = 7, useDingbats = FALSE)
pheatmap(vsd_grn_df_tfs[,1:59], cluster_rows=TRUE, show_rownames=TRUE,
                          cluster_cols=TRUE, show_colnames = FALSE, annotation_row = annoDF_tfs,
                          annotation_col=sample_ids_grn[,1:2], color = viridis(50))
dev.off()

#can try normalize the counts first within each study by zscores across all genes??

vsd_grn_df$ensembl = row.names(vsd_grn_df)
vsd_grn_df_long = gather(vsd_grn_df, key = "sample", value = "log2_Ncounts", -ensembl)
unique(vsd_grn_df_long$sample)
#add study type as a group to calculate zscores:
vsd_grn_df_long$study = sample_ids_grn$treatment[match(vsd_grn_df_long$sample, row.names(sample_ids_grn))]

#z-score by group (study column in this case):
vsd_grn_df_long$zscores = ave(vsd_grn_df_long$log2_Ncounts, vsd_grn_df_long$study, FUN = scale)
vsd_grn_df_zsc = spread(vsd_grn_df_long[,-c(3:4)], key = sample, value = zscores)

vsd_grn_df_zsc_tfs = vsd_grn_df_zsc[vsd_grn_df_zsc$ensembl %in% sel_TFs_noAge,]
vsd_grn_df_zsc_tfs$TF = GRNtable_sig_LSDp_3tissues_sel$TF[match(vsd_grn_df_zsc_tfs$ensembl, GRNtable_sig_LSDp_3tissues_sel$TFname)]

breaks_zsc = c(seq(min(vsd_grn_df_zsc_tfs[,2:60]),-0.199, length.out = 10), seq(-0.2,0.2, length.out = 5), seq(0.201,max(vsd_grn_df_zsc_tfs[,2:60]), length.out = 10))

pdf(file = file.path(dir1, "bunina/LSD1/GRN/output/plots/newGRN2020", "PheatmapSelTFsExprGRNsamples.pdf"), 
    width = 7, useDingbats = FALSE)
pheatmap(vsd_grn_df_tfs[,1:59], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, show_colnames = FALSE, labels_row = vsd_grn_df_tfs[,"TF"],
         annotation_col=sample_ids_grn[,1:3], color = viridis(50))
#now z-scored values:
pheatmap(vsd_grn_df_zsc_tfs[,2:60], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, show_colnames = FALSE, labels_row = vsd_grn_df_zsc_tfs[,"TF"], breaks = breaks_zsc,
         annotation_col=sample_ids_grn[,1:3], color = colorRampPalette(rev(brewer.pal(n = 5, name = "RdBu")))(25))
dev.off()

#nice, now the same for lineage expression data. all objects are in andrew's folder and can be loaded using the script 
# /g/scb2/zaugg/oheachte/LSD1/Tissue_comparison/lineageMarkers/lineageMarkers_dataCompilation.R

tissues_vsnCounts = list()

for ( tissue in c("endo", "card", "neur") ){
  dds <- readRDS(file = file.path(dir1, "oheachte/LSD1", paste0(tissue, "_RNAseq"), "Saved_objects/dds.rds"))
  
  #assign(paste0("counts_", tissue), counts(dds, normalized=TRUE))  # Count matrix for each tissue
  #assign(paste0("dds_", tissue), DESeq2::results(dds))   # dds results object for each tissue
  
  #make normalized coutns:
  vsd_tis <- varianceStabilizingTransformation(dds, blind=FALSE)
  vsd_tis_df = as.data.frame(assay(vsd_tis))
  vsd_tis_df$ensembl = row.names(vsd_tis_df)
  vsd_tis_df$tissue = tissue
  tissues_vsnCounts[[tissue]] <- vsd_tis_df
}
tissues_anno = as.data.frame(colData(dds))
tissues_anno$sample = gsub("_", "", tissues_anno$id)

head(tissues_vsnCounts[[1]])
#append and spread the tissues df:
colnames(tissues_vsnCounts[[3]]) <- colnames(tissues_vsnCounts[[2]]) #colnames are with underscores in 3rd object, but all the same - replace those
tissues_vsnCounts_df = merge(tissues_vsnCounts[[1]][,-20], tissues_vsnCounts[[2]][,-20], by = "ensembl", 
                             suffixes = c(".endo", ".card"), all = TRUE)
tissues_vsnCounts_df = merge(tissues_vsnCounts_df, tissues_vsnCounts[[3]][,-20], by = "ensembl", 
                             suffixes = c("", ".neur"), all = TRUE)
# colnames(tissues_vsnCounts_df[,38:55]) <- paste0(colnames(tissues_vsnCounts_df[,38:55]), ".neur") #doesn't work somehow!

{
  dds <- readRDS(file = file.path(dir1, "bunina/LSD1/RNAseq/iPSCs/output/Robjects/dds_iPSCsDadsKidsRNAseq_noPro3.Rds"))
  rownames(dds) <- gsub("\\..*","", rownames(dds))
  dds <- dds[rowSums(counts(dds)) >= 10]
  vsd_tis <- varianceStabilizingTransformation(dds, blind=FALSE)
  vsd_tis_df = as.data.frame(assay(vsd_tis))
  vsd_tis_df$ensembl = row.names(vsd_tis_df)
  tissues_vsnCounts[["iPSCs"]] <- vsd_tis_df
} ; rm(dds)

#now merge with the tissues wide df:
head(tissues_vsnCounts[["iPSCs"]])
tissues_vsnCounts_df = merge(tissues_vsnCounts_df, tissues_vsnCounts[["iPSCs"]], by = "ensembl", 
                             all = TRUE)
colnames(tissues_vsnCounts_df)

#select TFs:
tissues_vsnCounts_df_tfs = tissues_vsnCounts_df[tissues_vsnCounts_df$ensembl %in% sel_TFs_noAge,]
tissues_vsnCounts_df_tfs$TF = GRNtable_sig_LSDp_3tissues_sel$TF[match(tissues_vsnCounts_df_tfs$ensembl, GRNtable_sig_LSDp_3tissues_sel$TFname)]

lineages_anno = data.frame(samples = colnames(tissues_vsnCounts_df_tfs[,2:64]), tissue = NA, type = NA)
lineages_anno$tissue = c(rep("endo", 18), rep("card", 18), rep("neur", 18), rep("iPSCs", 9))
lineages_anno$type = substr(lineages_anno$samples, start = 1, stop = 3)
lineages_anno$condition = ifelse(lineages_anno$type %in% c("Dad", "KSF"), "Dad", "Kid")
row.names(lineages_anno) <- lineages_anno$samples #has to be the same as colnames of the df for pheatmap!

#now plot:
row.names(tissues_vsnCounts_df_tfs) <- tissues_vsnCounts_df_tfs$TF

pdf(file = file.path(dir1, "bunina/LSD1/GRN/output/plots/newGRN2020", "PheatmapSelTFsExprLineagesOurData.pdf"), 
    width = 7, useDingbats = FALSE)
pheatmap(tissues_vsnCounts_df_tfs[,2:64], cluster_rows=TRUE, show_rownames=TRUE, 
         cluster_cols=TRUE, show_colnames = FALSE, annotation_row = annoDF_tfs,
         annotation_col=lineages_anno[,c(2,4)], color = viridis(50))
dev.off()


# 05.08.20 Endoderm + K27ac -----------------------------------------------

### add endoderm K27ac roadmap data to the network and split network into "pluripotency" vs. lineage specific, both by gene expression and peaks:

k27ac_all_dba_rAll = readRDS(file = file.path(dir1, "bunina/LSD1/GRN/output/Robjects", "k27ac_ESCvsDiff_dba_rAll.rds"))
resallVSiPSC = readRDS(file.path(dir1, "bunina/LSD1/differentiation", "res_allLineagesVSiPSCs.rds")) #contrast is set so lfc is diff/iPSCs
resallVSiPSC #contrast is set so lfc is diff/iPSCs
resallVSiPSC$ensembl = sapply(strsplit(row.names(resallVSiPSC), "[.]"), "[", 1)
resallVSiPSC = as.data.frame(resallVSiPSC)

#add peak and gene IDs for matching:
k27ac_all_dba_rAll = as.data.frame(k27ac_all_dba_rAll)
k27ac_all_dba_rAll$peakID = paste0(k27ac_all_dba_rAll$seqnames, ":",
                                  k27ac_all_dba_rAll$start, "-",
                                  k27ac_all_dba_rAll$end) #Fold is iPSCs/diff

#add to subnetwork:

# GRNtable_sig_LSDp_3tissues_sel_k27ac = GRNtable_sig_LSDp_3tissues_sel

#or to the full network:
GRNtable_sig_LSDp_3tissues_sel_k27ac = GRNtable_sig_LSDp_3tissues
GRNtable_sig_LSDp_3tissues_sel_k27ac$k27ac_Diff_vs_iPSC = - k27ac_all_dba_rAll$Fold[match(GRNtable_sig_LSDp_3tissues_sel_k27ac$peak, k27ac_all_dba_rAll$peakID)]
summary(GRNtable_sig_LSDp_3tissues_sel_k27ac$k27ac_Diff_vs_iPSC)
GRNtable_sig_LSDp_3tissues_sel_k27ac$k27ac_Diff_vs_iPSC_FDR = k27ac_all_dba_rAll$FDR[match(GRNtable_sig_LSDp_3tissues_sel_k27ac$peak, k27ac_all_dba_rAll$peakID)]
GRNtable_sig_LSDp_3tissues_sel_k27ac$k27ac_Diff_vs_iPSC_sig = ifelse(GRNtable_sig_LSDp_3tissues_sel_k27ac$k27ac_Diff_vs_iPSC_FDR < 0.1, TRUE, FALSE)



### May2021: all sig TFs heatmap new Figure 3 ====

classif_GRN_TFmeanDiffAll_n = merge(classif_list_GRN[,c(1:3,9,13:16)], df_red_compare_all_sig_w, by.x = "TF", by.y = "Var1")

### filter for TFs that are not diff. expressed in lineages:
classif_GRN_TFmeanDiffAll_n$TFensembl = hocoIDs$ENSEMBL[match(classif_GRN_TFmeanDiffAll_n$TF, hocoIDs$HOCOID)]
classif_GRN_TFmeanDiffAll_n$TF_DEG = ifelse(classif_GRN_TFmeanDiffAll_n$TFensembl %in% 
                                            diff_genes_3layers$gene[diff_genes_3layers$sig_any == TRUE], TRUE, FALSE) #if sig_any - 8 TFs more to filter out; significant.endo

#add log2FC of TF expr:
diff_genes_3layers$maxDE = apply(diff_genes_3layers[,c("log2FoldChange.endo", "log2FoldChange.card", "log2FoldChange.neur")],1,max, na.rm=TRUE )
classif_GRN_TFmeanDiffAll_n$maxDE_TF = diff_genes_3layers$maxDE[match(classif_GRN_TFmeanDiffAll_n$TFensembl, diff_genes_3layers$gene)]

table(classif_GRN_TFmeanDiffAll_n[,c("TF_DEG", "classification_q0.1_final")])
table(classif_GRN_TFmeanDiffAll_n$TF_DEG)
summary(classif_GRN_TFmeanDiffAll_n$maxDE_TF)


### how many TFs in each category:
classif_GRN_TFmeanDiffAll_m = classif_GRN_TFmeanDiffAll_n
classif_GRN_TFmeanDiffAll_m$fib_sign = ifelse(classif_GRN_TFmeanDiffAll_m$fib > 0, "pos", "neg")

#keep only TFs that are non-NAs in at least one of the lineages: (used to also filter for nonNA in fibro, not anymore!!! !is.na(classif_GRN_TFmeanDiffAll_m$fib) &)
classif_GRN_TFmeanDiffAll_m = classif_GRN_TFmeanDiffAll_m[ abs(rowSums(classif_GRN_TFmeanDiffAll_m[,11:13], na.rm = TRUE)) > 0,]
classif_GRN_TFmeanDiffAll_m$TF_type = paste(classif_GRN_TFmeanDiffAll_m$classification_q0.1_final, classif_GRN_TFmeanDiffAll_m$type, 
                                            classif_GRN_TFmeanDiffAll_m$sigDiffTF, classif_GRN_TFmeanDiffAll_m$fib_sign, sep = ".")

table(classif_GRN_TFmeanDiffAll_m$TF_type)

# replace NAs by 0 for clustering and keep only non-DEG TFs: & new: filter by max log2FC of TF in tissues:
classif_GRN_TFmeanDiffAll_m = subset(classif_GRN_TFmeanDiffAll_m, TF_DEG == FALSE & abs(maxDE_TF) < 0.58)
table(classif_GRN_TFmeanDiffAll_m$TF_DEG)

#in fibroblasts too?
classif_GRN_TFmeanDiffAll_m$TF_DEG_fibro_ipsc = ifelse(classif_GRN_TFmeanDiffAll_m$TFensembl %in% 
                                                         c(res_fibOrdered_df$ensembl[res_fibOrdered_df$padj > 0.1 & abs(res_fibOrdered_df$log2FoldChange) < 0.58], 
                                                           row.names(res_ipsc_mod_sig_df[abs(res_ipsc_mod_sig_df$log2FoldChange) < 0.58 & res_ipsc_mod_sig_df$padj > 0.1])), FALSE, TRUE)
table(classif_GRN_TFmeanDiffAll_m[,c("TF_DEG", "TF_DEG_fibro_ipsc")])

#cluster only on nonNA fibro, the NAs ones add as a separate cluster later:
classif_GRN_TFmeanDiffAll_ms = subset(classif_GRN_TFmeanDiffAll_m, !is.na(fib))
classif_GRN_TFmeanDiffAll_z = classif_GRN_TFmeanDiffAll_ms
classif_GRN_TFmeanDiffAll_z[is.na(classif_GRN_TFmeanDiffAll_z)] <- 0
km = kmeans(classif_GRN_TFmeanDiffAll_z[,c(9:13)], 2)

classif_GRN_TFmeanDiffAll_ms$clust = km$cluster
#optional: merge with fibro-NAs:
fibroNAs = subset(classif_GRN_TFmeanDiffAll_m, is.na(fib))
fibroNAs$clust = 3
classif_GRN_TFmeanDiffAll_ms = rbind(classif_GRN_TFmeanDiffAll_ms, fibroNAs)
classif_GRN_TFmeanDiffAll_ms = classif_GRN_TFmeanDiffAll_ms[order(classif_GRN_TFmeanDiffAll_ms$clust),] #, classif_GRN_TFmeanDiffAll_ms$classification_q0.1_final

gapRowsA = which(row.names(classif_GRN_TFmeanDiffAll_ms) %in% as.character(aggregate(row.names(classif_GRN_TFmeanDiffAll_ms), list(classif_GRN_TFmeanDiffAll_ms$clust), tail, 1)[,2]))

breaks_targets = c(seq(-0.15,-0.05, length.out = 12), 
                seq(-0.049,0.049, length.out = 6), seq(0.05,0.2, length.out = 12)) #max(classif_GRN_TFmeanDiffAll_ms[,c(9:13)], na.rm = TRUE)

write.csv(classif_GRN_TFmeanDiffAll_ms, file = file.path(dir1, "bunina/LSD1/GRN/output", "classif_GRN_TFmeanDiffAll_ms.csv"))

pdf(file = file.path(dir1, "bunina/LSD1/GRN/output/plots/newGRN2020", "heatmap_meanDiff_AllTFsNonDE_noLSDpFilter.pdf"), 
    useDingbats = FALSE, width = 4, height = 5)
pheatmap(classif_GRN_TFmeanDiffAll_ms[,c(9:13)], cluster_cols = FALSE, cluster_rows = FALSE, cellwidth = 20, 
         color = colorRampPalette(rev(brewer.pal(n = 5, name = "PuOr")))(30), na_col = "grey",
         annotation_colors = ann_colors, breaks = breaks_targets, gaps_row = gapRowsA, border_color = NA,
         labels_row = classif_GRN_TFmeanDiffAll_ms$TF, fontsize = 2, annotation_row =  classif_GRN_TFmeanDiffAll_ms[,c(4),drop=FALSE])
dev.off()

classif_GRN_TFmeanDiffAll_ms2 = subset(classif_GRN_TFmeanDiffAll_ms, abs(endo) > 0.1 | abs(neur) > 0.1 | abs(card) > 0.1 )
gapRowsA2 = which(row.names(classif_GRN_TFmeanDiffAll_ms2) %in% as.character(aggregate(row.names(classif_GRN_TFmeanDiffAll_ms2), 
                                                                                       list(classif_GRN_TFmeanDiffAll_ms2$clust), tail, 1)[,2]))

pheatmap(classif_GRN_TFmeanDiffAll_ms2[,c(9:13)], cluster_cols = FALSE, cluster_rows = FALSE, cellwidth = 20, 
         color = colorRampPalette(rev(brewer.pal(n = 5, name = "PuOr")))(30), na_col = "grey",
         annotation_colors = ann_colors, breaks = breaks_targets, gaps_row = gapRowsA2,
         labels_row = classif_GRN_TFmeanDiffAll_ms2$TF, fontsize = 4, annotation_row =  classif_GRN_TFmeanDiffAll_ms2[,c(4),drop=FALSE])

#short format:
pdf(file = file.path(dir1, "bunina/LSD1/GRN/output/plots/newGRN2020", "heatmap_meanDiff_AllTFsNonDE_short.pdf"), 
    useDingbats = FALSE, width = 4, height = 4)
#sort by activ/repr:
# pheatmap(classif_GRN_TFmeanDiffAll_ms1[,c(9:13)], cluster_cols = FALSE, cluster_rows = FALSE, cellwidth = 10, 
#          color = colorRampPalette(rev(brewer.pal(n = 5, name = "PuOr")))(30), na_col = "grey",
#          annotation_colors = ann_colors, breaks = breaks_targets, gaps_row = gapRowsA1,
#          labels_row = classif_GRN_TFmeanDiffAll_ms1$TF, fontsize = 2, annotation_row =  classif_GRN_TFmeanDiffAll_ms1[,c(4),drop=FALSE])
#sort first by clusters:
pheatmap(classif_GRN_TFmeanDiffAll_ms[,c(9:13)], cluster_cols = FALSE, cluster_rows = FALSE, cellwidth = 10, 
         color = colorRampPalette(rev(brewer.pal(n = 5, name = "PuOr")))(30), na_col = "grey",
         annotation_colors = ann_colors, breaks = breaks_targets, gaps_row = gapRowsA, border_color = NA,
         labels_row = classif_GRN_TFmeanDiffAll_ms$TF, fontsize = 2, annotation_row =  classif_GRN_TFmeanDiffAll_ms[,c(4),drop=FALSE])
pheatmap(classif_GRN_TFmeanDiffAll_ms2[,c(9:13)], cluster_cols = FALSE, cluster_rows = FALSE, cellwidth = 20, 
         color = colorRampPalette(rev(brewer.pal(n = 5, name = "PuOr")))(30), na_col = "grey",
         annotation_colors = ann_colors, breaks = breaks_targets, gaps_row = gapRowsA2, border_color = NA,
         labels_row = classif_GRN_TFmeanDiffAll_ms2$TF, fontsize = 4, annotation_row =  classif_GRN_TFmeanDiffAll_ms2[,c(4),drop=FALSE])
dev.off()

### top 10 affected TFs in endoderm:
classif_GRN_TFmeanDiffAll_ms_10 = classif_GRN_TFmeanDiffAll_ms[classif_GRN_TFmeanDiffAll_ms$TF %in% subnetw_30TFs$TF[1:10],]
pheatmap(classif_GRN_TFmeanDiffAll_ms_10[,c(9:13)], 
         cluster_cols = FALSE, cluster_rows = FALSE, cellwidth = 10, 
         color = colorRampPalette(rev(brewer.pal(n = 5, name = "RdBu")))(30), na_col = "grey",
         annotation_colors = ann_colors, breaks = breaks_targets, 
         labels_row = classif_GRN_TFmeanDiffAll_ms_10$TF, 
         fontsize = 8, annotation_row =  classif_GRN_TFmeanDiffAll_ms_10[,c(4),drop=FALSE])

write.csv(classif_GRN_TFmeanDiffAll_ms_10[,c("TF", "ENSEMBL")], file = file.path(dir1, "bunina/LSD1/GRN/output", "top10endoTFs.csv"))


### top endo meanDiff TFs subnetwork graph Figure 3 ====

classif_GRN_TFmeanDiffAll_ms$maxLineage = apply(classif_GRN_TFmeanDiffAll_ms[,c("endo", "card", "neur")],1,min, na.rm=TRUE )

subnetw_30TFs = classif_GRN_TFmeanDiffAll_ms %>%
  arrange(maxLineage)
subnetw_30TFs = subnetw_30TFs[1:30,] #if selecting w/LSD1peak - too few genes left! take more TFs
subnetw_30TFs_full = GRNtable_sig_LSDp_3tissues_sel_k27ac %>%
  subset(TF %in% subnetw_30TFs$TF & LSDpeak == TRUE) %>%
  droplevels()

head(as.data.frame(table(subnetw_30TFs_full$TF)))
subnetw_30TFs_full_sort = subnetw_30TFs_full %>%
  count(TF) %>%
  arrange(desc(n)) %>%
  mutate(TF = factor(TF, levels = unique(TF)))

subnetw_30TFs_full_sort %>%
  ggplot(aes(x = TF, y = n)) + 
  geom_bar(stat = "identity") + theme_classic() + theme(axis.text.x = element_text(angle=45, hjust=1, size = 8))

#how many genes are sig DE? Do for any lineage!!
subnetw_30TFs_full %>%
  subset(significant.any == TRUE) %>%
  count(TF) %>%
  arrange(desc(n)) %>%
  mutate(TF = factor(TF, levels = unique(TF))) %>%
  ggplot(aes(x = TF, y = n)) + 
  geom_bar(stat = "identity") + theme_classic() + theme(axis.text.x = element_text(angle=45, hjust=1, size = 8))


# as beeswarm:
library(ggbeeswarm)
subnetw_30TFs_full %>%
  mutate(TF = factor(TF, levels = levels(subnetw_30TFs_full_sort$TF) )) %>%
  ggplot(aes(x = TF, y = log2FoldChange.endo, color = classif)) + 
  geom_boxplot() + theme_classic() + geom_hline(yintercept = 0, linetype = 2) + 
  theme(axis.text.x = element_text(angle=45, hjust=1, size = 9))

### plot top5 TFs network: (low N of connections)
library(ggraph)
library(igraph)

# a pic with all TFs: -> takes forever to plot! use only 20TF subnetwork
df_for_graph = subnetw_30TFs_full %>%
  subset(TF %in% subnetw_30TFs$TF[1:10] & significant.endo == TRUE) #& significant.endo == TRUE

GRN_vertices = data.frame(gene = c(as.character(unique(df_for_graph$TF)), unique(as.character(df_for_graph$SYMBOL))), 
                          type = c(rep("TF", length(as.character(unique(df_for_graph$TF)))), 
                                   rep("target", length(unique(df_for_graph$SYMBOL)))))
GRN_vertices$LSD1mut_LFC = df_for_graph$log2FoldChange.endo[match(GRN_vertices$gene, df_for_graph$SYMBOL)]
#with all targets might be duplicated gene names - if a TF is also a target. no duplicate vertices allowed!
# GRN_vertices$gene[GRN_vertices$gene == "BPTF" & GRN_vertices$type == "TF"] <- "BPTF.TF"

graph6 <- graph_from_data_frame(df_for_graph[,c("TF", "SYMBOL", "peak", "r_TF_peak", "classif")], 
                                vertices = GRN_vertices, directed = TRUE)
gg6 = ggraph(graph6, layout = "dh") +  
  geom_node_point(aes(shape = type, color = LSD1mut_LFC), size = 3) + 
  geom_edge_link0(edge_width = 0.2, edge_alpha = 0.6) + 
  scale_colour_gradient2(low = "#5e3c99", high = "#e66101", mid = "#e0e0e0") +
  theme_bw() + geom_node_text(aes(label = ifelse(type == "TF", as.character(name), NA_character_)), repel = TRUE, size = 2.5)
#low = "#4dac26", high = "#a6611a" # label = name

# pdf(file = file.path(dir1, "bunina/LSD1/GRN/output/plots/newGRN2020", "Graph10TFsSubnetw.pdf"), 
#     useDingbats = FALSE, width = 8.5, height = 6)
pdf(file = file.path(dir1, "bunina/LSD1/GRN/output/plots/newGRN2020", "Graph10TFsSubnetw_OnlyTFlab.pdf"), 
    useDingbats = FALSE, width = 8.5, height = 6)
gg6
dev.off()
#scale_color_manual(values = c("#999999", "#0072B2")) + 
# geom_node_point(aes(colour = factor(classif)), size = 0.3) + 
# scale_edge_color_manual(values = c("brown", "green", "#999999")) +

###
## how many targets are co-regulated by 2 or more TFs?:
###
df_for_graph %>%
  distinct(TF, ENSEMBL, .keep_all = TRUE) %>%
  filter(duplicated(.[["SYMBOL"]]) ) %>%
  distinct(SYMBOL) #25 targets

###
### are sig targets of 10 top TFs enriched in sth?
###

msig_10TFsigDEtargets = enricher(unique(df_for_graph$SYMBOL), universe = unique(as.character(GRNtable_sig_LSDp_3tissues_sel_k27ac$SYMBOL)),
                             TERM2GENE = m_t2g_mod, qvalueCutoff  = 0.1)
dotplot(msig_10TFsigDEtargets)
msig_10TFsigDEtargets = enricher(unique(df_for_graph$SYMBOL), 
                                 TERM2GENE = m_t2g, qvalueCutoff  = 0.1)
dotplot(msig_10TFsigDEtargets)
msig_10TFsigDEtargets_df = msig_10TFsigDEtargets@result

### some neuronal etc genes, are peaks enriched in polycomb? k27me3 peaks overlap (in iPSCs):
# table(GRNtable_sig_LSDp_3tissues_sel_k27ac$k27me3peak_iPSC)
# table(subnetw_30TFs_full$k27me3peak_iPSC)
# table(subnetw_30TFs_full[,c("k27me3peak_iPSC", "significant.any")])
# 
# fisher.test(data.frame(a = c(107,52871), b = c(14640,1494087))) #nope, even under-represented!

#are these TFs expressed in lineages?
seAllvsiPSC_m = readRDS(file.path(dir1, "bunina/LSD1/differentiation", "dds_allLineagesVSiPSCs.rds"))

plotManyCounts = function(genesTable, withNames = TRUE, ddsObject = dds_ipsc, modRowNames = TRUE, logScale = FALSE, boxplot = FALSE) {
  if (modRowNames == TRUE) {
    row.names(ddsObject) = sapply(strsplit(rownames(ddsObject), "[.]"), "[", 1)
  }
  genesList = rownames(ddsObject)[rownames(ddsObject) %in% genesTable$ensembl_gene_id] #to only take expressed here genes!!
  d = list()
  for (i in 1:length(genesList)) { 
    d[[genesList[[i]]]] <- plotCounts(ddsObject, gene=rownames(ddsObject)[pmatch(genesList[i], rownames(ddsObject))], 
                                      intgroup="condition", returnData=TRUE) 
  }
  d_m = do.call(what = cbind, args = d)
  cols = seq(1,length(genesList)*2,2)
  d_m = d_m[,c(cols, length(genesList)*2)]
  colnames(d_m)[1:length(genesList)] <- sapply(strsplit(colnames(d_m)[1:length(genesList)], "[.]"), "[", 1)
  colnames(d_m)[length(colnames(d_m))] = "condition"
  if (withNames == TRUE) {
    for (i in 1:length(colnames(d_m)[1:length(genesList)])) {
      colnames(d_m)[i] = as.character(genesTable[genesTable$ensembl_gene_id %in% colnames(d_m)[i] ,]$external_gene_name)
    }
  }
  d_m$id = row.names(d_m)
  #optional:
  d_m$line = colData(ddsObject)$line[match(d_m$id, row.names(colData(ddsObject)))]
  d_m$type = colData(ddsObject)$type[match(d_m$id, row.names(colData(ddsObject)))]
  #end of optional
  # d_long = gather(d_m, gene, normExpr, -condition, -id, factor_key = TRUE)
  # d_long = gather(d_m, gene, normExpr, -condition, -id, -replicate, -experiment, factor_key = TRUE)
  d_long = gather(d_m, gene, normExpr, -condition, -id, -line, -type, factor_key = TRUE) #changed this line
  # ggplot(d_long, aes(x=condition, y=normExpr, fill = gene)) + ggtitle("Expression of genes") +
  #   geom_boxplot()
  if (logScale == TRUE) {
    if (boxplot == TRUE) {
      ggplot(d_long, aes(x=gene, y=log2(normExpr), color = type)) + ylab("Normalized gene counts, log10") + 
        geom_boxplot(outlier.shape = NA) + theme_bw() + 
        geom_point(position = position_jitterdodge(dodge.width=0.9), size = 0.4)
    }
    else {
      ggplot(d_long, aes(x=type, y=log2(normExpr), color = line, shape = condition)) + ylab("Normalized gene counts, log10") + 
        geom_point() + facet_wrap(~gene, scales = "free") + theme_bw()
    }
    
  }
  else {
    if (boxplot == TRUE) {
      ggplot(d_long, aes(x=gene, y=normExpr, color = type)) + ylab("Normalized gene counts") + 
        geom_boxplot(outlier.shape = NA) + theme_bw() + 
        geom_point(position = position_jitterdodge(dodge.width=0.9), size = 0.4)
    }
    else {
      ggplot(d_long, aes(x=type, y=normExpr, color = line, shape = condition)) + ylab("Normalized gene counts") + 
        geom_point() + facet_wrap(~gene, scales = "free") + theme_bw() #position=position_jitter(w=0.3,h=0)
    }
    
  }
}
TFs10Genes = egs_hg19_all[egs_hg19_all$ensembl_gene_id %in% subnetw_30TFs$ENSEMBL[1:10],]
plotManyCounts(genesTable = TFs10Genes, ddsObject = seAllvsiPSC_m, boxplot = TRUE, logScale = TRUE)
plotManyCounts(genesTable = TFs10Genes, ddsObject = seAllvsiPSC_m)


### top cardio and neuro TFs ====

subnetw_30TFs_neuro = classif_GRN_TFmeanDiffAll_ms %>%
  arrange(neur)
write.csv(subnetw_30TFs_neuro[1:10,c("TF", "ENSEMBL")], file = file.path(dir1, "bunina/LSD1/GRN/output", "top10_neuroTFs.csv"))

subnetw_30TFs_card = classif_GRN_TFmeanDiffAll_ms %>%
  arrange(card)
write.csv(subnetw_30TFs_card[1:10,c("TF", "ENSEMBL")], file = file.path(dir1, "bunina/LSD1/GRN/output", "top10_mesodermTFs.csv"))



### GO enrichments of sig TF targets ====

subnetwAllSig_ug = subnetwAllSig %>%
  distinct(ENSEMBL, TFclusterNew, .keep_all = TRUE) %>%
  mutate(TFclusterNew = factor(TFclusterNew, levels = 1:3))

### do this only on unique genes in each cluster! remove all duplicated:
subnetwAllSig_ug_sel = subnetwAllSig_ug[!duplicated(subnetwAllSig_ug$ENSEMBL),]
table(subnetwAllSig_ug_sel$TFclusterNew)

msig_allTFs_u = compareCluster(SYMBOL ~ TFclusterNew, 
                             data = subnetwAllSig_ug_sel, fun = 'enricher', TERM2GENE = m_t2g_mod, 
                             universe = unique(as.character(GRNtable_sig_LSDp_3tissues_sel_k27ac$SYMBOL)), 
                             qvalueCutoff  = 0.1 )
dotplot(msig_allTFs_u, showCategory = 15)

### using MSigDB:
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
#select set H - hallmark genes; other options: C2 - curated gene set, for the other options see this vignette: https://yulab-smu.github.io/clusterProfiler-book/chapter3.html#msigdb-analysis
m_t2g_H <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, gene_symbol)

m_t2g <- msigdbr(species = "Homo sapiens", category = "C2") %>% 
  dplyr::select(gs_name, gene_symbol)
head(m_t2g)
m_t2g$source = sapply(strsplit(m_t2g$gs_name, "_"), "[", 1)
m_t2g_mod = subset(m_t2g, source %in% c("REACTOME", "KEGG"))
m_t2g_mod$gs_name <- sub("REACTOME_", "", m_t2g_mod$gs_name)
m_t2g_mod$gs_name <- sub("KEGG_", "", m_t2g_mod$gs_name)

msig_allTFs = compareCluster(SYMBOL ~ TFclusterNew, 
                             data = subnetwAllSig_ug, fun = 'enricher', TERM2GENE = m_t2g_mod, 
                             universe = unique(as.character(GRNtable_sig_LSDp_3tissues_sel_k27ac$SYMBOL)), 
               qvalueCutoff  = 0.1 )
dotplot(msig_allTFs, showCategory = 7)
msig_allTFs_df = as.data.frame(msig_allTFs@compareClusterResult)
msig_allTFs_df %>%
  group_by(Cluster) %>%
  arrange(desc(Count)) %>%
  top_n(7) %>%
  ggplot(aes(x = Cluster, y = Description, color = p.adjust, size = Count)) + geom_point() + theme_classic()

pdf(file = file.path(dir1, "bunina/LSD1/GRN/output/plots/newGRN2020", "GO_AllSigTFtargets.pdf"), 
    useDingbats = FALSE, width = 6, height = 4)
dotplot(msig_allTFs, showCategory = 10, font.size = 7)
msig_allTFs_df %>%
  group_by(Cluster) %>%
  arrange(desc(Count)) %>%
  top_n(10) %>%
  ggplot(aes(x = Cluster, y = Description, color = p.adjust, size = Count)) + geom_point() + theme_classic()
dev.off()



## SuppFig2G ====

### modified using 144 TFs subnetwork:
pdf(file = file.path(dir1, "bunina/LSD1/GRN/output/plots/newGRN2020", "ViolinEndoMutWT_endoVSipsc144TFs.pdf"), width = 3.5, height = 3)
subnetwAllSig %>%
  distinct(ENSEMBL, .keep_all = TRUE) %>%
  subset(!is.na(RNA_endo_vsPSC_dir) & significant.endo == TRUE ) %>%
  ggplot(aes(x = RNA_endo_vsPSC_dir, y = log2FoldChange.endo, fill = RNA_endo_vsPSC_dir )) + 
  geom_violin(draw_quantiles = 0.5) + geom_hline(yintercept = 0, linetype = 2) + theme_classic()
dev.off()

subnetwAllSig %>%
  distinct(ENSEMBL, .keep_all = TRUE) %>%
  subset(!is.na(RNA_endo_vsPSC_dir) & significant.endo == TRUE ) %>%
  group_by(RNA_endo_vsPSC_dir) %>%
  tally()



# 07.10.20 Giuseppe clusters k27ac signal ---------------------------------

###check k27ac signal at enhancers linked to Giuseppe genes in each cluster:

### add k27ac signal:
k27ac_endo_dba_rAll = readRDS(file.path(dir1, "bunina/LSD1/GRN/output/Robjects", "k27ac_endo_dba_rAll.rds"))
k27ac_endo_dba_rAll_df = as.data.frame(k27ac_endo_dba_rAll)
k27ac_endo_dba_rAll_df$peak = paste(k27ac_endo_dba_rAll_df$seqnames, paste(k27ac_endo_dba_rAll_df$start, k27ac_endo_dba_rAll_df$end, sep = "-"), sep = ":")

GRNtable_sig_LSDp_3tissues_sel_k27ac$k27ac_conc_PSC = k27ac_endo_dba_rAll_df$Conc_ESC[match(GRNtable_sig_LSDp_3tissues_sel_k27ac$peak, k27ac_endo_dba_rAll_df$peak)]

# select unique peaks:
GRNtable_sig_LSDp_3tissues_sel_k27ac_p = GRNtable_sig_LSDp_3tissues_sel_k27ac[!duplicated(GRNtable_sig_LSDp_3tissues_sel_k27ac$peak),]
colnames(GRNtable_sig_LSDp_3tissues_sel_k27ac_p)
unique(GRNtable_sig_LSDp_3tissues_sel_k27ac_p$giuseppeClust1) #this has all clusters

GRNtable_sig_LSDp_3tissues_sel_k27ac_p$giuseppeClust1 = factor(GRNtable_sig_LSDp_3tissues_sel_k27ac_p$giuseppeClust1, levels = 1:7)

ggplot(data = GRNtable_sig_LSDp_3tissues_sel_k27ac_p, aes(x = LSDpeak, y = k27ac_conc_PSC)) + geom_violin(aes(fill = giuseppeClust1), draw_quantiles = 0.5) + 
  theme_classic() + geom_hline(yintercept = 0, linetype = 2) #+ facet_wrap(~LSDpeak)

### what are these genes in clusters 2 & 5? string? -> not much in common... clust 5 has KRAS and couple other cell cycle proteins

#the distributions are too spread! do densities:
library(ggridges)

class(GRNtable_sig_LSDp_3tissues_sel_k27ac_p$k27ac_conc_PSC)
ggplot(data = GRNtable_sig_LSDp_3tissues_sel_k27ac_p, aes(x = k27ac_conc_PSC, y = fct_rev(giuseppeClust1), fill = stat(x))) + 
  geom_density_ridges_gradient() + facet_wrap(~LSDpeak) +
  scale_fill_viridis_c(name = "Giuseppe cluster", option = "C") + theme_classic() + labs(title = "K27ac levels at enhancers in PSCs")

ggplot(data = GRNtable_sig_LSDp_3tissues_sel_k27ac_p, aes(x = k27ac_conc_PSC, y = fct_rev(giuseppeClust1), group = giuseppeClust1)) + 
  geom_density_ridges(aes(fill = giuseppeClust1)) + facet_wrap(~LSDpeak) +
  theme_classic()

pdf(file = file.path(dir1, "bunina/LSD1/GRN/output/plots/newGRN2020", "k27ac_giuseppeClusters.pdf"), 
    height = 3, width = 5, useDingbats = FALSE)
GRNtable_sig_LSDp_3tissues_sel_k27ac_p %>%
  subset(giuseppeClust1 %in% c(1:6)) %>%
  ggplot(aes(x = k27ac_conc_PSC, y = fct_rev(giuseppeClust1), group = giuseppeClust1)) + 
  geom_density_ridges(aes(fill = giuseppeClust1)) + facet_wrap(~LSDpeak) +
  theme_classic()
dev.off()



# Jan2021 ggraph hairball + Zeb1 --------------------------------------------------

library(ggraph)
library(igraph)

# a pic with all TFs: -> takes forever to plot! use only 20TF subnetwork
sel_TFs = c("CREB5", "GLIS1", "ZBT18", "GLIS2", "GLIS3", "ZEB1", "HTF4")
df_for_graph = subset(GRNtable_sig_LSDp_3tissues_sel_k27ac, LSDpeak == TRUE & TF %in% sel_TFs) #levels(GRNtable_sig_LSDp_3tissues_selSig$TF)
df_for_graph$Zeb1 = ifelse(df_for_graph$TF == "ZEB1", TRUE, FALSE)
df_for_graph$Zeb1 = factor(df_for_graph$Zeb1, levels = c(TRUE, FALSE))
GRN_vertices = data.frame(gene = c(as.character(unique(df_for_graph$TF)), unique(as.character(df_for_graph$ENSEMBL))), 
                                              type = c(rep("TF", length(as.character(unique(df_for_graph$TF)))), 
                                                       rep("target", length(unique(df_for_graph$ENSEMBL)))))
graph6 <- graph_from_data_frame(df_for_graph[,c("TF", "ENSEMBL", "peak", "p.adj_peak_gene", "Zeb1")], 
                                vertices = GRN_vertices, directed = TRUE)
gg6 = ggraph(graph6, layout = "circle") + 
  geom_edge_fan(edge_alpha = 0.05, edge_width = 0.02) + 
  scale_color_manual(values = c("#999999", "#0072B2")) + 
  geom_edge_link0(aes(colour = factor(Zeb1))) +
  scale_edge_color_manual(values = c("#999999", "red")) +
  geom_node_point(aes(colour = factor(type)), size = 0.3) + theme_bw()
# gg6 #, levels = c(TRUE, FALSE)

png(file = file.path(dir1, "bunina/LSD1/GRN/output/plots/newGRN2020", "NetworkCircle_Zeb1_color.png"), 
    width = 7, height = 6, units = "in", res = 300)
gg6
dev.off()

#and no coloring (mimic full network):
gg6 = ggraph(graph6, layout = "circle") + 
  geom_edge_fan(edge_alpha = 0.05, edge_width = 0.02) + 
  scale_color_manual(values = c("#999999", "#0072B2")) + 
  geom_edge_link0(aes(colour = factor(Zeb1))) + scale_edge_color_manual(values = c("#999999", "black")) +
  geom_node_point(aes(colour = factor(type)), size = 0.3) + theme_bw()

png(file = file.path(dir1, "bunina/LSD1/GRN/output/plots/newGRN2020", "NetworkCircle_Zeb1_NOcolor.png"), 
    width = 7, height = 6, units = "in", res = 300)
gg6
dev.off()
