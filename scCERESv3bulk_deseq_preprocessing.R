#!/usr/bin/env /data/gersbachlab/mt274/software/miniconda3/envs/maria_r/R/bin/Rscript

# this script runs deseq on the counts data and does some preprocessing
# this preprocessing includes recentering the data appropriately on 0, as outlined above
# and re-affiliating guide information with this library
# the results of this are printed into the inputs for the ALlen Method

# load in necessary libraries
library(DESeq2) 
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(glue))      # supposedly in dplyr but not really
suppressPackageStartupMessages(library(stringr))   # for str_detect
suppressPackageStartupMessages(library(latex2exp)) # for pretty plotting
suppressPackageStartupMessages(library(tidyr))     # for pivoting tables
suppressPackageStartupMessages(library(ggplot2))   # for plotting
suppressPackageStartupMessages(library(GGally))    # for ggpairs

# prep for plotting
theme_set(theme_bw(base_size = 10))
theme_update(axis.text   = element_text(colour = "black"),
             axis.title  = element_text(color = "black"),
             legend.text = element_text(colour = "black"),
             legend.background = element_rect(fill="transparent"),
             panel.grid = element_blank(),
             panel.background = element_rect(color = "black", linewidth = 1, fill = "transparent"),
             plot.background  = element_rect(fill='transparent', color=NA),
             panel.spacing=unit(0.2, "lines"), legend.key = element_blank())
# define variables for saving information
plot.path = "/data/gersbachlab/mt274/scCERES_v3/scCERES_v3_bulk/deseq/plots/"

# load in guideinfo for after Deseq
guide_info = read.csv(file = "/data/gersbachlab/mt274/scCERES_v3/sc_k562_data/samplesheets/scCERESv3_guideinfo.csv")


# load in raw gRNA counts from alignment
counts = read.table("/data/gersbachlab/mt274/scCERES_v3/scCERES_v3_bulk/counts/WT_KRAB.counts.txt", header=TRUE)
# immediately plot pairwise correlations of counts between replicates
ggsave(paste0(plot.path,"counts_correlation.pdf"),
       ggpairs(counts %>% select(-grna), lower = list(continuous = wrap("points", alpha = 0.3, size=0.2))) +
           theme(axis.text.x = element_text(angle = 90)), 
       width =12, height=12)

## Define the DESeq datasets
# countData is just a table of read counts (no guide names)
countData = counts[c("WT_rep1", "WT_rep2","WT_rep3", "WT_rep4", 
                     "WT_rep5", "WT_rep6", "WT_rep7", "WT_rep8",
                     "KRAB_rep1","KRAB_rep2", "KRAB_rep3","KRAB_rep4",
                     "KRAB_rep5", "KRAB_rep6","KRAB_rep7","KRAB_rep8")]
# we can supply the gRNAs as rownames to the countData table
row.names(countData) <- counts[,"grna"]
# colData is a metadata table with info about the count columns
# e.g. condition, replicate, treatment, time, etc.
colData <- data.frame(condition = factor(c("WT","WT","WT","WT",
                                           "WT","WT","WT","WT",
                                           "KRAB","KRAB","KRAB","KRAB",
                                           "KRAB","KRAB","KRAB","KRAB")),
                      replicate = factor(c(1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8)))
# we can provide row names to our colData to reiterate the above info
row.names(colData) <- c("WT_rep1", "WT_rep2","WT_rep3", "WT_rep4",         
                        "WT_rep5", "WT_rep6", "WT_rep7", "WT_rep8",
                        "KRAB_rep1","KRAB_rep2", "KRAB_rep3","KRAB_rep4",
                        "KRAB_rep5", "KRAB_rep6","KRAB_rep7","KRAB_rep8")
# Define matrix and reference for DESeq
# "design" takes info about the counts from the colData table
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)
# relevel resets the control for comparison (unsorted, low, WT, etc)
dds$condition <- relevel(dds$condition, "WT")
dds <- DESeq(dds)   # Run DESeq

# Extract results for specific comparisons 
res_scceresv3 <- results(dds, contrast = c("condition","KRAB","WT"), 
                        independentFiltering = FALSE, cooksCutoff = Inf)

# shrink lfc values. Do not change lfcSE or stat or p
# shrinkage is just for visualization purposes.
resLFS <- lfcShrink(dds,coef="condition_KRAB_vs_WT", res=res_scceresv3, type="ashr")
res_scceresv3$log2FoldChange <- resLFS$log2FoldChange

# Plot Results
ggsave(filename = paste0(plot.path, "MAplot.pdf"), 
       plotMA(res_scceresv3, alpha = 0.1, main = "scCERES_v3 Bulk", ylim = c(-3,3)))
ggsave(filename = paste0(plot.path,"PCAplot.pdf"), 
       plotPCA( rlog(dds), intgroup = c( "condition", "replicate")))


# start preprocessing and labelling this data for plotting
res_scceresv3 = res_scceresv3 %>% as.data.frame() %>% tibble::rownames_to_column(var = "gRNAid") %>%
        dplyr::mutate(NT = str_detect(gRNAid, "SABETI_NT")) %>% 
        merge(guide_info, by = c("gRNAid", "NT"), all.x = TRUE) %>% 
        dplyr::rename(element = DHS)

        
# some functions for plotting
maria_pvalue_histogram <- function(my_data, p) {
  # use pvalue as an input so you can make histograms of transformed data
  # column NT serves for coloring
  ggplot(my_data %>% filter(!is.na({{p}})), aes(x = pvalue, fill = NT)) + 
    geom_histogram(binwidth = 0.01, boundary = 0, alpha =0.9, position = "identity") + 
    labs(x = "Unadjusted p-value", y = "Number of gRNAs",
         subtitle = sprintf("KS test p = %1.4f",
                    ks.test(my_data %>% filter(NT) %>% 
                            select({{p}}) %>% unlist(),"punif")$p.value)) +
    scale_fill_manual(values = c("grey", "red"), limits = c(TRUE,FALSE),
                      labels = c("Nontargeting", "Targeting"), name = "")
}

maria_violin_plot <- function(df,stat_col) {
    # z-test for difference of NT to mean
    sample = df %>% filter(str_detect(element,"nontargeting")) %>% select({{stat_col}}) %>% unlist()
    z_stat = mean(sample) / (sd(sample)/sqrt(length(sample)))
    p = pnorm(abs(mean(sample)), mean = 0, sd = sd(sample), lower.tail = FALSE)*2
    
    df %>% ggplot(aes(x = str_detect(element,"nontargeting"), y = {{stat_col}})) +
    ggbeeswarm::geom_quasirandom(alpha = 0.3, aes(color = str_detect(element,"nontargeting"))) +
    geom_hline(yintercept = mean(sample), linetype = "dashed") + 
    labs(x = "NT", y = "DESeq2 z-score", subtitle = sprintf("NT mean is 0, p = %0.3f",p)) + 
    theme(legend.position = "none")
}
# make some plots of log2FC and p-value to show the need for adjustment
ggsave(paste0(plot.path,"log2FC_original.pdf"), 
      maria_violin_plot(res_scceresv3,log2FoldChange))
ggsave(paste0(plot.path,"zscore_original.pdf"), 
      maria_violin_plot(res_scceresv3,stat))
ggsave(paste0(plot.path,"pvaluehist_original.pdf"),
      maria_pvalue_histogram(res_scceresv3,pvalue))
        
# adjust the log2FC and z-scores to be centered on 0      
mean_NT_log2FC = res_scceresv3 %>% filter(NT) %>% select(log2FoldChange) %>% unlist(use.names = FALSE) %>% mean()
mean_NT_z = res_scceresv3 %>% filter(NT) %>% select(stat) %>% unlist(use.names = FALSE) %>% mean()

# correct the effect sizes and p-values
res_scceresv3 = res_scceresv3 %>% 
  mutate(log2FoldChange = log2FoldChange - mean_NT_log2FC,
         stat = stat - mean_NT_z,
         pvalue = ifelse(stat < 0, pnorm(stat, lower.tail = T)*2,
                                   pnorm(stat, lower.tail = F)*2),
         padj = p.adjust(pvalue, method = "BH"))    

# remaking plots! (now adjusted)
ggsave(paste0(plot.path,"log2FC_adjusted.pdf"), 
      maria_violin_plot(res_scceresv3,log2FoldChange))
ggsave(paste0(plot.path,"zscore_adjusted.pdf"), 
      maria_violin_plot(res_scceresv3,stat))
ggsave(paste0(plot.path,"pvaluehist_adjusted.pdf"),
      maria_pvalue_histogram(res_scceresv3,pvalue))

# functions for more preprocessing
# function to apply significances given vectors of p-values and fdr
maria_addsignificance <- function(df, p, prefix, fdr) {
  df %>% mutate(significance = ifelse(df %>% select({{p}}) %>% unlist()<= fdr, 
                                      glue("{prefix}_Significant"), 
                                      glue("{prefix}_NotSignificant")),
                direction = ifelse(df %>% select(log2FoldChange) %>% unlist() >= 0,
                                   "Enriched", "Depleted"))
}
maria_addNTsig <- function(df, p, fdr) {
  # for df with NT column of logical values and a distanceToTSS column
  # customized to this particular dataset to make a some extra categories 
  # for HBG gRNAs (negative controls for growth, positive for gene expression)
  df %>% mutate(NTsig = 
  ifelse(df %>% select({{p}}) %>% unlist() >  fdr & df$NT == TRUE, 'NT_NotSignificant',
  ifelse(df %>% select({{p}}) %>% unlist() <= fdr & df$NT == TRUE, 'NT_Significant',
  ifelse(df %>% select({{p}}) %>% unlist() >  fdr & element == "chr11.784", "HBG_NotSignificant",
  ifelse(df %>% select({{p}}) %>% unlist() <= fdr & element == "chr11.784", "HBG_Significant",
  ifelse(df %>% select({{p}}) %>% unlist() >  fdr & df$NT == FALSE & 
         abs(df$distanceToTSS) <= 3000, 'Promoter_NotSignificant',
  ifelse(df %>% select({{p}}) %>% unlist() <= fdr & df$NT == FALSE &
         abs(df$distanceToTSS) <= 3000, 'Promoter_Significant',
  ifelse(df %>% select({{p}}) %>% unlist() <= fdr & df$NT == FALSE,
         'Targeting_Significant',
  ifelse(df %>% select({{p}}) %>% unlist() > fdr & df$NT == FALSE, 
         'Targeting_NotSignificant', NA))))) ))))
}

# finish preprocessing
res_scceresv3 = maria_addsignificance(df = res_scceresv3,
                                       p = padj, fdr = 0.05, 
                                       prefix = "scCERESv3bulk")
res_scceresv3 = maria_addNTsig(res_scceresv3, padj, 0.05)

# Save as a .csv file for input to the Allen method
write.csv(res_scceresv3, 
          "/data/gersbachlab/mt274/scCERES_v3/scCERES_v3_bulk/AllenMethod/inputs/scCERESv3bulk_deseq.csv")
sessionInfo()
