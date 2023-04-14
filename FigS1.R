## FigS1

library(tidyverse)
library(ggplot2)
'%ni%' <- Negate('%in%')
library(DESeq2)
library(qusage)
library(GSVA)
library(org.Hs.eg.db)
library(survival)
library(survminer)

proj.dir = ".../VIGEX,DAPPER/DAPPER/"
data.dir = ".../VIGEX,DAPPER/LMS/data/"
gene_sets.dir = ".../gene_set_analysis/_gene_sets/"
tcga.dir = ".../data_TCGA/"
fig.dir = ".../VIGEX,DAPPER/LMS/figures/GSVA_TCGA_sarcoma/"


# get TCGA data
library(UCSCXenaTools)
data(XenaData)

XenaData %>%
  filter(XenaHostNames == "toilHub") %>%
  filter(XenaCohorts == "TCGA Pan-Cancer (PANCAN)") %>%
  filter(XenaDatasets == "tcga_gene_expected_count")
## this is RSEM expected_count

## get expected counts
## https://xenabrowser.net/datapages/?dataset=tcga_gene_expected_count&host=https%3A%2F%2Ftoil.xenahubs.net&removeHub=http%3A%2F%2F127.0.0.1%3A7222

tcga_gene_expected_count = read.table(paste0(tcga.dir, "tcga_gene_expected_count"), header = T, sep = "\t")

head_cohort =
  XenaData %>%
  filter(XenaHostNames == "tcgaHub") %>% # select TCGA Hub
  XenaScan("TCGA")

cohorts.sarcoma =
  head_cohort %>%
  filter(str_detect(XenaCohorts, "Sarcoma")) %>%
  .$XenaCohorts %>% unique()

clin.sarcoma =
  XenaData %>%
  filter(XenaHostNames == "tcgaHub") %>%
  XenaScan("TCGA") %>%
  filter(DataSubtype == "phenotype", XenaCohorts == cohorts.sarcoma) %>%
  XenaGenerate() %>%
  XenaQuery() %>%
  XenaDownload() %>%
  XenaPrepare()
clin.sarcoma$SARC_clinicalMatrix
clin.sarcoma$SARC_survival.txt

clin.sarcoma.clinMat =
  clin.sarcoma$SARC_clinicalMatrix %>%
  mutate(sampleID = gsub("-",".",sampleID)) %>%
  filter(sample_type == "Primary Tumor") %>%
  # dplyr::select(sampleID, histological_type) %>% .$histological_type %>% table()
  filter(histological_type == "Leiomyosarcoma (LMS)") %>%
  filter(!sampleID %in% setdiff(clin.sarcoma$SARC_clinicalMatrix$sampleID %>% gsub("-",".",.), tcga_gene_expected_count %>% colnames()))

clin.sarcoma.surv =
  clin.sarcoma$SARC_survival.txt %>%
  mutate(sample = gsub("-",".",sample)) %>%
  filter(sample %in% clin.sarcoma.clinMat$sampleID) %>%
  filter(!sample %in% setdiff(clin.sarcoma$SARC_survival.txt$sample %>% gsub("-",".",.), tcga_gene_expected_count %>% colnames()))

tcga_gene_expected_countlog2.sarcoma =
  tcga_gene_expected_count %>%
  dplyr::select(sample, clin.sarcoma.clinMat$sampleID) %>%
  rename("sample"="gene_id")

## convert back to count
reverseLog2 = function(x){
  x = ifelse(x == 0, 0, 2^x)
  return(x)
}
tcga_gene_expected_count.sarcoma.mat =
  tcga_gene_expected_countlog2.sarcoma %>%
  column_to_rownames(var="gene_id") %>%
  mutate(across(.cols = everything(), .fns = reverseLog2)) %>%
  as.matrix() %>% round()

# DESeq2 vst normalization
clin.sarcoma.surv =
  clin.sarcoma.surv %>%
  mutate(os_bool = if_else(OS.time < clin.sarcoma.surv$OS.time %>% median(), "os_low",
                           if_else(OS.time >= clin.sarcoma.surv$OS.time %>% median(), "os_high", "ERROR"))) %>%
  mutate(pfi_bool = if_else(PFI.time < clin.sarcoma.surv$PFI.time %>% median(), "pfi_low",
                            if_else(PFI.time >= clin.sarcoma.surv$PFI.time %>% median(), "pfi_high", "ERROR")))

## get conditions
condition_os.sarcoma = factor(clin.sarcoma.surv$os_bool)
coldata_os.sarcoma <- data.frame(row.names=colnames(tcga_gene_expected_count.sarcoma.mat), cond = rep("sarcoma",length(colnames(tcga_gene_expected_count.sarcoma.mat))))

condition_pfs.sarcoma = factor(clin.sarcoma.surv$pfs_bool)
coldata_pfs.sarcoma <- data.frame(row.names=colnames(tcga_gene_expected_count.sarcoma.mat), cond = rep("sarcoma",length(colnames(tcga_gene_expected_count.sarcoma.mat))))

dds.sarcoma =
  DESeqDataSetFromMatrix(
    countData = tcga_gene_expected_count.sarcoma.mat,
    colData = coldata_os.sarcoma, # Data frame with annotation for our samples
    design = ~1 # Here we are not specifying a model
  )

dds_norm.sarcoma.vstBlind2 = vst(dds.sarcoma, blind = T)
vst_df.sarcoma.vstBlind2 = assay(dds_norm.sarcoma.vstBlind2) %>%
  as.data.frame() %>% # Make into a data frame
  tibble::rownames_to_column("Symbol") # Make Gene IDs into their own column

# convert gene_id to gene_name
## gencode.v23.annotation.gene.probemap
gencode_v23_geneID_SYMBOL_type = read.table(".../references/gencode_v23_geneID_SYMBOL_type.tsv", header = T, sep = "\t")
# gencode_v23_geneID_SYMBOL_type

vst_df.sarcoma.symbol.mat =
  vst_df.sarcoma.vstBlind2 %>%
  # mutate(Symbol = gsub("\\.[0-9]","",Symbol)) %>%
  left_join(., gencode_v23_geneID_SYMBOL_type, by=c("Symbol"="gene_id")) %>%
  dplyr::select(Symbol, gene_name, everything()) %>% # 60,498 rows
  filter(!is.na(gene_name)) %>%
  dplyr::select(-c(Symbol, gene_type)) %>%
  distinct(gene_name, .keep_all = TRUE) %>%
  column_to_rownames(var="gene_name") %>%
  as.matrix()
# vst_df.sarcoma.symbol.mat %>% dim() # [1] 58581   104

PredictIO_gene_list = read.table(".../gene_set_analysis/IO_signatures-main/PredictIO_gene_list.txt", header=T, sep="\t")
setdiff(PredictIO_gene_list$GeneID, vst_df.sarcoma.symbol.mat %>% rownames()) %>% length() # [1] 0


# --------------------------------------------

# GSVA

## calculate IO

source_dir = ".../gene_set_analysis/IO_signatures-main/"
source(".../gene_set_analysis/IO_signatures-main/PredictIO_signature_functions.R")
options("encoding" = "UTF-8")

IO_Sigs.sarcoma <- Calc_ALL_IO_Signatures(expr = vst_df.sarcoma.symbol.mat, tumor = "")
# IO_Sigs.sarcoma %>% write.table(., paste0(data.dir, "IO_signatures/IO_Sigs.sarcoma.vstBlind.primaryOnly.tsv"), quote = F, sep = "\t", row.names = T)
# IO_Sigs.sarcoma = read.table(paste0(data.dir, "IO_signatures/IO_Sigs.sarcoma.vstBlind.primaryOnly.tsv"), header = T, sep = "\t")

## VEGFRi
VEGFRi_signatures.list = list()
gmt.flist = list.files(paste0(gene_sets.dir, "VEGFRi_Signature/_genesets/"), ".gmt")
for (gmt in gmt.flist){
  sig = qusage::read.gmt(paste0(gene_sets.dir, "VEGFRi_Signature/_genesets/", gmt))
  VEGFRi_signatures.list = c(VEGFRi_signatures.list, sig)
}
TCGA_VEGFRi_sarcoma =
  gsva(
    vst_df.sarcoma.symbol.mat,
    VEGFRi_signatures.list,
    method = "gsva",
    kcdf = "Gaussian",   # Appropriate for vst transformed data
    min.sz = 2,          # Minimum gene set size
    max.sz = 500,        # Maximum gene set size
    mx.diff = TRUE,      # Compute Gaussian-distributed scores
    verbose = FALSE      # Don't print out the progress bar
  )
# TCGA_VEGFRi_sarcoma %>% write.table(., paste0(data.dir, "IO_signatures/VEGFRi_Sigs.sarcoma.vst_blind.primaryOnly.tsv"), quote = F, sep = "\t", row.names = T)
# TCGA_VEGFRi_sarcoma = read.table(paste0(data.dir, "IO_signatures/VEGFRi_Sigs.sarcoma.vst_blind.primaryOnly.tsv"), header = T, sep = "\t")

## t() and merge results
TCGA_GSVA_IO_VEGFRi_Sigs.sarcoma =
  full_join(IO_Sigs.sarcoma %>% t() %>% as.data.frame() %>% rownames_to_column(var="Samples"),
            TCGA_VEGFRi_sarcoma %>% t() %>% as.data.frame() %>% rownames_to_column(var="Samples"),
            by=c("Samples"="Samples"))

# --------------------------------------------

# correlate with survival
pathways.vec = TCGA_GSVA_IO_VEGFRi_Sigs.sarcoma %>% column_to_rownames(var="Samples") %>% colnames()

fig.subdir.os.sarcoma = ".../VIGEX,DAPPER/LMS/figures/GSVA_TCGA_sarcoma/OS/"; dir.create(fig.subdir.os.sarcoma)
fig.subdir.pfs.sarcoma = ".../VIGEX,DAPPER/LMS/figures/GSVA_TCGA_sarcoma/PFS/"; dir.create(fig.subdir.pfs.sarcoma)

p_cutoff = 0.05
for (pathway in pathways.vec){

  # pathway = pathways.vec[1]; pathway
  # pathway = "TIDE"

  pathway.df = TCGA_GSVA_IO_VEGFRi_Sigs.sarcoma %>% dplyr::select(Samples, all_of(pathway))

  if (pathway.df %>% filter(is.na(!!as.symbol(pathway))) %>% nrow() == 0){
    # print(paste0("Pathway: ", pathway))

    colnames(pathway.df) = c("Sample", "Score")
    pathway_median = pathway.df$Score %>% median; pathway_median

    ## OS ############################################
    pathway_clinInfo.os =
      left_join(pathway.df, clin.sarcoma.surv, by=c("Sample"="sample")) %>%
      dplyr::select(Sample, Score, OS.time, OS) %>%
      dplyr::mutate(Score_bool = if_else(Score > pathway_median, "Score_high",
                                         if_else(Score <= pathway_median, "Score_low", "ERROR")))
    # pathway_clinInfo.os

    ## get pval and HR
    pathway_clinInfo.os.coxphfit = coxph(Surv(OS.time, OS) ~ Score_bool, data = pathway_clinInfo.os)
    pathway_clinInfo.os.coxphfit.pval = summary(pathway_clinInfo.os.coxphfit)$sctest['pvalue'] %>% unname(); pathway_clinInfo.os.coxphfit.pval # [1] 0.01193941
    pathway_clinInfo.os.coxphfit.hr = exp(pathway_clinInfo.os.coxphfit$coefficients) %>% unname(); pathway_clinInfo.os.coxphfit.hr # [1] 2.468573

    ## plot
    if (pathway_clinInfo.os.coxphfit.pval < p_cutoff){
      pathway_clinInfo.os.fit = survfit(Surv(OS.time, OS) ~ Score_bool, data = pathway_clinInfo.os)
      ggsurv.pathway_clinInfo.os.fit =
        ggsurvplot(pathway_clinInfo.os.fit,
                   data = pathway_clinInfo.os,
                   fun = "pct",
                   # pval = TRUE,
                   risk.table = TRUE)

      ggsurv.pathway_clinInfo.os.fit$plot +
        ggplot2::annotate("text",
                          x = Inf, y = Inf,
                          vjust = 1, hjust = 1,
                          label = paste0("HR = ", round(pathway_clinInfo.os.coxphfit.hr, 2), "\n",
                                         "p = ", round(pathway_clinInfo.os.coxphfit.pval, 5)),
                          size = 5
        ) +
        labs(title = paste0("TCGA sarcoma GSVA OS\npathway = ", pathway))

      ggsave(paste0(fig.subdir.os.sarcoma, "TCGA_sarcoma_GSVA_OS_All_", pathway, ".pdf"),
             device = "pdf", width = 5, height = 4, units = "in", useDingbats=F)
    }

    ## PFS ############################################
    pathway_clinInfo.pfs =
      left_join(pathway.df, clin.sarcoma.surv, by=c("Sample"="sample")) %>%
      dplyr::select(Sample, Score, PFI.time, PFI) %>%
      dplyr::mutate(Score_bool = if_else(Score > pathway_median, "Score_high",
                                         if_else(Score <= pathway_median, "Score_low", "ERROR")))
    # pathway_clinInfo.pfs

    ## get pval and HR
    pathway_clinInfo.pfs.coxphfit = coxph(Surv(PFI.time, PFI) ~ Score_bool, data = pathway_clinInfo.pfs)
    pathway_clinInfo.pfs.coxphfit.pval = summary(pathway_clinInfo.pfs.coxphfit)$sctest['pvalue'] %>% unname(); pathway_clinInfo.pfs.coxphfit.pval # [1] 0.01193941
    pathway_clinInfo.pfs.coxphfit.hr = exp(pathway_clinInfo.pfs.coxphfit$coefficients) %>% unname(); pathway_clinInfo.pfs.coxphfit.hr # [1] 2.468573

    ## plot
    if (pathway_clinInfo.pfs.coxphfit.pval < p_cutoff){
      pathway_clinInfo.pfs.fit = survfit(Surv(PFI.time, PFI) ~ Score_bool, data = pathway_clinInfo.pfs)
      ggsurv.pathway_clinInfo.pfs.fit =
        ggsurvplot(pathway_clinInfo.pfs.fit,
                   data = pathway_clinInfo.pfs,
                   fun = "pct",
                   # pval = TRUE,
                   risk.table = TRUE)

      ggsurv.pathway_clinInfo.pfs.fit$plot +
        ggplot2::annotate("text",
                          x = Inf, y = Inf,
                          vjust = 1, hjust = 1,
                          label = paste0("HR = ", round(pathway_clinInfo.pfs.coxphfit.hr, 2), "\n",
                                         "p = ", round(pathway_clinInfo.pfs.coxphfit.pval, 5)),
                          size = 5
        ) +
        labs(title = paste0("TCGA LMS GSVA PFS\npathway = ", pathway))

      ggsave(paste0(fig.subdir.pfs.sarcoma, "TCGA_sarcoma_PFS_All_", pathway, ".pdf"),
             device = "pdf", width = 5, height = 4, units = "in", useDingbats=F)
    }
  }
}


# --------------------------------------------

# heatmap

## get annot_df
annot_df.sarcoma =
  clin.sarcoma.surv %>%
  dplyr::select(sample, os_bool, pfi_bool) %>%
  arrange(os_bool) %>%
  column_to_rownames(var="sample")
# annot_df.sarcoma

TCGA_GSVA_IO_VEGFRi_Sigs.sarcoma.Heatmap =
  TCGA_GSVA_IO_VEGFRi_Sigs.sarcoma %>%
  arrange(match(Samples, rownames(annot_df.sarcoma))) %>%
  column_to_rownames(var="Samples") %>%
  as.matrix() %>% t()

fig.sarcoma.dir = ".../VIGEX,DAPPER/LMS/figures/GSVA_TCGA_sarcoma"
pdf(paste0(fig.sarcoma.dir, "TCGA_GSVA_IO_VEGFRi_Sigs.sarcoma.Heatmap.pdf"), width=10, height=3.5, useDingbats = F)
pheatmap::pheatmap(TCGA_GSVA_IO_VEGFRi_Sigs.sarcoma.Heatmap,
                   annotation_col = annot_df.sarcoma,
                   scale = "none",
                   show_colnames = FALSE,
                   angle_col = 90,
                   cluster_rows = FALSE,
                   cluster_cols = FALSE,
                   fontsize_row = 6)
dev.off()

# EOF