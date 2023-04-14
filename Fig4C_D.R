## Fig4C and D

library(tidyverse)
library(ggplot2)
'%ni%' <- Negate('%in%')
library(DESeq2)
library(qusage)
library(GSVA)
library(org.Hs.eg.db)
library(survival)
library(survminer)

proj.dir = ".../VIGEX,DAPPER/LMS/"
clin.dir = ".../VIGEX,DAPPER/LMS/_clin_info/"
data.dir = ".../VIGEX,DAPPER/DAPPER/data/"
fig.dir = ".../VIGEX,DAPPER/LMS/figures/"
data.lms.dir = ".../VIGEX,DAPPER/LMS/data/"
gene_sets.dir = ".../gene_set_analysis/_gene_sets/"

# read in data
clin_info_DAPPER =
  readxl::read_xlsx(".../VIGEX,DAPPER/DAPPER/clin_info/DAPPER_sampleinfo_2022_10_21_1242.xlsx", sheet = "lean_LMS", na = "NA") %>%
  dplyr::select(RNAseqID, DAP_ID)
# clin_info_DAPPER

clin_info =
  readxl::read_xlsx(paste0(clin.dir, "Copy of DAPPER_data_Leiomyosarcoma_2022_10_21.xlsx"), sheet = "Response", na = "NA") %>%
  mutate(TreatmentArm = if_else(`Assigned Treatment Arm` == "Cohort A", "D+Olaparib",
                                if_else(`Assigned Treatment Arm` == "Cohort B", "D+Cediranib", "ERROR"))) %>%
  right_join(., clin_info_DAPPER, by=c("Subject"="DAP_ID")) %>%
  rename("Best_Response_RECISTv1"="RECIST_bool") %>%
  dplyr::select(RNAseqID, Subject, TreatmentArm, Diagnosis, Age, Gender, Race, os_time, os_status, pfs_time, pfs_status, RECIST_bool)
clin_info$os_time %>% median() # [1] 14.58727
clin_info$pfs_time %>% median() # [1] 2.792608
# clin_info

clin_info2 =
  clin_info %>%
  mutate(os_bool = case_when(os_time >= 14 ~ "os_ge14m",
                             os_time < 14  ~ "os_lt14m"),
         pfs_bool = case_when(pfs_time >= 3 ~ "pfs_ge3m",
                              pfs_time < 3 ~  "pfs_lt3m"))
# clin_info2

## gene expression matrix
gene_expr.df =
  read.table(paste0(data.dir, "GE_matrix/gene_expr.dupsSumd.counts.csv"), header=T, sep = ",") %>%
  column_to_rownames(var="Symbol") %>%
  dplyr::select(clin_info$RNAseqID)


## DESeq2 normalization

## get countdata
countdata = gene_expr.df %>% round() %>% dplyr::select(clin_info2$RNAseqID)
countdata %>% dim() # [1] 59050    17

## get conditions
condition_os = factor(clin_info2$os_bool)
coldata_os <- data.frame(row.names=colnames(countdata), condition_os)

condition_pfs = factor(clin_info2$pfs_bool)
coldata_pfs <- data.frame(row.names=colnames(countdata), condition_pfs)

dds = DESeqDataSetFromMatrix(
  countData = countdata,
  colData = coldata_os,
  design = ~1
)

dds_norm = vst(dds, blind = TRUE)
vst_df = assay(dds_norm) %>%
  as.data.frame()
  tibble::rownames_to_column("Symbol")
# vst_df


# --------------------------------------------

# GSVA

## calculate IO
source_dir = ".../gene_set_analysis/IO_signatures-main/"
source(".../gene_set_analysis/IO_signatures-main/PredictIO_signature_functions.R")
options("encoding" = "UTF-8")

vst.lms.mat = vst_df %>%
  tibble::column_to_rownames("Symbol") %>%
  as.matrix()
vst.lms.mat %>% dim() # [1] 59050    17

DAPPER.IO_Sigs.lms = Calc_ALL_IO_Signatures(expr = vst.lms.mat, tumor = "")
# DAPPER.IO_Sigs.lms %>% write.table(., paste0(data.lms.dir, "IO_signatures/DAPPER.IO_Sigs.lms.tsv"), quote = F, sep = "\t", row.names = T)
# DAPPER.IO_Sigs.lms = read.table(paste0(data.lms.dir, "IO_signatures/DAPPER.IO_Sigs.lms.tsv"), header = T, sep = "\t")

## VEGFRi

## get VEGFRi signatures
VEGFRi_signatures.list = list()
gmt.flist = list.files(paste0(gene_sets.dir, "VEGFRi_Signature/_genesets/"), ".gmt")
# gmt.flist
for (gmt in gmt.flist){
  # gmt = gmt.flist[2];gmt
  sig = qusage::read.gmt(paste0(gene_sets.dir, "VEGFRi_Signature/_genesets/", gmt))
  # print(paste0(names(sig), ": ", length(sig[[1]])))
  VEGFRi_signatures.list = c(VEGFRi_signatures.list, sig)
}

DAPPER_VEGFRi_LMS =
  gsva(
    vst.lms.mat,
    VEGFRi_signatures.list,
    method = "gsva",
    kcdf = "Gaussian",   # Appropriate for vst transformed data
    min.sz = 2,          # Minimum gene set size
    max.sz = 500,        # Maximum gene set size
    mx.diff = TRUE,      # Compute Gaussian-distributed scores
    verbose = FALSE      # Don't print out the progress bar
  )
# DAPPER_VEGFRi_LMS %>% write.table(., paste0(data.lms.dir, "IO_signatures/DAPPER_VEGFRi_Sigs.lms.tsv"), quote = F, sep = "\t", row.names = T)
# DAPPER_VEGFRi_LMS = read.table(paste0(data.lms.dir, "IO_signatures/DAPPER_VEGFRi_Sigs.lms.tsv"), header = T, sep = "\t")

DAPPER_LU_TUMOR_LMS =
  DAPPER_VEGFRi_LMS %>%
  filter(str_detect(rownames(.), "ENDOTHELIAL_MARKERS"))

## t() and merge results
DAPPER_GSVA_IO_VEGFRi_Sigs.LMS =
  full_join(DAPPER.IO_Sigs.lms %>% t() %>% as.data.frame() %>% rownames_to_column(var="Samples"),
            DAPPER_VEGFRi_LMS %>% t() %>% as.data.frame() %>% rownames_to_column(var="Samples"),
            by=c("Samples"="Samples"))

DAPPER_GSVA_IO_LU_TUMOR_Sigs.LMS =
  full_join(DAPPER.IO_Sigs.lms %>% t() %>% as.data.frame() %>% rownames_to_column(var="Samples"),
            DAPPER_LU_TUMOR_LMS %>% t() %>% as.data.frame() %>% rownames_to_column(var="Samples"),
            by=c("Samples"="Samples"))


# --------------------------------------------

# heatmap all

## get annot_df
annot_df.lms =
  clin_info2 %>%
  dplyr::select(RNAseqID, os_bool, pfs_bool, TreatmentArm, Age, Gender, Race, RECIST_bool) %>%
  arrange(os_bool) %>%
  column_to_rownames(var="RNAseqID")
# annot_df.lms

DAPPER_GSVA_IO_VEGFRi_Sigs.LMS.Heatmap =
  DAPPER_GSVA_IO_VEGFRi_Sigs.LMS %>%
  arrange(match(Samples, rownames(annot_df.lms))) %>%
  column_to_rownames(var="Samples") %>%
  # dplyr::select(PredictIO, IFNG, KDM5A, STAT1, LU_TUMOR_ENDOTHELIAL_MARKERS_DN, LU_TUMOR_VASCULATURE_DN) %>%
  as.matrix() %>% t()
# DAPPER_GSVA_IO_VEGFRi_Sigs.LMS.Heatmap[1:3,1:3]

fig.gsva.dir = ".../VIGEX,DAPPER/LMS/figures/GSVA_IO/"
pdf(paste0(fig.gsva.dir, "DAPPER_GSVA_IO_VEGFRi_Sigs.LMS.Heatmap.scale_none.pdf"), width=8, height=6, useDingbats = F)
pheatmap::pheatmap(DAPPER_GSVA_IO_VEGFRi_Sigs.LMS.Heatmap,
                   annotation_col = annot_df.lms,
                   scale = "none",
                   show_colnames = FALSE,
                   angle_col = 90,
                   cluster_rows = FALSE,
                   cluster_cols = FALSE,
                   fontsize_row = 6)
dev.off()


# --------------------------------------------

# correlate with survival

pathways.vec = DAPPER_GSVA_IO_LU_TUMOR_Sigs.LMS %>% column_to_rownames(var="Samples") %>% colnames()

fig.gsva.os.dir = ".../VIGEX,DAPPER/LMS/figures/GSVA_IO/OS/"
dir.create(fig.gsva.os.dir)
fig.gsva.pfs.dir = ".../VIGEX,DAPPER/LMS/figures/GSVA_IO/PFS/"
dir.create(fig.gsva.pfs.dir)

p_cutoff = 1

OSpathway= c(); OSpathway.median = c(); OSpathway.hr = c(); OSpathway.hr.CI5 = c(); OSpathway.hr.CI95 = c(); OSpathway.pvals = c()
PFSpathway= c(); PFSpathway.median = c(); PFSpathway.hr = c(); PFSpathway.hr.CI5 = c(); PFSpathway.hr.CI95 = c(); PFSpathway.pvals = c()
for (pathway in pathways.vec){

  # pathway = pathways.vec[1]; pathway
  # pathway = "TIDE"

  pathway.df = DAPPER_GSVA_IO_VEGFRi_Sigs.LMS %>% dplyr::select(Samples, all_of(pathway))

  if (pathway.df %>% filter(is.na(!!as.symbol(pathway))) %>% nrow() == 0){
    print(paste0("Pathway: ", pathway))

    colnames(pathway.df) = c("Sample", "Score")
    pathway_median = pathway.df$Score %>% median; pathway_median

    ## OS ############################################
    pathway_clinInfo.os =
      left_join(pathway.df, clin_info, by=c("Sample"="RNAseqID")) %>%
      dplyr::select(Sample, Score, os_time, os_status) %>%
      dplyr::mutate(Score_bool = if_else(Score > pathway_median, "Score_high",
                                         if_else(Score <= pathway_median, "Score_low", "ERROR")))
    pathway_clinInfo.os

    ## get pval and HR
    pathway_clinInfo.os.coxphfit = coxph(Surv(os_time, os_status) ~ Score_bool, data = pathway_clinInfo.os)
    pathway_clinInfo.os.coxphfit.pval = summary(pathway_clinInfo.os.coxphfit)$sctest['pvalue'] %>% unname(); pathway_clinInfo.os.coxphfit.pval
    pathway_clinInfo.os.coxphfit.hr = exp(pathway_clinInfo.os.coxphfit$coefficients) %>% unname(); pathway_clinInfo.os.coxphfit.hr
    pathway_clinInfo.os.coxphfit.hr.CI5 = confint(pathway_clinInfo.os.coxphfit)[1] %>% unname(); pathway_clinInfo.os.coxphfit.hr.CI5
    pathway_clinInfo.os.coxphfit.hr.CI95 = confint(pathway_clinInfo.os.coxphfit)[2] %>% unname(); pathway_clinInfo.os.coxphfit.hr.CI95

    OSpathway.median = c(OSpathway.median, setNames(pathway_median, pathway)); OSpathway.median
    OSpathway.hr = c(OSpathway.hr, setNames(pathway_clinInfo.os.coxphfit.hr, pathway)); OSpathway.hr
    OSpathway.hr.CI5 = c(OSpathway.hr.CI5, setNames(pathway_clinInfo.os.coxphfit.hr.CI5, pathway)); OSpathway.hr.CI5
    OSpathway.hr.CI95 = c(OSpathway.hr.CI95, setNames(pathway_clinInfo.os.coxphfit.hr.CI95, pathway)); OSpathway.hr.CI95
    OSpathway.pvals = c(OSpathway.pvals, setNames(pathway_clinInfo.os.coxphfit.pval, pathway)); OSpathway.pvals

    ## plot
    if (pathway_clinInfo.os.coxphfit.pval < p_cutoff){
      OSpathway = c(OSpathway, pathway)
      pathway_clinInfo.os.fit = survfit(Surv(os_time, os_status) ~ Score_bool, data = pathway_clinInfo.os)
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
        labs(title = paste0("DAPPER LMS GSVA OS\npathway = ", pathway))

      # ggsave(paste0(fig.gsva.os.dir, "DAPPER_LMS_GSVA_OS_", pathway, ".pdf"),
      #        device = "pdf", width = 5, height = 4, units = "in", useDingbats=F)
    }

    ## PFS ############################################
    pathway_clinInfo.pfs =
      left_join(pathway.df, clin_info, by=c("Sample"="RNAseqID")) %>%
      dplyr::select(Sample, Score, pfs_time, pfs_status) %>%
      dplyr::mutate(Score_bool = if_else(Score > pathway_median, "Score_high",
                                         if_else(Score <= pathway_median, "Score_low", "ERROR")))
    pathway_clinInfo.pfs

    ## get pval and HR
    pathway_clinInfo.pfs.coxphfit = coxph(Surv(pfs_time, pfs_status) ~ Score_bool, data = pathway_clinInfo.pfs)
    pathway_clinInfo.pfs.coxphfit.pval = summary(pathway_clinInfo.pfs.coxphfit)$sctest['pvalue'] %>% unname(); pathway_clinInfo.pfs.coxphfit.pval
    pathway_clinInfo.pfs.coxphfit.hr = exp(pathway_clinInfo.pfs.coxphfit$coefficients) %>% unname(); pathway_clinInfo.pfs.coxphfit.hr
    pathway_clinInfo.pfs.coxphfit.hr.CI5 = confint(pathway_clinInfo.pfs.coxphfit)[1] %>% unname(); pathway_clinInfo.pfs.coxphfit.hr.CI5
    pathway_clinInfo.pfs.coxphfit.hr.CI95 = confint(pathway_clinInfo.pfs.coxphfit)[2] %>% unname(); pathway_clinInfo.pfs.coxphfit.hr.CI95

    PFSpathway.median = c(PFSpathway.median, setNames(pathway_median, pathway)); PFSpathway.median
    PFSpathway.hr = c(PFSpathway.hr, setNames(pathway_clinInfo.pfs.coxphfit.hr, pathway)); PFSpathway.hr
    PFSpathway.hr.CI5 = c(PFSpathway.hr.CI5, setNames(pathway_clinInfo.pfs.coxphfit.hr.CI5, pathway)); PFSpathway.hr.CI5
    PFSpathway.hr.CI95 = c(PFSpathway.hr.CI95, setNames(pathway_clinInfo.pfs.coxphfit.hr.CI95, pathway)); PFSpathway.hr.CI95
    PFSpathway.pvals = c(PFSpathway.pvals, setNames(pathway_clinInfo.pfs.coxphfit.pval, pathway)); PFSpathway.pvals

    ## plot
    if (pathway_clinInfo.pfs.coxphfit.pval < p_cutoff){
      PFSpathway = c(PFSpathway, pathway)
      pathway_clinInfo.pfs.fit = survfit(Surv(pfs_time, pfs_status) ~ Score_bool, data = pathway_clinInfo.pfs)
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
        labs(title = paste0("DAPPER LMS GSVA PFS\npathway = ", pathway))

      # ggsave(paste0(fig.gsva.pfs.dir, "DAPPER_GSVA_LMS_PFS_", pathway, ".pdf"),
      #        device = "pdf", width = 5, height = 4, units = "in", useDingbats=F)
    }
  }
}

# OSpathway %>% dput()
OSpathway = c("LRRC15_CAF", "C_ECM", "B_cell_Budczies", "CD8_Jiang", "TLS", "PredictIO", "LU_TUMOR_ENDOTHELIAL_MARKERS_DN")

# PFSpathway %>% dput()
PFSpathway = c("IPRES", "EMT_hallmark", "LE_EMT", "LU_TUMOR_ENDOTHELIAL_MARKERS_UP",
                   "LU_TUMOR_VASCULATURE_UP", "PID_VEGF_VEGFR_PATHWAY", "REACTOME_VEGF_LIGAND_RECEPTOR_INTERACTIONS",
                   "REACTOME_VEGFR2_MEDIATED_VASCULAR_PERMEABILITY")


# wrangle outputs into df
OSpathway.p.adj = p.adjust(OSpathway.pvals, method = "BH")
PFSpathway.p.adj = p.adjust(sort(PFSpathway.pvals), method = "BH")

OSpathway.all.df =
  plyr::join_all(list(data.frame(OSpathway.median) %>% rownames_to_column(var="pathway"),
                      data.frame(OSpathway.hr) %>% rownames_to_column(var="pathway"),
                      data.frame(OSpathway.hr.CI5) %>% rownames_to_column(var="pathway"),
                      data.frame(OSpathway.hr.CI95) %>% rownames_to_column(var="pathway"),
                      data.frame(OSpathway.pvals) %>% rownames_to_column(var="pathway"),
                      data.frame(OSpathway.p.adj) %>% rownames_to_column(var="pathway")),
                 by="pathway", type = "full") %>%
  arrange(OSpathway.pvals)
# OSpathway.all.df %>% write.table(., paste0(data.lms.dir, "OSpathway.all.df.csv"), quote = F, sep = ",", row.names = F)

PFSpathway.all.df =
  plyr::join_all(list(data.frame(PFSpathway.median) %>% rownames_to_column(var="pathway"),
                      data.frame(PFSpathway.hr) %>% rownames_to_column(var="pathway"),
                      data.frame(PFSpathway.hr.CI5) %>% rownames_to_column(var="pathway"),
                      data.frame(PFSpathway.hr.CI95) %>% rownames_to_column(var="pathway"),
                      data.frame(PFSpathway.pvals) %>% rownames_to_column(var="pathway"),
                      data.frame(PFSpathway.p.adj) %>% rownames_to_column(var="pathway")),
                 by="pathway", type = "full") %>%
  arrange(PFSpathway.pvals)
# PFSpathway.all.df %>% write.table(., paste0(data.lms.dir, "PFSpathway.all.df.csv"), quote = F, sep = ",", row.names = F)

# --------------------------------------------

# heatmap OS significant

## get annot_df
annot_df.lms =
  clin_info2 %>%
  dplyr::select(RNAseqID, os_bool, pfs_bool, TreatmentArm, Age, Gender, Race, RECIST_bool) %>%
  arrange(os_bool) %>%
  column_to_rownames(var="RNAseqID")
# annot_df.lms

DAPPER_GSVA_IO_VEGFRi_Sigs.LMS.Heatmap.os =
  DAPPER_GSVA_IO_VEGFRi_Sigs.LMS %>%
  arrange(match(Samples, rownames(annot_df.lms))) %>%
  column_to_rownames(var="Samples") %>%
  dplyr::select(all_of(OSp005pathway)) %>%
  as.matrix() %>% t()
# DAPPER_GSVA_IO_VEGFRi_Sigs.LMS.Heatmap[1:3,1:3]

fig.gsva.dir = ".../VIGEX,DAPPER/LMS/figures/GSVA_IO/"
pdf(paste0(fig.gsva.dir, "DAPPER_GSVA_IO_VEGFRi_Sigs.LMS.Heatmap.OSsigPaths.scale_row.pdf"),
    width=7, height=7, useDingbats = F)
pheatmap::pheatmap(DAPPER_GSVA_IO_VEGFRi_Sigs.LMS.Heatmap.os,
                   annotation_col = annot_df.lms,
                   scale = "none",
                   # scale = "row",
                   show_colnames = FALSE,
                   angle_col = 90,
                   cluster_rows = FALSE, cellwidth = 12, cellheight = 12, fontsize = 8,
                   cluster_cols = FALSE,
                   fontsize_row = 6)
dev.off()

# EOF