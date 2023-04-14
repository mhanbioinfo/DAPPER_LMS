## Fig4A and B

library("tidyverse")
library("ggplot2")
'%ni%' <- Negate('%in%')
library("survival")
library("survminer")

proj.dir = ".../VIGEX,DAPPER/LMS/"
clin.dir = ".../VIGEX,DAPPER/LMS/_clin_info/"
data.dir = ".../VIGEX,DAPPER/DAPPER/"
fig.dir = ".../VIGEX,DAPPER/LMS/figures/"

# helper functions
give.n <- function(x){
  return(c(y = mean(x)*1.15, label = length(x)))
}

# read in data
clin_info_DAPPER =
  readxl::read_xlsx(".../VIGEX,DAPPER/DAPPER/clin_info/DAPPER_sampleinfo_2022_10_21_1242.xlsx", sheet = "lean_LMS", na = "NA") %>%
  dplyr::select(RNAseqID, DAP_ID)

clin_info =
  readxl::read_xlsx(paste0(clin.dir, "Copy of DAPPER_data_Leiomyosarcoma_2022_10_21.xlsx"), sheet = "Response", na = "NA") %>%
  mutate(TreatmentArm = if_else(`Assigned Treatment Arm` == "Cohort A", "D+Olaparib",
                                if_else(`Assigned Treatment Arm` == "Cohort B", "D+Cediranib", "ERROR"))) %>%
  right_join(., clin_info_DAPPER, by=c("Subject"="DAP_ID")) %>%
  dplyr::rename("RECIST_bool"="Best_Response_RECISTv1") %>%
  dplyr::select(RNAseqID, Subject, TreatmentArm, Diagnosis, Age, Gender, Race, os_time, os_status, pfs_time, pfs_status, RECIST_bool)
# clin_info$os_time %>% median() # [1] 14.58727
# clin_info$pfs_time %>% median() # [1] 2.792608

clin_info2 =
  clin_info %>%
  mutate(os_bool = case_when(os_time >= 14 ~ "os_ge14m",
                             os_time < 14  ~ "os_lt14m"),
         pfs_bool = case_when(pfs_time >= 3 ~ "pfs_ge3m",
                              pfs_time < 3 ~  "pfs_lt3m"))

res_cibersort_abs =
  read.table(".../VIGEX,DAPPER/DAPPER/cibersortx_online/CIBERSORTx_Job12_Results.txt", header = T, sep = "\t") %>%
  dplyr::select(-c(P.value, Correlation, RMSE)) %>%
  rename_with(.cols = -Mixture, .fn = ~ gsub("^", "ciber_", .x))

res_cibersort_abs.clin_info =
  clin_info2 %>%
  left_join(., res_cibersort_abs, by = c("RNAseqID" = "Mixture")) %>%
  dplyr::select(RNAseqID, RECIST_bool, os_bool, pfs_bool, Diagnosis, TreatmentArm, everything()) %>%
  filter(!is.na(RNAseqID))

plot_surv = function(save.dir, outcome, group, TIL, pthreshold){
  # outcome = "PFS_DAYS"
  # group = "all_indications"
  # TIL = "Absolute.score..sig.score"
  # pthreshold = 0.05
  fit = survfit(Surv(DAYS, STATUS) ~ TIL_hi, data = cibersort_abs_survival_df2)

  if (surv_pvalue(fit)$pval < pthreshold) {

      # print(paste0("TIL is: ",TIL,"; p-value is: ",surv_pvalue(fit)$pval))
      survp =
        ggsurvplot(fit,
                   data = cibersort_abs_survival_df2,
                   conf.int = FALSE,
                   pval = TRUE,
                   fun = "pct",
                   size = 1.5,
                   font.x = c(24, "plain", "black"),
                   font.y = c(24, "plain", "black"),
                   font.tickslab = c(20, "plain", "black"),
                   linetype = "strata",
                   palette = c("#E7B800","#2E9FDF"),
                   legend = "top",
                   legend.title = "CIBERSORT\nscore",
                   legend.labs = c("Low","High"),
                   risk.table = TRUE,
                   risk.table.height = 0.35,
                   risk.table.title="")
      survp = survp + labs(title = paste0(outcome, "_", group, "\n", TIL))
      # survp

      survp$plot +
        theme(legend.title = element_text(size = 25, color = "black"),
              legend.text = element_text(size = 25, color = "black"),
              axis.text.x = element_text(size = 22, color = "black"),
              axis.text.y = element_text(size = 22, color = "black"),
              axis.title.x = element_text(size = 25, color = "black"),
              axis.title.y = element_text(size = 25, color = "black"), )

      survp$table$layers[[1]]$aes_params$size = 7; survp
      survp$table$theme =
        theme_classic() +
        theme(text = element_text(size = 25),
              axis.text.y = element_text(colour = c("#2E9FDF", "#E7B800")),
              plot.title = element_text(size = 1, color = "white"))
              ## labs(title) repeats in risk table, manually hide it

      # pdf(paste0(save.dir, "CIBERSORTabs_survival_", outcome, "_", group, "_", TIL, ".pdf"))
      # print(survp, newpage = FALSE)
      # dev.off()

      return(survp)
  }
}


# M1

group = "all"
pthreshold = 1

## OS

TIL = "ciber_Macrophages.M1"
cibersort_abs_survival_df =
  res_cibersort_abs.clin_info %>%
  select(RNAseqID, Diagnosis, TreatmentArm, !!TIL, os_time, os_status) %>%
  filter(!is.na(os_status)) %>%
  rename("DAYS"="os_time", "STATUS"="os_status")
TIL_median = cibersort_abs_survival_df[[TIL]] %>% median()
# print(paste0("TIL: ", TIL, ", median: ", TIL_median)) # [1] 0.0291587775429776

cibersort_abs_survival_df2 =
  cibersort_abs_survival_df %>%
  mutate(TIL_hi = if_else(!!as.symbol(TIL) >= TIL_median, 1, 0))
# cibersort_abs_survival_df2

outcome = "OS_DAYS"
surv_plt = plot_surv(fig.dir.os, outcome, group, TIL, pthreshold)
# surv_plt$plot
fig.dir.cibersort = "/Users/minghan/bioinfoproj/VIGEX,DAPPER/LMS/figures/CIBERSORT_clin/"
pdf(paste0(fig.dir.cibersort, "CIBERSORTabs_survival_", outcome, "_", group, "_", TIL, ".pdf"))
print(surv_plt, newpage = FALSE)
dev.off()

## PFS

TIL = "ciber_Macrophages.M1"
cibersort_abs_survival_df =
  res_cibersort_abs.clin_info %>%
  select(RNAseqID, Diagnosis, TreatmentArm, !!TIL, pfs_time, pfs_status) %>%
  filter(!is.na(pfs_status)) %>%
  rename("DAYS"="pfs_time", "STATUS"="pfs_status")
TIL_median = cibersort_abs_survival_df[[TIL]] %>% median()
print(paste0("TIL: ", TIL, ", median: ", TIL_median)) # [1] 0.0291587775429776

cibersort_abs_survival_df2 =
  cibersort_abs_survival_df %>%
  mutate(TIL_hi = if_else(!!as.symbol(TIL) >= TIL_median, 1, 0))
# cibersort_abs_survival_df2

outcome = "PFS_DAYS"
surv_plt = plot_surv(fig.dir.os, outcome, group, TIL, pthreshold)
# surv_plt$plot
fig.dir.cibersort = "/Users/minghan/bioinfoproj/VIGEX,DAPPER/LMS/figures/CIBERSORT_clin/"
pdf(paste0(fig.dir.cibersort, "CIBERSORTabs_survival_", outcome, "_", group, "_", TIL, ".pdf"))
print(surv_plt, newpage = FALSE)
dev.off()


# M1/M2

## OS

cibersort_abs_survival_df.M1_M2.os =
  res_cibersort_abs.clin_info %>%
  mutate(M1_M2 = ciber_Macrophages.M1 / ciber_Macrophages.M2) %>%
  dplyr::select(RNAseqID, Diagnosis, TreatmentArm, M1_M2, os_time, os_status) %>%
  filter(!is.na(os_status)) %>%
  dplyr::rename("DAYS"="os_time", "STATUS"="os_status")
TIL_median = cibersort_abs_survival_df.M1_M2.os$M1_M2 %>% median(); TIL_median # 0.06470082

cibersort_abs_survival_df.M1_M2.os =
  cibersort_abs_survival_df.M1_M2.os %>%
  mutate(TIL_hi = if_else(M1_M2 >= TIL_median, 1, 0))
cibersort_abs_survival_df.M1_M2.os

pthreshold = 1
outcome = "OS_DAYS"
group = "all"
surv_plt = plot_surv(fig.dir.cibersort, outcome, group, "M1_M2", pthreshold)
surv_plt$plot
ggsave(paste0(fig.dir.cibersort, "CIBERSORT_M1_M2_median.os_v2.pdf"),
       device = "pdf", width = 6, height = 5, units = "in", useDingbats=FALSE)

### PFS

cibersort_abs_survival_df.M1_M2.pfs =
  res_cibersort_abs.clin_info %>%
  mutate(M1_M2 = ciber_Macrophages.M1 / ciber_Macrophages.M2) %>%
  dplyr::select(RNAseqID, Diagnosis, TreatmentArm, M1_M2, pfs_time, pfs_status) %>%
  filter(!is.na(pfs_status)) %>%
  dplyr::rename("DAYS"="pfs_time", "STATUS"="pfs_status")
cibersort_abs_survival_df.M1_M2.pfs
TIL_median = cibersort_abs_survival_df.M1_M2.pfs$M1_M2 %>% median(); TIL_median # [1] 0.06470082

cibersort_abs_survival_df.M1_M2.pfs =
  cibersort_abs_survival_df.M1_M2.pfs %>%
  mutate(TIL_hi = if_else(M1_M2 >= TIL_median, 1, 0))
cibersort_abs_survival_df.M1_M2.pfs

pthreshold = 1
outcome = "PFS_DAYS"
group = "all"
surv_plt = plot_surv(fig.dir.cibersort, outcome, group, "M1_M2", pthreshold)
surv_plt$plot
ggsave(paste0(fig.dir.cibersort, "CIBERSORT_M1_M2_median.pfs_v2.pdf"),
       device = "pdf", width = 6, height = 5, units = "in", useDingbats=FALSE)

# EOF