# Loading package ---------------------------------------------------------
library(tidyverse)
library(TwoSampleMR)
library(RadialMR)
library(meta)
library(phenoscanner)

# RA → CTS MR ------------------------------------------------------------------
RA = "ebi-a-GCST90013534"
CTS_finn = "finn-b-G6_CARPTU"

#Exposure data&Outcome data
exposure_data = extract_instruments(outcomes = RA)
outcome_data = extract_outcome_data(snps = exposure_data$SNP, outcomes = CTS_finn)

#harmoning
H_data = harmonise_data(exposure_dat = exposure_data,
                        outcome_dat = outcome_data)

#Remove outlier
egger_radial(
  r_input = H_data,
  alpha = 0.05,
  weights = 1,
  summary = TRUE
)
ivw_radial(
  r_input = H_data,
  alpha = 0.05,
  weights = 1,
  summary = TRUE
)

ivw.model = ivw_radial(H_data, 0.05, 1)
egger.model = egger_radial(H_data, 0.05, 1)


a1 = cbind(ivw.model$outliers[1])
a2 = cbind(egger.model$outliers[1])
a1 = as.matrix(a1)
a2 = as.matrix(a2)
a3 = rbind(a1, a2)
a3 = as.data.frame(a3)
a3 = unique(a3)

#IVWplot
pdf("IVWplot.pdf", width = 8, height = 5)
IVWplot1 = plot_radial(ivw.model, T, F, F)
IVWplot1
dev.off()

#eggerplot
pdf("eggerplot.pdf", width = 8, height = 5)
IVWplot1 = plot_radial(egger.model, T, F, F)
IVWplot1
dev.off()

#remove outlier
exposure_data = exposure_data %>%
  filter(!SNP %in% a3$SNP)
H_data = harmonise_data(
  exposure_dat = exposure_data,
  outcome_dat = outcome_data) %>%
  filter(mr_keep == T)

#MR analysis
mr_results_finn = generate_odds_ratios(mr(H_data))

#write results
write.csv(mr_results_finn, 'RA_CTS_MR_result_finn.csv')
write.csv(H_data, "RA_CTS_MR_Hdata_finn.csv")

# Phenoscanner ------------------------------------------------------------
# SNP
res = phenoscanner(snpquery=c(H_data$SNP))
write.csv(res$results,"snp_pheno.csv")


# Sensitivity Analysis -------------------------------------------------------------------
#Horizontal pleiotropy
write.csv(mr_pleiotropy_test(H_data),"RA_CTS_ple.csv")

#Heterogeneity statistics
write.csv(mr_heterogeneity(H_data, method_list=c("mr_egger_regression", "mr_ivw")),"RA_CTS_hete.csv")

#MR-PRESSO
presso = run_mr_presso(H_data,NbDistribution = 10000)
write.csv(presso,"mr_presso.csv")

# Plotting ----------------------------------------------------------------
pdf("scatter_plot.pdf",width=8,height=5)
plot1 = mr_scatter_plot(mr_results, H_data)
plot1
dev.off()

res_single = mr_singlesnp(H_data)
pdf("mr_forest_plot.pdf",width=8,height=5)
plot2 = mr_forest_plot(res_single)
plot2
dev.off()

pdf("mr_leaveoneout_plot.pdf",width=8,height=7)
res_loo = mr_leaveoneout(H_data)
plot3 = mr_leaveoneout_plot(res_loo)
plot3
dev.off()

pdf("mr_funnel_plot.pdf",width=6,height=4)
plot4 = mr_funnel_plot(res_single)
plot4
dev.off()

# Calculate R2&F ----------------------------------------------------------
eaf2maf(H_data$eaf.exposure)
tmp = H_data %>% 
  mutate(MAF = eaf2maf(H_data$eaf.exposure)) %>% 
  mutate(R2 = 2*(1-MAF)*MAF*beta.exposure/(se.exposure*sqrt(samplesize.exposure))) %>% 
  mutate(F = (samplesize.exposure-2)*R2/(1-R2))
sum_R2 = sum(tmp$R2)
sam_size = max(tmp$samplesize.exposure)
F_sta = sum_R2/(1-sum_R2)*(sam_size-nrow(H_data)-1)/(nrow(H_data))
write.csv(F_sta,"RA_CTS_F.csv")

# MVMR correcting for BMI  -----------------------------------------------------
BMI = "ieu-a-2"
id_exposure = c(BMI,RA) 
id_outcome = CTS
exposure_data = mv_extract_exposures(id_exposure)
dim(exposure_data)
exposure_data = exposure_data %>% 
  filter(!SNP %in% a3$SNP)
dim(exposure_data)
outcome_data = extract_outcome_data(exposure_data$SNP, id_outcome)
H_data = harmonise_data(
  exposure_dat = exposure_data, 
  outcome_dat = outcome_data)
mvdat = mv_harmonise_data(exposure_data,  outcome_data) # 对数据进行harmonisation
res = mv_multiple(mvdat) 
res1 = generate_odds_ratios(res$result)
write.csv(res1,"RA_BMI_MVMR.csv")



# Mediation analysis ------------------------------------------------------
#via BMI
#total effect
mr_total_effect_beta = RA_CTS_beta = mr_results_finn$b[3]
mr_total_effect_se = RA_CTS_se = mr_results_finn$se[3]
mr_total_effect_z = RA_CTS_z = RA_CTS_beta/RA_CTS_se
mr_total_effect_p = 2*pnorm(-abs(RA_CTS_z))

#direct effect
mr_direct_effect_beta= res$result$b[1]
mr_direct_effect_se= res$result$se[1]
mr_direct_effect_z = mr_direct_effect_beta/mr_direct_effect_se
mr_direct_effect_p = 2*pnorm(-abs(mr_direct_effect_z))

#indirect effect
mr_indirect_effect_beta =mr_total_effect_beta - mr_direct_effect_beta
mr_indirect_effect_se = (mr_total_effect_se^2 + mr_direct_effect_se^2)^(1/2)
mr_indirect_effect_z = mr_indirect_effect_beta/mr_indirect_effect_se
mr_indirect_effect_p = 2*pnorm(-abs(mr_indirect_effect_z))

#CI
lo_ci = mr_indirect_effect_beta - 1.96 * mr_indirect_effect_se
up_ci = mr_indirect_effect_beta + 1.96 * mr_indirect_effect_se

or = exp(mr_indirect_effect_beta)
or_lci95 = exp(lo_ci)
or_uci95 = exp(up_ci)

#via RA
exposure_data = extract_instruments(outcomes = BMI)#替换成你想做的病
outcome_data = extract_outcome_data(snps = exposure_data$SNP, outcomes = CTS)#替换成你想做的病
H_data = harmonise_data(
  exposure_dat = exposure_data, 
  outcome_dat = outcome_data)

mr_results=mr(H_data)  
mr_results
#total effect
mr_total_effect_beta = bmi_CTS_beta = mr_results$b[3]
mr_total_effect_se = bmi_CTS_se = mr_results$se[3]
mr_total_effect_z = bmi_CTS_z = bmi_CTS_beta/bmi_CTS_se
mr_total_effect_p = 2*pnorm(-abs(bmi_CTS_z))

#direct effect
mr_direct_effect_beta= res$result$b[2]
mr_direct_effect_se= res$result$se[2]
mr_direct_effect_z = mr_direct_effect_beta/mr_direct_effect_se
mr_direct_effect_p = 2*pnorm(-abs(mr_direct_effect_z))

#indirect effect
mr_indirect_effect_beta =mr_total_effect_beta - mr_direct_effect_beta
mr_indirect_effect_se = (mr_total_effect_se^2 + mr_direct_effect_se^2)^(1/2)
mr_indirect_effect_z = mr_indirect_effect_beta/mr_indirect_effect_se
mr_indirect_effect_p = 2*pnorm(-abs(mr_indirect_effect_z))

#CI
lo_ci = mr_indirect_effect_beta - 1.96 * mr_indirect_effect_se
up_ci = mr_indirect_effect_beta + 1.96 * mr_indirect_effect_se
or = exp(mr_indirect_effect_beta)
or_lci95 = exp(lo_ci)
or_uci95 = exp(up_ci)


# CTS → RA MR ---------------------------------------------------------------
exposure_data = extract_instruments(outcomes = "finn-b-G6_CARPTU")
outcome_data = extract_outcome_data(snps = exposure_data$SNP, outcomes = "ebi-a-GCST90013534")
H_data = harmonise_data(
  exposure_dat = exposure_data, 
  outcome_dat = outcome_data)
#MR result
mr_results2=generate_odds_ratios(mr(H_data))  
mr_results2
write.csv(mr_results2,"CTS RA MR.csv")

# Meta analysis -----------------------------------------------------------
finn_H_data = read.csv("RA_CTS_MR_Hdata_finn.csv")

exposure_data = extract_instruments(outcomes = RA)
exposure_data = exposure_data %>% 
  filter(SNP %in% finn_H_data$SNP)
outcome_data = read_outcome_data(
  snps = exposure_data$SNP,
  filename = "outcome_data_ukb.csv",
  sep = ",",
  snp_col = "name",
  beta_col = "beta",
  pval_col = "pval",
  eaf_col = "minor_AF",
  se_col = "se",
  effect_allele_col = "minor_allele",
  other_allele_col = "A2")
H_data = harmonise_data(
  exposure_dat = exposure_data, 
  outcome_dat = outcome_data)

mr_results_ukb= generate_odds_ratios(mr(H_data))  
write.csv(mr_results_ukb,"mr_results_ukb.csv")


mete_data = rbind(mr_results_finn[3,12:14],mr_results_ukb[3,12:14])
lnor= log(mete_data[,"or"])
lnuci= log(mete_data[,"or_uci95"])
lnlci= log(mete_data[,"or_lci95"])
selnor= (lnuci-lnlci)/(2*1.96)
pfs=metagen(lnor,selnor,sm="OR",data=mete_data,studlab=c("finn","ukb"))
forest(pfs)
summary(pfs)
