

knitr::opts_knit$set(root.dir = normalizePath("/Users/joshuagabbay/Desktop/QBIO/qbio_490_joshuagabbay/final_project")) 

if (!require(survival)){ 
  install.packages("survival")
}
library(survival)

if (!require(survminer)){
  install.packages("survminer")
}
library(survminer)

if (!require(ggplot2)){ 
  install.packages("ggplot2")
}
library(ggplot2)

if (!require(TCGAbiolinks)){ 
  install.packages("TCGAbiolinks")
}
library(TCGAbiolinks)

if (!require(BiocManager)){
  install.packages("BiocManager")
}
library(BiocManager)


#GDC on clinical data
clinical_query <- GDCquery(project = "TCGA-GBM", data.category = "Clinical", file.type = "xml" )
GDCdownload(clinical_query)
clinical <- GDCprepare_clinic(clinical_query, clinical.info = "patient")
gender_na_mask <- ifelse(is.na(clinical$gender), F, T)
gender_cleaned_clinical <- clinical[gender_na_mask, ]

#making a survival time column for survival plots
gender_cleaned_clinical$survival_time <- ifelse(is.na(gender_cleaned_clinical$days_to_death), gender_cleaned_clinical$survival_time <- gender_cleaned_clinical$days_to_last_followup, gender_cleaned_clinical$survival_time <- gender_cleaned_clinical$days_to_death)

#making a death event (T/F) column for survival plots
gender_cleaned_clinical$death_event <- ifelse(gender_cleaned_clinical$vital_status == "Alive", gender_cleaned_clinical$death_event <- FALSE, gender_cleaned_clinical$death_event <- TRUE)
#initializing a survival object
surv_object_gender <- Surv(time = gender_cleaned_clinical$survival_time, event = gender_cleaned_clinical$death_event)

#creating a fit object HERE
gender_fit <- survfit(surv_object_gender ~ gender_cleaned_clinical$gender, data = gender_cleaned_clinical)

#formats and creates KM plot
survplot_gender = ggsurvplot(gender_fit, pval=TRUE, ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), legend = "right")

jpeg("/Users/joshuagabbay/Desktop/QBIO/qbio_490_joshuagabbay/final_project/gender_KM_plot.jpg")
KM_plot_gender = survplot_gender$plot + theme_bw() + theme(axis.title = element_text(size=20), axis.text = element_text(size=16), legend.title = element_text(size=14), legend.text = element_text(size=12))
KM_plot_gender
dev.off()

age_na_mask <- ifelse(is.na(clinical$age_at_initial_pathologic_diagnosis), F, T)
age_cleaned_clinical <- clinical[age_na_mask, ]
young_mask <- ifelse(age_cleaned_clinical$age_at_initial_pathologic_diagnosis <= 35, T, F)
middle_mask <- ifelse(age_cleaned_clinical$age_at_initial_pathologic_diagnosis >= 35 & age_cleaned_clinical$age_at_initial_pathologic_diagnosis <= 50, T, F)
old_mask <- ifelse(age_cleaned_clinical$age_at_initial_pathologic_diagnosis >= 35, T, F)
age_cleaned_clinical$age_status <- ifelse(young_mask, "Young", ifelse(middle_mask, "Middle", "Old"))

#making a survival time column for survival plots
age_cleaned_clinical$survival_time <- ifelse(is.na(age_cleaned_clinical$days_to_death), age_cleaned_clinical$survival_time <- age_cleaned_clinical$days_to_last_followup, age_cleaned_clinical$survival_time <- age_cleaned_clinical$days_to_death)

#making a death event (T/F) column for survival plots
age_cleaned_clinical$death_event <- ifelse(age_cleaned_clinical$vital_status == "Alive", age_cleaned_clinical$death_event <- FALSE, age_cleaned_clinical$death_event <- TRUE)
#initializing a survival object
surv_object_age <- Surv(time = age_cleaned_clinical$survival_time, event = age_cleaned_clinical$death_event)

#creating a fit object HERE
age_fit <- survfit(surv_object_age ~ age_cleaned_clinical$age_status, data = age_cleaned_clinical)

#formats and creates KM plot
survplot_age = ggsurvplot(age_fit, pval=TRUE, ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), legend = "right")

jpeg("/Users/joshuagabbay/Desktop/QBIO/qbio_490_joshuagabbay/final_project/age_KM_plot.jpg")
KM_plot_age = survplot_age$plot + theme_bw() + theme(axis.title = element_text(size=20), axis.text = element_text(size=16), legend.title = element_text(size=14), legend.text = element_text(size=12))
KM_plot_age
dev.off()

#cleaning race data
race_na_mask <- ifelse(clinical$race_list == "", F, T)
race_cleaned_clinical <- clinical[race_na_mask, ]
race_mask <- ifelse(race_cleaned_clinical$race_list == "BLACK OR AFRICAN AMERICAN", T, F)
sum(race_mask)
#making a survival time column for survival plots
race_cleaned_clinical$survival_time <- ifelse(is.na(race_cleaned_clinical$days_to_death), race_cleaned_clinical$survival_time <- race_cleaned_clinical$days_to_last_followup, race_cleaned_clinical$survival_time <- race_cleaned_clinical$days_to_death)

#making a death event (T/F) column for survival plots
race_cleaned_clinical$death_event <- ifelse(race_cleaned_clinical$vital_status == "Alive", race_cleaned_clinical$death_event <- FALSE, race_cleaned_clinical$death_event <- TRUE)

#initializing a survival object
surv_object_race <- Surv(time = race_cleaned_clinical$survival_time, event = race_cleaned_clinical$death_event)

#creating a fit object HERE
race_fit <- survfit(surv_object_race ~ race_cleaned_clinical$race_list, data = race_cleaned_clinical)

#formats and creates KM plot
survplot_race = ggsurvplot(race_fit, pval=TRUE, ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), legend = "right")

jpeg("/Users/joshuagabbay/Desktop/QBIO/qbio_490_joshuagabbay/final_project/race_KM_plot.jpg")
KM_plot_race = survplot_race$plot + theme_bw() + theme(axis.title = element_text(size=20), axis.text = element_text(size=16), legend.title = element_text(size=14), legend.text = element_text(size=8))
KM_plot_race
dev.off()



BiocManager::install("DESeq2")
library(DESeq2)
BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)
library(SummarizedExperiment)


#GDC on transcriptomic data
rna_query <- GDCquery(project = "TCGA-GBM",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "STAR - Counts")

GDCdownload(rna_query)
rna_se <- GDCprepare(rna_query)

#extracting clinical, genes, and rna counts dataframes and cleaning vital status information
mask <- is.na(rna_se@colData$vital_status)
rna_clinical <- rna_se@colData[!mask, ]
new_mask <- ifelse(rna_clinical$vital_status == "Not Reported", F, T)
rna_clinical <- rna_clinical[new_mask, ]
rna_clinical <- as.data.frame(rna_clinical)

rna_genes <- rna_se@rowRanges@elementMetadata
rna_genes <- as.data.frame(rna_genes)
rownames(rna_genes) = rna_genes$gene_id
rna_counts <- rna_se@assays@data$unstranded
rna_counts <- rna_counts[ , !mask]
rna_counts <- rna_counts[, new_mask]
rownames(rna_counts) = rna_genes$gene_id
colnames(rna_counts) = rownames(rna_clinical)

rna_clinical$vital_status <- factor(rna_clinical$vital_status)

sum(is.na(rna_clinical$vital_status)) #make sure = 0

#running DESeq2
row_sums <- rowSums(rna_counts)
low_counts_mask <- ifelse(row_sums < 10, F, T)

rna_counts <- rna_counts[low_counts_mask,]
ncol(rna_counts)
nrow(rna_clinical)
dds = DESeqDataSetFromMatrix(countData = rna_counts,
                             colData = rna_clinical,
                             design = ~ vital_status)
dds_obj = DESeq(dds)

head(rna_clinical$vital_status)

results = results(dds_obj, format = "DataFrame", contrast = c("vital_status", "Alive", "Dead"))

#Formatting results
gene_mask <- ifelse(rna_genes$gene_id %in% results@rownames, T, F)
gene_name <- rna_genes[gene_mask, 7]
head(rna_genes)
results <- data.frame(gene_name, results@rownames, results$log2FoldChange, results$pvalue, results$padj, -log10(results$padj))
colnames(results) <- c("Gene_Name", "Gene_ID", "Log2FoldChange", "pvalue", "padj", "-log10(padj)")
par(mar=c(1,1,1,1))
rownames(results) <- results$Gene_ID
jpeg("/Users/joshuagabbay/Desktop/QBIO/qbio_490_joshuagabbay/final_project/enhanced_volcano_survivers.jpg")
EnhancedVolcano(results, lab = '', x = 'Log2FoldChange', y = 'pvalue', shape = 8, colAlpha = 1, labSize = 3)
dev.off()

