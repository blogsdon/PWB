synapser::synLogin()
foo <- data.table::fread('ROSMAP_PiL_dep_cog.csv',data.table=F)
View(foo)

#get design matrix
dmobj <- synapser::synGet('syn8456640')
dm <- data.table::fread(dmobj$path,data.table=F)

#get log cqn
countobj <- synapser::synGet('syn8456638')
countMatrix <- data.table::fread(countobj$path,data.table=F)
rownames(countMatrix) <- countMatrix$ensembl_gene_id
countMatrix <- countMatrix[,-1]
countMatrix <- t(countMatrix)
countMatrix <- data.frame(countMatrix,stringsAsFactors = F)
countMatrix$sampleId <- rownames(countMatrix)

#map sample ids to projids
key <- data.table::fread("key.csv",data.table=F)
foo <- dplyr::left_join(foo,key)
foo <- dplyr::filter(foo,!is.na(Sampleid))
foo <- dplyr::select(foo,Sampleid,purpose_total,fu_year)
foo <- dplyr::left_join(foo,dm,by=c('Sampleid'='SampleID'))
foo <- dplyr::left_join(foo,countMatrix,by=c('Sampleid'='sampleId'))


run_mixed_effects <- function(geneId,foo){

  str <- paste0("purpose_total ~ ",geneId," + DiagnosisAD + DiagnosisCONTROL + DiagnosisOTHER + Batch0 + Batch1 + Batch2 + Batch3 + Batch4 + Batch5 + Batch6 + Batch7 + msex0 + RINcontinuous + PCT_CODING_BASES + PCT_INTERGENIC_BASES + pmi + age_death + fu_year + (1|Sampleid)")
  #print(str)
  str <- as.formula(str)
  return(lme4::lmer(str,foo))
}

summary(lme4::lmer(purpose_total ~ + DiagnosisAD + DiagnosisCONTROL + DiagnosisOTHER + Batch0 + Batch1 + Batch2 + Batch3 + Batch4 + Batch5 + Batch6 + Batch7 + msex0 + RINcontinuous + PCT_CODING_BASES + PCT_INTERGENIC_BASES + pmi + age_death + fu_year+ (1|Sampleid),foo))
#gene_ids <- colnames(countMatrix)
#gene_ids <- gene_ids[-which(gene_ids=='sampleId')]
system.time(mixed_effect_results <- lapply(gene_ids[1:100],run_mixed_effects,foo))

