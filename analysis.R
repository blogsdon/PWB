synapser::synLogin()
foo <- data.table::fread('ROSMAP_PiL_dep_cog.csv',data.table=F)
View(foo)

#get design matrix
dmobj <- synapser::synGet('syn8456640')
dm <- data.table::fread(dmobj$path,data.table=F)

#get covariates
covobj <- synapser::synGet('syn3191087')
cova <- data.table::fread(covobj$path,data.table=F)
cova <- cova[,-1]
cova <- dplyr::select(cova,projid,cogdx)

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
foo <- dplyr::select(foo,Sampleid,projid,purpose_total,fu_year)
foo <- dplyr::left_join(foo,dm,by=c('Sampleid'='SampleID'))
foo <- dplyr::left_join(foo,cova)
foo <- dplyr::left_join(foo,countMatrix,by=c('Sampleid'='sampleId'))


run_mixed_effects <- function(geneId,foo,pb){

  str <- paste0("purpose_total ~ ",geneId," + cogdx + Batch0 + Batch1 + Batch2 + Batch3 + Batch4 + Batch5 + Batch6 + Batch7 + msex0 + RINcontinuous + PCT_CODING_BASES + PCT_INTERGENIC_BASES + pmi + age_death + fu_year + (1|Sampleid)")
  #print(str)
  str <- as.formula(str)
  pb$tick()
  res<-lme4::lmer(str,foo)
  return(summary(res)$coef[2,3])
}

sum1 <- summary(lme4::lmer(purpose_total ~ + cogdx + Batch0 + Batch1 + Batch2 + Batch3 + Batch4 + Batch5 + Batch6 + Batch7 + msex0 + RINcontinuous + PCT_CODING_BASES + PCT_INTERGENIC_BASES + pmi + age_death + fu_year+ (1|Sampleid),foo))

summary(lme4::lmer(purpose_total ~ ENSG00000227232 +  cogdx + Batch0 + Batch1 + Batch2 + Batch3 + Batch4 + Batch5 + Batch6 + Batch7 + msex0 + RINcontinuous + PCT_CODING_BASES + PCT_INTERGENIC_BASES + pmi + age_death + fu_year+ (ENSG00000227232||Sampleid),foo))

gene_ids <- colnames(countMatrix)
gene_ids <- gene_ids[-which(gene_ids=='sampleId')]

no_cores <- parallel::detectCores()
cl <- parallel::makeCluster(no_cores)
pb <- progress::progress_bar$new(total=(length(gene_ids)))
system.time(mixed_effect_results <- sapply(gene_ids,run_mixed_effects,foo,pb))

extract_last_followup <- function(x,foo){
  focus <- dplyr::filter(foo,projid==x)
  print(x)
  return(focus[which.max(focus$fu_year),])
}
unique_projid <- unique(foo$projid)
foo23<-sapply(unique_projid,extract_last_followup,foo)

pvals <- pt(abs(mixed_effect_results),4724,lower.tail=F)*2
gap::qqunif(pvals)


run_mixed_effects2 <- function(geneId,foo,pb){

  str <- paste0("purpose_total ~ ",geneId," + Batch0 + Batch1 + Batch2 + Batch3 + Batch4 + Batch5 + Batch6 + Batch7 + msex0 + RINcontinuous + PCT_CODING_BASES + PCT_INTERGENIC_BASES + pmi + age_death + fu_year + (1|Sampleid)")
  #print(str)
  str <- as.formula(str)
  pb$tick()
  res<-lme4::lmer(str,foo)
  return(summary(res)$coef[2,3])
}
#parallel::stopCluster(cl)
#cl <- parallel::makeCluster(no_cores)
pb <- progress::progress_bar$new(total=(length(gene_ids)))
system.time(mixed_effect_results2 <- sapply(gene_ids,run_mixed_effects2,foo,pb))


run_mixed_effects3 <- function(geneId,foo,pb){

  str <- paste0("purpose_total ~ ",geneId," + cogdx + Batch0 + Batch1 + Batch2 + Batch3 + Batch4 + Batch5 + Batch6 + Batch7 + msex0 + RINcontinuous + PCT_CODING_BASES + PCT_INTERGENIC_BASES + pmi + age_death + (1|Sampleid)")
  #print(str)
  str <- as.formula(str)
  pb$tick()
  res<-lme4::lmer(str,foo)
  return(summary(res)$coef[2,3])
}
pb <- progress::progress_bar$new(total=(length(gene_ids)))
system.time(mixed_effect_results3 <- sapply(gene_ids,run_mixed_effects3,foo,pb))
