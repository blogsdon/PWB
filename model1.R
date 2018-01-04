#load in PiL data
foo <- data.table::fread('ROSMAP_PiL_dep_cog.csv',data.table=F)

#restrict to baseline
foo <- dplyr::filter(foo,fu_year==0)

#get technical covariates
dmobj <- synapser::synGet('syn8456640')
dm <- data.table::fread(dmobj$path,data.table=F)

covobj <- synapser::synGet('syn3191087')
cova <- data.table::fread(covobj$path,data.table=F)
cova <- cova[,-1]
cova <- dplyr::select(cova,projid,cogdx)

countobj <- synapser::synGet('syn8456638')
countMatrix <- data.table::fread(countobj$path,data.table=F)
rownames(countMatrix) <- countMatrix$ensembl_gene_id
countMatrix <- countMatrix[,-1]
countMatrix <- t(countMatrix)
countMatrix <- data.frame(countMatrix,stringsAsFactors = F)
countMatrix$sampleId <- rownames(countMatrix)

key <- data.table::fread("key.csv",data.table=F)
foo <- dplyr::left_join(foo,key)
foo <- dplyr::filter(foo,!is.na(Sampleid))
foo <- dplyr::select(foo,Sampleid,projid,purpose_total,fu_year,dcfdx,age_at_visit)
foo <- dplyr::left_join(foo,dm,by=c('Sampleid'='SampleID'))
foo <- dplyr::left_join(foo,cova)
foo <- dplyr::left_join(foo,countMatrix,by=c('Sampleid'='sampleId'))

run_simple_model <- function(pheno,geneId,foo,pb){
  str <- paste0(pheno," ~ ",geneId," + dcfdx + Batch0 + Batch1 + Batch2 + Batch3 + Batch4 + Batch5 + Batch6 + Batch7 + msex0 + RINcontinuous + PCT_CODING_BASES + PCT_INTERGENIC_BASES + pmi + age_at_visit")
  str <- as.formula(str)
  pb$tick()
  #res<-lme4::lmer(str,foo)
  res <- lm(str,foo)
  return(summary(res)$coef[2,4])
}

gene_ids <- colnames(countMatrix)
gene_ids <- gene_ids[-which(gene_ids=='sampleId')]

pb <- progress::progress_bar$new(total=(length(gene_ids)))
system.time(pil_lm_results <- sapply(gene_ids,run_simple_model,foo,pb))

gap::qqunif(lm_results)
