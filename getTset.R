#Create a ToxicoSet using raw data downloaded from E-MEXP-2458 on ArrayExpress
library(ToxicoGx)
library(utils)
library(affy)
library(Biobase)
library(biomaRt)
library(dplyr)
library(SummarizedExperiment)
#library(biocompute)

#call CDF

library("hgu133plus2hsensgcdf")
cdf <- "hgu133plus2hsensgcdf"

#To be run for eset normalization
celfn <- list.files(path = '/pfs/EXP2458Array/raw/', pattern = '*.CEL$', full.names = TRUE)

eset_70 <- just.rma(filenames = celfn, verbose = TRUE, cdfname = cdf)
saveRDS(eset_70, "/pfs/out/eset_70_pre.rds")
#eset_70 <- readRDS("../data/eset_70_pre.rds")

storageMode(eset_70) <- "environment"
#subsetting out internal control probes
eset_70 <- subset(eset_70, substr(rownames(Biobase::exprs(eset_70)), 0, 4) != "AFFX")

#renaming columns
colnames(eset_70@assayData$exprs)<-gsub(".CEL", "", (colnames(eset_70@assayData$exprs)))
colnames(eset_70@assayData$se.exprs)<-gsub(".CEL", "", (colnames(eset_70@assayData$se.exprs)))
rownames(eset_70@protocolData@data)<-gsub(".CEL", "", (rownames(eset_70@protocolData@data)))
rownames(eset_70@phenoData@data)<-gsub(".CEL", "", (rownames(eset_70@protocolData@data)))

#lock eset@assayData environment again
storageMode(eset_70)<-"lockedEnvironment"

annotation(eset_70) <-"rna"

######################################################################################################################################################
#Creating featureData object
ensembl_data <- "hsapiens_gene_ensembl"
CELgenes <- rownames(eset_70@assayData$exprs)

ensembl<-useMart("ensembl", dataset = ensembl_data, host="uswest.ensembl.org",ensemblRedirect = FALSE)
results <- getBM(attributes=c("external_gene_name","ensembl_gene_id","gene_biotype","entrezgene_id","external_transcript_name","ensembl_transcript_id"), filters = "ensembl_gene_id",values=gsub("_at","",CELgenes),mart=ensembl)
uniqueBiomaRt <- results[!duplicated(results$ensembl_gene_id),]
names(uniqueBiomaRt)<-c("gene_name", "gene_id", "gene_biotype", "EntrezGene.ID", "transcript_name", "transcript_id")
#saveRDS(uniqueBiomaRt, "../data/uniqueBiomaRt.rds")

finalFeature <- uniqueBiomaRt

names(finalFeature)[1] <- "Symbol"
finalFeature$BEST <- NA
finalFeature<-arrange(finalFeature,finalFeature$gene_id)
finalFeature$gene_id <- paste(finalFeature$gene_id,"_at", sep = "")
rownames(finalFeature)<-finalFeature$gene_id

#adding 10 extra gene ids from eset
na_mat <- matrix(NA, ncol = ncol(finalFeature), nrow = 10)
colnames(na_mat) <- colnames(finalFeature)
finalFeature_appended <- rbind(finalFeature, na_mat)
finalFeature_appended$gene_id[19930:19939] <- setdiff(rownames(exprs(eset_70)), rownames(finalFeature))
rownames(finalFeature_appended)[19930:19939] <- finalFeature_appended$gene_id[19930:19939]
######################################################################################################################################################
#Creating phenoData object

#sample to file relatiosnhip
pheno_sdrf <- read.delim("/pfs/EXP2458Array/raw/E-MEXP-2458.sdrf.txt", sep = "\t",stringsAsFactors = F)
#Rename columns
colnames(pheno_sdrf) <- c("Source Name", "Characteristics [Organism]", "Term Source REF", "Term Accession Number", "Cell Type", "dataset.cellid",
                          "Term Source REF", "Term.Accession.Number.1", "Protocol REF", "Protocol.REF.1", "Extract Name", "Material Type", "Protocol.REF.2",
                          "Labeled Extract Name", "Label", "Material.Type.1", "Protocol.REF.3", "Hybridization.Name", "Array.Design.REF", "Protocol.REF.4",
                          "duration", "Unit [TimeUnit]", "concentration.raw", "Unit [ConcentrationUnit]", "Abbreviated name", "Term.Source.REF.2",
                          "Term.Accession.Number.2", "Factor Value [CELL_LINE]", "Term.Source.REF.3", "Term.Accession.Number.3", "Scan.Name", "Array.Data.File", 
                          "Comment [ArrayExpress FTP file]", "Protocol.REF.5", "Derived.Array.Data.Matrix.File", "Comment [Derived ArrayExpress FTP file]")

#convert nanomolar to micromolar
pheno_sdrf$concentration <- pheno_sdrf$concentration.raw
pheno_sdrf$dataset.drugid <- NA
pheno_sdrf$samplename <- gsub(".CEL","",pheno_sdrf$Scan.Name)

#add curated cell names used in lab to cellid col
replace_cell <- function(cl){
  if(cl == "HepaRG"){
    "HepaRG"
  }
  else if(cl == "HepG2"){
    "Hep-G2"
  }
}

pheno_sdrf$cellid <- as.character(lapply(pheno_sdrf$dataset.cellid, replace_cell))

#convert nanomolar to micromolar
for(un in 1:nrow(pheno_sdrf)){
  
  if(pheno_sdrf$`Unit [ConcentrationUnit]`[un] == "nM"){
    pheno_sdrf$concentration[un] <- pheno_sdrf$concentration.raw[un]*0.001 
  }
  
}

#CHECK DMSO CONC####

#Add expanded drug names from the paper
exp_names <- matrix(NA, ncol = 2, nrow = 6)
colnames(exp_names) <- c("abbr.name", "exp.name")
exp_names[,"abbr.name"] <- unique(pheno_sdrf$`Abbreviated name`)
exp_names[,"exp.name"] <- c("benzo[a]pyrene", "aflatoxin B1", "2,3,7,8-tetrachlorodibenzo-para-dioxin", "dimethyl sulfoxide", "cyclosporin A", "17b-estradiol")

for(dn in 1:nrow(pheno_sdrf)){
  for (ex in 1:nrow(exp_names)) {
    if(pheno_sdrf$`Abbreviated name`[dn] == exp_names[ex,1]){
      pheno_sdrf$dataset.drugid[dn] <- exp_names[ex,2]
    }
  }
}
rownames(pheno_sdrf) <- pheno_sdrf$samplename


#add curated compound names used in lab to drugid col
replace_drug <- function(dr){
  if(dr == "17b-estradiol"){
    "Estradiol"
  }
  else if(dr == "2,3,7,8-tetrachlorodibenzo-para-dioxin"){
    "TCDD"
  }
  else if(dr == "aflatoxin B1"){
    "Aflatoxin B1"
  }
  else if(dr == "benzo[a]pyrene"){
    "Benzo[a]pyrene"
  }
  else if(dr == "cyclosporin A"){
    "Cyclosporin A"
  }
  else if(dr == "dimethyl sulfoxide"){
    "DMSO"
  }
}
#TCDD (also used in paper and Pubchem synonym) replaced with non-human readable Pubchem name - not present in lab list
#Aflatoxin B1 - not present in lab list

pheno_sdrf$drugid <- as.character(lapply(pheno_sdrf$dataset.drugid, replace_drug))

#add dose_level and xp_type
assign_dose <- function(ds){
  if(ds == "DMSO"){
    "Control"
  }
  else {
    "High"
  }
}

pheno_sdrf$dose_level <- as.character(lapply(pheno_sdrf$drugid, assign_dose))

assign_xp <- function(xp){
  if(xp == "DMSO"){
    "control"
  }
  else {
    "perturbation"
  }
}

pheno_sdrf$xptype <- as.character(lapply(pheno_sdrf$drugid, assign_xp))
  
  
pheno_sdrf$batchid <- "NA"

#add replicate col

pheno_sdrf$individual_id <- NA

for(rp in 1:nrow(pheno_sdrf)){
  spl <- strsplit(pheno_sdrf$`Source Name`[rp], split = "_")
  pheno_sdrf$individual_id[rp] <- as.integer(spl[[1]][4])
}

pheno_sdrf <- pheno_sdrf[,c("samplename","Source Name", "Characteristics [Organism]", "Cell Type", "cellid", "dataset.cellid", "individual_id", "duration", "Unit [TimeUnit]",
                "concentration.raw", "Unit [ConcentrationUnit]","concentration", "Abbreviated name","drugid","dataset.drugid", "dose_level","Scan.Name","Array.Data.File",
                "xptype","batchid","Term Source REF","Term Accession Number","Term Source REF","Term.Accession.Number.1","Protocol REF","Protocol.REF.2","Extract Name","Material Type","Protocol.REF.1","Labeled Extract Name",
                "Label",  "Material.Type.1", "Protocol.REF.3", "Hybridization.Name", "Array.Design.REF", "Protocol.REF.4","Term.Source.REF.2",
                "Term.Accession.Number.2", "Factor Value [CELL_LINE]", "Term.Source.REF.3", "Term.Accession.Number.3","Comment [ArrayExpress FTP file]", "Protocol.REF.5", "Derived.Array.Data.Matrix.File", "Comment [Derived ArrayExpress FTP file]")]


pheno_sdrf <- pheno_sdrf[rownames(pData(eset_70)),]
finalFeature_appended <- finalFeature_appended[rownames(fData(eset_70)),]

stopifnot(all(rownames(pheno_sdrf) == rownames(pData(eset_70))))
stopifnot(all(rownames(finalFeature_appended) == rownames(fData(eset_70))))

#ASSIGN PHENO DATA AND FEATURE DATA TO ESET
pData(eset_70) <- pheno_sdrf
fData(eset_70) <- finalFeature_appended

eset_70 <- eset_70[,order(pData(eset_70)$dataset.drugid)]
#sorting rownames to maintain feature data mapping that is otherwise shuffled after converting to SE
fData(eset_70) <- fData(eset_70)[sort(rownames(fData(eset_70))),]
stopifnot(all(rownames(fData(eset_70)) == rownames(exprs(eset_70))))
stopifnot(all(rownames(pData(eset_70)) == colnames(exprs(eset_70))))

new_SE_EMEXP <-SummarizedExperiment::makeSummarizedExperimentFromExpressionSet(eset_70)

stopifnot(all(rownames(colData(new_SE_EMEXP)) == rownames(pData(eset_70))))
stopifnot(all(rownames(row(new_SE_EMEXP)) == rownames(fData(eset_70))))

#new_SE_EMEXP <-  as(eset_70, Class = "SummarizedExperiment")
saveRDS(new_SE_EMEXP, "/pfs/out/new_SE_EMEXP.rds")
######################################################################################################################################################
#CREATING CURATION DRUG
#curate dataset drug names to lab names (if  present)

curationDrug <- matrix(NA, ncol = 2, nrow = 6)
colnames(curationDrug) <- c("unique.drugid", "dataset.drugid")
curationDrug[,2] <- unique(pheno_sdrf$dataset.drugid)

curationDrug[,1] <- c("Benzo[a]pyrene", "Aflatoxin B1", "Cyclosporin A", "TCDD", "Estradiol", "DMSO")

rownames(curationDrug) <- curationDrug[,1]
######################################################################################################################################################
#CREATING CURATION CELL
curationCell <- matrix(NA, ncol = 2, nrow = 2)
colnames(curationCell) <- c("unique.cellid", "dataset.cellid")
curationCell[,2] <- unique(pheno_sdrf$dataset.cellid)
#HepaRG - not present in lab list
curationCell[,1] <- unique(pheno_sdrf$cellid)
rownames(curationCell) <- curationCell[,1]
######################################################################################################################################################
#CREATING CURATION TISSUE
curationTissue <- matrix(NA, ncol = 2, nrow = 2)
colnames(curationTissue) <- c("unique.tissueid", "dataset.tissueid")
curationTissue[,2] <- "Liver"
#HepaRG - not present in lab list
curationTissue[,1] <- c("Liver")
rownames(curationTissue) <- rownames(curationCell)
######################################################################################################################################################
#CREATING CELL
cell <- subset(pheno_sdrf, select = c("Cell Type", "cellid", "dataset.cellid"), drop = F)
cell$tissueid <- "Liver"
cell$species <- "Human"
cell <- unique(cell)
rownames(cell) <- rownames(curationCell)

#CREATING DRUG
drug <- subset(pheno_sdrf, select = c("dataset.drugid", "Abbreviated name"))
drug <- unique(drug)
drug$drugid <- rownames(curationDrug)
rownames(drug) <- rownames(curationDrug)
######################################################################################################################################################

#PUTTING TOXICOSET TOGETHER
EMEXP2458 <- ToxicoSet("EMEXP2458",
                         molecularProfiles=list("rna"= new_SE_EMEXP),
                         cell=cell,
                         drug=drug,
                         sensitivityInfo=NULL,
                         sensitivityRaw=NULL,
                         sensitivityProfiles=NULL,
                         curationDrug=curationDrug,
                         curationCell=curationCell,
                         curationTissue=curationTissue,
                         datasetType="perturbation",
                         verify = TRUE)
saveRDS(EMEXP2458, "/pfs/out/EMEXP2458.rds")

  
  
# ####CREATE BIOCOMPUTE OBJECT####
# 
# ###########################
# #####Provenance Domain#####
# ###########################
# 
# #Created and modified dates
# #Sys.setenv(TZ = "EST")
# created <- as.POSIXct(Sys.time(), format = "%Y-%m-%dT%H:%M:%S", tz = "EST")
# modified <- as.POSIXct(Sys.time(), format = "%Y-%m-%dT%H:%M:%S", tz = "EST")
# 
# #Contributions
# contributors <- data.frame(
#   "name" = c("Anthony Mammoliti", "Petr Smirnov", "Benjamin Haibe-Kains"),
#   "affiliation" = c(rep("University Health Network", 3)),
#   "email" = c("anthony.mammoliti@uhnresearch.ca", "petr.smirnov@utoronto.ca", "Benjamin.Haibe-Kains@uhnresearch.ca"),
#   "contribution" = c("createdBy","createdBy","authoredBy"),
#   "orcid" = c(NA,NA,"https://orcid.org/0000-0002-7684-0079"),
#   stringsAsFactors = FALSE
# )
# 
# #License
# license <- "https://opensource.org/licenses/Apache-2.0"
# 
# #Name of biocompute object
# name <- "EMEXP2458"
# 
# #Version of biocompute object
# bio_version <- "1.0.0"
# 
# #Embargo (none)
# embargo <- c()
# 
# #Derived from and obsolete after (none)
# derived_from <- c()
# obsolete_after <- c()
# 
# #reviewers (none)
# review <- c()
# 
# #compile domain
# provenance <- compose_provenance_v1.3.0(
#   name, bio_version, review, derived_from, obsolete_after,
#   embargo, created, modified, contributors, license
# )
# provenance %>% convert_json()
# 
# 
# ############################
# #####Description Domain#####
# ############################
# times_rnaseq <- as.POSIXct("2020-01-22T8:17:01", format = "%Y-%m-%dT%H:%M:%S", tz = "EST")
# #Keywords and platform info
# keywords <- c("Biomedical", "Toxicogenomics")
# platform <- c("Pachyderm", "ORCESTRA (orcestra.ca)", "Linux/Ubuntu")
# 
# #Metadata for each pipeline step
# pipeline_meta <- data.frame(
#   "step_number" = c("1","2","3"),
#   "name" = c("Microarray processing",
#              "Curated Sample and treatment identifier compilation",
#              "Build data object"),
#   "description" = c("Normalization of microarray data", 
#                     "Download of appropriate sample and treatment identifiers from GitHub (curations performed by BHK Lab - http://bhklab.ca)",
#                     "Building of ORCESTRA data object"),
#   "version" = c(1.0,1.0,1.0),
#   stringsAsFactors = FALSE
# )
# 
# #Inputs for each pipeline step
# pipeline_input <- data.frame(
#   "step_number" = c("1","2","3"),
#   "filename" = c("Raw microarray data",
#                  "Sample and treatment annotation data",
#                  "Script for data object generation"),
#   "uri" = c(
#     "https://www.ebi.ac.uk/arrayexpress/experiments/E-MEXP-2458/",
#     "/pfs/EXP2458Array/raw/E-MEXP-2458.sdrf.txt",
#     "https://github.com/BHKLAB-Pachyderm/getEXP2458/blob/main/getTset.R"
#   ),
#   "access_time" = c(times_rnaseq,created,created),
#   stringsAsFactors = FALSE
# )
# 
# 
# #Outputs for each pipeline step
# pipeline_output <- data.frame(
#   "step_number" = c("1","2","3"),
#   "filename" = c("Processed microarray data", 
#                  "Downloaded sample and treatment annotations",
#                  "Data object"),
#   "uri" = c(
#     "/pfs/out/new_SE_EMEXP.rds",
#     "/pfs/EXP2458Array/raw/E-MEXP-2458.sdrf.txt",
#     "/pfs/out/EMEXP2458.rds"
#   ),
#   "access_time" = c(times_rnaseq,created,created),
#   stringsAsFactors = FALSE
# )
# 
# #xref (none)
# xref <- c()
# 
# #pipeline prereq (none)
# pipeline_prerequisite <- c()
# 
# #compile domain
# description <- compose_description_v1.3.0(
#   keywords, xref, platform,
#   pipeline_meta, pipeline_prerequisite, pipeline_input, pipeline_output
# )
# description %>% convert_json()
# 
# 
# ############################
# ######Execution Domain######
# ############################
# 
# script <- c()
# script_driver <- c()
# 
# #software/tools and its versions used for data object creation
# software_prerequisites <- data.frame(
#   "name" = c("Pachyderm", "Docker Image"),
#   "version" = c("1.9.3", "v3"),
#   "uri" = c(
#     "https://www.pachyderm.com", "https://hub.docker.com/r/bhklab/toxicogx"
#   ),
#   stringsAsFactors = FALSE
# )
# 
# software_prerequisites[,"access_time"] <- rep(NA, length(software_prerequisites$name))
# software_prerequisites[,"sha1_chksum"] <- rep(NA, length(software_prerequisites$name))
# 
# external_data_endpoints <- c()
# environment_variables <- c()
# 
# execution <- compose_execution_v1.3.0(
#   script, script_driver, software_prerequisites, external_data_endpoints, environment_variables
# )
# execution %>% convert_json()
# 
# 
# ############################
# ######Extension Domain######
# ############################
# 
# #repo of scripts/data used
# scm_repository <- data.frame("extension_schema"= c("https://github.com/BHKLAB-Pachyderm"))
# scm_type <- "git"
# scm_commit <- c()
# scm_path <- c()
# scm_preview <- c()
# 
# scm <- compose_scm(scm_repository, scm_type, scm_commit, scm_path, scm_preview)
# scm %>% convert_json()
# 
# extension <- compose_extension_v1.3.0(scm)
# extension %>% convert_json()
# 
# ############################
# ######Parametric Domain#####
# ############################
# 
# df_parametric <- data.frame(
#   "param" = c(),
#   "value" = c(),
#   "step" = c(),
#   stringsAsFactors = FALSE
# )
# 
# parametric <- compose_parametric_v1.3.0(df_parametric)
# parametric %>% convert_json()
# 
# 
# 
# ############################
# ######Usability Domain######
# ############################
# 
# #usability of our data objects
# text <- c(
#   "Pipeline for creating EMEXP2458 data object through ORCESTRA (orcestra.ca), a platform for the reproducible and transparent processing, sharing, and analysis of biomedical data."
# )
# 
# usability <- compose_usability_v1.3.0(text)
# usability %>% convert_json()
# 
# 
# ######################
# ######I/O Domain######
# ######################
# 
# input_subdomain <- data.frame(
#   c("Raw microarray data",
#     "Sample and treatment annotation data",
#     "Script for data object generation"),
#   "uri" = c(
#     "https://www.ebi.ac.uk/arrayexpress/experiments/E-MEXP-2458/",
#     "/pfs/EXP2458Array/raw/E-MEXP-2458.sdrf.txt",
#     "https://github.com/BHKLAB-Pachyderm/getEXP2458/blob/main/getTset.R"
#   ),
#   "access_time" = c(times_rnaseq,created,created),
#   stringsAsFactors = FALSE
# )
# 
# output_subdomain <- data.frame(
#   "mediatype" = c("RDS", "csv", "RDS", "RDS"),
#   "uri" = c(
#     "/pfs/out/new_SE_EMEXP.rds",
#     "/pfs/EXP2458Array/raw/E-MEXP-2458.sdrf.txt",
#     "/pfs/out/EMEXP2458.rds"
#   ),
#   "access_time" = c(created,created,created),
#   stringsAsFactors = FALSE
# )
# 
# io <- compose_io_v1.3.0(input_subdomain, output_subdomain)
# io %>% convert_json()
# 
# 
# ########################
# ######Error Domain######
# ########################
# 
# empirical <- c()
# algorithmic <- c()
# 
# error <- compose_error(empirical, algorithmic)
# error %>% convert_json()
# 
# 
# ####Retrieve Top Level Fields####
# tlf <- compose_tlf_v1.3.0(
#   provenance, usability, extension, description,
#   execution, parametric, io, error
# )
# tlf %>% convert_json()
# 
# 
# ####Complete BCO####
# 
# bco <- biocompute::compose_v1.3.0(
#   tlf, provenance, usability, extension, description,
#   execution, parametric, io, error
# )
# bco %>% convert_json() %>% export_json("/pfs/out/EMEXP2458_BCO.json") %>% validate_checksum()
# 
# 
# 
# 
