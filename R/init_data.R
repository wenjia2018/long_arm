# analyse_count = function(){
rm(list=ls()) 

########################################################
# PACKAGES
########################################################

# Packages
library(tidyverse)
# tidyverse_update()

########################################################
# BIOCONDUCTOR
########################################################
# source("https://bioconductor.org/biocLite.R")
# biocLite()
pkg = c(
  "DESeq2",
  "edgeR",
  "limma",
  "Biobase",
  "recount",
  "globaltest",
  "gplots",
  "recount", # recount data repository
  "GEOquery", # geo data repository
  "GEOmetadb", # query geo metadata
  "biobroom",
  "compositions"
)
pkg %>%
  # map(biocLite) %>%  
  map(library, character.only = T)

########################################################
# CRAN
########################################################

pkg_cran = c(
  "ggformula",
  "modelr",
  "DescTools"
)
pkg_cran %>% 
  # map(install.packages) %>%
  map(library, character.only = TRUE)

########################################################
# HOMEMADE
########################################################

# source("R/utility.R") # redundant

# model.matrix() FUSSY ABOUT NA
old_na_action <- options('na.action')
options(na.action='na.pass')
# options(na.action = old_na_action)

########################################################
# LOAD DATA 
########################################################

########################################################
# PHENOTYPE
# GET THE LATEST PHENOTYPE AND RNA DATA
########################################################

# PHENOTYPE
if(str_detect(Sys.info()["nodename"], "Chumbley.jc.uzh.ch|JCPYD-Chumbley.local")){
  waves  <- readRDS("/volumes/share/preprocess/data/waves.rds") # path relative to local
} else {
  waves  <- readRDS("/home/share/preprocess/data/waves.rds")  # path relative to server root
}
# TO AVOID LATER MERGE CONFLICT
waves$Plate = NULL
waves$AvgCorrelogram100 = NULL

# RNA
dat       = readRDS("/volumes/share/preprocess/data/dt.rds") # expression set 
dat_ref   = readRDS("/volumes/share/preprocess/data/dt_steve.rds") # reference gene normalization

# augment phenotype data with that derived in init_pheno.R
# TO AVOID MERGE CONFLICTS
# pData(dat) <- pData(dat) %>% select(VialID, AID) # REMOVE LAURENS VARIABLES
pData(dat)     <- pData(dat) %>% left_join(waves, by = "AID")
pData(dat_ref) <- pData(dat_ref) %>% left_join(waves, by = "AID")

# LOOSE ALL DARK MATTER
dark_matter = str_which(featureNames(dat), "^ENSG.*")
dat         = dat[-dark_matter, ] 
dat_ref     = dat_ref[-dark_matter, ] 

########################################################
# DEFINE CTRA GENE SET
########################################################

source("R/define_gene_sets.R")

########################################################
# QUALITY CHECKS
########################################################
if(0){ 
  source("R/quality.R")
}

########################################################
# EXTRACT PHENOTYPE DATA FROM EXPRESSION SET
########################################################

phen = dat %>% pData
# phen$male = phen$sex_interv_m # BEWARE NAME CHANGE!!!!!!
# phen$raceth = phen$race_interv # BEWARE NAME CHANGE!!!!!!
 # phen$total_assetsw5 = ifelse(is.na(phen$assets_household_net_w5))
