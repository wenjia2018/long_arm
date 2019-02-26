########################################################
# HOUSEKEEPERS
########################################################
# Eisenberg, E., & Levanon, E. Y. (2013). Human housekeeping genes, revisited. Trends in Genetics, 29(10), 569-574.
housekeepers11 = c("C1orf43",  "CHMP2A", "EMC7",  "GPI", "PSMB2", 
                   "PSMB4", "RAB7A", "REEP5",  "SNRPD3", "VCP", "VPS29")

########################################################
# DEFINE CTRA SET (FROM LITERATURE)
########################################################

ctra_original = 
  c("IL1A", "IL1B", "IL6", "IL8", "TNF", "PTGS1", "PTGS2", "FOS", "FOSB", "FOSL1", 
    "FOSL2", "JUN", 'JUNB', "JUND", "NFKB1", "NFKB2", "REL", "RELA", "RELB", 
    "GBP1", "IFI16", "IFI27", "IFI27L1", "IFI27L2", "IFI30", "IFI35", "IFI44", 
    "IFI44L", "IFI6", "IFIH1", "IFIT1", "IFIT2", "IFIT3", "IFIT5", "IFIT1L", "IFITM1",
    "IFITM2", "IFITM3", "IFITM4P", "IFITM5", "IFNB1", "IRF2", "IRF7", "IRF8",
    "MX1", "MX2", "OAS1", "OAS2", "OAS3", "OASL", "IGJ", "IGLL1", "IGLL3")

(ctra_available = intersect(featureNames(dat), ctra_original))

#######################################################
# DEFINE CTRA GENE SET
# (FROM MY PREVIOUS PROJECT, WHICH UPDATED SOME OUTDATED GENE NAMES)
########################################################

# by their HGNC names (HUGO gene nomenclature committee).
inflammatory = c("IL1A", "IL1B", "IL6", "IL8", "TNF", "PTGS1", "PTGS2", "FOS",
                 "FOSB", "FOSL1", "FOSL2", "JUN", "JUNB", "JUND", "NFKB1",
                 "NFKB2", "REL", "RELA", "RELB") 

interferonTypeI = c("GBP1", "IFI16", "IFI27", "IFI27L1", "IFI27L2", "IFI30",
                    "IFI35", "IFI44", "IFI44L", "IFI6", "IFIH1", "IFIT1",
                    "IFIT2", "IFIT3", "IFIT5", "IFIT1L", "IFITM1", "IFITM2",
                    "IFITM3", "IFITM4P", "IFITM5", "IFNB1", "IRF2", "IRF7",
                    "IRF8", "MX1", "MX2", "OAS1", "OAS2", "OAS3", "OASL") 

antibody = c("IGJ", "IGLL1", "IGLL3") 

ctra0 = c(inflammatory, interferonTypeI, antibody)

ctraCore0 = c("IRF7", "JUN", "IGJ", "IL8", "IL1B", "FOSB", "FOSL2", "IFIT3",
              "IFI35", "IFI44L", "MX1", "OAS2")

ctraCore = c("IRF7", "JUN", "JCHAIN", "CXCL8", "IL1B", "FOSB", "FOSL2", "IFIT3",
             "IFI35", "IFI44L", "MX1", "OAS2")

########################################################
# INTERSECTION OF ABOVE SETS WITH THE GENES IN THIS SAMPLE
########################################################
ctra_available = intersect(featureNames(dat), ctra0)
inflammatory  = inflammatory %>% intersect(ctra_available)
interferonTypeI = interferonTypeI %>% intersect(ctra_available)
antibody = antibody %>%  intersect(ctra_available)

########################################################
# CELL TYPE MARKERS
########################################################
# From Steve's email: the set of genes used as covariates to crudely control for leukocyte subset prevalence is

cell_types = 
  list(T_0 = c("CD3D", "CD3E"), # for T cells (only occasionally the "CD3E")
       T_CD4 = "CD4", # for the CD4+ subset of T cells
       T_CD8 = "CD8A", # for the CD8+ subset of T cells
       B = "CD19", # for B cells
       NK = c("FCGR3A", "NCAM1"), # for NK cells
       M = "CD14") # for monocytes

cell_types_one =  map_chr(cell_types, 1) # hack: for now just take one gene for each type 

########################################################
# CECILIA
########################################################

# by their HGNC names (HUGO gene nomenclature committee).
infarction = c("ABCA1", "ACE", "ADD1", "AGT", "AGTR1", "ALOX5AP", "APOA1", "APOA5", "APOE", "CCL11", "CCR2", "CCR5","CD14", "CETP", "COMT", "CX3CR1", 
               "CYP11B2", "CYP2C9", "ENPP1", "ESR1", "F12", "F13A1", "F2", "F5", "FGB", "FTO", "GJA4", "GP1BA", "NR3C1", "GSTT1", "HFE", "HSL", "HNF1", 
               "HTR2A","ICAM1", "IL1B", "IL6", "IL18", "ITAGA2", "ITGB3", "KIF6", "LDL R", "LIPC", "LPL", "LRP1", "MGP", "MMP3", "MTHFR", "MTP", "MTR",
               "OLR1", "p22-PHOX","PAI1","PECAM", "PON1+2","PPARG","PTGS2","RECQL2", "SELE", "SELP","TFPI", "THBD","TLR4","TNF", "TNFRSF1A", "UCP-2")

depression = c("MSRA", "FDFT1", "C8orl12","c8orl13", "MTMR9","BLK","MFHAS1")
loneliness = c("TCF4", "PHF2", "NMUR2", "EPB41L2", "LOC100499466", "OSTF1", "ARFGEF2", "GBE1", "ERBB4", "NMUR2", "OR1S1", "HIVEP1", "NDUFS3", 
               "RAB9BP1", "MBD5")
social = c("BARHL2", "DPYD", "SLC4A10", "NUP35", "GNRHR", "ADH1B", "MIR2113", "DPP6", "ZNF462", "PAX2", "C17orf112", "PMAIP1", "TNRC6B", "NFIA",
           "BARHL2", "LRRN2", "KLHL29", "OTX1", "LONRF2", "CNTNAP5", "MST1", "IFT57", "PPA2", "LINC00461", "ODZ2", "FBXL4", "TAC1", "LHFPL3", "TNKS", 
           "PCDH17", "TCF4", "KANK4", "THSD7B", "CAMKV", "CADM2", "MIR1275", "MIR147A" )
diabetes = c("CD101", "CEP68", "EHHADH", "RP11-10L12.4", "ANKH", "POC5", "RREB1", "MICB", "HLA-DQB1", "CENPW", "ARG1", "MED23", "TP53INP1", "RPL8",
             "CAMK1D", "CWF19L1", "SNORA12", "PLEKHA1", "SSSCA1", "ARAP1", "P2RX4", "CAMKK2", "C15orf38",  "RCCD1", "ANKFY1", "ATP5G1", "UBE2Z", "PABPC4", "PABPC4",
             "LTA", "CUTA", "ARG1", "TP53INP1", "NUDT5", "CAMK1D", "LOC283070", "CWF19L1", "PLEKHA1", "KLHDC5", "P2RX4", "CAMKK2", "CAMKK2", "C15orf38",
             "ATP5G1", "ATP5G1")
alzeheimer = c("CR1", "BIN1", "CD2AP", "EPHA1", "CLU", "MS4A6A", "PICALM", "ABCA7", "APOE", "HLA-DRB5", "HLA-DRB1", "PTK2B", "SORL1", "SLC24A4", "RIN3", "INPP5D",
               "MEF2C", "NME8", "ZCWPW1", "CELF1", "FERMT2", "CASS4")

# not all GENES IN THE SET OF INTEREST are here (name change?)
(heart_available = intersect(featureNames(dat), infarction))
(depression_available = intersect(featureNames(dat), depression))
(loneliness_available = intersect(featureNames(dat), loneliness))
(social_available = intersect(featureNames(dat), social))
(diabetes_available = intersect(featureNames(dat), diabetes))
(alzhei_available = intersect(featureNames(dat), alzeheimer))

