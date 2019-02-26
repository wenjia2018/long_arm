########################################################
# INITIALIZE ANALYSES
########################################################
########################################################
# FIRST INIT DATA
########################################################

source("R/init_data.R")

########################################################
# REGRESSION MODEL SPECIFICATION: USED TO SELECT A DE GENE SET
########################################################

of_in = "obesew5" # effect of_interest, i.e. one of the below variables
vars = c("lowbirthweight", "obesew1or2", "obesew3", "obesew4", "obesew5",
         "re_2", "re_3", "re_4", "re_5", "Plate_Year1Plate02", "Plate_Year1Plate03", "Plate_Year1Plate04",
         "Plate_Year1Plate05", "Plate_Year1Plate06", "Plate_Year1Plate07", "Plate_Year1Plate08",
         "Plate_Year1Plate09", "Plate_Year1Plate10", "Plate_Year1Plate11", 'Plate_Year1Plate12',
         "AvgCorrelogram100", "sex_interv_m", "BirthY", "edu_max_votec", "edu_max_post", "edu_max_college",
         "bingedrink", "currentsmoke", "W5REGION_MW", "W5REGION_S", "W5REGION_W")
rhs   = str_c(vars, collapse = " + ")
keep  = phen %>% select(vars) %>% complete.cases
nsubs = sum(keep)

