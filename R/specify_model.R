# FUNCTIONS
source("R/init_data.R")

########################################################
# DEFINE DEPENDENT VARIABLE SETS AND LOAD INTO THE GLOBAL ENV
########################################################
# create outcome variables, either do the ilr here or within the function call

outcomes = 
  list(
    cell_type   = cell_types_one,
    ctra        = ctra_available,
    inflame     = inflammatory ,
    interferon  = interferonTypeI
    # housekeepers11, 
    # antibody
    # # # CECILIAS
    # heart = heart_available,
    # depression = depression_available,
    # RANK TOO LOW FOR THE FOLLOWING TWO APARENTLY
    # loneliness_available = loneliness_available,
    # social_available  = social_available,
    # diabetes  = diabetes_available,
    # alzhiemers = alzhei_available
  )

# INTO WORKING ENV
outcomes %>% 
  map(
    ~ dat[.x] %>% exprs %>% `+`(1) %>% t %>% acomp %>% ilr 
  ) %>% 
  list2env(.GlobalEnv)

# ALSO SAVE ALL INTO ONE LIST
outme = 
  list(
    cell_type,
    ctra,
    inflame,
    interferon
  )

########################################################
# DEFINE TREATMENT AND COVARIABLES
########################################################

controls = c( " ~ male + raceth + Plate + AvgCorrelogram100 +" )

# VARIABLES
focal_treatment = 
  c("edu_p",
    "pbills",
    #"poccrf", 
    #"poccrm",
    #"pneighqual1", 
    #"cdisab",
    #"chealth1",
    "psmoking", 
    "palchohol", 
    "edu_max", 
    "bills",
    "w5occupgroup", 
    "assets_household_net_w5",
    #"ahealth2",
    #"ahealth3",
    "srh",
    "w5currentsmoke",
    "drug",
    "bingedrink_month",
    "income_pp1_log",
    "smoke_cf_p1",
    "SEI_max_p_w12"
  )

# SAME VARIABLES, NOW IN VARIABLE SETS
f0  = controls %>% str_remove("~") %>% str_sub(1, str_length(.)-1)
f1s = 
  c(
    child_health = "smoke_cf_p1 + psmoking + palchohol",  #drop cdisab + chealth1 +                           
    child_ses    = "income_pp1_log + edu_p + pbills + SEI_max_p_w12 ",#drop + pneighqual1                         
    adult_health = "srh + w5currentsmoke + drug + bingedrink_month",  #drop ahealth2 + ahealth3
    adult_ses    = "edu_max + bills + w5occupgroup + assets_household_net_w5",              
    child        = "income_pp1_log + psmoking + palchohol + edu_p + pbills + SEI_max_p_w12 + smoke_cf_p1 ", # all child
    adult        = "srh + w5currentsmoke + drug + bingedrink_month + edu_max + bills + w5occupgroup + assets_household_net_w5" # adult
  ) %>%
  map_chr(~str_c(f0, .x , sep = " + "))  


# mod_spec = 
#   function(f_one, f_zero, model_name){
#     crossing(formula1 = f_one, 
#              formula0 = map_chr(str_c("\\+ ", f_zero), 
#                                 ~str_remove_all(f_one, .x))) %>% 
#       mutate(variable = f_zero,
#              model = model_name)
#   }
#  
# mod_spec(f1s[[1]], c( "psmoking", "palchohol", "smoke_cf_p1") ,  "child_health")
# mod_spec(f1s[[2]], c("income_pp1_log", "edu_p", "pbills", "SEI_max_p_w12"),  "child_ses")
# mod_spec(f1s[[3]], c("srh", "w5currentsmoke", "drug", "bingedrink_month"),  "adult_health")
# mod_spec(f1s[[4]], c("edu_max", "bills", "w5occupgroup", "assets_household_net_w5"),  "adult_ses")
# mod_spec(f1s[[5]], c("smoke_cf_p1", "income_pp1_log","psmoking", "palchohol", "edu_p", "pbills", "SEI_max_p_w12"),  "child")
# mod_spec(f1s[[6]], c( "srh", "w5currentsmoke", "drug", "bingedrink_month", "edu_max", "bills", "w5occupgroup", "assets_household_net_w5"),  "adult")

# fully conditnal
hyp1 = 
  crossing(formula1 = f1s[[1]], 
           formula0 = map_chr(str_c("\\+ ",c( "psmoking", "palchohol", "smoke_cf_p1")), 
                              ~str_remove_all(f1s[[1]], .x))) %>% 
  mutate(variable =                        c( "psmoking", "palchohol", "smoke_cf_p1"),
         model = "child_health")
hyp2 =
  crossing(formula1 = f1s[[2]], 
           formula0 = map_chr(str_c("\\+ ",c("income_pp1_log", "edu_p", "pbills", "SEI_max_p_w12")), 
                              ~str_remove_all(f1s[[2]], .x))) %>% 
  mutate(variable =                          c("income_pp1_log", "edu_p", "pbills", "SEI_max_p_w12"),
         model = "child_ses")
hyp3 = 
  crossing(formula1 = f1s[[3]], 
           formula0 = map_chr(str_c("\\+ ", c("srh", "w5currentsmoke", "drug", "bingedrink_month")), 
                              ~str_remove_all(f1s[[3]], .x))) %>% 
  mutate(variable =                         c("srh", "w5currentsmoke", "drug", "bingedrink_month"),
         model = "adult_health")
hyp4 = 
  crossing(formula1 = f1s[[4]], 
           formula0 = map_chr(str_c("\\+ ",c("edu_max", "bills", "w5occupgroup", "assets_household_net_w5")), 
                              ~str_remove_all(f1s[[4]], .x))) %>% 
  mutate(variable =                        c("edu_max", "bills", "w5occupgroup", "assets_household_net_w5"),
         model = "adult_ses")
hyp5 = 
  crossing(formula1 = f1s[[5]], 
           formula0 = map_chr(str_c("\\+ ",
                                    c("smoke_cf_p1", "income_pp1_log","psmoking", "palchohol", "edu_p", "pbills", "SEI_max_p_w12")), 
                              ~str_remove_all(f1s[[5]], .x))) %>% 
  mutate(variable =                 c("smoke_cf_p1", "income_pp1_log","psmoking", "palchohol", "edu_p", "pbills", "SEI_max_p_w12"),
         model = "child")

hyp6 = 
  crossing(formula1 = f1s[[6]], 
           formula0 = map_chr(str_c("\\+ ",c( "srh", "w5currentsmoke", "drug", "bingedrink_month", "edu_max", "bills", "w5occupgroup", "assets_household_net_w5" )), 
                              ~str_remove_all(f1s[[6]], .x))) %>% 
  mutate(variable =                         c( "srh", "w5currentsmoke", "drug", "bingedrink_month", "edu_max", "bills", "w5occupgroup", "assets_household_net_w5" ),
         model = "adult")

hyp = rbind(hyp1, hyp2, hyp3, hyp4, hyp5, hyp6)
rm(hyp1, hyp2, hyp3, hyp4, hyp5, hyp6)

hyper =
  rbind(
    mutate(hyp, outcome = list(cell_type),   out_var = rep("cell_type",  dim(hyp)[1])) ,
    mutate(hyp, outcome = list(ctra),        out_var = rep("ctra",       dim(hyp)[1])) ,
    mutate(hyp, outcome = list(inflame),     out_var = rep("inflame",    dim(hyp)[1])) ,
    mutate(hyp, outcome = list(interferon),  out_var = rep("interferon", dim(hyp)[1]))) %>% 
  as_tibble


# ########################################################
# # CREATE A DF OF FORMULAE FOR THE DESIRED REGRESSIONS (simplest model)
# ########################################################

args  =
  expand.grid(
    outcomes = names(outcomes), # toggle jointly with the previous comment
    controls = controls,
    focal_treatment = focal_treatment,
    stringsAsFactors = F
  ) %>%
  unite(model_formulae, sep = "", remove = F)

