########################################################
# SPECIFY MODEL
########################################################

source("R/specify_model.R")
source("R/utility_coda.R")

########################################################
# LOAD VERSUS CREATE PERMUTATION RESULTS
########################################################

run_from_scratch= F # if true, run permutation analysis. only do this on the server.

########################################################
# SEPERATE, ONE-BY-ONE "UNIVARIATE" LINEAR MODELS (WITH CONTROLS)
# FOR DIFFERENT FOCAL TREATMENTS AND DIFFERENT OUTCOMES
########################################################

# minimal example
perms   <- permute(phen,  3, edu_p)
cod     <- function(x)  lm(args$model_formulae[[1]], data = x, na.action = na.omit) %>% tidy %>% filter(str_detect(term, "edu_p"))
(models <- map(perms$perm, cod))

########################################################
# PARAMETER TESTS (SIMPLEST MODEL)
########################################################

if(run_from_scratch){
  
  library(furrr)
  plan(multiprocess, workers = detectCores()/2 - 2)##use multicore evalation
  
  diff_measures  = 
    args %>% 
    as_tibble %>% 
    mutate(out = future_pmap(., safely(get_permutate_p)))#, n_perm = 10)) 
  
  # EXTRACT 
  diff_measures = 
    diff_measures %>% 
    mutate(result   = map(out, "result"), 
           p_perm   = map_dbl(result, "p_perm"),
           p_par    = map_dbl(result, "p_par"))
  
  # SAVE
  diff_measures %>% saveRDS("data/stats/diff_measuress.rds")
  
} else {
  
  diff_measures = readRDS("data/stats/diff_measuress.rds")
  
}



########################################################
# MULTIPARAMETER (MODEL) COMPARISON
########################################################

# minimal example 
model_compare_par("male + raceth + Plate + AvgCorrelogram100", 
                  "male + raceth + Plate + AvgCorrelogram100 + edu_p + pbills + poccrf + poccrm + pneighqual1",
                  ctra)

if(run_from_scratch){
  
  ########################################################
  # CONSIDER OMBNIBUS TESTS OF EACH VARIABLE SUBSET, {CHILD, ADULT} X {SES, HEALTH}
  ########################################################
  
  ########################################################
  # PARAMETRIC 
  ########################################################
  
  (model_difference_par = 
     expand.grid(formula0 = f0, formula1 = f1s, outcome = outme, stringsAsFactors = F) %>%  
     as_tibble %>% 
     mutate(out = pmap(., model_compare_par)) %>% 
     mutate(p_par  = map_dbl(out, ~tidy(.x) %>%  `[`(2,8) %>% unlist),
            pillai = map_dbl(out, ~ .x$Pillai[2])) %>% 
     mutate(outcome_var = rep(c("cell_type", "ctra", "inflame", 'interferon'), each = length(f1s))))
  
  ########################################################
  # PERMUTATION 
  ########################################################
  
  # EXTRACT P-VAL
  model_difference_perm  = 
    model_difference_par %>% as_tibble %>% select(formula0, formula1, outcome, outcome_var) %>% # select the same inputs as for parametric approach
    mutate(perm_stats    = pmap(., model_compare_perm))
  # JOIN PARA AND PERM 
  model_diff_measures    = model_difference_par %>% left_join(model_difference_perm , by = c("formula0", "formula1", "outcome_var"))
  
  model_diff_measures = 
    model_diff_measures %>% 
    mutate(p_perm = map2_dbl(perm_stats, pillai, ~mean(.x > .y))) %>% 
    mutate(set = rep(names(f1s), length(outme))) %>% # just concise names
    select(p_par, p_perm, outcome_var, set, formula1)
  model_diff_measures %>% saveRDS("data/stats/model_diff_measures.rds")
  
} else {
  
  model_diff_measures = readRDS("data/stats/model_diff_measures.rds")
  
}

########################################################
# MODEL-CONDITIONAL PARAMETER INFERENCE
# FULL CONDITIONAL OF EACH PARAM WITHIN EACH PREDICTION SET
######################################################### 

if(run_from_scratch){
  (mods = 
     hyper %>% 
     mutate(out = pmap(list(formula0, formula1, outcome),
                       model_compare_par)) %>% 
     mutate(p_par  = map_dbl(out, ~ tidy(.x) %>%  `[`(2,8) %>% unlist),
            pillai = map_dbl(out, ~ .x$Pillai[2])) %>% 
     mutate(perm_stats = pmap(list(formula0,  formula1,  outcome,  out_var, 
                                   vars_only_in_the_alt_hyp = variable), 
                              model_compare_perm),
            p_perm = map2_dbl(perm_stats, pillai, ~mean(.x > .y))))  
  
  mods %>% saveRDS("data/stats/conditional_diff_measures.rds")
  
} else {
  
  conditional_diff_measures = readRDS("data/stats/conditional_diff_measures.rds")
  
}

########################################################
# TFBM
########################################################

library(dbr)

# dat = an expression set object
my_obj = dat_ref # steve's normalization
my_obj = dat

########################################################
# A MINIMAL EXAMPLE (SEE ALSO VIGNETTES)
########################################################

out = 
  get_tfbm_p_vals(rhs          =  conditional_diff_measures$formula1[1], 
                  of_in        = "palchohol", 
                  which_matrix = "utr",
                  dat          = my_obj, 
                  n_sim        = 1000)
out%>% extracter

########################################################
# SLIGHTLY MORE GENERALLY
########################################################

# Hack: The current functions cannot generally handle the case where the
# variable of interest (of_in) refers to multiple columns of the design!

my_obj$bingedrink_month = my_obj$bingedrink_month != "nonbinge"
my_obj$edu_max          = my_obj$edu_max != "high" 
my_obj$edu_p            = my_obj$edu_p != "high"
my_obj$w5occupgroup     = my_obj$w5occupgroup %in% c("construction" )

########################################################
# DB INFERENCE
########################################################

tfbm_res = 
  conditional_diff_measures %>% 
  filter(out_var == "cell_type") %>%  # eliminate things redundant to this analysis 
  select(rhs = formula1,
         of_in = variable, 
         model) %>% 
  # filter(of_in == "psmoking") %>% # temp
  mutate(tfbm_results_utr        = pmap(list(rhs, of_in), get_tfbm_p_vals, dat = my_obj, n_sim = 10000, which_matrix = "utr"),
         tfbm_results_exonic     = pmap(list(rhs, of_in), get_tfbm_p_vals, dat = my_obj, n_sim = 10000, which_matrix = "exonic"),
         tfbm_results_exonic_utr = pmap(list(rhs, of_in), get_tfbm_p_vals, dat = my_obj, n_sim = 10000, which_matrix = "exonic_utr"))

########################################################
# EXTRACT RESULTS
########################################################

tfbm_res %>% 
  mutate(tfbm_results_utr_pvals =
           tfbm_results_utr %>% 
           map(extracter)) %>%
  unnest(tfbm_results_utr_pvals)

tfbm_res %>% 
  mutate(tfbm_results_utr_pvals =
           tfbm_results_utr %>% 
           map(extracter)) %>%
  unnest(tfbm_results_utr_pvals) %>% 
  select(p_par, p_npar, p_uni, p_cov) %>% 
  plot(xlim = c(0,1), ylim = c(0,1))

########################################################
# OTHER
########################################################

# NULL STRESS TEST. DO I CONTROL TYPE 1?
pData(my_obj) = pData(my_obj)[sample(1:dim(pData(my_obj))[1]), ]  # SHUFFLE DATA!!!
tfbm_res_null = 
  conditional_diff_measures %>%
  `[`(1:12,) %>% 
  filter(out_var == "cell_type") %>%  # eliminate things redundant to this analysis
  mutate(tfbm_results = pmap(select(., rhs = formula1, 
                                    of_in = variable), 
                             get_tfbm_p_vals, dat = my_obj)) 


