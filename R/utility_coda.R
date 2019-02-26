
my_lm  = function(model_formulae){
  
  model_formulae %>% 
    lm(data  = phen, na.action = na.omit) %>% 
    manova %>% 
    tidy
  
} 

my_lm2  = function(model_formulae){
  
  model_formulae %>% 
    lm(data  = phen, na.action = na.omit) %>% 
    car::Manova()  
  # model_formulae %>% ffmanova::ffmanova(data  = phen, na.action = na.omit)  # no missing value handling, buggy
  
}

shifter =
  function(x, n = 1) {
    # ALL CYCLIC REORDERINGS OF A VECTOR
    if (n == 0) x else c(tail(x, -n), head(x, n))
  }


get_permutate_p = 
  function(model_formulae, outcomes, controls, focal_treatment,  n_perm = 100, pretty=F){
    
    ########################################################
    # PARAMETRIC CLASSICAL P VALUE
    ########################################################
    
    m        <- lm(model_formulae, data = phen, na.action = na.omit)  
    p_par    <- m %>% anova %>% tidy %>% filter(str_detect(term, focal_treatment)) %>% select(p.value) %>% as.numeric 
    
    ########################################################
    # NON PARAMETRIC PERMUTATION P-VAL FOR A SET OF CANDIT
    ########################################################
    
    perms    <- permute(phen, n_perm, focal_treatment)
    null_dn  <- perms$perm %>% map_dbl(get_statistic, model_formulae, focal_treatment)
    observed <- get_statistic(phen, model_formulae, focal_treatment)
    p_perm   <- mean(null_dn > observed)
    
    y = list(p_perm = p_perm,  p_par = p_par )
    
    if(pretty){
      tibble(null_dn) %>% 
        gf_histogram(~null_dn) %>% 
        gf_vline(xintercept = observed) %>% 
        gf_labs(title = str_c(outcomes, "~", focal_treatment, ": pval = ", p, collapse = T))
      ggsave(str_c("figs/univariate_", outcomes, "~", focal_treatment, ".pdf", collapse = "")) # don't save 
      graphics.off()
    }
    
    return(y)
    
  }

get_statistic =
  function(x, model_formulae, focal_treatment) {
    # calculate squared norm (sse)  
    lm(model_formulae, data = x, na.action = na.omit) %>%
      tidy %>% 
      filter(str_detect(term, focal_treatment))  %>% 
      select(estimate) %>% unlist %>% `^`(2) %>% sum
  }

########################################################
# Parametric omnibus/model tests
########################################################

model_compare_par = 
  function(formula0, formula1, outcome){
    
    y = outcome
    
    formula0 = str_c("y ~", formula0) %>% as.formula
    formula1 = str_c("y ~", formula1) %>% as.formula
    
    design = model.matrix(formula1, data = phen)
    keep   = complete.cases(design)
    y      = y[keep,]
    alt    = lm(formula1, data = phen[keep, ])
    null   = lm(formula0, data = phen[keep, ])
    anova(alt,  null)
    
  }

########################################################
# NON parametric omnibus/model tests
########################################################

model_compare_perm = 
  function(formula0, formula1, outcome, outcome_var, n_perm=100, vars_only_in_the_alt_hyp = NULL){
    
    y = outcome 
    
    formula0 = str_c("y ~", formula0) %>% as.formula
    formula1 = str_c("y ~", formula1) %>% as.formula
    
    design = model.matrix(formula1, data = phen)
    keep   = complete.cases(design)
    y      = y[keep,]
    phen_not_na = phen[keep, ]
    
    # if focal term is named use it, else take the final term in the equ
    vars_only_in_the_alt_hyp = ifelse(!is.null(vars_only_in_the_alt_hyp), 
                                      vars_only_in_the_alt_hyp, 
                                      formula1 %>% str_sub(str_length(formula0) + 1, str_length(formula1)) %>%
                                        str_split("\\+") %>%  unlist %>% str_trim() %>% str_subset("") )
    
    # vars_only_in_the_alt_hyp = formula1 %>% str_sub(str_length(formula0) + 1, str_length(formula1)) %>% str_split("\\+") %>%  unlist %>% str_trim() %>% str_subset("")
    
    print(vars_only_in_the_alt_hyp)
    
    my_perms = phen_not_na %>% permute(n_perm, vars_only_in_the_alt_hyp) 
    
    # browser()
    library(furrr)
    plan(multiprocess, workers = detectCores()/2 - 2)
    perm_stat = 
      map_dbl(
        my_perms$perm,
        function(x){
          alt    = lm(formula1, data = x)
          null   = lm(formula0, data = x)
          anova(alt,  null)$Pillai[2]
        }
      )
    return(perm_stat = perm_stat)
  }
