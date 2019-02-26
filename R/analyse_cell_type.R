########################################################
# SPECIFY RHS OF PREDICTIVE MODEL
########################################################

source("R/utility_coda.R")

########################################################
# CELL TYPE ANALYSIS
########################################################

n_sims = 1000

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
# THE BASIC RESULT
########################################################

(m =  
  lm(str_c("cell_type[keep, ] ~ ", rhs), 
   data = phen[keep, ]))
m %>%  car::Manova() 


########################################################
# ELABORATIONS
########################################################

# partial effect
b = coef(lm(str_c("cell_type[keep, ] ~ ", rhs), data = phen[keep, ]))[of_in, ] %>% ilrInv
names(b) = cell_types_one
barplot(b)
barplot(5*b)

# bootstrap uncertainty
my_boot = function(){
  y = cell_type[keep, ]
  x = phen[keep, ] 
  resampled = sample(nsubs, nsubs, replace = T) # with replacement
  y = y[resampled, ]
  x = x[resampled, ]
  coef(lm(as.formula(str_c("y ~ " , rhs, collapse = " + ")), data = x))["obesew5", ] %>% ilrInv
  # rhs_n   = str_c(c(vars[-5], vars[5]), collapse = " + ")
  # coef(lm(str_c("y ~ " , rhs_n, collapse = " + "), data = x))["obesew5", ] %>%  ilrInv
  # coef(lm(str_c("y ~ obesew5"), data = x))["obesew5", ] %>%  ilrInv  # without control
}
# sample with replacement
boots =  rerun(n_sims, my_boot())
#plot
boots %>% purrr::reduce(rbind) %>% acomp() %>% `*`(10) %>% plot
plot(acomp(c(1,1,1,1,1,1)), add = T, col = "red")


# permutation "test"
my_perm = function(){
  y = cell_type[keep, ]
  x = phen[keep, ] 
  permutations = sample(nsubs, nsubs, replace = F) # without replacement
  x$obesew5 = x$obesew5[permutations] # scramble the predictor
  coef(lm(as.formula(str_c("y ~ " , rhs, collapse = " + ")), data = x))["obesew5", ] %>%  ilrInv 
  # coef(lm(str_c("y ~ obesew5"), data = x))["obesew5", ] %>%  ilrInv  # without control
}
perms = rerun(n_sims, my_perm())
# plot scaled
perms %>% purrr::reduce(rbind) %>% acomp() %>% `*`(10) %>% plot
plot(10*b, add = T, col = "red")
# plot unscaled
perms %>% purrr::reduce(rbind) %>% acomp() %>% plot
plot(b, add = T, col = "red")
 