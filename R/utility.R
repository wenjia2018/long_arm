########################################################
# convenient
#########################################################

eData = exprs

# P-value distribution
my_p_val_head  = function(counts, design, n = 10){
  keep <- filterByExpr(counts, design) # For good results, the counts matrix should be filtered to remove remove rows with very low counts before running voom(). The filterByExpr function in the edgeR package can be used for that purpose.
  voom(counts = counts[keep, ],
       design = design) %>%
    # arrayWeights %>%
    lmFit %>%
    eBayes %>%
    topTable(number = n)
}# P-value distribution

my_p_val_plot = function(counts, design){
  keep <- filterByExpr(counts, design) # For good results, the counts matrix should be filtered to remove remove rows with very low counts before running voom(). The filterByExpr function in the edgeR package can be used for that purpose.
  voom(counts = counts[keep, ],
       design = design) %>%
    # arrayWeights %>%
    lmFit %>%
    eBayes %>%
    topTable(number = Inf) %>%
    gather(k,v, c("P.Value", "adj.P.Val")) %>%
    gf_histogram(~v) %>%
    gf_facet_wrap(~k)
}

my_p_val_pimped = function(counts, design, n_subsamps  = dim(design)[1]){
  
  if(0){
    ########################################################
    # REPRODUCE VOOM PREPROC?
    ########################################################
    
    (examine_mean_variance_relation =
       counts %>%  
       as_tibble %>%
       mutate_all(~log(.5 + .)) %>% # toggle
       mutate(mean = apply(.,1,mean),
              sd   = apply(.,1,sd)))  %>%
      gf_point(sd~mean) %>%
      gf_labs(title = "log cpm")
    
    ########################################################
    # VOOM
    ########################################################
    
    # Estimates relative quality weights for each array in a multi-array experiment.
    v =
      voom(counts = counts,
           design = design,
           plot=TRUE)
  }
  
  ########################################################
  # VOOM
  ########################################################
  
  # As above, but plot = FALSE
  v = voom(counts = counts,
           design = design)
  sample_w = v %>% arrayWeightsSimple
  
  # DOES NOT SCALE WELL! (HENCE SUBSAMPLING)
  subsamps = sample(dim(design)[1], n_subsamps)
  out = 
    voomWithQualityWeights(counts = counts[, subsamps],
                           design = design[subsamps, ]) %>%
    lmFit(weights = sample_w[subsamps]) %>%
    eBayes %>%
    topTable(number = Inf) 
  
}


my_cross_validation_logistic = function(phenotype, exprs_set_object){
  
  ########################################################
  # Simple N-fold cross validation for binary classification
  ########################################################
  D = 
    exprs(exprs_set_object) %>%
    cpm(log = TRUE) %>% 
    t %>%
    as_tibble %>%
    cbind(outcome = phenotype) %>% # DANGER
    as_tibble
  
  output = 
    D %>% 
    na.omit %>% 
    group_by(outcome) %>%  sample_n(300) %>%  # to avoid class imbalance??!
    crossv_kfold(10) %>% 
    mutate(model = map(train, ~ glm(outcome  ~ ., data = ., family = binomial))) %>% 
    unnest(map2(model, test, ~ augment(.x, newdata = .y, type.predict = "response"))) %>% 
    mutate(hard_pred = .fitted > .5)
  
  output %>% 
    gf_histogram(~.fitted)
  
  print( output %>% 
           select(hard_pred, outcome) %>% 
           table)
  
  break_it_down = FALSE
  if(break_it_down){
    D %>% 
      # modelr::crossv_kfold function takes a data frame and randomly partitions
      # it’s rows into k roughly equal groups. We’ve partitioned the row numbers
      # into k = 5 groups. The results are returned as a tibble (data frame) like
      # the one above. Each cell in the test column contains a resample object,
      # which is an efficient way of referencing a subset of rows in a data frame.
      # (Storing the complete resample each time would be very inefficient so this
      # class instead stores a "pointer" to the original dataset, and a vector of
      # row indexes. To turn this into a regular data frame, call as.data.frame, to
      # extract the indices, use as.integer.
      crossv_kfold(10) %>% 
      # Fit model to each training partition (adds a model column)
      mutate(model = map(train, ~ glm(sex_interv ~ ., family = binomial, data = .))) %>% 
      # General way to add prediction, then unnest
      mutate(predicted = map2(model, test, predict, type = "response")) %>%
      unnest(predicted) %>%
      gf_histogram(~predicted)
    # # OR equivalent shortcut for some models
    # unnest(map2(model, test, ~ augment(.x, newdata = .y, type.predict = "response"))) %>%
    # gf_histogram(~.fitted)
  }
  if(0){ 
    output %>% 
      group_by(.id) %>% 
      # 0-1 loss, based on a hard - most probable - predicted sex > .5
      mutate(zero_one_loss = (hard_prediction & (sex_interv == "f")) | (!hard_prediction & (sex_interv == "m"))) %>% 
      summarize(within_fold_classification_success = mean(zero_one_loss)) %>% # within fold classification success
      summarize(mean_within_fold_classification_success = mean(within_fold_classification_success))
  }
}
