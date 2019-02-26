source("R/analyse_count.R")

if(0){
  ########################################################
  # MASS UNIVARIATE SCREEN
  ########################################################
  
  source("R/my_coda_regression_DV_general.R")
  
  ########################################################
  # MASS UNIVARIATE SCREEN
  ########################################################
  
  source("R/my_coda_regression_DV_nonparametric.R")
  
  
  ########################################################
  # TIDY 
  ########################################################
  
  dat[c(ctra_available, housekeepers11), ] %>% 
    tidy  %>% 
    spread(gene, value) 
  
  ########################################################
  # MDS
  ########################################################
  
  source("R/my_mds.R")
  
  ########################################################
  # REPRODUCE STEVE's PROPOSED NORMALIZATION 
  ########################################################
  
  if(0){
    normalize_to_ctra_total = FALSE
    if(normalize_to_ctra_total){ 
      # normalize to CTRA total
      dat_preproc       = readRDS("../data/dat.rds") # expression set
      tcounts_ctra = log2(pmax(cpm(exprs(dat[ctra_available,])), 1)) 
      tcounts_steve_ctra = exprs(dat_preproc[ctra_available,]) 
    } else { 
      # normalize to total read count
      dat_preproc       = readRDS("../data/dat.rds") # expression set
      tcounts = log2(pmax(cpm(exprs(dat)), 1)) # me, trying to reproduce steve
      tcounts_steve = exprs(dat_preproc) # steve
      # just convenience: subset to CTRA
      tcounts_ctra       = tcounts[ctra_available,] # me 
      tcounts_steve_ctra = tcounts_steve[ctra_available,] # steve
      # CTRA expressions agree pretty well
      all.equal(tcounts_ctra, tcounts_steve_ctra) 
      c(tcounts_ctra - tcounts_steve_ctra) %>% 
        hist(main = "CTRA differences are very small: numeric precision") 
    }
  }
  
  ########################################################
  # NORMALIZATION, STEVE's WAY
  # THE EFFECT OF THE GENE SET (AS OPPOSSED TO THE FORMULA)
  ########################################################
  
  # Reproducing Steve's preprocessing (see above)
  phen  <- 
    pData(dat) %>%
    `rownames<-`(.$VialID) 
  dat_normalize_to_all <- 
    ExpressionSet(assayData = log2(pmax(cpm(exprs(dat)), 1)), # steve's preprocessing
                  phenoData = new("AnnotatedDataFrame", data = phen))  
  
  # Same but normalizing to total CTRA read count (not total overall)
  # First create an expressionset object
  phen  <- 
    pData( dat[ctra_available,]) %>%
    `rownames<-`(.$VialID) 
  dat_normalize_to_ctra <- 
    ExpressionSet(assayData = log2(pmax(cpm(exprs(dat[ctra_available,])), 1)), # steve's preprocessing, but on the CTRA subset
                  phenoData = new("AnnotatedDataFrame", data = phen))  
  
  # Same but normalizing to total CTRA read count (not total overall)
  # First create an expressionset object
  phen  <-
    pData( dat ) %>%
    `rownames<-`(.$VialID)
  dat_normalize_limma <- 
    ExpressionSet(assayData = dat %>% exprs %>% normalizeBetweenArrays, # steve's preprocessing, but on the CTRA subset
                  phenoData = new("AnnotatedDataFrame", data = phen))
  
  ########################################################
  # combine and compare different normalizations
  ########################################################
  # CHECK filterByExpr documentation
  badgenes = filterByExpr(dat[ctra_available,])[filterByExpr(dat[ctra_available,]) == FALSE] %>% names()
  tidy(dat_normalize_to_all[ctra_available,]) %>% 
    left_join(tidy(dat_normalize_to_ctra), by = c("gene", "sample")) %>% 
    left_join(tidy(dat_normalize_limma[ctra_available]), by = c("gene", "sample")) %>% 
    dplyr::rename(to_all = "value.x", to_ctra = "value.y", limma = "value")  %>% 
    select(-limma) %>% # limma is totally different!
    dplyr::filter(!(gene %in% badgenes)) %>% 
    gather(k,v, -c(gene, sample)) %>% 
    gf_histogram(~v, color =~ k, bins = 100, alpha = .2, position = position_dodge()) %>% 
    gf_facet_wrap(~gene) %>% 
    gf_theme(theme_minimal())
  
  # combine and compare both normalizations
  dat_normalize_to_all[ctra_available,] %>%
    tidy %>% 
    left_join(dat_normalize_to_ctra %>% tidy, by = c("gene", "sample")) %>% 
    gf_point(value.x~value.y) %>% 
    gf_facet_wrap(~gene)
  
  ########################################################
  # PICK ONE NORMALIZATION AND EXAMINE BIVARIATE MARGINS
  ########################################################
  
  # take 5 at a time and examine their bivariate margins
  map(list(1:3, 4:10), 
      ~dat_normalize_to_ctra[.x,] %>% 
        exprs %>%
        t %>% 
        pairs)
  
  ########################################################
  # TIDY VISUALIZATIONS
  ########################################################
  
  # All this only necessary because addPheno argement in tidy(., addPheno = TRUE)
  # is crap. It overwrites all phenodata with NA.
  # genes = inflammatory
  genes = ctra_available
  phen_tidy =
    dat_normalize_to_ctra[genes, ] %>% 
    pData %>% 
    mutate(sample = VialID) # to match the below
  dat_tidy  =
    dat_normalize_to_ctra[genes, ] %>% 
    tidy %>% 
    left_join(phen_tidy, by = "sample")
  dat_tidy %>% 
    gf_violin(value~edu_p) %>% 
    gf_facet_wrap(~gene)
  
  if(0){
    
    ########################################################
    #  OUTLIER AND HETEROSKEDASTISITY
    ########################################################
    
    # PCA: ON WHICH VARIABLES and FEATURES (ctra normalized, unnormalized, etc)
    # On CTRA set
    # keep only higher counts
    keep <- filterByExpr(dat[ctra_available,])  
    t_counts = t(dat[ctra_available[keep],] %>% exprs)
    # On some other input space
    # keep only higher counts AND and only interpretable genes "ENSG" genes
    keep = -str_which(featureNames(dat), "^ENSG.*") # more concise than below but not yet confirmed to not break something downstream. checked (hence overridden next)
    keep <- filterByExpr(dat) & (!(featureNames(dat) %in% featureNames(dat)[str_detect(featureNames(dat), "ENSG")]))
    t_counts = t(dat[keep,] %>% exprs)
    t_counts = t(dat_normalize_limma[keep,] %>% exprs)
    t_counts = t(dat_normalize_to_all[keep,] %>% exprs)
    # Normal and robust PCA
    normal_pca = t_counts %>% prcomp %>% augment(data = cbind(t_counts, pData(dat)))
    robust_pca = t_counts %>% pcaPP::PCAproj(k = 2)
    # normal pca plot
    normal_pca %>%
      gf_point(.fittedPC1 ~ .fittedPC2, color =~factor(male)) %>%
      gf_density2d(.fittedPC1 ~ .fittedPC2 )  %>%  
      # gf_text(label = ~.rownames, vjust = 1, hjust = 1)
      ggExtra::ggMarginal(type = "histogram")
    # robust pca
    robust_pca$scores %>% 
      as_tibble %>%
      rownames_to_column %>%
      cbind(pData(dat)) %>%  # danger
      gf_point(Comp.1~Comp.2, color =~factor(male)) %>%
      gf_density2d(Comp.1 ~ Comp.2)  %>%  
      # gf_text(label = ~VialID, vjust = 1, hjust = 1) %>% 
      ggExtra::ggMarginal(type = "histogram")
    
    ########################################################
    # MDS
    ########################################################
    
    # Plot samples on a two-dimensional scatterplot so that distances on the plot
    # approximate the typical log2 fold changes between the samples. This function
    # is a variation on the usual multdimensional scaling (or principle coordinate)
    # plot, in that a distance measure particularly appropriate for the microarray
    # context is used. The distance between each pair of samples (columns) is the
    # root-mean-square deviation (Euclidean distance) for the top top genes.
    # Distances on the plot can be interpreted as leading log2-fold-change, meaning
    # the typical (root-mean-square) log2-fold-change between the samples for the
    # genes that distinguish those samples.
    
    t_counts %>% t %>% plotMDS # nb this expects counts, not t(counts)
    
    
    ########################################################
    # HEATMAP SLOW
    ########################################################
    
    # descriptive clustering
    counts[keep,] %>% heatmap
    
  }
  
  ########################################################
  # DECODING METHODS
  ########################################################
  ########################################################
  # CV ACCURACY OF PHENOTYPE PREDICTION FROM A PRIORI GENE SET
  ########################################################
  
  # sex from sex genes
  my_cross_validation_logistic(dat$sex_interv, dat[rownames(sex_genes), ])
  # parental education from inflamation
  my_cross_validation_logistic(dat$edu_p != "high", dat[ctra_available, ])
  my_cross_validation_logistic(dat$edu_p != "high", dat[inflammatory, ])
  
  ########################################################
  # GLOBAL TEST OF PREDICTION ACCURACY
  ########################################################
  
  outcome = dat$edu_p
  outcome_not_missing = !is.na(outcome)
  Y = outcome[outcome_not_missing] 
  X = dat[inflammatory,outcome_not_missing] %>% 
    exprs %>% 
    cpm(log = TRUE) %>% 
    t %>% 
    as_tibble
  mod = gt(Y, ~ ., ~ 1, data = X)
  (tab1 <- mod %>% summary)
  (fig1 <- mod %>% covariates)
  
  # But not specific to CTRA
  (junk = sample(featureNames(dat), length(ctra_available)))
  X = dat[junk,outcome_not_missing] %>% 
    exprs %>% 
    t %>% 
    as_tibble
  gt(Y, ~ ., ~1, data = X) %>% summary
  gt(Y, ~ ., ~1, data = X) %>% covariates 
  
  outcome = dat$sss_5 - dat$sss_4
  outcome = dat$edu_max
  
  ########################################################
  # ENCODING METHODS
  ########################################################
  ########################################################
  # WHOLE TRANSCRIPTOME 
  ########################################################
  # featureNames(dat)[str_detect(featureNames(dat), "ENSG")]) 
  
  # data admissible for a design looking at ...
  design = model.matrix(~ dat$edu_p ) 
  counts = dat[filterByExpr(dat), complete.cases(design)] %>% exprs 
  design = design[complete.cases(design), ]
  (out = my_p_val_head(counts, design))
  (out = my_p_val_pimped(counts, design))
  out %>% rownames_to_column %>% filter(adj.P.Val < 0.05)
  
  # data admissible for a design looking at ...
  design = model.matrix(~ dat$w5bmi) 
  counts = dat[filterByExpr(dat), complete.cases(design)] %>% exprs 
  design = design[complete.cases(design), ]
  (out = my_p_val_plot_pimped(counts, design))
  
  if(0){
    # whole genome: memory fails
    X = dat %>% 
      exprs %>% 
      t %>% 
      as_tibble
    gt(Y, ~ ., ~1, data = X) %>% summary
    gt(Y, ~ ., ~1, data = X) %>% covariates
  }
  
  if(0){
    ########################################################
    # fake data
    ########################################################
    
    (counts = matrix(rnbinom(10*10000,size=1,mu=10),10000,10))
    (design = model.matrix(~rep(1:2, each = 5) %>% factor))
    
    # counts per million
    (y <- matrix(rnbinom(20,size=1,mu=10),5,4))
    y[,1]/sum(y[,1]) * 1000000
    # shortcut
    cpm(y)
    # other normalization methods
    (counts =  counts %>% normalizeBetweenArrays)
    
    my_func(counts, design)
    
  }
  
  out = list(fig1 = fig1,
             tab1 = tab1)
  return(out = out)
}



