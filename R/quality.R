source("R/utility.R")
# load 11 duplicate participants and normalize 
dat_dupes = readRDS("data/dt_dupes.rds") # duplicated subjects only
exprs(dat_dupes) <- dat_dupes %>% exprs %>% cpm(log = TRUE) 

########################################################
# MDS: TECHNICAL REPLICATES ARE CLOSE OVERALL
######################################################## 

# A plot can be made in limma using the plotMDS function. The first
# dimension represents the leading-fold-change that best separates samples and
# explains the largest proportion of variation in the data, with subsequent
# dimensions having a smaller effect and being orthogonal to the ones before it.
# When experimental design involves multiple factors, it is recommended that
# each factor is examined over several dimensions. If samples cluster by a given
# factor in any of these dimensions, it suggests that the factor contributes to
# expression differences and is worth including in the linear modelling. On the
# other hand, factors that show little or no effect may be left out of
# downstream analysis.
# source: https://www.bioconductor.org/packages/devel/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html

# In this dataset, samples can be seen to cluster well within experimental
# groups over dimension 1 and 2, and then separate by sequencing lane (sample
# batch) over dimension 3 (shown in the plot below). Keeping in mind that the
# first dimension explains the largest proportion of variation in the data,
# notice that the range of values over the dimensions become smaller as we move
# to higher dimensions.


lcpm <- dat_dupes %>% exprs
library(RColorBrewer) 
brewer.pal.info
group = dat_dupes$AID %>% factor
par(mfrow=c(1,2))
col.group <- group
levels(col.group) <-  brewer.pal(nlevels(col.group), "PuOr")
col.group <- as.character(col.group)
plotMDS(lcpm, labels=group, col=col.group)
title(main="A. Duplicate pairs are similar")
group = dat_dupes$male %>% factor
col.group <- group
levels(col.group) <-  c("red", "green")
col.group <- as.character(col.group)
plotMDS(lcpm, labels=group, col=col.group)
title(main="B. Sex")

########################################################
# SEX CHECKS
########################################################

########################################################
# DUPLICATES HAVE THE SAME INFERED GENDER
########################################################

# Using only non-duplicates to identify the sex genes 
counts = dat[, !(dat$VialID %in% dat_dupes$VialID)] %>% exprs 
design = dat[, !(dat$VialID %in% dat_dupes$VialID)] %>% pData
design = model.matrix(~factor(design$male))
(sex_genes = my_p_val_head(counts, design)) 

########################################################
# USE THE ABOVE GENE SET TO 
########################################################
# Is the infered sex of the two supposed duplicates equal?

# extract their values on the top sex genes above
(sexy = dat_dupes[rownames(sex_genes),])
sexy %>%
  tidy(addPheno = TRUE) %>% 
  gf_point(value~factor(AID), 
           color = ~factor(AID),
           shape = ~factor(male)) %>% 
  gf_facet_wrap(~gene) %>% 
  gf_labs(title = "Duplicates are of the same infered biological gender")
