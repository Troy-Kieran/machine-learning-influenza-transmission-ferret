########################################################################################
########################################################################################
### Transmission spinoff
### Influenza-Ferret Predictive Analytics
### Identification of optimal thresholds and key parameters for predicting 
###   influenza A virus transmission events in ferrets using machine learning 
###   on virological data
### 11 March 2024 - 30 October 2024
### Troy J. Kieran
########################################################################################
########################################################################################
### Load packages
library(tidyverse)
library(tidylog) ## detailed log of dplyr/tidyr functions
library(funModeling) ## loads Hmisc
library(caret) ## modeling
library(caretEnsemble) ## ensemble/stacked models
library(patchwork)
library(gbm)
library(glmnet)

########################################################################################
### Import Data

## replace file.csv with name of file

## download data
# Pathogenesis Laboratory Team, Influenza Division, CDC. 
# An aggregated dataset of serially collected influenza A virus morbidity and 
# titer measurements from virus-infected ferrets.  
# https://data.cdc.gov/National-Center-for-Immunization-and-Respiratory-D/An-aggregated-dataset-of-serially-collected-influe/cr56-k9wj/about_data
#
# fullData <- read.csv("file1.csv", 
#                      header = TRUE, check.names = FALSE)
# 
# ## validation with published literature data
# Compiled from:
# Kieran TJ, et al. Machine learning approaches for influenza A virus 
#   risk assessment identifies predictive correlates using ferret model 
#   in vivo data. Commun Biol 7, 927 (2024).
#   
# pubLit <- read.csv("file2.csv", header = TRUE) 
# 
# ## validation with standardization exercise data
# Compiled from:
# Belser JA, et al. Robustness of the Ferret Model for Influenza Risk 
#   Assessment Studies: a Cross-Laboratory Exercise. mBio 13, e0117422 (2022).
#   
# standardized <- read.csv("file3.csv", header = TRUE) 

########################################################################################

trans_dataFilt <- fullData

janitor::tabyl(trans_dataFilt$RD_33)
janitor::tabyl(trans_dataFilt$RD_50)
janitor::tabyl(trans_dataFilt$RD_67)
janitor::tabyl(trans_dataFilt$DC_33)
janitor::tabyl(trans_dataFilt$DC_50)
janitor::tabyl(trans_dataFilt$DC_67)

##

## slope1-3 categorical

## explore slopes
slope <- unique(fullData$slope13_v)
intercept <- rep(0, length(slope))
x <- rep(1:5, length.out = length(slope))
y <- c(rep(seq(from = -0.7, to = 0.7, by = 0.1), length.out = length(slope)))
df <- data.frame(cbind(slope, intercept, x, y))

df <- df %>% 
  mutate(direction = ifelse(slope < 0, 'neg', 
                            ifelse(slope > 0, 'pos', 'zero')),
         direction2 = ifelse(slope < -0.005, 'neg', 
                             ifelse(slope > 0.005, 'pos', 'zero'))) %>%
  drop_na()

janitor::tabyl(df$direction)
janitor::tabyl(df$direction2)

ggplot(data = df) +
  geom_point(alpha = 0, aes(x, y)) +
  theme_bw() +
  geom_abline(data = df, aes(slope = slope, intercept = intercept, color = factor(direction2))) +
  geom_abline(aes(slope = 0, intercept = 0), size = 2, linetype = 'dotted', alpha = 0.7) +
  xlim (2.5, 3.5) +
  ylim(0.23, -0.23)

## add categorical to trans_dataFilt
trans_dataFilt <- trans_dataFilt %>% 
  mutate(slope13_vcat = ifelse(slope13_v < -0.005, 'neg', 
                             ifelse(slope13_v > 0.005, 'pos', 'zero')),
         slope13_fcat = ifelse(slope13_f < -0.005, 'neg', 
                              ifelse(slope13_f > 0.005, 'pos', 'zero')))
janitor::tabyl(trans_dataFilt$slope13_vcat)
janitor::tabyl(trans_dataFilt$slope13_fcat)

## add categorical to standardized data
standardized <- standardized %>%
  mutate(slope13_vcat = ifelse(slope13_v < -0.005, 'neg', 
                               ifelse(slope13_v > 0.005, 'pos', 'zero')),
         slope13_fcat = ifelse(slope13_f < -0.005, 'neg', 
                               ifelse(slope13_f > 0.005, 'pos', 'zero')))

###

## external validation data
pubLit[pubLit == ""] <- NA

pubLit <- pubLit %>%
  rename(Origin_orig = Origin) %>%
  mutate(Origin = if_else(Origin_orig == "avian" , "avian", "mammal")) %>%
  group_by(PMID, Virus) %>%
  ## convert yes/no to percentage by Virus
  mutate(RD_trans_yes_p = mean(RD_trans == 'yes', na.rm = TRUE),
         DC_trans_yes_p = mean(DC_trans == 'yes', na.rm = TRUE)) %>%
  fill(slope13_vcat, .direction = "down" ) %>%
  ungroup() %>%  
  ## Handle missing values in 'RD_trans_yes_p' and 'DC_trans_yes_p'
  ## Replace them with NA if the corresponding 'RD_trans_yes' or 'DC_trans_yes' value is NA
  mutate(RD_trans_yes_p = ifelse(is.na(RD_trans), NA_character_, RD_trans_yes_p) %>% as.numeric(),
         DC_trans_yes_p = ifelse(is.na(DC_trans), NA_character_, DC_trans_yes_p) %>% as.numeric()) 

pubLit <- pubLit %>%
  mutate(RD_33 = ifelse((RD_trans_yes_p < 0.33 | 
                           is.na(RD_trans_yes_p) & DC_trans_yes_p < 0.1), "no",
                        ifelse(RD_trans_yes_p > 0.33, "yes", "no")),
         RD_50 = ifelse((RD_trans_yes_p < 0.5 | 
                           is.na(RD_trans_yes_p) & DC_trans_yes_p < 0.1), "no",
                        ifelse(RD_trans_yes_p > 0.5, "yes", "no")),
         RD_67 = ifelse((RD_trans_yes_p < 0.67 | 
                           is.na(RD_trans_yes_p) & DC_trans_yes_p < 0.1), "no",
                        ifelse(RD_trans_yes_p > 0.67, "yes", "no")),
         DC_33 = ifelse(DC_trans_yes_p < 0.33 & !is.na(DC_trans_yes_p), "no",
                        ifelse(DC_trans_yes_p > 0.33 | 
                                 (is.na(DC_trans_yes_p) & RD_trans_yes_p == 1), "yes", "no")),
         DC_50 = ifelse(DC_trans_yes_p < 0.5 & !is.na(DC_trans_yes_p), "no",
                        ifelse(DC_trans_yes_p > 0.5 | 
                                 (is.na(DC_trans_yes_p) & RD_trans_yes_p == 1), "yes", "no")),
         DC_67 = ifelse(DC_trans_yes_p < 0.67 & !is.na(DC_trans_yes_p), "no",
                        ifelse(DC_trans_yes_p > 0.67 | 
                                 (is.na(DC_trans_yes_p) & RD_trans_yes_p == 1), "yes", "no"))) 

###

standardized[standardized == ""] <- NA

standardized <- standardized %>%
  drop_na(slope13_v) %>%
  group_by(Virus, Group) %>%
  ## convert yes/no to percentage by Virus
  mutate(RD_trans = ifelse(RD_titer == 'yes' &
                             RD_sero == 'yes', 'yes', 'no')) %>%
  mutate(RD_trans_yes_p = mean(RD_trans == 'yes', na.rm = TRUE)) %>%
  #fill(slope13_vcat, .direction = "down" ) %>%
  ungroup() %>%  
  ## Handle missing values in 'RD_trans_yes_p' and 'DC_trans_yes_p'
  ## Replace them with NA if the corresponding 'RD_trans_yes' or 'DC_trans_yes' value is NA
  mutate(RD_trans_yes_p = ifelse(is.na(RD_trans), NA_character_, RD_trans_yes_p) %>% as.numeric()) 

standardized <- standardized %>%
  mutate(RD_33 = ifelse(RD_trans_yes_p < 0.33, "no", "yes"),
         RD_50 = ifelse(RD_trans_yes_p <= 0.5, "no", "yes"),
         RD_67 = ifelse(RD_trans_yes_p < 0.67, "no", "yes"))


########################################################################################

## Matthew's Correlation Coefficient Function
## manual calculation/function
## T = true, F = false, P = positive, N = negative
## use data from confusion matrix
mcc_func <- function(TP, FP, FN, TN){
  ((TP*TN) - (FP*FN))/
    sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))}

########################################################################################

## Predictive Power Scores

## RD_33
trans_dataFilt %>%
  dplyr::select(RD_33, AUC_4_v, AUC_6_v, AUC_8_v, AUC_9_v, AUC_RD,
                AUC_4_f, AUC_6_f, AUC_8_f, AUC_9_f, peak_inoc) %>%
  ppsr::score_predictors(df = ., y = 'RD_33')

trans_dataFilt %>%
  dplyr::select(RD_33, slope13_v, slope24_v, slope35_v, slope46_v, slope57_v, slope68_v,
                slope13_f, slope24_f, slope35_f, slope46_f, slope57_f, slope68_f, 
                slope13_vcat, slope13_fcat) %>%
  ppsr::score_predictors(df = ., y = 'RD_33')

trans_dataFilt %>%
  dplyr::select(RD_33, Origin, HA, Subtype, RBS, PBS) %>%
  ppsr::score_predictors(df = ., y = 'RD_33')

# trans_dataFilt %>%
#   dplyr::select(-c(CDCID, Virus, ends_with('_euth'), 17:PB2_701)) %>%
#   ppsr::score_predictors(df = ., y = 'RD_33') %>% View()

###

## RD_50
trans_dataFilt %>%
  dplyr::select(RD_50, AUC_4_v, AUC_6_v, AUC_8_v, AUC_9_v, AUC_RD,
                AUC_4_f, AUC_6_f, AUC_8_f, AUC_9_f, peak_inoc) %>%
  ppsr::score_predictors(df = ., y = 'RD_50')

trans_dataFilt %>%
  dplyr::select(RD_50, slope13_v, slope24_v, slope35_v, slope46_v, slope57_v, slope68_v,
                slope13_f, slope24_f, slope35_f, slope46_f, slope57_f, slope68_f, 
                slope13_vcat, slope13_fcat) %>%
  ppsr::score_predictors(df = ., y = 'RD_50')

trans_dataFilt %>%
  dplyr::select(RD_50, Origin, HA, Subtype, RBS, PBS) %>%
  ppsr::score_predictors(df = ., y = 'RD_50')

###

## RD_67
trans_dataFilt %>%
  dplyr::select(RD_67, AUC_4_v, AUC_6_v, AUC_8_v, AUC_9_v, AUC_RD,
                AUC_4_f, AUC_6_f, AUC_8_f, AUC_9_f, peak_inoc) %>%
  ppsr::score_predictors(df = ., y = 'RD_67')

trans_dataFilt %>%
  dplyr::select(RD_67, slope13_v, slope24_v, slope35_v, slope46_v, slope57_v, slope68_v,
                slope13_f, slope24_f, slope35_f, slope46_f, slope57_f, slope68_f, 
                slope13_vcat, slope13_fcat) %>%
  ppsr::score_predictors(df = ., y = 'RD_67')

trans_dataFilt %>%
  dplyr::select(RD_67, Origin, HA, Subtype, RBS, PBS) %>%
  ppsr::score_predictors(df = ., y = 'RD_67')

### ### ### ### ### ###

## DC_33
trans_dataFilt %>%
  dplyr::select(DC_33, AUC_4_v, AUC_6_v, AUC_8_v, AUC_9_v, AUC_RD,
                AUC_4_f, AUC_6_f, AUC_8_f, AUC_9_f, peak_inoc) %>%
  ppsr::score_predictors(df = ., y = 'DC_33')

trans_dataFilt %>%
  dplyr::select(DC_33, slope13_v, slope24_v, slope35_v, slope46_v, slope57_v, slope68_v,
                slope13_f, slope24_f, slope35_f, slope46_f, slope57_f, slope68_f, 
                slope13_vcat, slope13_fcat) %>%
  ppsr::score_predictors(df = ., y = 'DC_33')

trans_dataFilt %>%
  dplyr::select(DC_33, Origin, HA, Subtype, RBS, PBS) %>%
  ppsr::score_predictors(df = ., y = 'DC_33')

# trans_dataFilt %>%
#   dplyr::select(-c(CDCID, Virus, ends_with('_euth'), 17:PB2_701)) %>%
#   ppsr::score_predictors(df = ., y = 'DC_33') %>% View()

###

## DC_50
trans_dataFilt %>%
  dplyr::select(DC_50, AUC_4_v, AUC_6_v, AUC_8_v, AUC_9_v, AUC_RD,
                AUC_4_f, AUC_6_f, AUC_8_f, AUC_9_f, peak_inoc) %>%
  ppsr::score_predictors(df = ., y = 'DC_50')

trans_dataFilt %>%
  dplyr::select(DC_50, slope13_v, slope24_v, slope35_v, slope46_v, slope57_v, slope68_v,
                slope13_f, slope24_f, slope35_f, slope46_f, slope57_f, slope68_f, 
                slope13_vcat, slope13_fcat) %>%
  ppsr::score_predictors(df = ., y = 'DC_50')

trans_dataFilt %>%
  dplyr::select(DC_50, Origin, HA, Subtype, RBS, PBS) %>%
  ppsr::score_predictors(df = ., y = 'DC_50')

###

## DC_67
trans_dataFilt %>%
  dplyr::select(DC_67, AUC_4_v, AUC_6_v, AUC_8_v, AUC_9_v, AUC_RD,
                AUC_4_f, AUC_6_f, AUC_8_f, AUC_9_f, peak_inoc) %>%
  ppsr::score_predictors(df = ., y = 'DC_67')

trans_dataFilt %>%
  dplyr::select(DC_67, slope13_v, slope24_v, slope35_v, slope46_v, slope57_v, slope68_v,
                slope13_f, slope24_f, slope35_f, slope46_f, slope57_f, slope68_f, 
                slope13_vcat, slope13_fcat) %>%
  ppsr::score_predictors(df = ., y = 'DC_67')

trans_dataFilt %>%
  dplyr::select(DC_67, Origin, HA, Subtype, RBS, PBS) %>%
  ppsr::score_predictors(df = ., y = 'DC_67')

########################################################################################

trans_dataFilt2 <- trans_dataFilt

## convert columns to factor, then to 0/1
columns <- c("RD_33", "RD_50", "RD_67", "DC_33", "DC_50", "DC_67")
for (col in columns) {
  trans_dataFilt2[[col]] <- factor(trans_dataFilt2[[col]], levels = c("yes", "no"))
  trans_dataFilt2[[col]] <- as.integer(trans_dataFilt2[[col]] == "no")}

## rotate through
TRANS_COL <- "RD_33"
TRANS_COL <- "RD_50"
TRANS_COL <- "RD_67"
TRANS_COL <- "DC_33"
TRANS_COL <- "DC_50"
TRANS_COL <- "DC_67"

trans_dataFilt3 <- trans_dataFilt2 %>%
  dplyr::select({{TRANS_COL}}, Origin, HPAI, HPAI_MBAA, HA, Subtype,
                slope13_v, slope35_v, AUC_6_v, AUC_8_v,
                d1_inoc, d3_inoc, peak_inoc) %>% drop_na()

trans_dataFilt4 <- trans_dataFilt3 %>%
  dplyr::select({{TRANS_COL}}, Origin, HPAI_MBAA, Subtype,
                slope13_v, AUC_8_v, d1_inoc) %>% drop_na()

trans_dataFilt5 <- trans_dataFilt3 %>%
  dplyr::select({{TRANS_COL}}, Origin, HPAI_MBAA, Subtype,
                slope13_v, d1_inoc) %>% drop_na()

trans_dataFilt6 <- trans_dataFilt3 %>%
  dplyr::select({{TRANS_COL}}, Origin, HPAI_MBAA, Subtype,
                slope13_v) %>% drop_na()

trans_dataFilt7 <- trans_dataFilt3 %>%
  dplyr::select({{TRANS_COL}}, Origin, HPAI_MBAA, HA,
                slope13_v) %>% drop_na()

trans_dataFilt8 <- trans_dataFilt3 %>%
  dplyr::select({{TRANS_COL}}, Origin, HPAI_MBAA, HA,
                slope13_v, AUC_8_v) %>% drop_na()

trans_dataFilt9 <- trans_dataFilt3 %>%
  dplyr::select({{TRANS_COL}}, Origin, HPAI_MBAA, HA,
                slope13_v, AUC_6_v) %>% drop_na()

trans_dataFilt10 <- trans_dataFilt3 %>%
  dplyr::select({{TRANS_COL}}, Origin, HPAI_MBAA, HA,
                slope13_v, peak_inoc) %>% drop_na()

trans_dataFilt11 <- trans_dataFilt2 %>%
  dplyr::select({{TRANS_COL}}, Origin, HPAI_MBAA, HA,
                slope13_vcat) %>% drop_na()

trans_dataFilt12 <- trans_dataFilt2 %>%
  dplyr::select({{TRANS_COL}}, Origin, HPAI_MBAA, HA,
                slope13_vcat, AUC_8_v) %>% drop_na()

trans_dataFilt13 <- trans_dataFilt2 %>%
  dplyr::select({{TRANS_COL}}, Origin, HPAI_MBAA, HA,
                slope13_fcat) %>% drop_na()

trans_dataFilt14 <- trans_dataFilt2 %>%
  dplyr::select({{TRANS_COL}}, Origin, HPAI_MBAA, HA,
                slope13_fcat, AUC_8_v) %>% drop_na()

trans_dataFilt15 <- trans_dataFilt2 %>%
  dplyr::select({{TRANS_COL}}, Origin, HPAI_MBAA, Subtype,
                slope13_v, slope35_v) %>% drop_na()

trans_dataFilt16 <- trans_dataFilt2 %>%
  dplyr::select({{TRANS_COL}}, Origin, HPAI_MBAA, Subtype,
                slope35_v) %>% drop_na()

trans_dataFilt17 <- trans_dataFilt2 %>%
  dplyr::select({{TRANS_COL}}, Origin, HPAI_MBAA, Subtype,
                slope13_v) %>% drop_na()

trans_dataFilt18 <- trans_dataFilt2 %>%
  dplyr::select({{TRANS_COL}}, Origin, HPAI_MBAA, HA,
                slope13_v) %>% drop_na()

trans_dataFilt19 <- trans_dataFilt2 %>%
  dplyr::select({{TRANS_COL}}, Origin, HPAI_MBAA, Subtype,
                slope13_f) %>% drop_na()

trans_dataFilt20 <- trans_dataFilt2 %>%
  dplyr::select({{TRANS_COL}}, Origin, Subtype,
                slope13_v) %>% drop_na()

trans_dataFilt21 <- trans_dataFilt2 %>%
  dplyr::select({{TRANS_COL}}, Origin, HPAI_MBAA, Subtype,
                slope13_vcat) %>% drop_na()

trans_dataFilt22 <- trans_dataFilt2 %>%
  dplyr::select({{TRANS_COL}}, Origin, HPAI_MBAA, HA,
                slope13_vcat) %>% drop_na()

###

## elastic net regression
DF_filt <- trans_dataFilt3
set.seed(9876)
X1 <- subset(DF_filt, select = -DF_filt[[TRANS_COL]])
XX <- model.matrix(~ ., data = X1)
yy <- DF_filt[[TRANS_COL]]

elastic_model <- cv.glmnet(x = XX, y = yy, alpha = 0.35, nfolds = 10, 
                           family = "binomial")
best_lambda <- elastic_model$lambda.min

best_model <- glmnet(XX, y = yy, alpha = 0.35, lambda = best_lambda, 
                     family = "binomial")
coef(best_model)

y_predicted <- predict(elastic_model, s = best_lambda, newx = XX)

## ROC curve
roc_curve <- roc.glmnet(best_model, newx = XX, newy = yy)
plot(roc_curve)

## Calculate RMSE - lower is better
#sqrt(mean((yy - y_predicted)^2)) # manual
Metrics::rmse(yy, y_predicted) # function
## other metrics
#Metrics::mae(yy, y_predicted) # lower is better
Metrics::auc(yy, y_predicted) # higher is better

########################################################################################

## Recursive Feature Elimination
rfeCtrl <- rfeControl(functions = rfFuncs,
                      method = "repeatedcv",
                      number = 10,
                      repeats = 2,
                      saveDetails = TRUE,
                      verbose = FALSE)

rfeProfile <- rfe(x = X1, ## from elastic net model above
                  y = yy, ## from elastic net model above
                  sizes = c(1:16),
                  rfeControl = rfeCtrl)
rfeProfile
predictors(rfeProfile)
## plot the results
plot(rfeProfile, type = c("g", "o"))

########################################################################################

## Data Distribution Checks
variables <- c("slope35_v", "slope13_v", "slope13_f", "slope13_vcat", 
               "slope13_fcat", "d1_inoc", "d3_inoc", "d3_inoc_avg",
               "peak_inoc", "AUC_8_v", "AUC_6_v", "AUC_4_v")

TRANSMISSION <- 'RD_33'
TRANSMISSION <- 'RD_50'
TRANSMISSION <- 'RD_67'
TRANSMISSION <- 'DC_33'
TRANSMISSION <- 'DC_50'
TRANSMISSION <- 'DC_67'

## loop through variables to plot
plots <- list()
for (variable in variables) {
  plot <- trans_dataFilt %>%
    drop_na({{TRANSMISSION}}) %>%
    ggplot(aes(x = .data[[variable]], fill = .data[[TRANSMISSION]])) + 
    geom_density(alpha = 0.5) + 
    scale_fill_viridis_d(direction = -1, end = 0.8) +
    theme_bw() +
    theme(legend.position = 'none')
  plots[[variable]] <- plot}
## combine all looped plots
combined_plot <- wrap_plots(plots, ncol = 3)
combined_plot + theme(legend.position = 'bottom')

## single plots for facet_wrap
trans_dataFilt %>%
  drop_na({{TRANSMISSION}}) %>%
  ggplot(aes(x = slope13_vcat, fill = .data[[TRANSMISSION]])) + 
  geom_density(alpha = 0.5) +
  scale_fill_viridis_d(direction = -1, end = 0.8) +
  theme_bw() +
  facet_wrap(~slope13_vcat)

trans_dataFilt %>%
  drop_na({{TRANSMISSION}}) %>%
  ggplot(aes(x = slope13_fcat, fill = .data[[TRANSMISSION]])) + 
  geom_density(alpha = 0.5) +
  scale_fill_viridis_d(direction = -1, end = 0.8) +
  theme_bw() +
  facet_wrap(~slope13_fcat)

########################################################################################


########################################################################################
########################################################################################
########################################################################################


########################################################################################

## Respiratory Droplet (RD) Transmission Model - Machine Learning

########################################################################################

RD_TRANSMISSION <- 'RD_33'
FORMULA <- RD_33 ~ .
# RD_TRANSMISSION <- 'RD_50'
# FORMULA <- RD_50 ~ .
# RD_TRANSMISSION <- 'RD_67'
# FORMULA <- RD_67 ~ .

fitControl <- trainControl(method = "repeatedcv",   
                           number = 10,  # number of folds
                           repeats = 2,  # repeated two times = 20 folds
                           savePredictions = 'final',
                           classProbs = TRUE,
                           summaryFunction = multiClassSummary)

###

RD_test1 <- trans_dataFilt %>% 
  drop_na(RD_TRANSMISSION) %>%
  dplyr::select(c(RD_TRANSMISSION, Origin, RBS, PBS, slope13_v, HA)) 

RD_test2 <- trans_dataFilt %>% 
  drop_na(RD_TRANSMISSION) %>%
  dplyr::select(c(RD_TRANSMISSION, Origin, RBS, PBS, slope13_v, Subtype)) 

RD_test3 <- trans_dataFilt %>% 
  drop_na(RD_TRANSMISSION) %>%
  dplyr::select(c(RD_TRANSMISSION, Origin, RBS, PBS, slope13_vcat, HA))

RD_test4 <- trans_dataFilt %>% 
  drop_na(RD_TRANSMISSION) %>%
  dplyr::select(c(RD_TRANSMISSION, Origin, RBS, PBS, slope13_fcat, HA))

RD_test5 <- trans_dataFilt %>% 
  drop_na(RD_TRANSMISSION) %>%
  dplyr::select(c(RD_TRANSMISSION, Origin, RBS, PBS, slope13_vcat, HA, AUC_8_v))

RD_test6 <- trans_dataFilt %>% 
  drop_na(RD_TRANSMISSION) %>%
  dplyr::select(c(RD_TRANSMISSION, Origin, RBS, PBS, slope13_vcat, HA, HPAI_MBAA))

###

RD_test1_dummy <- fastDummies::dummy_cols(
  RD_test1, select_columns =
    c('Origin', 'RBS', 'PBS', 'HA'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

RD_test2_dummy <- fastDummies::dummy_cols(
  RD_test2, select_columns =
    c('Origin', 'RBS', 'PBS', 'Subtype'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

RD_test3_dummy <- fastDummies::dummy_cols(
  RD_test3, select_columns =
    c('Origin', 'RBS', 'PBS', 'HA', 'slope13_vcat'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

RD_test4_dummy <- fastDummies::dummy_cols(
  RD_test4, select_columns =
    c('Origin', 'RBS', 'PBS', 'HA', 'slope13_fcat'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

RD_test5_dummy <- fastDummies::dummy_cols(
  RD_test5, select_columns =
    c('Origin', 'RBS', 'PBS', 'HA', 'slope13_vcat'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

RD_test6_dummy <- fastDummies::dummy_cols(
  RD_test6, select_columns =
    c('Origin', 'RBS', 'PBS', 'HA', 'slope13_vcat', 'HPAI_MBAA'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)


###

RD_test1_dummy$RD_33 <- as.factor(RD_test1_dummy$RD_33)
RD_test2_dummy$RD_33 <- as.factor(RD_test2_dummy$RD_33)
RD_test3_dummy$RD_33 <- as.factor(RD_test3_dummy$RD_33)
RD_test4_dummy$RD_33 <- as.factor(RD_test4_dummy$RD_33)
RD_test5_dummy$RD_33 <- as.factor(RD_test5_dummy$RD_33)
RD_test6_dummy$RD_33 <- as.factor(RD_test6_dummy$RD_33)

RD_test1_dummy$RD_50 <- as.factor(RD_test1_dummy$RD_50)
RD_test2_dummy$RD_50 <- as.factor(RD_test2_dummy$RD_50)
RD_test3_dummy$RD_50 <- as.factor(RD_test3_dummy$RD_50)
RD_test4_dummy$RD_50 <- as.factor(RD_test4_dummy$RD_50)
RD_test5_dummy$RD_50 <- as.factor(RD_test5_dummy$RD_50)
RD_test6_dummy$RD_50 <- as.factor(RD_test6_dummy$RD_50)

RD_test1_dummy$RD_67 <- as.factor(RD_test1_dummy$RD_67)
RD_test2_dummy$RD_67 <- as.factor(RD_test2_dummy$RD_67)
RD_test3_dummy$RD_67 <- as.factor(RD_test3_dummy$RD_67)
RD_test4_dummy$RD_67 <- as.factor(RD_test4_dummy$RD_67)
RD_test5_dummy$RD_67 <- as.factor(RD_test5_dummy$RD_67)
RD_test6_dummy$RD_67 <- as.factor(RD_test6_dummy$RD_67)

## split data into train/test using the tidymodels/rsample package
set.seed(9595)
#RD_test_dumSplit <- rsample::initial_split(RD_test1_dummy, prop = 0.70)
#RD_test_dumSplit <- rsample::initial_split(RD_test2_dummy, prop = 0.70)
RD_test_dumSplit <- rsample::initial_split(RD_test3_dummy, prop = 0.70)
#RD_test_dumSplit <- rsample::initial_split(RD_test4_dummy, prop = 0.70)
#RD_test_dumSplit <- rsample::initial_split(RD_test5_dummy, prop = 0.70)
#RD_test_dumSplit <- rsample::initial_split(RD_test6_dummy, prop = 0.70)

trainData_RD <- rsample::training(RD_test_dumSplit)
testData_RD <- rsample::testing(RD_test_dumSplit)

###

model_list <- caretList(
  FORMULA, data = trainData_RD,
  trControl = fitControl,
  na.action = 'na.exclude',
  preProcess = c("nzv", "scale", "center"),
  methodList = c("nnet", "ranger", "rf", "gbm", "svmRadial"))

xyplot(resamples(model_list))
modelCor(resamples(model_list))

## resample metrics
resamps_RD <- resamples(model_list)

summary(resamps_RD)
bwplot(resamps_RD)
dotplot(resamps_RD)
modelCor(resamps_RD)
splom(resamps_RD)

varImp(model_list$nnet)
varImp(model_list$ranger)
varImp(model_list$rf)
varImp(model_list$gbm)
varImp(model_list$svm)

###

## check a range of probability values
thresholds <- seq(0.5, 0.95, by = 0.05)

## empty list to populate
prob_results <- list()
prob_results_df <- list()
## function
## Note: [[RD_TRANSMISSION]] used for list, {{RD_TRANSMISSION}} used for tidyverse
for(threshold in thresholds){
  probTest <- as.data.frame(predict(model_list, newdata = testData_RD))
  probTest <- probTest %>%
    mutate(
      nnet = factor(ifelse(nnet >= threshold, 'no', 'yes')),
      ranger = factor(ifelse(ranger >= threshold, 'no', 'yes')),
      rf = factor(ifelse(rf >= threshold, 'no', 'yes')),
      gbm = factor(ifelse(gbm >= threshold, 'no', 'yes')),
      svmRadial = factor(ifelse(svmRadial >= threshold, 'no', 'yes')))
  
  probTruth <- testData_RD %>%
    na.omit() %>%
    dplyr::select({{RD_TRANSMISSION}}) %>%
    cbind(., probTest)
  
  for(col_name in names(probTest)){
    confusion_matrix <- confusionMatrix(probTest[[col_name]], 
                                        probTruth[[RD_TRANSMISSION]], 
                                        mode = "everything")
    
    sensitivity <- confusion_matrix$byClass['Sensitivity']
    specificity <- confusion_matrix$byClass['Specificity']
    balanced_accuracy <- (sensitivity + specificity) / 2
    precision <- confusion_matrix$byClass['Precision']
    recall <- confusion_matrix$byClass['Recall']
    F1 <- confusion_matrix$byClass['F1']
    
    result_entry <- data.frame(
      threshold = threshold,
      model = col_name,
      sensitivity = sensitivity,
      specificity = specificity,
      balanced_accuracy = balanced_accuracy,
      precision = precision,
      recall = recall,
      F1 = F1)
    
    prob_results_df <- rbind(prob_results_df, result_entry)}
  
  for(col_name in names(probTest)){  # Exclude the first column from the loop
    confusion_matrix <- confusionMatrix(probTest[[col_name]], 
                                        probTruth[[RD_TRANSMISSION]], 
                                        mode = "everything")

    prob_results[[as.character(threshold)]][[col_name]] <- list(
      threshold = threshold,
      confusion_matrix = confusion_matrix)}}

## check outputs
prob_results_df

## Access the results for each threshold, e.g., results[['0.5']]
prob_results[['0.5']]
prob_results[['0.55']]
prob_results[['0.6']]
prob_results[['0.65']]
prob_results[['0.7']]
prob_results[['0.75']]
prob_results[['0.8']]
prob_results[['0.85']]
prob_results[['0.9']]
prob_results[['0.95']]

########################################################################################

## Matthew's Correlation Coefficients

## TP, FP, FN, TN

## RD_33
## RD_test1
mcc_func(40, 3, 0, 41) ## nnet
mcc_func(41, 2, 3, 38) ## ranger
mcc_func(43, 0, 0, 41) ## rf
mcc_func(42, 1, 0, 41) ## gbm
mcc_func(38, 5, 3, 38) ## svm
## RD_test2
mcc_func(37, 6, 1, 40) ## svm
## RD_test3
mcc_func(33, 10, 1, 40) ## gbm
## RD_test4
mcc_func(39, 6, 3, 38) ## nnet
mcc_func(40, 5, 1, 40) ## rf
mcc_func(40, 5, 2, 39) ## gbm
mcc_func(40, 5, 3, 38) ## svm
## RD_test5
mcc_func(37, 3, 0, 41) ## nnet
mcc_func(38, 2, 0, 41) ## rf
mcc_func(34, 6, 1, 40) ## svm
## RD_test6
mcc_func(36, 7, 0, 41) ## nnet
mcc_func(30, 13, 0, 41) ## ranger
mcc_func(32, 11, 0, 41) ## gbm

###

## RD_50
## RD_test1
mcc_func(47, 6, 3, 28)
mcc_func(53, 0, 3, 28)
mcc_func(52, 1, 3, 28)
mcc_func(46, 7, 10, 21)
## RD_test2
mcc_func(50, 3, 3, 28)
mcc_func(49, 4, 3, 28)
mcc_func(53, 0, 3, 28)
mcc_func(52, 1, 3, 28)
mcc_func(46, 7, 7, 24)
## RD_test3
mcc_func(49, 4, 3, 28)
mcc_func(45, 8, 3, 28)
mcc_func(51, 2, 10, 21)
## RD_test4
mcc_func(47, 8, 4, 27)
mcc_func(47, 8, 3, 28)
mcc_func(47, 8, 6, 25)
mcc_func(48, 7, 10, 21)
## RD_test5
mcc_func(48, 2, 3, 28)
mcc_func(50, 0, 3, 28)
mcc_func(46, 4, 3, 28)
mcc_func(48, 2, 6, 25)

###

## RD_67
## RD_test1
mcc_func(56, 4, 3, 21)
mcc_func(60, 0, 3, 21)
mcc_func(59, 1, 3, 21)
mcc_func(54, 6, 9, 15)
## RD_test2
mcc_func(57, 3, 3, 21)
mcc_func(57, 3, 8, 16)
## RD_test3
mcc_func(51, 9, 3, 21)
mcc_func(54, 6, 7, 17)
## RD_test4
mcc_func(52, 10, 5, 19)
mcc_func(51, 11, 4, 20)
mcc_func(51, 11, 7, 17)
mcc_func(55, 7, 9, 15)
## RD_test5
mcc_func(52, 5, 8, 16)
mcc_func(57, 0, 3, 21)
mcc_func(56, 1, 3, 21)
mcc_func(51, 6, 7, 17)

########################################################################################

## hyperparameter tuning
modelLookup(model = 'rf') ## mtry
RD_test_grid1 <- expand.grid(mtry = seq(from = 2, to = 20, by = 1))

RD_test1_final <- expand.grid(mtry = 13)
RD_test3_final <- expand.grid(mtry = 2)

set.seed(9595)
RD_test1_tune <- train(RD_33 ~ ., data = trainData_RD,
                  method = "rf",
                  na.action = na.exclude,
                  tuneGrid = RD_test1_final,
                  metric = "Balanced_Accuracy",
                  preProcess = c("nzv", "scale", "center"),
                  trControl = fitControl)
RD_test1_tune

# saveRDS(RD_test1_tune, "./final_models/RD_test1_tune_slope13_num_final_model.rds")
RD_test1_tune <- read_rds("./final_models/RD_test1_tune_slope13_num_final_model.rds")

varImp(RD_test1_tune)

set.seed(9595)
RD_test3_tune <- train(RD_33 ~ ., data = trainData_RD,
                       method = "rf",
                       na.action = na.exclude,
                       tuneGrid = RD_test3_final,
                       metric = "Balanced_Accuracy",
                       preProcess = c("nzv", "scale", "center"),
                       trControl = fitControl)
RD_test3_tune

# saveRDS(RD_test3_tune, "./final_models/RD_test3_tune_slope13_vcat_final_model.rds")
RD_test3_tune <- read_rds("./final_models/RD_test3_tune_slope13_vcat_final_model.rds")

varImp(RD_test3_tune)

########################################################################################

### external validation data from standardization exercise

## set up train/test data
RD_33_standardized <- standardized %>%
  dplyr::select(c(RD_33, Origin, HA, RBS, PBS, slope13_v)) %>%
  drop_na()

RD_33_standardized <- standardized %>%
  dplyr::select(c(RD_33, Origin, HA, RBS, PBS, slope13_vcat)) %>%
  drop_na()

## one hot code molecular data
RD_33_standardized_dummy <- fastDummies::dummy_cols(
  RD_33_standardized, select_columns =
    c('Origin', 'HA', 'RBS', 'PBS'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

RD_33_standardized_dummy <- fastDummies::dummy_cols(
  RD_33_standardized, select_columns =
    c('Origin', 'HA', 'RBS', 'PBS', 'slope13_vcat'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

RD_33_standardized_dummy$RD_33 <- as.factor(RD_33_standardized_dummy$RD_33)

# set.seed(9595)
# RD_standardized_test_dumSplit <- rsample::initial_split(RD_33_standardized_dummy, prop = 0.02)
# 
# trainData_standardized_RD <- rsample::training(RD_standardized_test_dumSplit)
# testData_standardized_RD <- rsample::testing(RD_standardized_test_dumSplit)
# 
# rbind(trainData_standardized_RD, testData_standardized_RD)
# ## add missing dummy var columns
# testData_standardized_RD <- testData_standardized_RD %>% 
#   mutate(HA_H2 = 0, HA_H3 = 0, HA_H5 = 0, HA_H7 = 0, HA_H9 = 0)

## use data without splitting
testData_standardized_RD <- RD_33_standardized_dummy %>% 
  mutate(HA_H2 = 0, HA_H3 = 0, HA_H5 = 0, HA_H7 = 0, HA_H9 = 0,
         RBS_A = 0)

testData_standardized_RD <- RD_33_standardized_dummy %>% 
  mutate(HA_H2 = 0, HA_H3 = 0, HA_H5 = 0, HA_H7 = 0, HA_H9 = 0,
         RBS_A = 0, slope13_vcat_pos = 0)

###

### external validation data from published literature

## set up train/test data
RD_33_pubLit <- pubLit %>%
  dplyr::select(c(RD_33, Origin, HA, RBS, PBS, slope13_vcat)) %>%
  drop_na()

## one hot code molecular data
RD_33_pubLit_dummy <- fastDummies::dummy_cols(
  RD_33_pubLit, select_columns =
    c('Origin', 'HA', 'RBS', 'PBS', 'slope13_vcat'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

RD_33_pubLit_dummy$RD_33 <- as.factor(RD_33_pubLit_dummy$RD_33)

# set.seed(9595)
# RD_pubLit_test_dumSplit <- rsample::initial_split(RD_33_pubLit_dummy, prop = 0.01)
# 
# trainData_pubLit_RD <- rsample::training(RD_pubLit_test_dumSplit)
# testData_pubLit_RD <- rsample::testing(RD_pubLit_test_dumSplit)
# ## add missing dummy var column
# testData_pubLit_RD <- testData_pubLit_RD %>% mutate(HA_H2 = 0)

## use data without splitting
testData_pubLit_RD <- RD_33_pubLit_dummy %>% mutate(HA_H2 = 0)

###

## run the model and evaluate
thresholds <- seq(0.5, 0.95, by = 0.05)
tuned_model <- RD_test1_tune ## testData_RD & testData_standardized_RD
#tuned_model <- RD_test3_tune ## testData_RD & testData_pubLit_RD

#TESTDATA <- testData_standardized_RD ## use validation (standardization) data w/ tuned model 1
#TESTDATA <- testData_pubLit_RD ## use validation (pubLit) data w/ tuned model 3
TESTDATA <- testData_RD ## use test data w/ both tuned model

prob_results <- list()

for (threshold in thresholds) {
  probTest <- predict(tuned_model, TESTDATA, type = 'prob')
  probTest <- factor(ifelse(probTest$no >= threshold, 'no', 'yes')) %>%
    as.data.frame()
  probTruth <- TESTDATA %>%
    na.omit() %>%
    dplyr::select(RD_33) %>%
    cbind(., probTest) %>%
    drop_na(RD_33)
  confusionMatrixResult <- confusionMatrix(probTruth$., probTruth$RD_33, mode = "everything")
  prob_results[[as.character(threshold)]] <- confusionMatrixResult}

prob_results[['0.5']]
prob_results[['0.55']]
prob_results[['0.6']]
prob_results[['0.65']]
prob_results[['0.7']]
prob_results[['0.75']]
prob_results[['0.8']]
prob_results[['0.85']]
prob_results[['0.9']]
prob_results[['0.95']]

mcc_func(43, 0, 0, 41) ## tune1 - testData
mcc_func(8, 28, 4, 32) ## tune1 - standardized

mcc_func(37, 6, 1, 40) ## tune3 - testData
mcc_func(109, 16, 3, 40) ## tune3 - pubLit
mcc_func(32, 4, 4, 32) ## tune3 - standardized

varImp(prob_results)

########################################################################################


########################################################################################
########################################################################################
########################################################################################


########################################################################################

## Direct Contact (DC) Transmission Model - Machine Learning

########################################################################################

# DC_TRANSMISSION <- 'DC_33'
# FORMULA_DC <- DC_33 ~ .
# DC_TRANSMISSION <- 'DC_50'
# FORMULA_DC <- DC_50 ~ .
DC_TRANSMISSION <- 'DC_67'
FORMULA_DC <- DC_67 ~ .

fitControl <- trainControl(method = "repeatedcv",   
                           number = 10,  # number of folds
                           repeats = 2,  # repeated two times = 20 folds
                           savePredictions = 'final',
                           classProbs = TRUE,
                           summaryFunction = multiClassSummary)

###

DC_test1 <- trans_dataFilt %>% 
  drop_na(DC_TRANSMISSION) %>%
  dplyr::select(c(DC_TRANSMISSION, Origin, Subtype, RBS, PBS)) 

DC_test2 <- trans_dataFilt %>% 
  drop_na(DC_TRANSMISSION) %>%
  dplyr::select(c(DC_TRANSMISSION, Origin, HA, RBS, PBS)) 

DC_test3 <- trans_dataFilt %>% 
  drop_na(DC_TRANSMISSION) %>%
  dplyr::select(c(DC_TRANSMISSION, Origin, HA, RBS, PBS, slope13_v))

DC_test4 <- trans_dataFilt %>% 
  drop_na(DC_TRANSMISSION) %>%
  dplyr::select(c(DC_TRANSMISSION, Origin, HA, RBS, PBS, peak_inoc))

DC_test5 <- trans_dataFilt %>% 
  drop_na(DC_TRANSMISSION) %>%
  dplyr::select(c(DC_TRANSMISSION, Origin, HA, RBS, PBS, AUC_8_v))

DC_test6 <- trans_dataFilt %>% 
  drop_na(DC_TRANSMISSION) %>%
  dplyr::select(c(DC_TRANSMISSION, Origin, HA, RBS, PBS, slope13_v, HPAI_MBAA))

DC_test7 <- trans_dataFilt %>% 
  drop_na(DC_TRANSMISSION) %>%
  dplyr::select(c(DC_TRANSMISSION, Origin, HA, RBS, PBS, HPAI_MBAA))

DC_test8 <- trans_dataFilt %>% 
  drop_na(DC_TRANSMISSION) %>%
  dplyr::select(c(DC_TRANSMISSION, Origin, HA, RBS, PBS, slope13_vcat, HPAI_MBAA))

DC_test9 <- trans_dataFilt %>% 
  drop_na(DC_TRANSMISSION) %>%
  dplyr::select(c(DC_TRANSMISSION, Origin, HA, RBS, PBS, slope13_vcat))

DC_test1_dummy <- fastDummies::dummy_cols(
  DC_test1, select_columns =
    c('Origin', 'Subtype', 'RBS', 'PBS'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

DC_test2_dummy <- fastDummies::dummy_cols(
  DC_test2, select_columns =
    c('Origin', 'HA', 'RBS', 'PBS'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

DC_test3_dummy <- fastDummies::dummy_cols(
  DC_test3, select_columns =
    c('Origin', 'HA', 'RBS', 'PBS'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

DC_test4_dummy <- fastDummies::dummy_cols(
  DC_test4, select_columns =
    c('Origin', 'HA', 'RBS', 'PBS'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

DC_test5_dummy <- fastDummies::dummy_cols(
  DC_test5, select_columns =
    c('Origin', 'HA', 'RBS', 'PBS'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

DC_test6_dummy <- fastDummies::dummy_cols(
  DC_test6, select_columns =
    c('Origin', 'HA', 'RBS', 'PBS', 'HPAI_MBAA'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

DC_test7_dummy <- fastDummies::dummy_cols(
  DC_test7, select_columns =
    c('Origin', 'HA', 'RBS', 'PBS', 'HPAI_MBAA'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

DC_test8_dummy <- fastDummies::dummy_cols(
  DC_test8, select_columns =
    c('Origin', 'HA', 'RBS', 'PBS', 'HPAI_MBAA', 'slope13_vcat'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

DC_test9_dummy <- fastDummies::dummy_cols(
  DC_test9, select_columns =
    c('Origin', 'HA', 'RBS', 'PBS', 'slope13_vcat'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

DC_test1_dummy$DC_33 <- as.factor(DC_test1_dummy$DC_33)
DC_test2_dummy$DC_33 <- as.factor(DC_test2_dummy$DC_33)
DC_test3_dummy$DC_33 <- as.factor(DC_test3_dummy$DC_33)
DC_test4_dummy$DC_33 <- as.factor(DC_test4_dummy$DC_33)
DC_test5_dummy$DC_33 <- as.factor(DC_test5_dummy$DC_33)
DC_test6_dummy$DC_33 <- as.factor(DC_test6_dummy$DC_33)
DC_test7_dummy$DC_33 <- as.factor(DC_test7_dummy$DC_33)
DC_test8_dummy$DC_33 <- as.factor(DC_test8_dummy$DC_33)
DC_test9_dummy$DC_33 <- as.factor(DC_test9_dummy$DC_33)

DC_test1_dummy$DC_50 <- as.factor(DC_test1_dummy$DC_50)
DC_test2_dummy$DC_50 <- as.factor(DC_test2_dummy$DC_50)
DC_test3_dummy$DC_50 <- as.factor(DC_test3_dummy$DC_50)
DC_test4_dummy$DC_50 <- as.factor(DC_test4_dummy$DC_50)
DC_test5_dummy$DC_50 <- as.factor(DC_test5_dummy$DC_50)
DC_test6_dummy$DC_50 <- as.factor(DC_test6_dummy$DC_50)
DC_test7_dummy$DC_50 <- as.factor(DC_test7_dummy$DC_50)
DC_test8_dummy$DC_50 <- as.factor(DC_test8_dummy$DC_50)
DC_test9_dummy$DC_50 <- as.factor(DC_test9_dummy$DC_50)

DC_test1_dummy$DC_67 <- as.factor(DC_test1_dummy$DC_67)
DC_test2_dummy$DC_67 <- as.factor(DC_test2_dummy$DC_67)
DC_test3_dummy$DC_67 <- as.factor(DC_test3_dummy$DC_67)
DC_test4_dummy$DC_67 <- as.factor(DC_test4_dummy$DC_67)
DC_test5_dummy$DC_67 <- as.factor(DC_test5_dummy$DC_67)
DC_test6_dummy$DC_67 <- as.factor(DC_test6_dummy$DC_67)
DC_test7_dummy$DC_67 <- as.factor(DC_test7_dummy$DC_67)
DC_test8_dummy$DC_67 <- as.factor(DC_test8_dummy$DC_67)
DC_test9_dummy$DC_67 <- as.factor(DC_test9_dummy$DC_67)

###

set.seed(9595)
#DC_test_dumSplit <- rsample::initial_split(DC_test1_dummy, prop = 0.70)
#DC_test_dumSplit <- rsample::initial_split(DC_test2_dummy, prop = 0.70)
#DC_test_dumSplit <- rsample::initial_split(DC_test3_dummy, prop = 0.70)
#DC_test_dumSplit <- rsample::initial_split(DC_test4_dummy, prop = 0.70)
#DC_test_dumSplit <- rsample::initial_split(DC_test5_dummy, prop = 0.70)
#DC_test_dumSplit <- rsample::initial_split(DC_test6_dummy, prop = 0.70)
#DC_test_dumSplit <- rsample::initial_split(DC_test7_dummy, prop = 0.70)
#DC_test_dumSplit <- rsample::initial_split(DC_test8_dummy, prop = 0.70)
DC_test_dumSplit <- rsample::initial_split(DC_test9_dummy, prop = 0.70)

trainData_DC <- rsample::training(DC_test_dumSplit)
testData_DC <- rsample::testing(DC_test_dumSplit)

###

model_listDC <- caretList(
  FORMULA_DC, data = trainData_DC,
  trControl = fitControl,
  na.action = 'na.exclude',
  preProcess = c("nzv", "scale", "center"),
  methodList = c("nnet", "ranger", "rf", "gbm", "svmRadial"))

xyplot(resamples(model_listDC))
modelCor(resamples(model_listDC))

## resample metrics
resamps_DC <- resamples(model_listDC)

summary(resamps_DC)
bwplot(resamps_DC)
dotplot(resamps_DC)
modelCor(resamps_DC)
splom(resamps_DC)

#library(gbm)
varImp(model_listDC$nnet)
varImp(model_listDC$ranger)
varImp(model_listDC$rf)
varImp(model_listDC$gbm)
varImp(model_listDC$svm)

###

## check a range of probability values
thresholds <- seq(0.5, 0.95, by = 0.05)

## empty list to populate
prob_results <- list()
prob_results_df <- list()
## function
## Note: [[DC_TRANSMISSION]] used for list, {{DC_TRANSMISSION}} used for tidyverse
for(threshold in thresholds){
  probTest <- as.data.frame(predict(model_listDC, newdata = testData_DC))
  probTest <- probTest %>%
    mutate(
      nnet = factor(ifelse(nnet >= threshold, 'no', 'yes')),
      ranger = factor(ifelse(ranger >= threshold, 'no', 'yes')),
      rf = factor(ifelse(rf >= threshold, 'no', 'yes')),
      gbm = factor(ifelse(gbm >= threshold, 'no', 'yes')),
      svmRadial = factor(ifelse(svmRadial >= threshold, 'no', 'yes')))
  
  probTruth <- testData_DC %>%
    na.omit() %>%
    dplyr::select({{DC_TRANSMISSION}}) %>%
    cbind(., probTest)
  
  for(col_name in names(probTest)){
    confusion_matrix <- confusionMatrix(probTest[[col_name]], 
                                        probTruth[[DC_TRANSMISSION]], 
                                        mode = "everything")
    
    sensitivity <- confusion_matrix$byClass['Sensitivity']
    specificity <- confusion_matrix$byClass['Specificity']
    balanced_accuracy <- (sensitivity + specificity) / 2
    precision <- confusion_matrix$byClass['Precision']
    recall <- confusion_matrix$byClass['Recall']
    F1 <- confusion_matrix$byClass['F1']
    
    result_entry <- data.frame(
      threshold = threshold,
      model = col_name,
      sensitivity = sensitivity,
      specificity = specificity,
      balanced_accuracy = balanced_accuracy,
      precision = precision,
      recall = recall,
      F1 = F1)
    
    prob_results_df <- rbind(prob_results_df, result_entry)}
  
  for(col_name in names(probTest)){  # Exclude the first column from the loop
    confusion_matrix <- confusionMatrix(probTest[[col_name]], 
                                        probTruth[[DC_TRANSMISSION]], 
                                        mode = "everything")
    
    prob_results[[as.character(threshold)]][[col_name]] <- list(
      threshold = threshold,
      confusion_matrix = confusion_matrix)}}

## check outputs
prob_results_df

## Access the results for each threshold, e.g., results[['0.5']]
prob_results[['0.5']]
prob_results[['0.55']]
prob_results[['0.6']]
prob_results[['0.65']]
prob_results[['0.7']]
prob_results[['0.75']]
prob_results[['0.8']]
prob_results[['0.85']]
prob_results[['0.9']]
prob_results[['0.95']]

########################################################################################

## Matthew's Correlation Coefficients

## TP, FP, FN, TN

## DC_33
## DC_test1
mcc_func(9, 11, 0, 71)
mcc_func(0, 20, 0, 71)
mcc_func(8, 12, 0, 71)
## DC_test2
mcc_func(12, 8, 2, 69)
mcc_func(11, 9, 0, 71)
mcc_func(13, 7, 2, 69)
## DC_test3
mcc_func(17, 2, 1, 69)
mcc_func(18, 1, 0, 70)
mcc_func(13, 6, 1, 69)
## DC_test4
mcc_func(14, 6, 0, 71)
mcc_func(12, 8, 1, 70)
mcc_func(13, 7, 2, 69)
mcc_func(11, 9, 2, 69)
## DC_test5
mcc_func(11, 7, 0, 70)
mcc_func(15, 3, 0, 70)
mcc_func(12, 6, 0, 70)
mcc_func(9, 9, 0, 70)
## DC_test6
mcc_func(15, 4, 0, 70)
mcc_func(12, 7, 1, 69)
## DC_test6
mcc_func(14, 5, 1, 69)
mcc_func(13, 6, 1, 69)
## DC_test9
mcc_func(12, 7, 1, 69)

## DC_50
## DC_test1
mcc_func(26, 4, 4, 57)
mcc_func(28, 2, 4, 57)
mcc_func(25, 5, 4, 57)
## DC_test2
mcc_func(27, 3, 4, 57)
mcc_func(20, 10, 1, 60)
mcc_func(22, 8, 2, 59)
## DC_test3
mcc_func(29, 0, 1, 59)
mcc_func(28, 1, 0, 60)
mcc_func(28, 1, 2, 58)
## DC_test4
mcc_func(24, 6, 3, 58)
mcc_func(23, 7, 3, 58)
mcc_func(21, 9, 3, 58)
mcc_func(19, 11, 1, 60)
## DC_test5
mcc_func(20, 7, 2, 59)
mcc_func(24, 3, 0, 61)
mcc_func(23, 4, 0, 61)
mcc_func(20, 7, 0, 61)
## DC_test6
mcc_func(24, 5, 0, 60)
mcc_func(22, 4, 0, 60)
mcc_func(27, 2, 0, 60)
mcc_func(23, 6, 0, 60)
## DC_test7
mcc_func(25, 5, 1, 60)
mcc_func(25, 5, 3, 58)
## DC_test8
mcc_func(25, 4, 2, 58)
## DC_test9
mcc_func(25, 4, 0, 60)
mcc_func(20, 9, 0, 60)
mcc_func(21, 8, 0, 60)

## DC_67
## DC_test1
mcc_func(35, 0, 4, 52)
mcc_func(33, 2, 4, 52)
## DC_test2
mcc_func(35, 0, 2, 54)
mcc_func(33, 2, 2, 54)
## DC_test3
mcc_func(32, 1, 0, 56)
mcc_func(28, 5, 0, 56)
## DC_test4
mcc_func(33, 2, 1, 55)
mcc_func(35, 0, 3, 53)
mcc_func(30, 5, 1, 55)
mcc_func(28, 7, 1, 55)
mcc_func(32, 3, 2, 54)
## DC_test5
mcc_func(26, 6, 0, 56)
mcc_func(32, 0, 2, 54)
## DC_test6
mcc_func(31, 2, 0, 56)
mcc_func(30, 3, 0, 56)
## DC_test8
mcc_func(33, 0, 4, 52)
mcc_func(32, 1, 2, 54)
mcc_func(30, 3, 0, 56)
mcc_func(33, 0, 2, 54)
## DC_test9
mcc_func(33, 0, 2, 54)
mcc_func(25, 8, 0, 56)

########################################################################################

## hyperparameter tuning
modelLookup(model = 'rf') ## mtry
DC_test_grid1 <- expand.grid(mtry = seq(from = 2, to = 20, by = 1))

DC_test3_final <- expand.grid(mtry = 10)
DC_test9_final <- expand.grid(mtry = 2)

set.seed(9595)
DC_test3_tune <- train(DC_67 ~ ., data = trainData_DC,
                       method = "rf",
                       na.action = na.exclude,
                       tuneGrid = DC_test3_final,
                       metric = "Balanced_Accuracy",
                       preProcess = c("nzv", "scale", "center"),
                       trControl = fitControl)
DC_test3_tune

# saveRDS(DC_test3_tune, "./final_models/DC_test3_tune_slope13_num_final_model.rds")
DC_test3_tune <- read_rds("./final_models/DC_test3_tune_slope13_num_final_model.rds")

varImp(DC_test3_tune)

set.seed(9595)
DC_test9_tune <- train(DC_67 ~ ., data = trainData_DC,
                       method = "rf",
                       na.action = na.exclude,
                       tuneGrid = DC_test9_final,
                       metric = "Balanced_Accuracy",
                       preProcess = c("nzv", "scale", "center"),
                       trControl = fitControl)
DC_test9_tune

# saveRDS(DC_test9_tune, "./final_models/DC_test9_tune_slope13_vcat_final_model.rds")
DC_test9_tune <- read_rds("./final_models/DC_test9_tune_slope13_vcat_final_model.rds")

varImp(DC_test9_tune)  

########################################################################################

### external validation data from published literature

## set up train/test data
DC_67_pubLit <- pubLit %>%
  dplyr::select(c(DC_67, Origin, HA, RBS, PBS, slope13_vcat)) %>%
  drop_na()

## one hot code molecular data
DC_67_pubLit_dummy <- fastDummies::dummy_cols(
  DC_67_pubLit, select_columns =
    c('Origin', 'HA', 'RBS', 'PBS', 'slope13_vcat'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

DC_67_pubLit_dummy$DC_67 <- as.factor(DC_67_pubLit_dummy$DC_67)

## use data without splitting
testData_pubLit_DC <- DC_67_pubLit_dummy %>% mutate(HA_H2 = 0)

###

## run the model and evaluate
thresholds <- seq(0.5, 0.95, by = 0.05)
tuned_model <- DC_test3_tune ## testData_DC
#tuned_model <- DC_test9_tune ## testData_DC & testData_pubLit_DC

TESTDATA <- testData_DC ## use test data w/ both tuned model
#TESTDATA <- testData_pubLit_DC ## use validation (pubLit) data w/ tuned model 9

prob_results <- list()

for (threshold in thresholds) {
  probTest <- predict(tuned_model, TESTDATA, type = 'prob')
  probTest <- factor(ifelse(probTest$no >= threshold, 'no', 'yes')) %>%
    as.data.frame()
  probTruth <- TESTDATA %>%
    na.omit() %>%
    dplyr::select(DC_67) %>%
    cbind(., probTest) %>%
    drop_na(DC_67)
  confusionMatrixResult <- confusionMatrix(probTruth$., probTruth$DC_67, mode = "everything")
  prob_results[[as.character(threshold)]] <- confusionMatrixResult}

prob_results[['0.5']]
prob_results[['0.55']]
prob_results[['0.6']]
prob_results[['0.65']]
prob_results[['0.7']]
prob_results[['0.75']]
prob_results[['0.8']]
prob_results[['0.85']]
prob_results[['0.9']]
prob_results[['0.95']]

mcc_func(33, 0, 0, 56) ## tune3 - testData

mcc_func(33, 0, 2, 54) ## tune9 - testData
mcc_func(47, 3, 9, 45) ## tune9 - pubLit


########################################################################################
########################################################################################

### Section below combines data from RD & DC models for plot making

########################################################################################
########################################################################################

### import summary data

features <- read.csv("ResultsSummary_Rinputs/Features4Heatmap.csv", header = TRUE)
modelMetrics <- read.csv("ResultsSummary_Rinputs/ModelMetrics4Heatmap.csv", header = TRUE)
sample_size_data <- read.csv("ResultsSummary_Rinputs/SampleSizeData.csv", header = TRUE)
SelectedImportance_RD <- read.csv("ResultsSummary_Rinputs/SelectedImportance_RD.csv", header = TRUE)
SelectedImportance_DC <- read.csv("ResultsSummary_Rinputs/SelectedImportance_DC.csv", header = TRUE)
pps_trans <- read.csv("ResultsSummary_Rinputs/PPS_transmission.csv", header = TRUE)

########################################################################################

### Sample Sizes, Model Metrics, and Features Heatmaps Figure

## Select relevant columns for creating the presence-absence matrix
presence_matrix <- features[, 5:10] != ""
## Convert to binary matrix format (True for presence, False for absence)
binary_matrix <- as.matrix(presence_matrix)
## drop text Feature columns
features <- features[, 1:4]
## combine binary to features
features_binary <- cbind(features, binary_matrix)   
## convert TRUE/FALSE to 1/0
features_binary <- features_binary %>%
  mutate(across(where(is.logical), ~ if_else(.x, 1, 0)))

gathered_data <- features_binary %>%
  gather(key = "Feature", value = "Presence", 5:10)

## should separate out the Types and then patchwork together rather than facet.
gathered_data$Test <- as.factor(gathered_data$Test)

gathered_data$Type <- factor(gathered_data$Type, 
                             levels=c('test', 'stan', 'pub'))

gathered_data$Feature <- factor(gathered_data$Feature,
                                levels=c('Origin', 'HA', 'RBS', 'PBS',
                                         'slope13_v', 'slope13_vcat'))
## subset data
RD_data <- gathered_data %>% filter(Classification == 'RDT')
DC_data <- gathered_data %>% filter(Classification == 'DCT')

## make subplots
RD_plot <- 
  ggplot(data = RD_data) +
  geom_tile(aes(x = Feature, y = Test, fill = Presence)) +
  scale_fill_viridis_c(direction = -1, end = 0.7) +
  theme_minimal() +
  facet_grid(~ Classification, scales = 'free', space = 'free') +
  ggtitle('Features') +
  theme(legend.position = 'none', axis.title.y = element_blank(),
        strip.text.x = element_blank(), axis.text.y = element_blank(),
        plot.title = element_text(hjust = 0.5), axis.title.x = element_blank()) +
  scale_x_discrete(position = "top")

DC_plot <- 
  ggplot(data = DC_data) +
  geom_tile(aes(x = Feature, y = Test, fill = Presence)) +
  scale_fill_viridis_c(direction = -1, end = 0.7) +
  theme_minimal() +
  facet_grid(~ Classification, scales = 'free', space = 'free') +
  theme(legend.position = 'none', axis.title.y = element_blank(),
        strip.text.x = element_blank(), axis.text.y = element_blank(),
        axis.title.x = element_blank()) +
  scale_x_discrete(position = "top")

###

## final models metric heatmap
modelMetrics <- modelMetrics %>% 
  gather(RDT.num.test:DCT.cat.pub, key = Model, value = Value, factor_key = TRUE) %>%
  dplyr::filter(Metric != 'test') %>%
  drop_na()

modelMetrics$Metric <- factor(modelMetrics$Metric, 
                              levels=c('BA', 'F1', 'MCC'))

modelMetrics_RD <- modelMetrics %>% filter(str_detect(Model, "^RDT"))
modelMetrics_DC <- modelMetrics %>% filter(str_detect(Model, "^DCT"))

modelMetrics_RD$Model <- factor(modelMetrics_RD$Model, 
                                levels=c('RDT.num.test', 'RDT.num.stan', 'RDT.cat.test', 
                                         'RDT.cat.stan', 'RDT.cat.pub'))

modelMetrics_DC$Model <- factor(modelMetrics_DC$Model, 
                                levels=c('DCT.num.test', 'DCT.cat.test', 'DCT.cat.pub'))

modelMetrics_RD$Value <- as.numeric(modelMetrics_RD$Value)
modelMetrics_DC$Value <- as.numeric(modelMetrics_DC$Value)

viridis_scale <- scale_fill_viridis_c(limits = c(0.1, 1), direction = -1, end = 0.9)

modelMetrics_RDplot <- 
  ggplot(modelMetrics_RD, aes(Metric, Model, fill = Value)) + 
  geom_tile() +
  geom_text(aes(label = round(Value, digits = 3)), color = "white", size = 3) +
  viridis_scale + 
  theme_minimal() +
  ggtitle('Metrics') +
  theme(legend.position = 'none', plot.title = element_text(hjust = 0.5), 
        axis.title.x = element_blank(), axis.title.y = element_blank()) +
  scale_x_discrete(position = "top")

modelMetrics_DCplot <- 
  ggplot(modelMetrics_DC, aes(Metric, Model, fill = Value)) + 
  geom_tile() +
  geom_text(aes(label = round(Value, digits = 3)), color = "white", size = 3) +
  viridis_scale + 
  #scale_fill_viridis_d(direction = -1, begin = 0.1, end = 0.85) +
  theme_minimal() +
  theme(legend.position = 'none', axis.title.x = element_blank(), 
        axis.title.y = element_blank()) +
  scale_x_discrete(position = "top")

###

sample_size_data_RD <- sample_size_data %>% filter(str_detect(Model, "^RDT"))
sample_size_data_DC <- sample_size_data %>% filter(str_detect(Model, "^DCT"))

sample_size_data_RD$Model <- factor(sample_size_data_RD$Model, 
                                    levels=c('RDT-num-test', 'RDT-num-stan', 'RDT-cat-test', 
                                             'RDT-cat-stan', 'RDT-cat-pub'))

sample_size_data_DC$Model <- factor(sample_size_data_DC$Model, 
                                    levels=c('DCT-num-test', 'DCT-cat-test', 'DCT-cat-pub'))

sample_size_RDplot <- ggplot(sample_size_data_RD, aes(x = x, y = Model)) +
  geom_tile(fill = 'dodgerblue4') +
  geom_text(aes(label = Virus), color = 'white', size = 3) +
  theme_minimal() +
  ggtitle('Virus') +
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        plot.title = element_text(size = 9, hjust = 0.5))

sample_size_RDplot2 <- ggplot(sample_size_data_RD, aes(x = x, y = Model)) +
  geom_tile(fill = 'dodgerblue3') +
  geom_text(aes(label = Obs_yes), color = 'white', size = 2.8) +
  theme_minimal() +
  ggtitle('Obs(yes)') +
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        plot.title = element_text(size = 9, hjust = 0.5))

sample_size_DCplot <- ggplot(sample_size_data_DC, aes(x = x, y = Model)) +
  geom_tile(fill = 'dodgerblue4') +
  geom_text(aes(label = Virus), color = 'white', size = 3) +
  theme_minimal() +
  ggtitle('Virus') +
  theme(axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size = 9, hjust = 0.5))

sample_size_DCplot2 <- ggplot(sample_size_data_DC, aes(x = x, y = Model)) +
  geom_tile(fill = 'dodgerblue3') +
  geom_text(aes(label = Obs_yes), color = 'white', size = 2.8) +
  theme_minimal() +
  ggtitle('Obs(yes)') +
  theme(axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size = 9, hjust = 0.5))

layout <- ((sample_size_RDplot | sample_size_RDplot2) / 
             (sample_size_DCplot | sample_size_DCplot2) |
             (modelMetrics_RDplot / modelMetrics_DCplot) | 
             (RD_plot / DC_plot))

layout + patchwork::plot_layout(widths = c(0.55, 0.8, 1.4)) 


########################################################################################

### Feature Importance Figure

## barplots with selected importance variables

## RD
SelectedImportance_RD_num <- SelectedImportance_RD %>%
  filter(Model == 'RDT-num-test') #%>%
#mutate(Feature = fct_reorder(factor(Feature), Importance, .desc = FALSE))

SelectedImportance_RD_cat <- SelectedImportance_RD %>%
  filter(Model != 'RDT-num-test')

SelectedImportance_RD_num_plot <- 
  ggplot(SelectedImportance_RD_num, 
         aes(x = Feature, y = Importance, fill = Feature)) +
  geom_col() +
  facet_wrap(~ Model, scales = "free_y") +
  coord_flip() +
  ylab("Relative Ranked Importance") +
  scale_fill_viridis_d(option = 'B', begin = 0.1, end = 0.8) +
  theme_minimal() +
  theme(legend.position = 'none')

SelectedImportance_RD_cat_plot <- 
  ggplot(SelectedImportance_RD_cat, 
         aes(x = Feature, y = Importance, fill = Feature)) +
  geom_col() +
  facet_wrap(~ Model) +
  coord_flip() +
  ylab("Relative Ranked Importance") +
  scale_fill_viridis_d(option = 'G', end = 0.8) +
  theme_minimal() +
  theme(legend.position = 'none') #+
  #scale_y_continuous(expand = expansion(add = c(0, 25)))

###

## DC
SelectedImportance_DC_num <- SelectedImportance_DC %>%
  filter(Model == 'DCT-num-test') #%>%
#mutate(Feature = fct_reorder(factor(Feature), Importance, .desc = FALSE))

SelectedImportance_DC_cat <- SelectedImportance_DC %>%
  filter(Model != 'DCT-num-test')

SelectedImportance_DC_num_plot <- 
  ggplot(SelectedImportance_DC_num, 
         aes(x = Feature, y = Importance, fill = Feature)) +
  geom_col() +
  facet_wrap(~ Model, scales = "free_y") +
  coord_flip() +
  ylab("Relative Ranked Importance") +
  scale_fill_viridis_d(option = 'B', begin = 0.1, end = 0.8) +
  theme_minimal() +
  theme(legend.position = 'none')

SelectedImportance_DC_cat_plot <- 
  ggplot(SelectedImportance_DC_cat, 
         aes(x = Feature, y = Importance, fill = Feature)) +
  geom_col() +
  facet_wrap(~ Model) +
  coord_flip() +
  ylab("Relative Ranked Importance") +
  scale_fill_viridis_d(option = 'G', end = 0.8) +
  theme_minimal() +
  theme(legend.position = 'none')

## combine lethality & morbidity plots
(SelectedImportance_RD_num_plot | SelectedImportance_RD_cat_plot) / 
  (SelectedImportance_DC_num_plot | SelectedImportance_DC_cat_plot) + 
  patchwork::plot_annotation(tag_levels = 'A')

########################################################################################

## Predictive Power Scores (PPS) & MCC RD v DC for each threshold

pps_trans$Transmission <- factor(pps_trans$Transmission, 
                                 levels = c("RDT", "DCT"))
## remove categorical values
pps_trans <- pps_trans %>% filter(Feature != 'slope13_vcat' & 
                                    Feature != 'MCC_cat')

pps_trans$Feature <- factor(pps_trans$Feature, 
                                 levels = c('MCC', 'Origin', 'PBS',
                                            'RBS', 'slope13', 'AUC_8_v', 'peak_inoc'))

ggplot(pps_trans, aes(factor(Threshold), Feature, fill = PPS)) + 
  geom_tile() +
  geom_text(aes(label = round(PPS, digits = 3)), color = "white", size = 3) +
  scale_fill_viridis_c(direction = -1, end = 0.8) + 
  facet_wrap(~Transmission) +
  theme_minimal() +
  ggtitle('Predictive Power Scores (PPS) & Matthew\'s Correlation Coeffecient (MCC)') +
  theme(legend.position = 'none', plot.title = element_text(hjust = 0.5, size = 12), 
        axis.title.x = element_blank(), axis.title.y = element_blank()) +
  scale_x_discrete(position = "top")


########################################################################################
########################################################################################
### End of Code
########################################################################################
########################################################################################
