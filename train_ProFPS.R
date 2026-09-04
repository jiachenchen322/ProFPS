################################################################################
## ProFPS: a plasma proteomic signature of frailty pace
##
## Training pipeline in the FHS Offspring cohort (N = 1,553):
##   Step 1  Olink NPX QC and restriction to the 2,799 proteins shared across
##           FHS Offspring, FHS Gen 3 and UK Biobank
##   Step 2  kNN imputation + rank-based inverse-normal transformation
##   Step 3  Merge repeated physical measures (Exams 8/9/10) and standardize
##   Step 4  Participant-specific slopes for body weight, gait and grip (LME)
##   Step 5  PCA over the three component slopes -> frailty pace (PC1)
##   Step 6  Elastic net regression of frailty pace on protein levels
##   Step 7  Export the final 112-protein model
##
## Individual-level FHS data are not distributable and must be obtained through
## dbGaP; see README.md. This script is provided so that the derivation of the
## frailty pace phenotype and the ProFPS model can be inspected and reproduced
## by investigators with approved access to the source data.
################################################################################

rm(list = ls())

library(caret)     # knnImpute via preProcess
library(nlme)      # linear mixed-effects models
library(glmnet)    # elastic net
library(haven)     # read_sas
library(tidyr)     # pivot_longer
library(dplyr)

## ---------------------------------------------------------------------------
## Configuration: point these at your own copies of the source data.
## ---------------------------------------------------------------------------
data_dir <- "path/to/data"          # <-- EDIT
out_dir  <- "."

f_frailty   <- file.path(data_dir, "FHS_FRAILTY_PHENOTYPE_ex8910_timepointsinfoweight.csv")
f_shared    <- file.path(data_dir, "impute.rda")        # defines the 2,799 shared proteins
f_olink     <- file.path(data_dir, "FHS_proteomics_organize.rda")   # dataGen2 (long NPX)
f_olinkwide <- file.path(data_dir, "dataGen2_wide.rda")            # dataGen2_wide (wide NPX)
f_linker    <- file.path(data_dir, "FHS_OlinkHT_Mapping_09032024_framid.csv")
f_apoe      <- file.path(data_dir, "apoe_genotype.sas7bdat")   # FHS Offspring APOE genotype
f_wkthru    <- file.path(data_dir, "exam_covariates.sas7bdat") # FHS Offspring covariates

set.seed(12345678)

## Rank-based inverse-normal transformation (Blom)
inormal <- function(x) {
  qnorm((rank(x, na.last = "keep") - 3/8) / (sum(!is.na(x)) - 1/4))
}

## ===========================================================================
## Step 1. Olink NPX: QC filtering and the 2,799 shared proteins
## ===========================================================================
frailty8910 <- read.csv(f_frailty)

## `complete` holds the harmonized protein set common to FHS and UKB.
load(f_shared)
complete <- rbind(complete_train, complete_test)
ukb_protein_list <- colnames(complete)[1:2799]

load(f_olink)                       # -> dataGen2
IDmatch <- read.csv(f_linker)
dataGen2$SampleID <- as.numeric(dataGen2$SampleID)

## Drop samples failing Olink sample QC; keep assay (non-control) measurements.
dataGen2_qc <- dataGen2[dataGen2$SampleQC != "FAIL", ]
dataGen2_qc_framid <- merge(dataGen2_qc, IDmatch[, c(1, 4)], by = "SampleID")
dataGen2_qc_framid_assay <-
  dataGen2_qc_framid[dataGen2_qc_framid$AssayType == "assay", ]

load(f_olinkwide)                   # -> dataGen2_wide (framid x OlinkID)
colnames(dataGen2_wide)[2:5348] <- gsub("NPX\\.", "", colnames(dataGen2_wide)[2:5348])

unique_olink_uniprot <-
  dataGen2_qc_framid_assay[!duplicated(dataGen2_qc_framid_assay$OlinkID),
                           c("OlinkID", "UniProt")]

## Two UniProt accessions are measured by more than one Olink assay. For each,
## retain the single assay carried forward to the shared protein set
## (OID43743 for P32455 and OID43956 for Q02750); the redundant assays are
## dropped so that protein identifiers map one-to-one.
drop_oid <- c("OID42498", "OID44673", "OID44768", "OID42734")
useOLINKID <- unique_olink_uniprot[!unique_olink_uniprot$OlinkID %in% drop_oid, ]

## Relabel columns by UniProt accession and keep only the shared proteins.
colnames(dataGen2_wide)[-1] <-
  useOLINKID$UniProt[match(colnames(dataGen2_wide)[-1], useOLINKID$OlinkID)]
dataGen2_wide <- dataGen2_wide[, !is.na(colnames(dataGen2_wide))]
stopifnot(sum(colnames(dataGen2_wide) %in% ukb_protein_list) == 2799)

dataGen2_wide_use <- dataGen2_wide[, colnames(dataGen2_wide) %in% ukb_protein_list]
dataGen2_wide_use$framid <- dataGen2_wide$framid   # framid becomes column 2800

## ===========================================================================
## Step 2. Imputation and inverse-normal transformation
## ===========================================================================
## k = 45 ~ sqrt(n) for this sample.
impute_train <- preProcess(dataGen2_wide_use[, -2800], method = "knnImpute", k = 45)

complete_FHS <- predict(impute_train, dataGen2_wide_use[, -2800])
complete_FHS <- as.data.frame(apply(complete_FHS, 2, inormal))
complete_FHS$framid <- dataGen2_wide_use$framid

## ===========================================================================
## Step 3. Merge Exam 8 covariates and repeated physical measures
## ===========================================================================
offspringAPOE <- read_sas(f_apoe)
offspringAPOE <- offspringAPOE[offspringAPOE$idtype == 1, ]
colnames(offspringAPOE)[2] <- "ID"

complete_FHS_APOE <- merge(complete_FHS, offspringAPOE[, c(1, 4)],
                           by = "framid", all.x = TRUE)

offspringfhs <- read_sas(f_wkthru)
offspringfhs_exam8 <- offspringfhs[, c(2, 3, which(grepl("8", colnames(offspringfhs)) == 1))]
offspringfhs_exam8$framid <- offspringfhs_exam8$ID + 80000
E8_complete_FHS_APOE <- merge(complete_FHS_APOE, offspringfhs_exam8, by = "framid")

frailty_comb <- merge(E8_complete_FHS_APOE, frailty8910, by = "framid")
frailty_comb_use <- frailty_comb[, c(1:2806, 2831:2833, 2840:2842, 2853:2854,
                                     2857, 2867:2868, 2871, 2881:2882, 2885)]

## NOTE: `walk_min` is walking TIME (minutes over a fixed course), so a larger
## value means slower gait. The column is named `walkspeed*` for historical
## reasons only; the sign convention in Step 4 accounts for this.
frailty_comb_use$walkspeed_ex8  <- frailty_comb_use$walk_min_ex8
frailty_comb_use$walkspeed_ex9  <- frailty_comb_use$walk_min_ex9
frailty_comb_use$walkspeed_ex10 <- frailty_comb_use$walk_min_ex10

colnames(frailty_comb_use)[2807:2824] <-
  c("age8", "age9", "age10", "bmi8", "bmi9", "bmi10",
    "grip8", "walkmin8", "weight8", "grip9", "walkmin9", "weight9",
    "grip10", "walkmin10", "weight10", "walkspeed8", "walkspeed9", "walkspeed10")

## Standardize every visit to the Exam 8 (baseline) mean and SD, so that all
## three visits share one scale and within-person change is preserved.
std_to_baseline <- function(df, vars, newnames) {
  m <- mean(df[[vars[1]]], na.rm = TRUE)
  s <- sd(df[[vars[1]]], na.rm = TRUE)
  for (i in seq_along(vars)) df[[newnames[i]]] <- (df[[vars[i]]] - m) / s
  df
}
frailty_comb_use <- std_to_baseline(frailty_comb_use,
  c("weight8", "weight9", "weight10"), c("stdweight8", "stdweight9", "stdweight10"))
frailty_comb_use <- std_to_baseline(frailty_comb_use,
  c("bmi8", "bmi9", "bmi10"), c("stdbmi8", "stdbmi9", "stdbmi10"))
frailty_comb_use <- std_to_baseline(frailty_comb_use,
  c("grip8", "grip9", "grip10"), c("stdgrip8", "stdgrip9", "stdgrip10"))
frailty_comb_use <- std_to_baseline(frailty_comb_use,
  c("walkspeed8", "walkspeed9", "walkspeed10"), c("stdwalk8", "stdwalk9", "stdwalk10"))

## Reshape to one row per participant-visit.
long_frail <- frailty_comb_use[c(2807:2836)] %>%
  pivot_longer(everything(), cols_vary = "slowest",
               names_to = c(".value", "set"), names_pattern = "([a-zA-Z]+)(\\d+)")
colnames(long_frail)[1] <- "visit"
long_frail$framid <- rep(frailty_comb_use$framid, 3)
long_frail$SEX    <- rep(frailty_comb_use$SEX, 3)

long_frail_com <- long_frail[complete.cases(long_frail), ]

## Slopes require at least two observed visits per participant.
keep_id <- names(table(long_frail_com$framid))[table(long_frail_com$framid) > 1]
long_frail_3 <- long_frail[long_frail$framid %in% keep_id, ]

## ===========================================================================
## Step 4. Participant-specific rates of change (random-slope LME)
## ===========================================================================
## Total slope = fixed age effect + participant random age effect. Each slope is
## then divided by the cohort mean slope, so that 1 = the average rate of change
## in the cohort. `flip` orients each component so that HIGHER = FASTER decline:
##   weight  declines with frailty          -> flip
##   grip    declines with frailty          -> flip
##   walktime increases with frailty        -> no flip
component_slope <- function(formula, data, name, flip) {
  fit <- lme(formula, random = ~ age | framid, data = data, na.action = na.omit)
  re  <- ranef(fit)
  total <- fixef(fit)["age"] + re$age
  scaled <- (if (flip) -total else total) / mean(total, na.rm = TRUE)
  out <- data.frame(framid = unique(fit$groups$framid), scaled)
  colnames(out)[2] <- name
  out
}

slope_wt   <- component_slope(stdweight ~ age + SEX,       long_frail_3, "scaled_slope_wt",   flip = TRUE)
slope_walk <- component_slope(stdwalk   ~ age + SEX,       long_frail_3, "scaled_slope_walk", flip = FALSE)
## Grip strength is additionally adjusted for BMI.
slope_grip <- component_slope(stdgrip   ~ age + SEX + bmi, long_frail_3, "scaled_slope_grip", flip = TRUE)

E8_complete_FHS_APOE_grip <- E8_complete_FHS_APOE %>%
  merge(slope_wt,   by = "framid") %>%
  merge(slope_walk, by = "framid") %>%
  merge(slope_grip, by = "framid")

## ===========================================================================
## Step 5. Frailty pace = first principal component of the three slopes
## ===========================================================================
pace_vars <- c("scaled_slope_wt", "scaled_slope_walk", "scaled_slope_grip")

pc <- prcomp(E8_complete_FHS_APOE_grip[, pace_vars], center = TRUE, scale. = TRUE)
summary(pc)
print(pc)

E8_complete_FHS_APOE_grip$PC1 <-
  predict(pc, E8_complete_FHS_APOE_grip[, pace_vars])[, 1]

## Scale frailty pace to unit SD.
E8_complete_FHS_APOE_grip$PC1 <-
  E8_complete_FHS_APOE_grip$PC1 / sd(E8_complete_FHS_APOE_grip$PC1)

cor.test(E8_complete_FHS_APOE_grip$AGE8, E8_complete_FHS_APOE_grip$PC1)

## ===========================================================================
## Step 6. Elastic net regression of frailty pace on protein levels
## ===========================================================================
## Both the elastic net mixing parameter (alpha) and the regularization
## parameter (lambda) were selected by grid search with five-fold
## cross-validation in the training sample, minimizing mean squared error. The
## final model was then refit in the full training sample at the selected
## (alpha, lambda).
set.seed(1234)

X_train <- as.matrix(E8_complete_FHS_APOE_grip[, 2:2800])   # 2,799 proteins
y_train <- E8_complete_FHS_APOE_grip$PC1
dim(X_train)                                                # expected: 1553 x 2799

alpha_grid <- seq(0.05, 1, by = 0.05)   # alpha = 1 is lasso; pure ridge (0) gives no sparsity
nfolds <- 5

## A single fold assignment shared by every candidate alpha, so that the
## cross-validated errors are compared on identical partitions of the sample.
foldid <- sample(rep(seq_len(nfolds), length.out = nrow(X_train)))

cv_list <- lapply(alpha_grid, function(a) {
  cv.glmnet(x = X_train, y = y_train,
            family = "gaussian", type.measure = "mse",
            alpha = a, foldid = foldid, lambda.min.ratio = 0.01,
            intercept = TRUE, standardize = TRUE)
})

## Cross-validated MSE at each alpha's best lambda
cv_grid <- do.call(rbind, lapply(seq_along(alpha_grid), function(i) {
  fit <- cv_list[[i]]
  j <- which(fit$lambda == fit$lambda.min)
  data.frame(alpha      = alpha_grid[i],
             lambda.min = fit$lambda.min,
             cvm        = fit$cvm[j],
             cvsd       = fit$cvsd[j],
             nzero      = as.integer(fit$nzero[j]))
}))
cv_grid

## Optimal tuning parameters over the joint (alpha, lambda) grid
best_i      <- which.min(cv_grid$cvm)
best_alpha  <- cv_grid$alpha[best_i]     # selected in the published model: 0.5
best_lambda <- cv_grid$lambda.min[best_i] # selected in the published model: ~0.08415
c(alpha = best_alpha, lambda = best_lambda,
  cv_mse = cv_grid$cvm[best_i], nzero = cv_grid$nzero[best_i])

age.fit <- cv_list[[best_i]]   # cv.glmnet object at the selected alpha

## CV MSE profile across the alpha grid (+/- 1 CV SE)
plot(cv_grid$alpha, cv_grid$cvm, type = "b", pch = 16,
     ylim = range(c(cv_grid$cvm - cv_grid$cvsd, cv_grid$cvm + cv_grid$cvsd)),
     xlab = "Elastic net mixing parameter (alpha)",
     ylab = "Five-fold CV mean squared error")
arrows(cv_grid$alpha, cv_grid$cvm - cv_grid$cvsd,
       cv_grid$alpha, cv_grid$cvm + cv_grid$cvsd,
       length = 0.03, angle = 90, code = 3, col = "grey60")
abline(v = best_alpha, col = "red", lty = 2)

## Lambda path at the selected alpha
plot(age.fit)

## Final model: refit in the full training sample at the selected parameters
age.fit.min <- glmnet(x = X_train, y = y_train,
                      family = "gaussian", alpha = best_alpha,
                      intercept = TRUE, standardize = TRUE,
                      lambda = best_lambda)

## Expected in the published model: 112 non-zero coefficients, intercept 0.05080
sum(as.matrix(age.fit.min$beta) != 0)
age.fit.min$a0

## ===========================================================================
## Step 7. In-sample fit and export of the final model
## ===========================================================================
sig <- X_train %*% as.matrix(age.fit.min$beta) + age.fit.min$a0

sqrt(mean((y_train - sig)^2))                               # RMSE
cor.test(sig, y_train)
1 - sum((y_train - sig)^2) / sum((y_train - mean(y_train))^2)   # R-squared

plot(y_train, sig, main = "FHS Offspring",
     xlab = "Frailty pace (PC1)", ylab = "ProFPS")
abline(lm(sig ~ y_train), col = "red")

beta <- as.matrix(age.fit.min$beta)[, 1]
ProFPS_coef <- data.frame(uniprot     = names(beta)[beta != 0],
                          coefficient = as.numeric(beta[beta != 0]))
ProFPS_coef <- ProFPS_coef[order(ProFPS_coef$coefficient), ]

## Only the protein coefficients are exported. ProFPS scores are standardized
## within cohort before analysis, so the model intercept shifts every score by
## the same constant and does not affect any downstream result.
write.csv(ProFPS_coef, file.path(out_dir, "ProFPS_112_coefficients.csv"),
          row.names = FALSE)

## The fixed model is applied without refitting to FHS Offspring (N = 2,078),
## FHS Gen 3 (N = 3,130) and UK Biobank (N = 44,949):
##   ProFPS = X %*% beta
## where X holds inverse-normal-transformed NPX values for the same proteins,
## column-matched by UniProt accession. ProFPS is then z-scored within cohort.
