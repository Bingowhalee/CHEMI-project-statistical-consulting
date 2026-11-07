# =============================================================
# CHEMI project - gut microbiota vs plasma metabolites and gut health markers
#
# 0. Setup
# 1. Prepare data          (filter, transform, covariates, readable names)
# 2. Describe data         (key numbers, prevalence, distributions)
# 3. Exploratory analysis  (PCoA, stage, metabolite arrows)
# 4. Research question 1   (genus-by-genus regression + heatmap)
# 5. Research question 2   (global test ranking)
# 6. Supplementary         (elastic net + LOOCV, sparse PLS)
# 7. Sensitivity analyses  (stage only, 30% threshold, FDR per outcome)
# 8. Diagnostics           (key associations)
#
# Run one section at a time. Packages to install once:
#   install.packages(c("vegan", "ggplot2", "tidyr", "pheatmap", "glmnet", "patchwork"))
#   BiocManager::install(c("mixOmics", "globaltest"))
# =============================================================


# 0. Setup 
rm(list = ls())
library(mixOmics)     # sparse PLS
library(vegan)        # Hellinger, PERMANOVA, envfit
library(globaltest)   # global test
library(glmnet)       # elastic net
library(ggplot2); library(tidyr); library(pheatmap); library(patchwork)

SEED <- 4562550       # one seed for the whole analysis (reported in the methods)

data <- read.csv("Data_CHEMI_statistics_LU_20260929.csv", check.names = FALSE)


#1. Prepare data 
meta     <- data[, 1:18]                        # stage, creatinine, 16 outcomes
genus    <- as.matrix(data[, 19:ncol(data)])    # relative abundances (rows sum to 1)
outcomes <- setdiff(names(meta), c("stage", "Creat_B"))

genus_nz <- genus[, colSums(genus) > 0]         # drop genera absent in all patients
prev     <- colMeans(genus_nz > 0)              # prevalence of each genus

hellinger <- function(cutoff) decostand(genus_nz[, prev >= cutoff], method = "hellinger")
X   <- hellinger(0.10)                          # main analysis: genera in >= 10% of patients
X30 <- hellinger(0.30)                          # sensitivity analysis: >= 30%

Y        <- scale(log2(meta[, outcomes]))       # outcomes: log2, then mean 0 / SD 1
stage    <- factor(meta$stage)
log_crea <- log2(meta$Creat_B)

# Readable genus names
short_taxon <- function(x) {
  order  <- sub(".*o__([^_]+(?:_[^_]+)*)_f__.*", "\\1", x)
  family <- sub(".*f__([^_]+(?:_[^_]+)*)_g__.*", "\\1", x)
  genus  <- sub(".*g__", "", x)
  genus  <- gsub("^_+|_+$", "", genus)
  if (nchar(genus) > 0 && !genus %in% c("uncultured", "uncultured_bacterium")) return(genus)
  if (family != x && !family %in% c("uncultured", "uncultured_bacterium")) return(paste0(family, "_uncultured"))
  paste0(order, "_unclassified")
}


nm <- setNames(make.unique(sapply(colnames(X), short_taxon), sep = "_"), colnames(X))
                        
lab <- c(Ace_B = "Acetate", Pro_B = "Propionate", But_B = "Butyrate", Vlr_B = "Valerate",
         iVlr_B = "Isovalerate", iBut_B = "Isobutyrate", aMB_B = "a-Methylbutyrate",
         IS_B = "Indoxyl sulfate", IAM_B = "Indole-3-acetamide", IAA_B = "Indole-3-acetate",
         ILA_B = "Indole-3-lactate", IAld_B = "Indole-3-aldehyde", IPA_B = "Indole-3-propionate",
         TMAO_B = "TMAO", LBP_B = "LBP", Calpr_B = "Calprotectin")


# 2. Describe data
c(patients = nrow(data), genera = ncol(genus), non_empty = ncol(genus_nz),
  genera_10pct = ncol(X), genera_30pct = ncol(X30),
  zeros_pct = round(100 * mean(genus_nz == 0)))
table(stage)
# no missing or zero values expected
colSums(is.na(meta)); colSums(meta[, -1] <= 0)  

# save the rounded zeros distribution plot
png("Final/prevalence_histogram.png", width = 6, height = 5, units = "in", res = 300)
hist(prev, breaks = seq(0, 1, by = 0.05), col = "grey80", border = "black", main = "",
     xlab = "Proportion of patients in which the genus is present", ylab = "Number of genera",
     ylim = c(0, max(table(cut(prev, seq(0, 1, 0.05), include.lowest = TRUE))) + 5))
abline(v = 0.10, lty = 2, col = "grey40")
abline(v = 0.30, lty = 2, col = "red")
text(0.11, par("usr")[4] * 0.97, "10%", col = "grey40", adj = 0)
text(0.31, par("usr")[4] * 0.97, "30%", col = "red",     adj = 0)
dev.off()

met_long <- pivot_longer(meta[, outcomes], everything(), names_to = "metabolite")

#distribution of raw metabolite/markers concentration
ggplot(met_long, aes(value)) +           
  geom_density(fill = "red", alpha = 0.4) +
  facet_wrap(~ metabolite, scales = "free") + labs(x = "raw concentration") + theme_minimal()

#Distribution of log2 transformed metabolite/makers concentration
ggplot(met_long, aes(log2(value))) +           
  geom_density(fill = "steelblue", alpha = 0.4) +
  facet_wrap(~ metabolite, scales = "free") + labs(x = "log2 concentration") + theme_minimal()


# 3. Exploratory analysis 
pcoa <- cmdscale(dist(X), k = 2, eig = TRUE)
ve   <- round(100 * pcoa$eig[1:2] / sum(pcoa$eig[pcoa$eig > 0]), 1)   # % variation per axis

set.seed(SEED)
adonis2(dist(X) ~ stage, permutations = 9999)   # stage: difference in composition
anova(betadisper(dist(X), stage))               # stage: difference in spread

Y_env   <- Y; colnames(Y_env) <- lab[colnames(Y)]
fit_env <- envfit(pcoa$points, Y_env, permutations = 9999)
fit_env

#save pcoa plot with envfit arrows
png("Final/pcoa_arrows.png", width = 8, height = 7, units = "in", res = 300)
plot(pcoa$points, pch = 19, col = "grey60",
     xlab = paste0("PCoA 1 (", ve[1], "%)"), ylab = paste0("PCoA 2 (", ve[2], "%)"))
plot(fit_env, p.max = 0.05, col = "red")
dev.off()


# 4. Research question 1: genus-by-genus regression 
# One model per genus x outcome: outcome ~ genus + stage (+ creatinine)
run_regression <- function(X, adjust_crea = TRUE) {
  Xz  <- scale(X)
  res <- expand.grid(genus = colnames(Xz), outcome = outcomes, stringsAsFactors = FALSE)
  est <- t(mapply(function(g, m) {
    fit <- if (adjust_crea) lm(Y[, m] ~ Xz[, g] + stage + log_crea) else lm(Y[, m] ~ Xz[, g] + stage)
    c(coef(summary(fit))[2, c(1, 4)], confint(fit)[2, ])
  }, res$genus, res$outcome))
  res[, c("beta", "p", "lower", "upper")] <- est
  res$q     <- p.adjust(res$p, method = "BH")                                       # FDR across all tests
  res$q_out <- ave(res$p, res$outcome, FUN = function(p) p.adjust(p, method = "BH")) # FDR within outcome
  res$genus_lab   <- nm[res$genus]
  res$outcome_lab <- lab[res$outcome]
  res[order(res$p), ]
}

reg <- run_regression(X)
nrow(reg)                                                             # number of tests
show_cols <- c("genus_lab", "outcome_lab", "beta", "lower", "upper", "p", "q", "q_out")
head(reg[, show_cols], 20)                                            # strongest associations
reg[reg$q < 0.05, show_cols]                                          # FDR-significant (all tests)
reg[reg$q_out < 0.05, show_cols]                                      # FDR-significant (per outcome)
write.csv(reg, "Final/regression_results.csv", row.names = FALSE)

# Heatmap: genera with at least one p < 0.01; # = FDR q < 0.05
plot_heatmap <- function(res, file) {
  r    <- res[res$genus %in% res$genus[res$p < 0.01], ]
  beta <- tapply(r$beta, list(r$genus_lab, r$outcome_lab), mean)
  sig  <- tapply(r$q,    list(r$genus_lab, r$outcome_lab), mean) < 0.05
  lim  <- max(abs(beta))
  pheatmap(beta, color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
           breaks = seq(-lim, lim, length.out = 101), display_numbers = ifelse(sig, "#", ""),
           filename = file, width = 8, height = 10)
}
plot_heatmap(reg, "Final/heatmap_regression.png")


# 5. Research question 2: global test
# One test per outcome: do all genera together predict the outcome?
run_gt <- function(X, adjust_crea = TRUE) {
  null_f <- if (adjust_crea) ~ stage + log_crea else ~ stage
  p <- sapply(outcomes, function(m)
    p.value(gt(Y[, m], alternative = X, null = null_f, data = data.frame(stage, log_crea))))
  res <- data.frame(outcome = outcomes, outcome_lab = lab[outcomes], p = p,
                    q = p.adjust(p, method = "BH"))
  res[order(res$p), ]
}

set.seed(SEED)
gt_main <- run_gt(X)
gt_main

# Figure for the slides: strength of evidence per outcome (longer bar = stronger evidence)
plot_gt <- function(res, file = "Final/hierarchy_globaltest.png") {
  res$sig   <- res$q < 0.05
  res$score <- -log10(res$q)                                   # evidence on the FDR-adjusted scale
  res$label <- ifelse(res$q < 0.001, "q < 0.001", sprintf("q = %.3f", res$q))
  g <- ggplot(res, aes(score, reorder(outcome_lab, score), fill = sig)) +
    geom_col(width = 0.7) +
    geom_vline(xintercept = -log10(0.05), linetype = 2, colour = "grey30") +
    annotate("text", x = -log10(0.05), y = nrow(res) + 0.9, label = "significance threshold (q = 0.05)", hjust = -0.05, size = 4, colour = "grey30") +
    geom_label(aes(label = label), hjust = -0.05, size = 4, fill = "white", label.size = 0) +   # white box keeps the dashed line from crossing the labels
    scale_y_discrete(expand = expansion(add = c(0.6, 1.2))) +
    scale_fill_manual(values = c(`TRUE` = "#0B1B5E", `FALSE` = "grey75"), guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.25))) +
    labs(x = "Strength of evidence  (longer bar = stronger)", y = NULL) +
    theme_minimal(base_size = 16) +
    theme(panel.grid = element_blank(), axis.text.x = element_blank(),      # the numbers on the axis (-log10 q) mean little to a client,
          axis.title.x = element_text(margin = margin(t = 10)))                # so only the q-values next to the bars are shown
  ggsave(file, g, width = 9, height = 5.8, dpi = 300, bg = "white")
  g
}
plot_gt(gt_main)

# ---- 6a. Supplementary: elastic net + leave-one-out cross-validation ----
# One elastic net per outcome with ALL genera in the same model (stage and creatinine never penalised).
# alpha = 0.5 (between ridge and lasso); the penalty strength is chosen by 10-fold CV.
covs <- model.matrix(~ stage + log_crea)[, -1]            # stage II, stage III, creatinine
Xen  <- cbind(covs, scale(X))
pen  <- c(rep(0, ncol(covs)), rep(1, ncol(X)))            # 0 = covariates always kept
fit_enet <- function(x, y) cv.glmnet(x, y, alpha = 0.5, penalty.factor = pen, nfolds = 10)

# For one outcome: leave each patient out, fit on the other 59, predict the left-out patient.
run_loocv <- function(m) {
  y <- Y[, m]; n <- length(y)
  pred <- base <- numeric(n)
  sel  <- matrix(0, n, ncol(X), dimnames = list(NULL, colnames(X)))
  for (i in 1:n) {
    cv       <- fit_enet(Xen[-i, ], y[-i])                                     # penalty chosen without patient i
    pred[i]  <- predict(cv, Xen[i, , drop = FALSE], s = "lambda.min")
    sel[i, ] <- coef(cv, s = "lambda.min")[colnames(X), 1] != 0
    base[i]  <- sum(coef(lm(y[-i] ~ covs[-i, ])) * c(1, covs[i, ]))           # stage + creatinine only
  }
  r2 <- function(p) 1 - sum((y - p)^2) / sum((y - mean(y))^2)
  list(R2_enet = r2(pred), R2_base = r2(base),
       stab = colMeans(sel),                                                    # selection frequency per genus
       beta = coef(fit_enet(Xen, y), s = "lambda.min")[colnames(X), 1])         # direction (all 60 patients)
}

set.seed(SEED)
loo <- setNames(lapply(outcomes, run_loocv), outcomes)                          # takes a few minutes

# Does the microbiota improve the prediction beyond stage + creatinine? (R2 < 0 = worse than the average)
res_enet <- data.frame(outcome = lab[outcomes],
                       R2_base = sapply(loo, `[[`, "R2_base"),
                       R2_enet = sapply(loo, `[[`, "R2_enet"))
res_enet$gain <- res_enet$R2_enet - res_enet$R2_base
res_enet[order(-res_enet$gain), ]

# Genera selected in >= freq of the 60 leave-one-out models
stable <- function(m, freq = 0.5) {
  s <- loo[[m]]$stab; keep <- names(s)[s >= freq]
  out <- data.frame(genus = nm[keep], freq = round(s[keep], 2), beta = round(loo[[m]]$beta[keep], 3))
  out <- out[order(-out$freq), ]; rownames(out) <- NULL; out
}
stable("IPA_B"); stable("Calpr_B")

# Figure: stable genera (selected in >= 90% of the models, |coefficient| >= 0.05)
plot_stable <- function(m) {
  s <- subset(stable(m, 0.9), abs(beta) >= 0.05)
  ggplot(s, aes(reorder(genus, beta), beta, fill = beta > 0)) +
    geom_col() + coord_flip(ylim = c(-0.3, 0.3)) +          # was ylim(-0.2, 0.25): it dropped bars
    scale_fill_manual(values = c("TRUE" = "#B2182B", "FALSE" = "#2166AC"), guide = "none") +
    labs(x = NULL, y = "Elastic net coefficient", title = lab[m]) + theme_minimal(base_size = 13)
}
plot_stable("IPA_B") + plot_stable("Calpr_B")
ggsave("Final/elastic_net_stable.png", width = 12, height = 5, dpi = 300, bg = "white")


# Small table for the slides: prediction with and without the microbiota
tab_enet <- data.frame(Outcome = res_enet$outcome,
                       `Stage + creatinine only` = round(res_enet$R2_base, 2),
                       `+ gut bacteria`          = round(res_enet$R2_enet, 2),
                       Gain                      = round(res_enet$gain, 2), check.names = FALSE)
tab_enet[order(-tab_enet$Gain), ]

# ---- Elastic net with a chosen set of covariates (sensitivity: stage only) ----
run_loocv_cov <- function(m, covs) {
  Xen2 <- cbind(covs, scale(X))
  pen2 <- c(rep(0, ncol(covs)), rep(1, ncol(X)))                       # covariates never penalised
  fit2 <- function(x, y) cv.glmnet(x, y, alpha = 0.5, penalty.factor = pen2, nfolds = 10)
  y <- Y[, m]; n <- length(y)
  pred <- base <- numeric(n)
  sel  <- matrix(0, n, ncol(X), dimnames = list(NULL, colnames(X)))
  for (i in 1:n) {
    cv       <- fit2(Xen2[-i, ], y[-i])
    pred[i]  <- predict(cv, Xen2[i, , drop = FALSE], s = "lambda.min")
    sel[i, ] <- coef(cv, s = "lambda.min")[colnames(X), 1] != 0
    base[i]  <- sum(coef(lm(y[-i] ~ covs[-i, , drop = FALSE])) * c(1, covs[i, ]))
  }
  r2 <- function(p) 1 - sum((y - p)^2) / sum((y - mean(y))^2)
  list(R2_enet = r2(pred), R2_base = r2(base), stab = colMeans(sel),
       beta = coef(fit2(Xen2, y), s = "lambda.min")[colnames(X), 1])
}

covs_stage <- model.matrix(~ stage)[, -1, drop = FALSE]                  # stage only (no creatinine)
key <- c("Calpr_B", "LBP_B")
set.seed(SEED)
loo_stage <- setNames(lapply(key, run_loocv_cov, covs = covs_stage), key)

# 1. Prediction: main (stage + creatinine) versus stage only
cmp_enet <- data.frame(outcome         = lab[key],
                       base_main       = sapply(loo[key],       `[[`, "R2_base"),
                       enet_main       = sapply(loo[key],       `[[`, "R2_enet"),
                       base_stage_only = sapply(loo_stage[key], `[[`, "R2_base"),
                       enet_stage_only = sapply(loo_stage[key], `[[`, "R2_enet"))
cmp_enet[, -1] <- round(cmp_enet[, -1], 2)
cmp_enet

# 2. Stable genera (selected in >= 50% of the 60 models), stage only
stable_stage <- function(m, freq = 0.5) {
  s <- loo_stage[[m]]$stab; keep <- names(s)[s >= freq]
  out <- data.frame(genus = nm[keep], freq = round(s[keep], 2), beta = round(loo_stage[[m]]$beta[keep], 3))
  out <- out[order(-out$freq), ]; rownames(out) <- NULL; out
}
stable_stage("Calpr_B"); stable_stage("LBP_B")
plot_stable("Calpr_B")

# 3. Do the same genera come out? (main versus stage only)
for (m in key) {
  a <- stable(m)$genus; b <- stable_stage(m)$genus
  cat(lab[m], ": main", length(a), "genera | stage only", length(b), "| in both", length(intersect(a, b)), "\n")
}

# 6b. Supplementary: sparse PLS
X_spls <- X; colnames(X_spls) <- nm[colnames(X)]
Y_spls <- Y; colnames(Y_spls) <- lab[colnames(Y)]
nzv <- nearZeroVar(X_spls)$Position
if (length(nzv) > 0) X_spls <- X_spls[, -nzv]          # near-zero variance genera break the CV

set.seed(SEED)
tune <- tune.spls(X_spls, Y_spls, ncomp = 2,
                  test.keepX = c(5, 10, 20, 30, 50), test.keepY = c(3, 5, 8, 12, 16),
                  validation = "Mfold", folds = 5, nrepeat = 10, measure = "cor",
                  BPPARAM = BiocParallel::SerialParam())
fit_spls <- spls(X_spls, Y_spls, ncomp = 2, keepX = tune$choice.keepX, keepY = tune$choice.keepY)

png("Final/spls_heatmap.png", width = 9, height = 9, units = "in", res = 300)
cim(fit_spls, comp = 1:2, margins = c(10, 15))
dev.off()
#Sensitivity Analysis:
# 30% prevalence
X_spls30 <- X30; colnames(X_spls30) <- nm[colnames(X30)]
nzv30 <- nearZeroVar(X_spls30)$Position
if (length(nzv30) > 0) X_spls30 <- X_spls30[, -nzv30]          # near-zero variance genera break the CV

set.seed(SEED)
tune30 <- tune.spls(X_spls30, Y_spls, ncomp = 2,
                  test.keepX = c(5, 10, 20, 30, 50), test.keepY = c(3, 5, 8, 12, 16),
                  validation = "Mfold", folds = 5, nrepeat = 10, measure = "cor",
                  BPPARAM = BiocParallel::SerialParam())
fit_spls30 <- spls(X_spls30, Y_spls, ncomp = 2, keepX = tune30$choice.keepX, keepY = tune30$choice.keepY)

png("Final/spls_heatmap30.png", width = 9, height = 9, units = "in", res = 300)
cim(fit_spls30, comp = 1:2, margins = c(10, 15))
dev.off()



# ---- 7. Sensitivity analyses ---------------------------------------
# (i) adjusted for stage only, (ii) genera in >= 30% of patients
reg_stage <- run_regression(X, adjust_crea = FALSE)
reg_30    <- run_regression(X30)
gt_stage  <- run_gt(X, adjust_crea = FALSE)
gt_30     <- run_gt(X30)


# Heatmap: genera with at least one p < 0.01; # = FDR q < 0.05
plot_heatmap(reg_30,"Final/heatmap_regression_30.png")
plot_heatmap(reg_stage,"Final/heatmap_regression_stage.png")
# 7a. One table for the client: do the key findings hold?
key_pair <- function(res, genus_label, m) {
  r <- res[res$genus_lab == genus_label & res$outcome == m, ]
  if (nrow(r) == 0) return("not tested")                  # genus removed by the filter
  sprintf("%.2f (q = %.3f)", r$beta, r$q)
}
settings <- list(`Main analysis`       = list(reg,       gt_main),
                 `Stage only`          = list(reg_stage, gt_stage),
                 `Genera in >= 30%`    = list(reg_30,    gt_30))
sens_table <- sapply(settings, function(s) {
  r <- s[[1]]; g <- s[[2]]                                # g is sorted by p-value, so position = rank
  c(`Coprococcus 2 - IPA: beta (q)`               = key_pair(r, "Coprococcus_2", "IPA_B"),
    `[R.] torques group - calprotectin: beta (q)` = key_pair(r, "[Ruminococcus]_torques_group", "Calpr_B"),
    `Rank of IPA (global test, out of 16)`        = which(g$outcome == "IPA_B"),
    `Rank of calprotectin (global test)`          = which(g$outcome == "Calpr_B"),
    `Significant genus-metabolite pairs`          = sum(r$q < 0.05),
    `Number of tests`                             = nrow(r))
})
sens_table
write.csv(sens_table, "Final/sensitivity_table.csv")

# 7b. One figure: each dot is a genus-metabolite pair; dots on the red line = same result
pair_beta <- function(sens, label) {
  m <- merge(reg[, c("genus", "outcome", "beta")], sens[, c("genus", "outcome", "beta")],
             by = c("genus", "outcome"), suffixes = c("_main", "_sens"))
  m$setting <- label
  m
}
pb <- rbind(pair_beta(reg_stage, "Adjusted for stage only"),
            pair_beta(reg_30,    "Genera in >= 30% of patients"))
ggplot(pb, aes(beta_main, beta_sens)) +
  geom_abline(linetype = 2, colour = "red") +
  geom_point(alpha = 0.4, size = 1) +
  facet_wrap(~ setting) + coord_equal() +
  labs(x = "Coefficient, main analysis", y = "Coefficient, sensitivity analysis") +
  theme_minimal(base_size = 12)
ggsave("Final/sensitivity_scatter.png", width = 9, height = 4.5, dpi = 300, bg = "white")


# 7c. New or lost associations compared with the main analysis
# status: "new" = significant only in the sensitivity analysis, "lost" = only in the main analysis
changes <- function(sens, label) {
  m <- merge(reg[,  c("genus_lab", "outcome_lab", "beta", "p", "q")],
             sens[, c("genus_lab", "outcome_lab", "beta", "p", "q")],
             by = c("genus_lab", "outcome_lab"), suffixes = c("_main", "_sens"))
  m$status <- ifelse(m$q_sens < 0.05 & m$q_main >= 0.05, "new",
                     ifelse(m$q_main < 0.05 & m$q_sens >= 0.05, "lost",
                            ifelse(m$q_main < 0.05 & m$q_sens < 0.05, "both", NA)))
  m <- m[!is.na(m$status), ]
  if (nrow(m) == 0) { cat(label, ": no significant associations in either analysis\n"); return(invisible(m)) }
  m$setting <- label
  m[order(m$status, m$p_sens), c("setting", "status", "genus_lab", "outcome_lab",
                                 "beta_main", "beta_sens", "p_main", "p_sens", "q_main", "q_sens")]
}
changes(reg_stage, "Stage only")
changes(reg_30,    "Genera in >= 30%")

# Global test: q-values of all outcomes in each setting (new = significant only in a sensitivity analysis)
gt_q <- data.frame(outcome    = lab[outcomes],
                   main       = gt_main$q[match(outcomes, gt_main$outcome)],
                   stage_only = gt_stage$q[match(outcomes, gt_stage$outcome)],
                   prev_30    = gt_30$q[match(outcomes, gt_30$outcome)])
gt_q[order(gt_q$main), ]


# ---- 8. Diagnostics for key associations ---------------------------
# Six panels: residuals vs fitted, normal Q-Q, scale-location, residuals vs leverage,
# the scatter plot with fitted line, and Cook's distance per patient.
# Set file = "name.png" to save the figure; leave file = NULL to draw it on screen.
check_assoc <- function(genus_label, m, file = NULL) {
  g <- names(nm)[nm == genus_label]
  if (length(g) != 1)
    stop("No unique match for '", genus_label, "'. Similar names: ",
         paste(nm[grepl(substr(genus_label, 1, 6), nm, fixed = TRUE)], collapse = ", "))
  x <- scale(X[, g])[, 1]; y <- Y[, m]; n <- length(y)
  fit <- lm(y ~ x + stage + log_crea)
  cd  <- cooks.distance(fit)
  infl <- which(cd > 4 / n)
  
  if (!is.null(file)) png(file, width = 11, height = 7, units = "in", res = 300)
  op <- par(mfrow = c(2, 3), mar = c(4.2, 4.2, 3, 1), oma = c(0, 0, 2, 0))
  plot(fit, which = c(1, 2, 3, 5), sub.caption = "")                                  # the four standard residual plots
  plot(X[, g], y, pch = 19, col = ifelse(X[, g] > 0, "steelblue", "grey60"),
       xlab = paste(genus_label, "(Hellinger)"), ylab = paste(lab[m], "(standardised log2)"),
       main = "Data and fitted line")
  abline(lm(y ~ X[, g]), col = "red")
  plot(cd, type = "h", xlab = "Patient", ylab = "Cook's distance", main = "Influential patients")
  abline(h = 4 / n, lty = 2, col = "red")
  mtext(paste(genus_label, "-", lab[m]), outer = TRUE, side = 3, line = 0, cex = 1.1, font = 2)
  par(op)
  if (!is.null(file)) dev.off()
  
  fit2 <- if (length(infl) > 0) lm(y[-infl] ~ x[-infl] + stage[-infl] + log_crea[-infl]) else fit
  sp   <- cor.test(X[, g], y, method = "spearman", exact = FALSE)
  list(patients_with_genus = sum(X[, g] > 0), influential = infl,
       all_patients        = coef(summary(fit))[2, c(1, 4)],
       without_influential = coef(summary(fit2))[2, c(1, 4)],
       spearman            = c(rho = unname(sp$estimate), p = sp$p.value))
}

graphics.off()                                   # close any device left open (e.g. an unfinished png())
check_assoc("Coprococcus_2", "IPA_B")                                                  # plots on screen
check_assoc("Coprococcus_2", "IPA_B", file = "diagnostics_Coprococcus2_IPA.png")       # saved version
check_assoc("[Ruminococcus]_torques_group", "Calpr_B", file = "diagnostics_torques_calprotectin.png")


# 9. Acetate and LBP 
# Results of each method are reported separately, as for the other outcomes.

# 9a. Are they related to IPA and calprotectin themselves?
cor.test(log2(meta$Ace_B), log2(meta$IPA_B),   method = "spearman", exact = FALSE)
cor.test(log2(meta$LBP_B), log2(meta$Calpr_B), method = "spearman", exact = FALSE)

# 9b. Genus-by-genus regression (main analysis): ten strongest associations per outcome
top_reg <- function(m, top = 10) head(reg[reg$outcome == m, show_cols], top)
top_reg("Ace_B")
top_reg("LBP_B")

# 9c. Elastic net: genera selected in >= 50% of the 60 leave-one-out models
stable("Ace_B")
stable("LBP_B")
res_enet[c("Ace_B", "LBP_B"), ]                  # predictive R2 with and without the microbiota

# 9d. Sparse PLS: see the heat map (spls_heatmap.png); acetate and LBP are columns of the plot

# 9e. Figure: ten strongest genera with 95% CI (red = FDR q < 0.05)
forest <- function(m, top = 10) {
  r <- head(reg[reg$outcome == m, ], top)
  ggplot(r, aes(beta, reorder(genus_lab, beta))) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.2, orientation = "y") +
    geom_point(aes(colour = q < 0.05), size = 2.5) +
    scale_colour_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "grey40"), guide = "none") +
    labs(x = "Standardised coefficient (95% CI)", y = NULL, title = lab[m]) +
    theme_minimal(base_size = 12)
}
forest("Ace_B") + forest("LBP_B")
ggsave("Final/acetate_lbp_genera.png", width = 12, height = 5, dpi = 300, bg = "white")

# ---- 10. Comparison: genus-by-genus regression vs elastic net ------
# For one outcome: the k strongest genera of each method side by side.
# Regression: ranked by p-value. Elastic net: ranked by selection frequency, then |coefficient|.
compare_methods <- function(m, k = 5) {
  rr <- reg[reg$outcome == m, ]                          # sorted by p-value
  s  <- loo[[m]]$stab; b <- loo[[m]]$beta
  top_r <- head(rr$genus, k)
  top_e <- names(s)[order(-s, -abs(b))][1:k]
  top_e <- top_e[s[top_e] > 0]                          # ignore genera never selected
  g  <- union(top_r, top_e)
  rr <- rr[match(g, rr$genus), ]
  data.frame(outcome   = lab[m],
             genus     = nm[g],
             top_in    = ifelse(g %in% top_r & g %in% top_e, "both",
                                ifelse(g %in% top_r, "regression", "elastic net")),
             reg_beta  = round(rr$beta, 2),
             reg_q     = signif(rr$q, 2),
             enet_freq = round(s[g], 2),
             enet_beta = round(b[g], 3),
             same_sign = ifelse(b[g] == 0, NA, sign(rr$beta) == sign(b[g])),
             row.names = NULL)
}

# Key outcomes in one table
cmp_key <- do.call(rbind, lapply(c("IPA_B", "Calpr_B", "Ace_B", "LBP_B"), compare_methods))
cmp_key
write.csv(cmp_key, "comparison_regression_elasticnet.csv", row.names = FALSE)

# Summary for all 16 outcomes: overlap of the top 5 and the predictive gain of the elastic net
cmp_summary <- do.call(rbind, lapply(outcomes, function(m) {
  x <- compare_methods(m)
  data.frame(outcome           = lab[m],
             sig_regression    = sum(reg$outcome == m & reg$q < 0.05),   # FDR-significant genera
             top5_in_both      = sum(x$top_in == "both"),
             same_sign_in_both = sum(x$same_sign[x$top_in == "both"], na.rm = TRUE),
             enet_R2_gain      = round(res_enet[m, "gain"], 2))
}))
cmp_summary[order(-cmp_summary$top5_in_both, -cmp_summary$enet_R2_gain), ]
