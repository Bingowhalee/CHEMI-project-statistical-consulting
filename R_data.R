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

# Readable genus names (renaming function provided by the client)
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
any(duplicated(nm))                           
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


# ---- 3. Exploratory analysis --------------------------------------
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

gt_main <- run_gt(X)
gt_main

gt_main$sig <- ifelse(gt_main$q < 0.05, "FDR q < 0.05", "Not significant")
ggplot(gt_main, aes(reorder(outcome_lab, -p), -log10(p), fill = sig)) +
  geom_col() + coord_flip() +
  scale_fill_manual(values = c("FDR q < 0.05" = "#1F4E79", "Not significant" = "grey75"), name = NULL) +
  labs(x = NULL, y = expression(-log[10]~"(p-value, global test)")) +
  theme_minimal(base_size = 13) + theme(legend.position = "bottom")
ggsave("Final/hierarchy_globaltest.png", width = 8, height = 6, dpi = 300, bg = "white")


# 6a. Supplementary: elastic net + LOOCV
covs <- model.matrix(~ stage + log_crea)[, -1]
Xen  <- cbind(covs, scale(X))
pen  <- c(rep(0, ncol(covs)), rep(1, ncol(X)))       # covariates are never penalised
fit_enet <- function(x, y) cv.glmnet(x, y, alpha = 0.5, penalty.factor = pen, nfolds = 10)

run_loocv <- function(m) {
  y <- Y[, m]; n <- length(y)
  pred <- base <- numeric(n)
  sel  <- matrix(0, n, ncol(X), dimnames = list(NULL, colnames(X)))
  for (i in 1:n) {
    cv       <- fit_enet(Xen[-i, ], y[-i])                                 # penalty chosen without patient i
    pred[i]  <- predict(cv, Xen[i, , drop = FALSE], s = "lambda.min")
    sel[i, ] <- coef(cv, s = "lambda.min")[colnames(X), 1] != 0
    base[i]  <- sum(coef(lm(y[-i] ~ covs[-i, ])) * c(1, covs[i, ]))       # covariates only
  }
  r2 <- function(p) 1 - sum((y - p)^2) / sum((y - mean(y))^2)
  list(R2_enet = r2(pred), R2_base = r2(base), stab = colMeans(sel),
       beta = coef(fit_enet(Xen, y), s = "lambda.min")[colnames(X), 1])
}

set.seed(SEED)
loo <- setNames(lapply(outcomes, run_loocv), outcomes)                   # takes a few minutes

res_enet <- data.frame(outcome = lab[outcomes],
                       R2_base = sapply(loo, `[[`, "R2_base"),
                       R2_enet = sapply(loo, `[[`, "R2_enet"))
res_enet$gain <- res_enet$R2_enet - res_enet$R2_base
res_enet[order(-res_enet$gain), ]

stable <- function(m, freq = 0.5) {                 # genera selected in >= freq of the 60 models
  s <- loo[[m]]$stab; keep <- names(s)[s >= freq]
  out <- data.frame(genus = nm[keep], freq = round(s[keep], 2), beta = round(loo[[m]]$beta[keep], 3))
  out[order(-out$freq), ]
}
plot_stable <- function(m) {
  s <- subset(stable(m, 0.9), abs(beta) >= 0.05)
  ggplot(s, aes(reorder(genus, beta), beta, fill = beta > 0)) +
    geom_col() + coord_flip() + ylim(-0.2, 0.25) +
    scale_fill_manual(values = c("TRUE" = "#B2182B", "FALSE" = "#2166AC"), guide = "none") +
    labs(x = NULL, y = "Elastic net coefficient", title = lab[m]) + theme_minimal(base_size = 13)
}
plot_stable("IPA_B") + plot_stable("Calpr_B")
#acetate and LBP
plot_stable("Ace_B")
plot_stable("LBP")
ggsave("Final/elastic_net_stable.png", width = 12, height = 5, dpi = 300, bg = "white")


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


#7. Sensitivity analyses
# (i) adjusted for stage only, (ii) genera in >= 30% of patients
reg_stage <- run_regression(X, adjust_crea = FALSE)
reg_30    <- run_regression(X30)
gt_stage  <- run_gt(X, adjust_crea = FALSE)
gt_30     <- run_gt(X30)

compare_reg <- function(sens, label) {          # regression: main vs sensitivity
  m <- merge(reg[, c("genus", "outcome", "beta")], sens[, c("genus", "outcome", "beta")],
             by = c("genus", "outcome"), suffixes = c("_main", "_sens"))
  cat(label, "| beta correlation:", round(cor(m$beta_main, m$beta_sens), 3),
      "| significant (q < 0.05): main", sum(reg$q < 0.05), "- sensitivity", sum(sens$q < 0.05), "\n")
}
compare_reg(reg_stage, "Stage only")
compare_reg(reg_30,    "30% threshold")
reg_30[reg_30$q < 0.05, show_cols]              # expected: more significant pairs (fewer tests)

ranks <- data.frame(outcome = lab[outcomes],    # global test ranking in each setting
                    main       = rank(gt_main$p)[match(outcomes, gt_main$outcome)],
                    stage_only = rank(gt_stage$p)[match(outcomes, gt_stage$outcome)],
                    prev_30    = rank(gt_30$p)[match(outcomes, gt_30$outcome)])
ranks[order(ranks$main), ]
cor(ranks[, -1], method = "spearman")


#8. Diagnostics for key associations
check_assoc <- function(genus_label, m) {
  g   <- names(nm)[nm == genus_label]
  x   <- scale(X[, g])[, 1]; y <- Y[, m]
  fit <- lm(y ~ x + stage + log_crea)
  par(mfrow = c(2, 2)); plot(fit, main = paste(genus_label, "-", lab[m])); par(mfrow = c(1, 1))
  infl <- which(cooks.distance(fit) > 4 / length(y))
  fit2 <- if (length(infl) > 0) lm(y[-infl] ~ x[-infl] + stage[-infl] + log_crea[-infl]) else fit
  sp   <- cor.test(X[, g], y, method = "spearman", exact = FALSE)
  list(patients_with_genus = sum(X[, g] > 0), influential = infl,
       all_patients        = coef(summary(fit))[2, c(1, 4)],
       without_influential = coef(summary(fit2))[2, c(1, 4)],
       spearman            = c(rho = unname(sp$estimate), p = sp$p.value))
}
check_assoc("Coprococcus 2", "IPA_B")
check_assoc("[Ruminococcus] torques group", "Calpr_B")
