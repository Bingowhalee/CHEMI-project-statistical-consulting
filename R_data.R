# =============================================================
# CHEMI project - gut microbiota vs plasma metabolites and gut health markers
#
# 0. Setup            1. Prepare data        2. Describe data
# 3. Exploration      4. RQ1: genus-by-genus regression
# 5. RQ2: global test 6. Supplementary (elastic net, sparse PLS)
# 7. Sensitivity      8. Diagnostics
#
# Packages (once):
#   install.packages(c("vegan", "ggplot2", "tidyr", "pheatmap", "glmnet", "patchwork"))
#   BiocManager::install(c("mixOmics", "globaltest"))
# =============================================================

# 0. Setup ----
rm(list = ls())
library(mixOmics); library(vegan); library(globaltest); library(glmnet)
library(ggplot2); library(tidyr); library(pheatmap); library(patchwork)

SEED <- 4562550                                   # one seed for the whole analysis
dir.create("Final", showWarnings = FALSE)
data <- read.csv("Data_CHEMI_statistics_LU_20260929.csv", check.names = FALSE)

# 1. Prepare data ----
meta     <- data[, 1:18]                          # stage, creatinine, 16 outcomes
genus    <- as.matrix(data[, 19:ncol(data)])      # relative abundances
outcomes <- setdiff(names(meta), c("stage", "Creat_B"))

genus_nz <- genus[, colSums(genus) > 0]
prev     <- colMeans(genus_nz > 0)
hellinger <- function(cutoff) decostand(genus_nz[, prev >= cutoff], method = "hellinger")
X   <- hellinger(0.10)                            # main: genera in >= 10% of patients
X30 <- hellinger(0.30)                            # sensitivity: >= 30%

Y        <- scale(log2(meta[, outcomes]))         # log2, then mean 0 / SD 1
stage    <- factor(meta$stage)
log_crea <- log2(meta$Creat_B)

short_taxon <- function(x) {
  order  <- sub(".*o__([^_]+(?:_[^_]+)*)_f__.*", "\\1", x)
  family <- sub(".*f__([^_]+(?:_[^_]+)*)_g__.*", "\\1", x)
  g      <- gsub("^_+|_+$", "", sub(".*g__", "", x))
  bad    <- c("uncultured", "uncultured_bacterium")
  if (nchar(g) > 0 && !g %in% bad) return(g)
  if (family != x && !family %in% bad) return(paste0(family, "_uncultured"))
  paste0(order, "_unclassified")
}
nm  <- setNames(make.unique(sapply(colnames(X), short_taxon), sep = "_"), colnames(X))
lab <- c(Ace_B = "Acetate", Pro_B = "Propionate", But_B = "Butyrate", Vlr_B = "Valerate",
         iVlr_B = "Isovalerate", iBut_B = "Isobutyrate", aMB_B = "a-Methylbutyrate",
         IS_B = "Indoxyl sulfate", IAM_B = "Indole-3-acetamide", IAA_B = "Indole-3-acetate",
         ILA_B = "Indole-3-lactate", IAld_B = "Indole-3-aldehyde", IPA_B = "Indole-3-propionate",
         TMAO_B = "TMAO", LBP_B = "LBP", Calpr_B = "Calprotectin")
Y_lab <- Y; colnames(Y_lab) <- lab[colnames(Y)]   # same outcomes, readable names
red_blue <- c("TRUE" = "#B2182B", "FALSE" = "#2166AC")

# 2. Describe data ----
c(patients = nrow(data), genera = ncol(genus), non_empty = ncol(genus_nz),
  genera_10pct = ncol(X), genera_30pct = ncol(X30), zeros_pct = round(100 * mean(genus_nz == 0)))
table(stage)
colSums(is.na(meta)); colSums(meta[, -1] <= 0)    # expect no missing / zero values

png("Final/prevalence_histogram.png", width = 6, height = 5, units = "in", res = 300)
hist(prev, breaks = seq(0, 1, 0.05), col = "grey80", main = "",
     xlab = "Proportion of patients in which the genus is present", ylab = "Number of genera",
     ylim = c(0, max(table(cut(prev, seq(0, 1, 0.05), include.lowest = TRUE))) + 5))
abline(v = c(0.10, 0.30), lty = 2, col = c("grey40", "red"))
text(c(0.11, 0.31), par("usr")[4] * 0.97, c("10%", "30%"), col = c("grey40", "red"), adj = 0)
dev.off()

met_long <- pivot_longer(meta[, outcomes], everything(), names_to = "metabolite")
ggplot(met_long, aes(log2(value))) + geom_density(fill = "steelblue", alpha = 0.4) +
  facet_wrap(~ metabolite, scales = "free") + labs(x = "log2 concentration") + theme_minimal()

# 3. Exploration: PCoA, stage, metabolite arrows ----
pcoa <- cmdscale(dist(X), k = 2, eig = TRUE)
ve   <- round(100 * pcoa$eig[1:2] / sum(pcoa$eig[pcoa$eig > 0]), 1)

set.seed(SEED)
adonis2(dist(X) ~ stage, permutations = 9999)     # composition differs by stage?
anova(betadisper(dist(X), stage))                 # spread differs by stage?
fit_env <- envfit(pcoa$points, Y_lab, permutations = 9999); fit_env

png("Final/pcoa_arrows.png", width = 8, height = 7, units = "in", res = 300)
plot(pcoa$points, pch = 19, col = "grey60",
     xlab = paste0("PCoA 1 (", ve[1], "%)"), ylab = paste0("PCoA 2 (", ve[2], "%)"))
plot(fit_env, p.max = 0.05, col = "red")
dev.off()

# 4. RQ1: genus-by-genus regression ----
# outcome ~ genus + stage (+ log2 creatinine); BH-FDR over all tests and within outcome
run_regression <- function(X, adjust_crea = TRUE) {
  Xz  <- scale(X)
  res <- expand.grid(genus = colnames(Xz), outcome = outcomes, stringsAsFactors = FALSE)
  res[, c("beta", "p", "lower", "upper")] <- t(mapply(function(g, m) {
    fit <- if (adjust_crea) lm(Y[, m] ~ Xz[, g] + stage + log_crea) else lm(Y[, m] ~ Xz[, g] + stage)
    c(coef(summary(fit))[2, c(1, 4)], confint(fit)[2, ])
  }, res$genus, res$outcome))
  res$q     <- p.adjust(res$p, "BH")
  res$q_out <- ave(res$p, res$outcome, FUN = function(p) p.adjust(p, "BH"))
  res$genus_lab <- nm[res$genus]; res$outcome_lab <- lab[res$outcome]
  res[order(res$p), ]
}
plot_heatmap <- function(res, file) {             # genera with >= 1 p < 0.01; # = q < 0.05
  r    <- res[res$genus %in% res$genus[res$p < 0.01], ]
  beta <- tapply(r$beta, list(r$genus_lab, r$outcome_lab), mean)
  sig  <- tapply(r$q,    list(r$genus_lab, r$outcome_lab), mean) < 0.05
  lim  <- max(abs(beta))
  pheatmap(beta, color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
           breaks = seq(-lim, lim, length.out = 101), display_numbers = ifelse(sig, "#", ""),
           filename = file, width = 8, height = 10)
}

reg <- run_regression(X)
show_cols <- c("genus_lab", "outcome_lab", "beta", "lower", "upper", "p", "q", "q_out")
head(reg[, show_cols], 20)
reg[reg$q < 0.05, show_cols]
write.csv(reg, "Final/regression_results.csv", row.names = FALSE)
plot_heatmap(reg, "Final/heatmap_regression.png")

# 5. RQ2: global test (all genera together, per outcome) ----
run_gt <- function(X, adjust_crea = TRUE) {
  null_f <- if (adjust_crea) ~ stage + log_crea else ~ stage
  p <- sapply(outcomes, function(m)
    p.value(gt(Y[, m], alternative = X, null = null_f, data = data.frame(stage, log_crea))))
  res <- data.frame(outcome = outcomes, outcome_lab = lab[outcomes], p = p, q = p.adjust(p, "BH"))
  res[order(res$p), ]
}
set.seed(SEED)
gt_main <- run_gt(X); gt_main

res <- transform(gt_main, sig = q < 0.05, score = -log10(q),
                 label = ifelse(q < 0.001, "q < 0.001", sprintf("q = %.3f", q)))
ggplot(res, aes(score, reorder(outcome_lab, score), fill = sig)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = -log10(0.05), linetype = 2, colour = "grey30") +
  annotate("text", x = -log10(0.05), y = nrow(res) + 0.9, label = "significance threshold (q = 0.05)",
           hjust = -0.05, size = 4, colour = "grey30") +
  geom_label(aes(label = label), hjust = -0.05, size = 4, fill = "white", label.size = 0) +
  scale_y_discrete(expand = expansion(add = c(0.6, 1.2))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.25))) +
  scale_fill_manual(values = c(`TRUE` = "#0B1B5E", `FALSE` = "grey75"), guide = "none") +
  labs(x = "Strength of evidence  (longer bar = stronger)", y = NULL) +
  theme_minimal(base_size = 16) +
  theme(panel.grid = element_blank(), axis.text.x = element_blank(),
        axis.title.x = element_text(margin = margin(t = 10)))
ggsave("Final/hierarchy_globaltest.png", width = 9, height = 5.8, dpi = 300, bg = "white")

# 6a. Elastic net + leave-one-out CV ----
# All genera in one model per outcome; covariates never penalised; alpha = 0.5, lambda by 10-fold CV
covs_main  <- model.matrix(~ stage + log_crea)[, -1]
covs_stage <- model.matrix(~ stage)[, -1, drop = FALSE]

run_loocv <- function(m, covs = covs_main) {
  Xen <- cbind(covs, scale(X)); pen <- c(rep(0, ncol(covs)), rep(1, ncol(X)))
  fit <- function(x, y) cv.glmnet(x, y, alpha = 0.5, penalty.factor = pen, nfolds = 10)
  y <- Y[, m]; n <- length(y); pred <- base <- numeric(n)
  sel <- matrix(0, n, ncol(X), dimnames = list(NULL, colnames(X)))
  for (i in 1:n) {
    cv       <- fit(Xen[-i, ], y[-i])
    pred[i]  <- predict(cv, Xen[i, , drop = FALSE], s = "lambda.min")
    sel[i, ] <- coef(cv, s = "lambda.min")[colnames(X), 1] != 0
    base[i]  <- sum(coef(lm(y[-i] ~ covs[-i, , drop = FALSE])) * c(1, covs[i, ]))
  }
  r2 <- function(p) 1 - sum((y - p)^2) / sum((y - mean(y))^2)
  list(R2_enet = r2(pred), R2_base = r2(base), stab = colMeans(sel),
       beta = coef(fit(Xen, y), s = "lambda.min")[colnames(X), 1])
}
set.seed(SEED)
loo <- setNames(lapply(outcomes, run_loocv), outcomes)          # takes a few minutes

res_enet <- data.frame(outcome = lab[outcomes],
                       R2_base = sapply(loo, `[[`, "R2_base"), R2_enet = sapply(loo, `[[`, "R2_enet"))
res_enet$gain <- res_enet$R2_enet - res_enet$R2_base
res_enet[order(-res_enet$gain), ]                               # R2 < 0 = worse than the mean

stable <- function(m, freq = 0.5, L = loo) {                    # genera selected in >= freq of 60 models
  s <- L[[m]]$stab; keep <- names(s)[s >= freq]
  out <- data.frame(genus = nm[keep], freq = round(s[keep], 2), beta = round(L[[m]]$beta[keep], 3))
  out[order(-out$freq), ]
}
coef_bars <- function(d, title) {
  ggplot(d, aes(reorder(genus, beta), beta, fill = beta > 0)) + geom_col() +
    coord_flip(ylim = c(-0.3, 0.3)) + scale_fill_manual(values = red_blue, guide = "none") +
    labs(x = NULL, y = "Elastic net coefficient", title = title) + theme_minimal(base_size = 13)
}
plot_stable <- function(m) coef_bars(subset(stable(m, 0.9), abs(beta) >= 0.05), lab[m])
stable("IPA_B"); stable("Calpr_B")
plot_stable("IPA_B") + plot_stable("Calpr_B")
ggsave("Final/elastic_net_stable.png", width = 12, height = 5, dpi = 300, bg = "white")

# Final model on all 60 patients (lambda by leave-one-out), 12 largest coefficients
final_beta <- function(m) {
  Xen <- cbind(covs_main, scale(X)); pen <- c(rep(0, ncol(covs_main)), rep(1, ncol(X)))
  cv  <- cv.glmnet(Xen, Y[, m], alpha = 0.5, penalty.factor = pen, nfolds = nrow(Xen), grouped = FALSE)
  b   <- coef(cv, s = "lambda.min")[colnames(X), 1]; b[b != 0]
}
b_ipa <- final_beta("IPA_B")
b <- head(b_ipa[order(-abs(b_ipa))], 12)
coef_bars(data.frame(genus = gsub("_", " ", nm[names(b)]), beta = b), "Indole-3-propionate")
ggsave("Final/elastic_net_final.png", width = 8, height = 5, dpi = 300, bg = "white")

# Sensitivity: stage only as covariate (calprotectin, LBP)
key <- c("Calpr_B", "LBP_B")
set.seed(SEED)
loo_stage <- setNames(lapply(key, run_loocv, covs = covs_stage), key)
round(data.frame(row.names = lab[key],
                 base_main  = sapply(loo[key], `[[`, "R2_base"),       enet_main  = sapply(loo[key], `[[`, "R2_enet"),
                 base_stage = sapply(loo_stage, `[[`, "R2_base"),      enet_stage = sapply(loo_stage, `[[`, "R2_enet")), 2)
for (m in key) {
  a <- stable(m)$genus; s <- stable(m, L = loo_stage)$genus
  cat(lab[m], ": main", length(a), "| stage only", length(s), "| in both", length(intersect(a, s)), "\n")
}

# 6b. Sparse PLS: fixed sparsity + stability selection ----
# tune.spls() with n = 60 picked a different keepX for every seed/filter, so the heat map
# changed between runs. Now: keepX fixed a priori; each genus gets a selection frequency over
# B subsamples of 80% of patients (stability selection, Meinshausen & Buehlmann 2010);
# the heat map shows only genera selected in >= 70% of subsamples.
KEEP_X <- c(20, 20); KEEP_Y <- c(16, 16); B <- 200; THR <- 0.7

prep_spls <- function(Xm) {
  colnames(Xm) <- nm[colnames(Xm)]
  nzv <- nearZeroVar(Xm)$Position
  if (length(nzv) > 0) Xm[, -nzv] else Xm
}
spls_stability <- function(Xs) {                 # comp 1 and 2 pooled (they can swap order)
  set.seed(SEED)
  count <- setNames(numeric(ncol(Xs)), colnames(Xs))
  for (b in 1:B) {
    idx <- sample(nrow(Xs), round(0.8 * nrow(Xs)))
    Xb  <- Xs[idx, ]; Xb <- Xb[, apply(Xb, 2, sd) > 0]
    fit <- spls(Xb, Y_lab[idx, ], ncomp = 2, keepX = KEEP_X, keepY = KEEP_Y)
    sel <- unique(unlist(lapply(1:2, function(k) selectVar(fit, comp = k)$X$name)))
    count[sel] <- count[sel] + 1
  }
  sort(count / B, decreasing = TRUE)
}
spls_heatmap <- function(Xs, freq, file) {       # refit on the stable genera only, then cim()
  g   <- names(freq)[freq >= THR]
  fit <- spls(Xs[, g], Y_lab, ncomp = 2, keepX = rep(length(g), 2), keepY = KEEP_Y)
  png(file, width = 9, height = 9, units = "in", res = 300)
  cim(fit, comp = 1:2, margins = c(10, 15))
  dev.off()
}

Xs10 <- prep_spls(X); Xs30 <- prep_spls(X30)
stab10 <- spls_stability(Xs10)                   # ~1 min each
stab30 <- spls_stability(Xs30)
spls_heatmap(Xs10, stab10, "Final/spls_heatmap.png")
spls_heatmap(Xs30, stab30, "Final/spls_heatmap30.png")

stab_tab <- merge(data.frame(genus = names(stab10), freq10 = stab10),
                  data.frame(genus = names(stab30), freq30 = stab30), all = TRUE)  # NA = not in 30% set
stab_tab$stable_both <- with(stab_tab, freq10 >= THR & !is.na(freq30) & freq30 >= THR)
stab_tab <- stab_tab[order(-stab_tab$freq10), ]
head(stab_tab, 30)
write.csv(stab_tab, "Final/spls_stability.csv", row.names = FALSE)

top <- subset(stab_tab, pmax(freq10, freq30, na.rm = TRUE) >= 0.5)
pivot_longer(top, c(freq10, freq30), names_to = "filter", values_to = "freq") |>
  ggplot(aes(freq, reorder(genus, freq), colour = filter)) +
  geom_vline(xintercept = THR, linetype = 2, colour = "grey50") + geom_point(size = 2.5, na.rm = TRUE) +
  scale_colour_manual(values = c(freq10 = "grey30", freq30 = "#B2182B"),
                      labels = c(">= 10% prevalence", ">= 30% prevalence")) +
  labs(x = sprintf("sPLS selection frequency (%d subsamples, 80%% of patients)", B), y = NULL, colour = NULL) +
  theme_minimal(base_size = 12) + theme(legend.position = "top")
ggsave("Final/spls_stability.png", width = 8, height = 7, dpi = 300, bg = "white")

# 7. Sensitivity: stage only / genera in >= 30% ----
reg_stage <- run_regression(X, adjust_crea = FALSE); reg_30 <- run_regression(X30)
gt_stage  <- run_gt(X, adjust_crea = FALSE);         gt_30  <- run_gt(X30)
plot_heatmap(reg_30,    "Final/heatmap_regression_30.png")
plot_heatmap(reg_stage, "Final/heatmap_regression_stage.png")

key_pair <- function(r, g, m) {
  r <- r[r$genus_lab == g & r$outcome == m, ]
  if (nrow(r) == 0) "not tested" else sprintf("%.2f (q = %.3f)", r$beta, r$q)
}
settings <- list(`Main analysis` = list(reg, gt_main), `Stage only` = list(reg_stage, gt_stage),
                 `Genera in >= 30%` = list(reg_30, gt_30))
sens_table <- sapply(settings, function(s) {
  r <- s[[1]]; g <- s[[2]]
  c(`Coprococcus 2 - IPA: beta (q)`               = key_pair(r, "Coprococcus_2", "IPA_B"),
    `[R.] torques group - calprotectin: beta (q)` = key_pair(r, "[Ruminococcus]_torques_group", "Calpr_B"),
    `Rank of IPA (global test, out of 16)`        = which(g$outcome == "IPA_B"),
    `Rank of calprotectin (global test)`          = which(g$outcome == "Calpr_B"),
    `Significant genus-metabolite pairs`          = sum(r$q < 0.05),
    `Number of tests`                             = nrow(r))
})
sens_table
write.csv(sens_table, "Final/sensitivity_table.csv")

pair_beta <- function(sens, label)
  transform(merge(reg[, c("genus", "outcome", "beta")], sens[, c("genus", "outcome", "beta")],
                  by = c("genus", "outcome"), suffixes = c("_main", "_sens")), setting = label)
ggplot(rbind(pair_beta(reg_stage, "Adjusted for stage only"), pair_beta(reg_30, "Genera in >= 30% of patients")),
       aes(beta_main, beta_sens)) +
  geom_abline(linetype = 2, colour = "red") + geom_point(alpha = 0.4, size = 1) +
  facet_wrap(~ setting) + coord_equal() +
  labs(x = "Coefficient, main analysis", y = "Coefficient, sensitivity analysis") + theme_minimal(base_size = 12)
ggsave("Final/sensitivity_scatter.png", width = 9, height = 4.5, dpi = 300, bg = "white")

# New / lost FDR-significant pairs versus the main analysis
changes <- function(sens) {
  cols <- c("genus_lab", "outcome_lab", "beta", "p", "q")
  m <- merge(reg[, cols], sens[, cols], by = c("genus_lab", "outcome_lab"), suffixes = c("_main", "_sens"))
  m$status <- with(m, ifelse(q_main < 0.05 & q_sens < 0.05, "both",
                      ifelse(q_sens < 0.05, "new", ifelse(q_main < 0.05, "lost", NA))))
  m <- m[!is.na(m$status), ]
  m[order(m$status, m$p_sens), ]
}
changes(reg_stage); changes(reg_30)
data.frame(outcome = lab[outcomes],
           main = gt_main$q[match(outcomes, gt_main$outcome)],
           stage_only = gt_stage$q[match(outcomes, gt_stage$outcome)],
           prev_30 = gt_30$q[match(outcomes, gt_30$outcome)])

# 8. Diagnostics for key associations ----
check_assoc <- function(genus_label, m, file) {
  g <- names(nm)[nm == genus_label]; stopifnot(length(g) == 1)
  x <- scale(X[, g])[, 1]; y <- Y[, m]; n <- length(y)
  fit <- lm(y ~ x + stage + log_crea); cd <- cooks.distance(fit); infl <- which(cd > 4 / n)
  png(file, width = 11, height = 7, units = "in", res = 300)
  par(mfrow = c(2, 3), mar = c(4.2, 4.2, 3, 1), oma = c(0, 0, 2, 0))
  plot(fit, which = c(1, 2, 3, 5), sub.caption = "")
  plot(X[, g], y, pch = 19, col = ifelse(X[, g] > 0, "steelblue", "grey60"), main = "Data and fitted line",
       xlab = paste(genus_label, "(Hellinger)"), ylab = paste(lab[m], "(standardised log2)"))
  abline(lm(y ~ X[, g]), col = "red")
  plot(cd, type = "h", xlab = "Patient", ylab = "Cook's distance", main = "Influential patients")
  abline(h = 4 / n, lty = 2, col = "red")
  mtext(paste(genus_label, "-", lab[m]), outer = TRUE, line = 0, cex = 1.1, font = 2)
  dev.off()
  fit2 <- if (length(infl)) lm(y[-infl] ~ x[-infl] + stage[-infl] + log_crea[-infl]) else fit
  sp   <- cor.test(X[, g], y, method = "spearman", exact = FALSE)
  list(patients_with_genus = sum(X[, g] > 0), influential = infl,
       all_patients = coef(summary(fit))[2, c(1, 4)], without_influential = coef(summary(fit2))[2, c(1, 4)],
       spearman = c(rho = unname(sp$estimate), p = sp$p.value))
}
graphics.off()
check_assoc("Coprococcus_2", "IPA_B", "diagnostics_Coprococcus2_IPA.png")
check_assoc("[Ruminococcus]_torques_group", "Calpr_B", "diagnostics_torques_calprotectin.png")

# Acetate and LBP: ten strongest genera with 95% CI (red = q < 0.05)
cor.test(log2(meta$Ace_B), log2(meta$IPA_B),   method = "spearman", exact = FALSE)
cor.test(log2(meta$LBP_B), log2(meta$Calpr_B), method = "spearman", exact = FALSE)
forest <- function(m, top = 10) {
  ggplot(head(reg[reg$outcome == m, ], top), aes(beta, reorder(genus_lab, beta))) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.2, orientation = "y") +
    geom_point(aes(colour = q < 0.05), size = 2.5) +
    scale_colour_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "grey40"), guide = "none") +
    labs(x = "Standardised coefficient (95% CI)", y = NULL, title = lab[m]) + theme_minimal(base_size = 12)
}
forest("Ace_B") + forest("LBP_B")
ggsave("Final/acetate_lbp_genera.png", width = 12, height = 5, dpi = 300, bg = "white")

# Regression vs elastic net: top-k genera of each method side by side
compare_methods <- function(m, k = 5) {
  rr <- reg[reg$outcome == m, ]; s <- loo[[m]]$stab; b <- loo[[m]]$beta
  top_r <- head(rr$genus, k)
  top_e <- names(s)[order(-s, -abs(b))][1:k]; top_e <- top_e[s[top_e] > 0]
  g <- union(top_r, top_e); rr <- rr[match(g, rr$genus), ]
  data.frame(outcome = lab[m], genus = nm[g],
             top_in = ifelse(g %in% top_r & g %in% top_e, "both", ifelse(g %in% top_r, "regression", "elastic net")),
             reg_beta = round(rr$beta, 2), reg_q = signif(rr$q, 2),
             enet_freq = round(s[g], 2), enet_beta = round(b[g], 3),
             same_sign = ifelse(b[g] == 0, NA, sign(rr$beta) == sign(b[g])), row.names = NULL)
}
cmp_key <- do.call(rbind, lapply(c("IPA_B", "Calpr_B", "Ace_B", "LBP_B"), compare_methods))
write.csv(cmp_key, "comparison_regression_elasticnet.csv", row.names = FALSE)
cmp_summary <- do.call(rbind, lapply(outcomes, function(m) {
  x <- compare_methods(m)
  data.frame(outcome = lab[m], sig_regression = sum(reg$outcome == m & reg$q < 0.05),
             top5_in_both = sum(x$top_in == "both"),
             same_sign_in_both = sum(x$same_sign[x$top_in == "both"], na.rm = TRUE),
             enet_R2_gain = round(res_enet[m, "gain"], 2))
}))
cmp_summary[order(-cmp_summary$top5_in_both, -cmp_summary$enet_R2_gain), ]
