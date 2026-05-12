library(vars); library(svars); library(fredr); library(zoo)
library(ggplot2); library(dplyr); library(scales)

fredr_set_key("be238a2a428573d02dc7925d8742c9df")

start_date <- as.Date("1999-01-01")
end_date   <- as.Date("2019-12-01")

pull_fred <- function(id) {
  fredr(series_id = id, observation_start = start_date,
        observation_end = end_date, frequency = "m")$value
}

hicp   <- pull_fred("ITACPIALLMINMEI")
ip     <- pull_fred("ITAPROINDMISMEI")
mro    <- pull_fred("ECBMRRFR")
btp10  <- pull_fred("IRLTLT01ITM156N")
bund10 <- pull_fred("IRLTLT01DEM156N")

credit_q <- fredr(series_id = "CRDQITBPABIS",
                  observation_start = as.Date("1998-10-01"),
                  observation_end   = end_date, frequency = "q")
dates_m  <- seq(start_date, end_date, by = "month")
credit   <- approx(as.numeric(credit_q$date), credit_q$value,
                   xout = as.numeric(dates_m), method = "linear")$y

min_len <- min(lengths(list(hicp, ip, mro, btp10, bund10, credit, dates_m)))
hicp   <- hicp[1:min_len];   ip     <- ip[1:min_len]
btp10  <- btp10[1:min_len];  bund10 <- bund10[1:min_len]
credit <- credit[1:min_len]

mro_filled <- na.locf(mro[1:min_len], na.rm = FALSE)
mro_filled <- na.locf(mro_filled, fromLast = TRUE)

Y5 <- na.trim(ts(data.frame(
  HICP   = 100 * log(hicp),
  IP     = 100 * log(ip),
  MRO    = mro_filled,
  SPREAD = btp10 - bund10,
  CREDIT = 100 * log(credit)
), start = c(1999, 1), frequency = 12))

Y3 <- Y5[, c("HICP", "IP", "MRO")]
Y4 <- Y5[, c("HICP", "IP", "MRO", "SPREAD")]

p_lags <- 6
SB_cal <- as.Date("2012-07-01")
SB_idx <- as.integer(difftime(SB_cal, start_date, units = "days") / 30.44) + 1
SB_res <- SB_idx - p_lags

var3 <- VAR(Y3, p = p_lags, type = "const")
var4 <- VAR(Y4, p = p_lags, type = "const")
var5 <- VAR(Y5, p = p_lags, type = "const")

cv3 <- id.cv(var3, SB = SB_res)
cv4 <- id.cv(var4, SB = SB_res)
cv5 <- id.cv(var5, SB = SB_res)

cat("AIC — M1:", AIC(var3), "  M2:", AIC(var4), "  M3:", AIC(var5), "\n")


extract_cv <- function(cv_obj, p) {
  U  <- residuals(cv_obj$VAR)
  T_ <- nrow(U)
  sb <- cv_obj$SB
  S1 <- t(U[1:(sb-1), ]) %*% U[1:(sb-1), ] / (sb - 1)
  S2 <- t(U[sb:T_,    ]) %*% U[sb:T_,    ] / (T_ - sb + 1)
  C  <- t(chol(S1))
  W  <- solve(C) %*% S2 %*% t(solve(C))
  eg <- eigen(W, symmetric = TRUE)
  ord <- order(eg$values, decreasing = TRUE)
  mu  <- eg$values[ord]
  Q   <- eg$vectors[, ord]
  list(U = U, S1 = S1, S2 = S2, C = C, W = W, Q = Q, B = C %*% Q, mu = mu, sb = sb)
}

spectral_gap <- function(mu) {
  n <- length(mu)
  sapply(seq_len(n), function(j) {
    min(if (j > 1) mu[j-1] - mu[j] else Inf,
        if (j < n) mu[j] - mu[j+1] else Inf)
  })
}

op_norm <- function(M) max(svd(M)$d)

align_B <- function(B_ref, B_est) {
  n <- ncol(B_ref)
  out  <- B_est
  used <- integer(0)
  perm <- integer(n)
  for (j in seq_len(n)) {
    corrs <- sapply(seq_len(n), function(k) {
      if (k %in% used) return(-Inf)
      abs(sum(B_ref[, j] * B_est[, k])) /
        (sqrt(sum(B_ref[, j]^2)) * sqrt(sum(B_est[, k]^2)) + 1e-12)
    })
    best <- which.max(corrs)
    perm[j] <- best
    used <- c(used, best)
    if (sum(B_ref[, j] * B_est[, best]) < 0) out[, best] <- -B_est[, best]
  }
  out[, perm]
}

align_B_full <- function(B_ref_3x3, B_full) {
  n_ref <- ncol(B_ref_3x3)
  K     <- nrow(B_full)
  out   <- matrix(0, K, n_ref)
  used  <- integer(0)
  for (j in seq_len(n_ref)) {
    scores <- sapply(seq_len(ncol(B_full)), function(k) {
      if (k %in% used) return(-Inf)
      b_ref <- B_ref_3x3[, j]
      b_can <- B_full[1:n_ref, k]
      abs(sum(b_ref * b_can)) /
        (sqrt(sum(b_ref^2)) * sqrt(sum(b_can^2)) + 1e-12)
    })
    best <- which.max(scores)
    used <- c(used, best)
    col  <- B_full[, best]
    if (sum(B_ref_3x3[, j] * col[1:n_ref]) < 0) col <- -col
    out[, j] <- col
  }
  out
}


dk_diagnostics <- function(label, cv_small, cv_large, n_shared = NULL) {
  e_s <- extract_cv(cv_small, p_lags)
  e_l <- extract_cv(cv_large, p_lags)
  K_s <- ncol(e_s$U); K_l <- ncol(e_l$U)
  if (is.null(n_shared)) n_shared <- K_s
  
  cv_sub <- function(e, n) {
    S1 <- e$S1[1:n, 1:n]; S2 <- e$S2[1:n, 1:n]
    C  <- t(chol(S1))
    W  <- solve(C) %*% S2 %*% t(solve(C))
    eg <- eigen(W, symmetric = TRUE)
    ord <- order(eg$values, decreasing = TRUE)
    mu  <- eg$values[ord]; Q <- eg$vectors[, ord]
    list(C = C, W = W, Q = Q, mu = mu, B = C %*% Q)
  }
  l <- cv_sub(e_l, n_shared)
  s <- cv_sub(e_s, n_shared)
  
  B_s     <- align_B(l$B, s$B)
  Q_s_aln <- align_B(l$Q, s$Q)
  E_W     <- s$W - l$W
  E_W_op  <- op_norm(E_W)
  delta   <- spectral_gap(l$mu)
  dk_bnd  <- ifelse(delta > 0, 2 * E_W_op / delta, Inf)
  
  angles_deg <- sapply(seq_len(n_shared), function(j) {
    u <- l$B[, j] / sqrt(sum(l$B[, j]^2))
    v <- B_s[, j] / sqrt(sum(B_s[, j]^2))
    acos(min(abs(sum(u * v)), 1)) * 180 / pi
  })
  
  dC <- s$C - l$C; dQ <- Q_s_aln - l$Q
  T1 <- dC %*% l$Q; T2 <- l$C %*% dQ; T3 <- dC %*% dQ
  
  cat("\n", strrep("=", 55), "\n", label, "\n", strrep("=", 55), "\n")
  cat(sprintf("  ||E_W||_op : %.5f\n", E_W_op))
  cat(sprintf("  mu_j       : %s\n", paste(round(l$mu, 4), collapse = ", ")))
  cat(sprintf("  delta_j    : %s\n", paste(round(delta, 5), collapse = ", ")))
  cat(sprintf("  DK bound   : %s\n", paste(round(dk_bnd, 4), collapse = ", ")))
  cat(sprintf("  angle (deg): %s\n", paste(round(angles_deg, 2), collapse = ", ")))
  cat(sprintf("  ||dB||     : %.5f\n", op_norm(T1 + T2 + T3)))
  cat(sprintf("  ||T1||     : %.5f  ||T2||: %.5f  ||T3||: %.5f\n",
              op_norm(T1), op_norm(T2), op_norm(T3)))
  
  invisible(list(
    E_W = E_W, E_W_op = E_W_op, delta = delta, dk_bnd = dk_bnd,
    angles = angles_deg, B_hat = B_s, B_tilde = l$B,
    T1 = T1, T2 = T2, T3 = T3,
    mu_large = l$mu, mu_small = s$mu,
    W_small = s$W, W_large = l$W
  ))
}

res_3v4 <- dk_diagnostics("M1 vs M2 (omitted: Spread)",         cv3, cv4, n_shared = 3)
res_3v5 <- dk_diagnostics("M1 vs M3 (omitted: Spread, Credit)", cv3, cv5, n_shared = 3)
res_4v5 <- dk_diagnostics("M2 vs M3 (omitted: Credit)",         cv4, cv5, n_shared = 3)


get_whitened_angles <- function(res) {
  eig_sorted <- function(W) {
    eg <- eigen(W, symmetric = TRUE)
    eg$vectors[, order(eg$values, decreasing = TRUE)]
  }
  Q_l <- eig_sorted(res$W_large)
  Q_s <- eig_sorted(res$W_small)
  Q_s_aln <- align_B(Q_l, Q_s)
  sapply(seq_len(ncol(Q_l)), function(j) {
    acos(min(abs(sum(Q_l[, j] * Q_s_aln[, j])), 1)) * 180 / pi
  })
}

aw_3v4 <- get_whitened_angles(res_3v4)
aw_3v5 <- get_whitened_angles(res_3v5)
aw_4v5 <- get_whitened_angles(res_4v5)

cat("--- Chow test ---\n")
chow5 <- chow.test(var5, SB = SB_res, nboot = 500)
print(summary(chow5))


e3 <- extract_cv(cv3, p_lags)
e4 <- extract_cv(cv4, p_lags)
e5 <- extract_cv(cv5, p_lags)

Phi3 <- Phi(var3, nstep = 48)
Phi4 <- Phi(var4, nstep = 48)
Phi5 <- Phi(var5, nstep = 48)

B_M1_full <- e3$B
B_M2_full <- align_B_full(e3$B, e4$B)
B_M3_full <- align_B_full(e3$B, e5$B)

for (j in 1:3) {
  b1 <- e3$B[, j] / sqrt(sum(e3$B[, j]^2))
  b2 <- B_M2_full[1:3, j]; b2 <- b2 / sqrt(sum(b2^2))
  b3 <- B_M3_full[1:3, j]; b3 <- b3 / sqrt(sum(b3^2))
  cat(sprintf("  sh%d: M1-M2 cos=%.4f   M1-M3 cos=%.4f\n", j,
              abs(sum(b1 * b2)), abs(sum(b1 * b3))))
}

make_irf_df <- function(Phi_arr, B_full, vnames_resp, slabels, mlabel, n_ahead = 48) {
  nh <- n_ahead + 1L
  stopifnot(dim(Phi_arr)[2] == nrow(B_full))
  out <- vector("list", length(vnames_resp) * ncol(B_full))
  r <- 0L
  for (j in seq_len(ncol(B_full))) {
    bj <- B_full[, j]
    for (i in seq_along(vnames_resp)) {
      r <- r + 1L
      out[[r]] <- data.frame(
        horizon  = 0L:n_ahead,
        irf      = vapply(seq_len(nh), function(h) sum(Phi_arr[i, , h] * bj), numeric(1L)),
        shock    = slabels[j],
        response = vnames_resp[i],
        model    = mlabel,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, out)
}

vnames <- c("HICP", "IP", "MRO")

slabels_base <- sprintf("Sh%d  \u03bc=%.2f", 1:3, e3$mu[1:3])
slabels_comp <- sprintf("Sh%d  (\u03bc=%.2f)", 1:3, e3$mu[1:3])

resp_labels <- c(HICP = "CPI", IP = "IP", MRO = "MRO")

theme_publication <- function(base_size = 11) {
  theme_classic(base_size = base_size, base_family = "serif") +
    theme(
      panel.border      = element_rect(colour = "black", fill = NA, linewidth = 0.4),
      panel.grid        = element_blank(),
      panel.background  = element_rect(fill = "white"),
      strip.background  = element_rect(fill = "white", colour = "black", linewidth = 0.4),
      strip.text        = element_text(size = base_size * 0.85),
      axis.line         = element_blank(),
      axis.ticks        = element_line(linewidth = 0.3),
      axis.ticks.length = unit(2, "pt"),
      axis.text         = element_text(size = base_size * 0.8),
      axis.title        = element_text(size = base_size * 0.9),
      legend.background = element_rect(fill = "white", colour = "black", linewidth = 0.3),
      legend.key        = element_rect(fill = "white"),
      legend.key.width  = unit(1.8, "cm"),
      legend.text       = element_text(size = base_size * 0.85),
      legend.title      = element_blank(),
      legend.position   = "bottom",
      plot.title        = element_text(size = base_size, face = "plain", hjust = 0),
      plot.subtitle     = element_text(size = base_size * 0.85, colour = "grey30", hjust = 0),
      plot.margin       = margin(8, 8, 6, 6)
    )
}

grid_theme <- theme(
  panel.grid.major = element_line(linetype = "dashed", linewidth = 0.3, colour = "grey70"),
  panel.grid.minor = element_line(linetype = "dotted", linewidth = 0.2, colour = "grey85")
)

lty_bw <- c("M1: 3-var" = "solid",    "M2: 4-var" = "longdash", "M3: 5-var" = "dotted")
lwd_bw <- c("M1: 3-var" = 0.9,        "M2: 4-var" = 0.7,        "M3: 5-var" = 0.6)

df_base <- make_irf_df(Phi3, B_M1_full, vnames, slabels_base, "M1")

p_base <- ggplot(df_base, aes(x = horizon, y = irf)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "black") +
  geom_line(linewidth = 0.85, colour = "black") +
  facet_grid(response ~ shock, scales = "free_y",
             labeller = labeller(response = resp_labels)) +
  labs(title    = "Impulse responses \u2014 Baseline model M1 (CPI, IP, MRO)",
       subtitle = "CV-SVAR, Italy 1999\u20132019. Break: July 2012 (Draghi).",
       x = "Months after shock", y = "Response") +
  theme_publication(base_size = 11) + grid_theme

ggsave("irf_M1_baseline.pdf", p_base, width = 16, height = 10,
       units = "cm", device = cairo_pdf)

df_comp <- rbind(
  make_irf_df(Phi3, B_M1_full, vnames, slabels_comp, "M1: 3-var"),
  make_irf_df(Phi4, B_M2_full, vnames, slabels_comp, "M2: 4-var"),
  make_irf_df(Phi5, B_M3_full, vnames, slabels_comp, "M3: 5-var")
)
df_comp$model <- factor(df_comp$model, levels = c("M1: 3-var","M2: 4-var","M3: 5-var"))

p_comp <- ggplot(df_comp, aes(x = horizon, y = irf, linetype = model, linewidth = model)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "black") +
  geom_line(colour = "black") +
  facet_grid(response ~ shock, scales = "free_y",
             labeller = labeller(response = resp_labels)) +
  scale_linetype_manual(values = lty_bw) +
  scale_linewidth_manual(values = lwd_bw) +
  guides(linetype = guide_legend(nrow = 1, override.aes = list(linewidth = c(0.9, 0.7, 0.6))),
         linewidth = "none") +
  labs(title    = "Impulse responses \u2014 M1 / M2 / M3 comparison",
       subtitle = "Shocks aligned to M1 baseline (full-dimensional IRF).",
       x = "Months after shock", y = "Response") +
  theme_publication(base_size = 11) + grid_theme

ggsave("irf_confronto_M1M2M3.pdf", p_comp, width = 20, height = 12,
       units = "cm", device = cairo_pdf)


N_BOOT  <- 5000
N_AHEAD <- 48
CONF    <- 0.68

bootstrap_irf <- function(var_obj, SB_res, p_lags, B_ref_3x3, align_fn,
                          n_resp = 3, n_sh = 3,
                          n_boot = N_BOOT, n_ahead = N_AHEAD, conf = CONF) {
  U_hat  <- residuals(var_obj)
  Y_orig <- as.matrix(var_obj$y)
  K      <- ncol(U_hat)
  T_res  <- nrow(U_hat)
  T_tot  <- nrow(Y_orig)
  n1     <- SB_res - 1
  n2     <- T_res - SB_res + 1
  alpha  <- (1 - conf) / 2
  A_list <- Acoef(var_obj)
  nu     <- sapply(coef(var_obj), function(eq) eq["const", 1])
  irf_acc <- array(NA_real_, dim = c(n_resp, n_sh, n_ahead + 1L, n_boot))
  n_ok <- 0L
  for (b in seq_len(n_boot)) {
    idx1   <- sample(1:n1,         n1, replace = TRUE)
    idx2   <- sample(SB_res:T_res, n2, replace = TRUE)
    U_star <- U_hat[c(idx1, idx2), ]
    Y_star <- matrix(0, T_tot, K)
    Y_star[1:p_lags, ] <- Y_orig[1:p_lags, ]
    for (t in (p_lags + 1):T_tot) {
      y_t <- nu
      for (lag in seq_len(p_lags))
        y_t <- y_t + A_list[[lag]] %*% Y_star[t - lag, ]
      Y_star[t, ] <- y_t + U_star[t - p_lags, ]
    }
    Y_star <- ts(Y_star, start = start(var_obj$y), frequency = frequency(var_obj$y))
    var_b <- tryCatch(VAR(Y_star, p = p_lags, type = "const"), error = function(e) NULL)
    if (is.null(var_b)) next
    cv_b  <- tryCatch(id.cv(var_b, SB = SB_res), error = function(e) NULL)
    if (is.null(cv_b)) next
    e_b  <- extract_cv(cv_b, p_lags)
    B_b  <- align_fn(B_ref_3x3, e_b$B)
    Phi_b <- Phi(var_b, nstep = n_ahead)
    for (j in seq_len(n_sh)) {
      bj <- B_b[, j]
      for (i in seq_len(n_resp))
        irf_acc[i, j, , b] <- vapply(seq_len(n_ahead + 1L),
                                     function(h) sum(Phi_b[i, , h] * bj), numeric(1L))
    }
    n_ok <- n_ok + 1L
  }
  cat(sprintf("  Valid replications: %d / %d\n", n_ok, n_boot))
  list(
    lo = apply(irf_acc, 1:3, quantile, probs = alpha,     na.rm = TRUE),
    hi = apply(irf_acc, 1:3, quantile, probs = 1 - alpha, na.rm = TRUE)
  )
}

ci_to_df <- function(ci_obj, vnames, slabels, mlabel, n_ahead = N_AHEAD) {
  rows_lo <- rows_hi <- vector("list", length(vnames) * length(slabels))
  r <- 0L
  for (j in seq_along(slabels)) {
    for (i in seq_along(vnames)) {
      r <- r + 1L
      base <- data.frame(horizon = 0:n_ahead, shock = slabels[j],
                         response = vnames[i], model = mlabel,
                         stringsAsFactors = FALSE)
      rows_lo[[r]] <- cbind(base, irf = ci_obj$lo[i, j, ])
      rows_hi[[r]] <- cbind(base, irf = ci_obj$hi[i, j, ])
    }
  }
  list(lo = do.call(rbind, rows_lo), hi = do.call(rbind, rows_hi))
}

cat("\n[1/3] Bootstrap M1...\n")
ci_M1 <- bootstrap_irf(var3, SB_res, p_lags, B_ref_3x3 = e3$B, align_fn = align_B)
cat("[2/3] Bootstrap M2...\n")
ci_M2 <- bootstrap_irf(var4, SB_res, p_lags, B_ref_3x3 = e3$B, align_fn = align_B_full)
cat("[3/3] Bootstrap M3...\n")
ci_M3 <- bootstrap_irf(var5, SB_res, p_lags, B_ref_3x3 = e3$B, align_fn = align_B_full)

saveRDS(list(ci_M1 = ci_M1, ci_M2 = ci_M2, ci_M3 = ci_M3), "bootstrap_ci.rds")

df_ci_M1_base <- ci_to_df(ci_M1, vnames, slabels_base, "M1")

p_base_ci <- ggplot(df_base, aes(x = horizon, y = irf)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "black") +
  geom_line(data = df_ci_M1_base$lo, aes(x = horizon, y = irf),
            linetype = "dashed", linewidth = 0.45, colour = "black") +
  geom_line(data = df_ci_M1_base$hi, aes(x = horizon, y = irf),
            linetype = "dashed", linewidth = 0.45, colour = "black") +
  geom_line(linewidth = 0.85, colour = "black") +
  facet_grid(response ~ shock, scales = "free_y",
             labeller = labeller(response = resp_labels)) +
  labs(title    = "Impulse responses \u2014 Baseline model M1 (CPI, IP, MRO)",
       subtitle = "CV-SVAR, Italy 1999\u20132019. Break: July 2012. Dashed: 68% recursive bootstrap CI.",
       x = "Months after shock", y = "Response") +
  theme_publication(base_size = 11) + grid_theme

ggsave("irf_M1_baseline_ci.pdf", p_base_ci, width = 16, height = 10,
       units = "cm", device = cairo_pdf)

df_ci_M3_comp <- ci_to_df(ci_M3, vnames, slabels_comp, "M3: 5-var")
df_ribbon_M3  <- merge(
  df_ci_M3_comp$lo[, c("horizon","shock","response","irf")],
  df_ci_M3_comp$hi[, c("horizon","shock","response","irf")],
  by = c("horizon","shock","response"), suffixes = c("_lo","_hi")
)

p_comp_ci <- ggplot(df_comp, aes(x = horizon, y = irf, linetype = model, linewidth = model)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "black") +
  geom_ribbon(data = df_ribbon_M3,
              aes(x = horizon, ymin = irf_lo, ymax = irf_hi),
              inherit.aes = FALSE, fill = "grey80", alpha = 0.55) +
  geom_line(colour = "black") +
  facet_grid(response ~ shock, scales = "free_y",
             labeller = labeller(response = resp_labels)) +
  scale_linetype_manual(values = lty_bw) +
  scale_linewidth_manual(values = lwd_bw) +
  guides(linetype = guide_legend(nrow = 1, override.aes = list(linewidth = c(0.9, 0.7, 0.6))),
         linewidth = "none") +
  labs(title    = "Impulse responses \u2014 M1 / M2 / M3 comparison",
       subtitle = "Shocks aligned to M1 baseline. Shaded: 68% recursive bootstrap CI of benchmark M3.",
       x = "Months after shock", y = "Response") +
  theme_publication(base_size = 11) + grid_theme

ggsave("irf_confronto_M1M2M3_ci.pdf", p_comp_ci, width = 20, height = 12,
       units = "cm", device = cairo_pdf)

for (sh in unique(df_comp$shock)) {
  df_sh     <- df_comp[df_comp$shock == sh, ]
  df_rib_sh <- df_ribbon_M3[df_ribbon_M3$shock == sh, ]
  
  p_sh <- ggplot(df_sh, aes(x = horizon, y = irf, linetype = model, linewidth = model)) +
    geom_hline(yintercept = 0, linewidth = 0.25, colour = "black") +
    geom_ribbon(data = df_rib_sh,
                aes(x = horizon, ymin = irf_lo, ymax = irf_hi),
                inherit.aes = FALSE, fill = "grey80", alpha = 0.55) +
    geom_line(colour = "black") +
    facet_wrap(~response, scales = "free_y", ncol = 3,
               labeller = labeller(response = resp_labels)) +
    scale_linetype_manual(values = lty_bw) +
    scale_linewidth_manual(values = lwd_bw) +
    scale_y_continuous(breaks = pretty_breaks(n = 3)) +
    scale_x_continuous(breaks = c(0, 12, 24, 36, 48)) +
    guides(linetype = guide_legend(nrow = 1, override.aes = list(linewidth = c(0.9, 0.7, 0.6))),
           linewidth = "none") +
    labs(title    = paste0("Impulse responses to ", sh),
         subtitle = "CV-SVAR, Italy 1999\u20132019. Shaded: 68% recursive bootstrap CI of benchmark M3.",
         x = "Months after shock", y = "Response") +
    theme_publication(base_size = 10) +
    theme(
      panel.grid.major  = element_line(linetype = "dotted", linewidth = 0.2, colour = "grey85"),
      panel.grid.minor  = element_line(linetype = "dotted", linewidth = 0.2, colour = "grey85"),
      legend.key.width  = unit(1.2, "cm"),
      legend.key.height = unit(0.35, "cm"),
      legend.margin     = margin(3, 6, 3, 6),
      legend.box.margin = margin(0, 0, 0, 0)
    )
  
  ggsave(paste0("irf_bw_ci_", gsub("[^A-Za-z0-9]", "_", sh), ".pdf"),
         p_sh, width = 16, height = 9, units = "cm", device = cairo_pdf)
}