library(vars); library(svars); library(clue)
library(dplyr); library(tidyr); library(ggplot2)
library(doParallel); library(foreach); library(scales)

#############################################
# DGP
#############################################

A1 <- matrix(c(
  0.40, 0.08,-0.05, 0.10,
  0.05, 0.45, 0.08, 0.07,
  0.00, 0.15, 0.50, 0.05,
  0.10, 0.05, 0.05, 0.50
), 4, 4, byrow=TRUE)

A2 <- matrix(c(
  0.05, 0.00, 0.00, 0.02,
  0.00, 0.05, 0.00, 0.02,
  0.00, 0.00, 0.05, 0.01,
  0.02, 0.01, 0.01, 0.10
), 4, 4, byrow=TRUE)

B_true_4x4 <- matrix(c(
  1.3,  1.0, -0.7,  0.8,
  0.9, -0.8, -0.6,  0.7,
  0.7,  0.6,  1.3,  0.6,
  -0.6,  0.7, -0.5,  1.5
), 4, 4, byrow=TRUE)

B_true_3x3 <- B_true_4x4[1:3, 1:3]
LAMBDA_FIXED <- 0.8

LAMBDA4_GRID <- c(2.0, 4.0, 8.0)
LAMBDA4_LABS <- c("Low contamination (lambda4=2)",
                  "Med contamination (lambda4=4)",
                  "High contamination (lambda4=8)")
LAMBDA4_IRR  <- 1.0 #VARIANCE OF Y4 DOES NOT CHANGE
lev_B <- c(LAMBDA4_LABS, "Irrelevant omission")

make_scenario <- function(scale_col4=1.0, scale_lag4=1.0) {
  B   <- B_true_4x4; B[,4]      <- B[,4]      * scale_col4
  A1s <- A1;         A1s[1:3,4] <- A1[1:3,4]  * scale_lag4
  A2s <- A2;         A2s[1:3,4] <- A2[1:3,4]  * scale_lag4
  list(B=B, A1=A1s, A2=A2s)
}

sc_relevant   <- make_scenario(1.0, 1.0)
sc_irrelevant <- make_scenario(0.1, 0.1)
burn <- 200
Nrep <- 500

#############################################
#  B* — THEORETICAL PLATEAU (Proposition 1)
#############################################

compute_Bstar_frob <- function(sc, lambda4) {
  # 1. CONTAMINATED POPULATION COVARIANCE MATRICES (contemporaneous channel only)
  Bt    <- B_true_3x3
  c_vec <- sc$B[1:3, 4]
  lam   <- LAMBDA_FIXED
  P1    <- 1.0                  * outer(c_vec, c_vec)
  P2    <- (1 + lam*(lambda4-1)) * outer(c_vec, c_vec)
  Lam2  <- (1-lam)*diag(3) + lam*diag(c(2,3,4))
  S1    <- Bt %*% t(Bt) + P1
  S2    <- Bt %*% Lam2 %*% t(Bt) + P2
  # 2. CV IDENTIFICATION APPLIED TO CONTAMINATED COVARIANCES

  C_    <- t(chol(S1))
  W_    <- solve(C_) %*% S2 %*% t(solve(C_)); W_ <- (W_+t(W_))/2
  Bs    <- C_ %*% eigen(W_, symmetric=TRUE)$vectors
  # 3. OPTIMAL ALIGNMENT (permutation and sign)
  nc    <- function(x) x/sqrt(sum(x^2))
  S_    <- abs(t(apply(Bt,2,nc)) %*% apply(Bs,2,nc))
  perm  <- as.integer(clue::solve_LSAP(1-S_))
  Bp    <- Bs[,perm]
  for (j in 1:3) if (sum(Bp[,j]*Bt[,j])<0) Bp[,j] <- -Bp[,j]
  # 4. FROBENIUS DISTANCE 
  sqrt(sum((Bp - Bt)^2))
}

plateau_df <- data.frame(
  lambda4     = LAMBDA4_GRID,
  lambda4_lab = factor(LAMBDA4_LABS, levels=LAMBDA4_LABS),
  plateau     = sapply(LAMBDA4_GRID, function(l4)
    compute_Bstar_frob(sc_relevant, l4))
)
cat("=== Theoretical plateaus ===\n"); print(plateau_df)

#############################################
#  TRUE EIGENVECTORS 
#############################################

compute_Q_tilde <- function(lam_obs=LAMBDA_FIXED) {
  Lam2 <- (1-lam_obs)*diag(3) + lam_obs*diag(c(2,3,4))
  C_   <- t(chol(B_true_3x3 %*% t(B_true_3x3)))
  W_   <- solve(C_) %*% (B_true_3x3 %*% Lam2 %*% t(B_true_3x3)) %*% t(solve(C_))
  eigen((W_+t(W_))/2, symmetric=TRUE)$vectors
}
Q_TILDE <- compute_Q_tilde()

#############################################
#  j* — MOST-AFFECTED EIGENVECTOR
#############################################
{
  c_vec_pert <- sc_relevant$B[1:3, 4]
  lam        <- LAMBDA_FIXED
  P1_pert    <- 1.0                             * outer(c_vec_pert, c_vec_pert)  # Var(w4,r1)*cc'
  P2_pert    <- (1 + lam*(max(LAMBDA4_GRID)-1)) * outer(c_vec_pert, c_vec_pert)  # Var(w4,r2)*cc'
  Lam2_pert  <- (1-lam)*diag(3) + lam*diag(c(2,3,4))
  S1_pert    <- B_true_3x3 %*% t(B_true_3x3) + P1_pert                    # contaminated Sigma^(1)
  S2_pert    <- B_true_3x3 %*% Lam2_pert %*% t(B_true_3x3) + P2_pert     # contaminated Sigma^(2)
  C_pert     <- t(chol(S1_pert))
  W_pert     <- solve(C_pert) %*% S2_pert %*% t(solve(C_pert))
  W_pert     <- (W_pert + t(W_pert))/2
  Q_star     <- eigen(W_pert, symmetric=TRUE)$vectors   # eigenvectors of contaminated W (population)
  JSTAR_GLOBAL <- which.min(sapply(1:3, function(j)     # most misaligned = most rotated
    max(abs(t(Q_star) %*% Q_TILDE[,j]))))
}
cat(sprintf("=== j* (exact, population rotation) = %d ===\n", JSTAR_GLOBAL))
#############################################
#  HELPERS (alingment and frob error)
#############################################

align_B <- function(B_true, B_hat) {
  nc   <- function(x) x/sqrt(sum(x^2))
  S    <- abs(t(apply(B_true,2,nc)) %*% apply(B_hat,2,nc))
  perm <- clue::solve_LSAP(1-S)
  Bp   <- B_hat[,perm]
  for (j in 1:ncol(Bp)) if (sum(Bp[,j]*B_true[,j])<0) Bp[,j] <- -Bp[,j]
  Bp
}
frob_err <- function(A,B) sqrt(sum((A-B)^2))
?clue::solve_LSAP()
#############################################
#  DGP FUNCTIONS
#############################################

generate_shocks <- function(lambda4, T_tot, n=4L) {
  e    <- matrix(rnorm(n*T_tot), n, T_tot)
  Sig2 <- (1-LAMBDA_FIXED)*diag(n) + LAMBDA_FIXED*diag(c(2,3,4,lambda4))
  half <- floor((T_tot-burn)/2) + burn
  e[,(half+1):T_tot] <- chol(Sig2) %*% e[,(half+1):T_tot]
  e
}

simulate_VAR <- function(eps, T_tot, sc) {
  y <- matrix(0, 4L, T_tot)
  y[,1] <- rnorm(4); y[,2] <- rnorm(4)
  for (t in 3:T_tot)
    y[,t] <- sc$A1%*%y[,t-1] + sc$A2%*%y[,t-2] + sc$B%*%eps[,t]
  t(y[1:3, (burn+1):T_tot])
}

#############################################
#  ONE REP CV
#############################################

one_rep_cv <- function(y3, Q_tilde, jstar_global) {
  out <- list(frob=NA_real_, sintheta=NA_real_, sintheta_jstar=NA_real_, ok=FALSE)
  
  #--- 1. VAR ESTIMATION AND CV IDENTIFICATION ---
  fit <- try(VAR(y3, p=2, type="none"), silent=TRUE)
  if (inherits(fit,"try-error") || !all(Mod(roots(fit))<1)) return(out)  
  sv  <- try(id.cv(fit, SB=floor(nrow(y3)/2)), silent=TRUE)              
  if (inherits(sv,"try-error")) return(out)
  
  #--- 2. FROBENIUS ERROR ON B ---
  out$ok   <- TRUE
  out$frob <- frob_err(align_B(B_true_3x3, sv$B), B_true_3x3)     
  
  #--- 3. ANGULAR BIAS IN WHITENED SPACE ---
  U  <- residuals(fit); T_ <- nrow(U); T2 <- floor(T_/2)
  Ch <- try(t(chol(crossprod(U[1:T2,])/T2)), silent=TRUE)               
  if (!inherits(Ch,"try-error")) {
    Wh  <- solve(Ch) %*% (crossprod(U[(T2+1):T_,])/(T_-T2)) %*% t(solve(Ch)) 
    Qh  <- eigen((Wh+t(Wh))/2, symmetric=TRUE)$vectors                          
    
    out$sintheta <- mean(sapply(1:3, function(j) {                        # mean sin(theta) over all j
      best <- max(abs(t(Qh) %*% Q_tilde[,j]))                            # best-match cosine (sign-invariant)
      sqrt(max(0, 1-min(best,1)^2))                                       # sin = sqrt(1 - cos^2)
    }))
    
    best_jstar         <- max(abs(t(Qh) %*% Q_tilde[, jstar_global]))    # sin(theta) for j* only
    out$sintheta_jstar <- sqrt(max(0, 1-min(best_jstar,1)^2))            # j* is DGP property, not sample estimate
  }
  out
}

#############################################
#  CLUSTER
#############################################

cl <- makeCluster(max(1L, detectCores()-1L)); registerDoParallel(cl)
clusterExport(cl, c("B_true_3x3","burn","LAMBDA_FIXED","Q_TILDE",
                    "JSTAR_GLOBAL",
                    "generate_shocks","simulate_VAR","align_B",
                    "frob_err","one_rep_cv"))
clusterEvalQ(cl, { library(vars); library(svars); library(clue) })

#############################################
#  EXPERIMENT B — asymptotic behavior
#############################################

T_grid_B <- c(300, 1000, 3000, 10000, 30000)

scenarios_B <- c(
  lapply(seq_along(LAMBDA4_GRID), function(li)
    list(sc=sc_relevant,   lam4=LAMBDA4_GRID[li], lab=LAMBDA4_LABS[li], irr=FALSE)),
  list(list(sc=sc_irrelevant, lam4=LAMBDA4_IRR, lab="Irrelevant omission", irr=TRUE))
)

#SIMULATION LOOP: iterate over scenarios and sample sizes
rows_B <- list(); idx <- 1L
for (scn in scenarios_B) {
  for (Tk in T_grid_B) {
    Nrep_k <- max(30L, round(Nrep*min(1, 1000/Tk))) # fewer reps at large T to save time
    cat(sprintf("B | %s | T=%d\n", scn$lab, Tk))
    
    #PARALLEL MONTE CARLO: one rep = simulate DGP + run CV 
    res <- foreach(r=seq_len(Nrep_k), .combine=rbind,
                   .export=c("B_true_3x3","burn","LAMBDA_FIXED","Q_TILDE",
                             "JSTAR_GLOBAL",
                             "generate_shocks","simulate_VAR","align_B",
                             "frob_err","one_rep_cv")) %dopar% {
                               set.seed(r + 2e4*match(Tk,T_grid_B) +
                                          ifelse(scn$irr, 9e5, 1e6*match(scn$lam4,LAMBDA4_GRID)))  # unique seed per (r, T, scenario)
                               eps <- generate_shocks(scn$lam4, Tk)
                               y3  <- simulate_VAR(eps, Tk, scn$sc)
                               out <- one_rep_cv(y3, Q_TILDE, JSTAR_GLOBAL)
                               data.frame(frob           = out$frob,
                                          sintheta       = out$sintheta,
                                          sintheta_jstar = out$sintheta_jstar,
                                          ok             = out$ok)
                             }
    
    #RESULTS: mean and 90% CI for each metric
    ok  <- res$ok
    z90 <- qnorm(0.95)
    se  <- function(x) sd(x, na.rm=TRUE) / sqrt(sum(!is.na(x)))
    fv  <- res$frob[ok]; sv <- res$sintheta[ok]; sj <- res$sintheta_jstar[ok]
    rows_B[[idx]] <- data.frame(
      T=Tk, lambda4=scn$lam4, lambda4_lab=scn$lab,
      scenario   = ifelse(scn$irr, "Irrelevant omission", "Relevant omission"),
      frob_mean  = mean(fv), frob_lo = mean(fv)-z90*se(fv), frob_hi = mean(fv)+z90*se(fv),
      sin_mean   = mean(sv), sin_lo  = mean(sv)-z90*se(sv), sin_hi  = mean(sv)+z90*se(sv),
      sinj_mean  = mean(sj), sinj_lo = mean(sj)-z90*se(sj), sinj_hi = mean(sj)+z90*se(sj))
    idx <- idx+1L
  }
}

#SAVING
df_B <- bind_rows(rows_B)
df_B$lambda4_lab <- factor(df_B$lambda4_lab, levels=lev_B)
df_B$scenario    <- factor(df_B$scenario, levels=c("Relevant omission","Irrelevant omission"))
saveRDS(df_B, "expB_cv.rds")

#EXPERIMENT A3
# ══════════════════════════════════════════════════════════════

#VARYING CONTAMINATION AT FIXED T: cost of omission as lambda4 grows
T_A3            <- 10000
LAMBDA4_GRID_A3 <- c(1.5, 2.0, 3.0, 4.0, 6.0, 8.0)

rows_A3 <- list()
for (li in seq_along(LAMBDA4_GRID_A3)) {
  cat(sprintf("A3 | lambda4=%.1f | T=%d\n", LAMBDA4_GRID_A3[li], T_A3))
  
 #PARALLEL MC: relevant omission at each contamination level
  res <- foreach(r=seq_len(Nrep), .combine=rbind,
                 .export=c("B_true_3x3","burn","LAMBDA_FIXED","Q_TILDE",
                           "JSTAR_GLOBAL","generate_shocks","simulate_VAR",
                           "align_B","frob_err","one_rep_cv")) %dopar% {
                             set.seed(r + 3e5*li)
                             eps <- generate_shocks(LAMBDA4_GRID_A3[li], T_A3)
                             y3  <- simulate_VAR(eps, T_A3, sc_relevant)
                             out <- one_rep_cv(y3, Q_TILDE, JSTAR_GLOBAL)
                             data.frame(frob=out$frob, ok=out$ok)
                           }
  fv <- res$frob[res$ok]
  rows_A3[[li]] <- data.frame(lambda4   = LAMBDA4_GRID_A3[li],
                              frob_mean = mean(fv),
                              frob_se   = sd(fv)/sqrt(length(fv)))
}

#IRRELEVANT BASELINE: frob error when omission has no effect (P1=P2=0)
cat(sprintf("A3 | Irrelevant | T=%d\n", T_A3))
res_irr <- foreach(r=seq_len(Nrep), .combine=rbind,
                   .export=c("B_true_3x3","burn","LAMBDA_FIXED","Q_TILDE",
                             "JSTAR_GLOBAL","generate_shocks","simulate_VAR",
                             "align_B","frob_err","one_rep_cv")) %dopar% {
                               set.seed(r + 9e5)
                               eps <- generate_shocks(LAMBDA4_IRR, T_A3)
                               y3  <- simulate_VAR(eps, T_A3, sc_irrelevant)
                               out <- one_rep_cv(y3, Q_TILDE, JSTAR_GLOBAL)
                               data.frame(frob=out$frob, ok=out$ok)
                             }
frob_irr <- mean(res_irr$frob[res_irr$ok])                         # baseline: sampling error only

#OMISSION COST = relevant frob - irrelevant baseline
df_A3 <- bind_rows(rows_A3) %>%
  mutate(cost    = frob_mean - frob_irr,                           # net bias due to omission
         cost_lo = cost - qnorm(0.95)*frob_se,
         cost_hi = cost + qnorm(0.95)*frob_se)
saveRDS(df_A3, "expA3_cv.rds")

stopCluster(cl)

#############################################
#  THEME
#############################################

CONT_PAL <- c("Low contamination (lambda4=2)"  = "#2ca02c",
              "Med contamination (lambda4=4)"  = "#ff7f0e",
              "High contamination (lambda4=8)" = "#d62728",
              "Irrelevant omission"            = "#aaaaaa")

THEME_PAPER <- theme_bw(base_size=11, base_family="serif") +
  theme(strip.background = element_rect(fill="grey92", colour="grey70"),
        strip.text       = element_text(family="serif", size=10),
        legend.position  = "bottom",
        plot.title       = element_text(family="serif", size=11, face="bold"),
        plot.subtitle    = element_text(family="serif", size=9, colour="grey30"),
        axis.title       = element_text(family="serif", size=10),
        axis.text        = element_text(family="serif", size=9),
        panel.grid.minor = element_blank())

# ── Fig A3 ────────────────────────────────────────────────────
p_a3 <- ggplot(df_A3, aes(x=lambda4, y=cost)) +
  geom_ribbon(aes(ymin=cost_lo, ymax=cost_hi), alpha=0.15, fill="#1f77b4") +
  geom_line(linewidth=1.3, colour="#1f77b4") +
  geom_point(size=3, colour="#1f77b4") +
  geom_hline(yintercept=0, linetype="dashed", colour="grey50") +
  scale_x_continuous(breaks=LAMBDA4_GRID_A3) +
  labs(title    = sprintf("CV: omission cost vs contamination  (T = %d)", T_A3),
       subtitle = paste(
         "Cost = Frob(relevant) \u2212 Frob(irrelevant).  Ribbon = 90% CI.",
         "DK prediction: cost increases with \u03bb\u2084 (larger \u2225E_W\u2225, fixed \u03b4_j).",
         sep="\n"),
       x = expression(lambda[4]~~"(contamination intensity)"),
       y = "Frobenius cost of omission") +
  THEME_PAPER + theme(legend.position="none")

ggsave("figA3_cost_vs_lambda4.pdf", p_a3, width=7, height=4.5)
cat("figA3 saved.\n")

# ── Fig B1 ───────────────────────────────────────────────────
p_b1 <- ggplot(df_B, aes(x=T, y=frob_mean, colour=lambda4_lab,
                         fill=lambda4_lab, group=lambda4_lab)) +
  geom_ribbon(aes(ymin=frob_lo, ymax=frob_hi), alpha=0.15, colour=NA) +
  geom_line(linewidth=1.2) + geom_point(size=2.5) +
  geom_hline(data = plateau_df %>%
               mutate(scenario=factor("Relevant omission",
                                      levels=c("Relevant omission","Irrelevant omission"))),
             aes(yintercept=plateau, colour=lambda4_lab),
             linetype="dashed", linewidth=0.9, alpha=0.8) +
  facet_wrap(~scenario) +
  scale_colour_manual(values=CONT_PAL) + scale_fill_manual(values=CONT_PAL) +
  scale_x_log10(breaks=T_grid_B, labels=comma) +
  labs(title    = "CV: Frobenius Error vs T",
       subtitle = paste("Relevant: plateaus at B* \u2014 inconsistency (Proposition 1).",
                        "Irrelevant: declines \u2014 consistent when P_1 = P_2.",
                        "Dashed = theoretical B*.", sep="\n"),
       x="Sample size T  (log scale)", y="Mean Frobenius Error",
       colour="Contamination", fill="Contamination") +
  THEME_PAPER
ggsave("figB1_asymptotic.pdf", p_b1, width=12, height=5)
cat("figB1 saved.\n")

# ── Fig B2 ───────────────────────────────────────────

p_b2 <- ggplot(df_B, aes(x=T, colour=lambda4_lab, fill=lambda4_lab,
                         group=lambda4_lab)) +
  # ribbon e solida: su tutti i dati
  geom_ribbon(aes(ymin=sin_lo, ymax=sin_hi), alpha=0.10, colour=NA) +
  geom_line(aes(y=sin_mean), linewidth=1.2) +
  geom_point(aes(y=sin_mean), size=2.5) +
  geom_line(data = df_B %>% filter(scenario == "Relevant omission"),
            aes(y=sinj_mean), linewidth=0.9, linetype="dashed") +
  facet_wrap(~scenario) +
  scale_colour_manual(values=CONT_PAL) +
  scale_fill_manual(values=CONT_PAL) +
  scale_x_log10(breaks=T_grid_B, labels=comma) +
  labs(title    = "CV: Mean sin(\u03b8) vs T",
       subtitle = paste(
         "Solid = average over j=1,2,3.  Dashed = most-affected eigenvector j* (from Q\u0303).",
         "Relevant: both plateau \u2014 persistent angular bias.",
         "Irrelevant: both \u2192 0 \u2014 no angular bias.",
         sep="\n"),
       x      = "Sample size T  (log scale)",
       y      = expression(sin(theta)),
       colour = "Contamination",
       fill   = "Contamination") +
  THEME_PAPER
ggsave("figB2_sintheta_asymptotic.pdf", p_b2, width=12, height=5)
cat("figB2 saved.\n")

##############################################################
#  EXPERIMENT B3 — spectral gap as resilience index
##############################################################

LAMBDA4_B3      <- 2.0
LAMBDA_OBS_GRID <- c(0.1, 0.3, 0.5, 0.7, 0.9)
T_grid_B3       <- c(300, 1000, 3000, 10000, 30000)

VAR4_R2_B3 <- (1 - LAMBDA_FIXED) + LAMBDA_FIXED * LAMBDA4_B3  # = 1.8

# Perturbation direction: based on population P1, P2 for omitted shock
{
  c_vec_b3 <- sc_relevant$B[1:3, 4]
  P1_b3    <- 1.0         * outer(c_vec_b3, c_vec_b3)   # Var(w4, r1) = 1
  P2_b3    <- VAR4_R2_B3  * outer(c_vec_b3, c_vec_b3)   # Var(w4, r2) = 1.8
  S1_b3    <- B_true_3x3 %*% t(B_true_3x3) + P1_b3
  C_b3     <- t(chol(S1_b3))
  EW_b3    <- solve(C_b3) %*% (P2_b3 - P1_b3) %*% t(solve(C_b3))
  V_B3     <- svd(EW_b3)$u[, 1]   # perturbation direction
}


Q_tilde_b3 <- lapply(LAMBDA_OBS_GRID, compute_Q_tilde)
jstar_b3   <- sapply(Q_tilde_b3, function(Qt) which.max(abs(t(Qt) %*% V_B3)))
cat("=== B3: j* per lam_obs ===\n")
cat(sprintf("  lam_obs=%s -> j*=%s\n",
            paste(LAMBDA_OBS_GRID, collapse=","),
            paste(jstar_b3, collapse=",")))

simulate_var_b3 <- function(sc, lam_obs, lambda4, T_tot) {

  sqL2 <- sqrt(c((1-lam_obs) + lam_obs*c(2,3,4),
                 (1-LAMBDA_FIXED) + LAMBDA_FIXED*lambda4))
  y    <- matrix(0, 4L, T_tot)
  y[,1] <- rnorm(4); y[,2] <- rnorm(4)
  half  <- floor((T_tot-burn)/2) + burn
  for (t in 3:T_tot) {
    e     <- rnorm(4) * if (t > half) sqL2 else rep(1, 4)
    y[,t] <- sc$A1%*%y[,t-1] + sc$A2%*%y[,t-2] + sc$B%*%e
  }
  t(y[1:3, (burn+1):T_tot])
}

one_rep_b3 <- function(sc, lam_obs, lambda4, T_tot, Q_tilde, jstar) {
  y3  <- simulate_var_b3(sc, lam_obs, lambda4, T_tot)
  fit <- try(VAR(y3, p=2, type="none"), silent=TRUE)
  if (inherits(fit,"try-error") || !all(Mod(roots(fit))<1)) return(NA_real_)
  sv  <- try(id.cv(fit, SB=floor(nrow(y3)/2)), silent=TRUE)
  if (inherits(sv,"try-error")) return(NA_real_)
  U  <- residuals(fit); T_ <- nrow(U); T2 <- floor(T_/2)
  Ch <- try(t(chol(crossprod(U[1:T2,])/T2)), silent=TRUE)
  if (inherits(Ch,"try-error")) return(NA_real_)
  Wh    <- solve(Ch)%*%(crossprod(U[(T2+1):T_,])/(T_-T2))%*%t(solve(Ch))
  Q_hat <- eigen((Wh+t(Wh))/2, symmetric=TRUE)$vectors
  best <- max(abs(t(Q_hat) %*% Q_tilde[, jstar]))
  sqrt(max(0, 1-min(best,1)^2))
}

cl2 <- makeCluster(max(1L, detectCores()-1L)); registerDoParallel(cl2)
clusterExport(cl2, c("B_true_3x3","sc_relevant","burn","LAMBDA_FIXED",
                     "simulate_var_b3","one_rep_b3","LAMBDA4_B3","VAR4_R2_B3"))
clusterEvalQ(cl2, { library(vars); library(svars) })

z90 <- qnorm(0.95); rows_B3 <- list()
for (i in seq_along(LAMBDA_OBS_GRID)) {
  Qt <- Q_tilde_b3[[i]]; js <- jstar_b3[i]
  clusterExport(cl2, c("Qt","js"))
  for (k in seq_along(T_grid_B3)) {
    Tk     <- T_grid_B3[k]
    Nrep_k <- max(30L, round(Nrep*min(1, 1000/Tk)))
    cat(sprintf("B3 | lam_obs=%.1f | T=%d\n", LAMBDA_OBS_GRID[i], Tk))
    sins <- foreach(r=seq_len(Nrep_k), .combine=c,
                    .export=c("sc_relevant","B_true_3x3","burn","LAMBDA_FIXED",
                              "simulate_var_b3","one_rep_b3",
                              "LAMBDA4_B3","VAR4_R2_B3","Qt","js"),
                    .packages=c("vars","svars")) %dopar% {
                      set.seed(r + 1000*k + 10000*i)
                      one_rep_b3(sc_relevant, LAMBDA_OBS_GRID[i], LAMBDA4_B3, Tk, Qt, js)
                    }
    sins <- sins[!is.na(sins)]
    rows_B3[[length(rows_B3)+1]] <- data.frame(
      lam_obs  = LAMBDA_OBS_GRID[i],
      T        = Tk,
      sin_mean = mean(sins),
      sin_lo   = mean(sins) - z90*sd(sins)/sqrt(length(sins)),
      sin_hi   = mean(sins) + z90*sd(sins)/sqrt(length(sins)))
  }
}

df_B3 <- bind_rows(rows_B3)
saveRDS(df_B3, "expB3_cv.rds")

T_B3_fixed   <- 30000
df_B3_fixed  <- df_B3 %>% filter(T == T_B3_fixed)

p_b3 <- ggplot(df_B3_fixed, aes(x=lam_obs, y=sin_mean)) +
  geom_ribbon(aes(ymin=sin_lo, ymax=sin_hi), alpha=0.15, fill="#1f77b4") +
  geom_line(linewidth=1.3, colour="#1f77b4") +
  geom_point(size=3.5, colour="#1f77b4") +
  scale_x_continuous(breaks=LAMBDA_OBS_GRID) +
  scale_y_continuous(limits=c(0, NA)) +
  labs(title    = sprintf(
    "CV: spectral gap as resilience index  (\u03bb\u2084 = %.1f,  T = %d)",
    LAMBDA4_B3, T_B3_fixed),
    subtitle = paste(
      sprintf("Omitted shock: Var(r2)/Var(r1) = %.1f (consistent with main DGP).", VAR4_R2_B3),
      "sin(\u03b8) on eigenvector j* most aligned with perturbation direction.",
      "Larger spectral gap \u03b4_{j*} \u2192 lower angular bias (Davis\u2013Kahan denominator).",
      sep="\n"),
    x = expression(lambda[obs]~~"(spectral gap driver)"),
    y = expression(paste(sin, theta[j^"*"],
                         "  (most-affected eigenvector)"))) +
  THEME_PAPER
ggsave("figB3_gap_resilience.pdf", p_b3, width=7, height=4.5)
cat("figB3 saved.\n")

##############################################################
#  DECOMPOSITION FIGURE
##############################################################

one_rep_decomp <- function(y3) {
  
  #--- 1. TRUE POPULATION QUANTITIES (Proposition 2.2 benchmark) ---
  Bt      <- B_true_3x3
  C_tilde <- t(chol(Bt%*%t(Bt)))          # true Cholesky: C_tilde C_tilde' = B_tilde B_tilde'
  Q_tilde <- solve(C_tilde, Bt)            # true whitened directions: C_tilde Q_tilde = B_tilde exactly
  
  #--- 2. VAR ESTIMATION AND CV IDENTIFICATION ---
  fit <- try(VAR(y3, p=2, type="none"), silent=TRUE)
  if (inherits(fit,"try-error") || !all(Mod(roots(fit))<1))
    return(data.frame(chan1=NA_real_, chan2=NA_real_, cross=NA_real_))
  sv <- try(id.cv(fit, SB=floor(nrow(y3)/2)), silent=TRUE)
  if (inherits(sv,"try-error"))
    return(data.frame(chan1=NA_real_, chan2=NA_real_, cross=NA_real_))
  
  #--- 3. CONTAMINATED CHOLESKY AND WHITENED MATRIX ---
  U     <- residuals(fit); T_ <- nrow(U); T2 <- floor(T_/2)
  C_hat <- try(t(chol(crossprod(U[1:T2,])/T2)), silent=TRUE)             # contaminated Cholesky (Channel 1)
  if (inherits(C_hat,"try-error"))
    return(data.frame(chan1=NA_real_, chan2=NA_real_, cross=NA_real_))
  Wh        <- solve(C_hat)%*%(crossprod(U[(T2+1):T_,])/(T_-T2))%*%t(solve(C_hat))  # W_hat
  
  #--- 4. ESTIMATED B AND ALIGNMENT ---
  B_hat_raw <- C_hat %*% eigen((Wh+t(Wh))/2, symmetric=TRUE)$vectors     # B_hat = C_hat * Q_hat
  nc        <- function(x) x/sqrt(sum(x^2))
  perm      <- as.integer(clue::solve_LSAP(
    1 - abs(t(apply(Bt,2,nc)) %*% apply(B_hat_raw,2,nc))))               # optimal column permutation
  B_hat     <- B_hat_raw[,perm]
  for (j in 1:3) if (sum(B_hat[,j]*Bt[,j])<0) B_hat[,j] <- -B_hat[,j]  # sign alignment
  Q_hat <- solve(C_hat, B_hat)                                             # Q_hat = C_hat^{-1} B_hat
  
  #--- 5. THREE-TERM DECOMPOSITION (Proposition 2.2) ---
  # B_hat - B_tilde = (C_hat-C_tilde)Q_tilde        [Channel 1: whitening error]
  #                 + C_tilde(Q_hat-Q_tilde)         [Channel 2: eigenvector rotation]
  #                 + (C_hat-C_tilde)(Q_hat-Q_tilde) [Cross term: interaction]
  dC <- C_hat-C_tilde; dQ <- Q_hat-Q_tilde
  data.frame(chan1=sqrt(sum((dC%*%Q_tilde)^2)),                           # ||Channel 1||_F
             chan2=sqrt(sum((C_tilde%*%dQ)^2)),                           # ||Channel 2||_F
             cross=sqrt(sum((dC%*%dQ)^2)))                                # ||Cross term||_F
}

T_grid_D <- c(300, 1000, 3000, 10000)
clusterExport(cl2, c("one_rep_decomp","align_B","frob_err"))
clusterEvalQ(cl2, library(clue))

rows_D <- list()
for (li in seq_along(LAMBDA4_GRID)) {
  for (k in seq_along(T_grid_D)) {
    Tk     <- T_grid_D[k]
    Nrep_k <- max(30L, round(Nrep*min(1, 1000/Tk)))
    cat(sprintf("Decomp | %s | T=%d\n", LAMBDA4_LABS[li], Tk))
    res <- foreach(r=seq_len(Nrep_k), .combine=rbind,
                   .export=c("B_true_3x3","burn","LAMBDA_FIXED","generate_shocks",
                             "simulate_VAR","sc_relevant","one_rep_decomp"),
                   .packages=c("vars","svars","clue")) %dopar% {
                     set.seed(r + 7e4*k + 7e5*li)
                     one_rep_decomp(simulate_VAR(generate_shocks(LAMBDA4_GRID[li], Tk),
                                                 Tk, sc_relevant))
                   }
    res <- na.omit(res)
    rows_D[[length(rows_D)+1]] <- data.frame(
      T=Tk, lambda4=LAMBDA4_GRID[li], lambda4_lab=LAMBDA4_LABS[li],
      chan1=mean(res$chan1), chan2=mean(res$chan2), cross=mean(res$cross))
  }
}

stopCluster(cl2)

df_decomp <- bind_rows(rows_D) %>%
  mutate(lambda4_lab=factor(lambda4_lab, levels=LAMBDA4_LABS)) %>%
  pivot_longer(c(chan1,chan2,cross), names_to="term", values_to="norm") %>%
  mutate(term=factor(term,
                     levels=c("chan1","chan2","cross"),
                     labels=c("Channel 1: Cholesky",
                              "Channel 2: eigenvectors",
                              "Cross term")))

saveRDS(df_decomp, "decomp_cv.rds")

p_decomp <- ggplot(df_decomp, aes(x=factor(T), y=norm, fill=term)) +
  geom_col(position=position_dodge(width=0.75), width=0.7) +
  facet_wrap(~lambda4_lab, nrow=1) +
  scale_fill_manual(values=c("Channel 1: Cholesky"    = "#1f77b4",
                             "Channel 2: eigenvectors" = "#d62728",
                             "Cross term"              = "#aaaaaa"),
                    name=NULL) +
  labs(title    = "Error decomposition: two channels of identification bias",
       subtitle = paste(
         "B\u0302\u2212B\u0303 = (C\u0302\u2212C\u0303)Q\u0303 + C\u0303(Q\u0302\u2212Q\u0303) + (C\u0302\u2212C\u0303)(Q\u0302\u2212Q\u0303).  Identity exact.",
         "Channel 1 (Cholesky): constant across \u03bb\u2084.  Channel 2 (eigenvectors): grows with \u03bb\u2084.",
         sep="\n"),
       x="Sample size T", y="Mean Frobenius norm of each term") +
  THEME_PAPER + theme(legend.key.size=unit(0.4,"cm"))
ggsave("figD1_decomposition.pdf", p_decomp, width=12, height=5)
cat("figD1 saved.\n")

