library(fredr)
library(tidyverse)
library(vars)
library(svars)
library(tseries)

# ----------- LOADING DATASETS ---------------------

fredr_set_key("be238a2a428573d02dc7925d8742c9df")

industrial_prod <- fredr(
  series_id = "INDPRO",
  observation_start = as.Date("1965-01-01"),
  observation_end = as.Date("1996-12-01")
)

cpi <- fredr(
  series_id = "CPIAUCSL",
  observation_start = as.Date("1965-01-01"),
  observation_end = as.Date("1996-12-01")
)

csv <- read.csv("commodity.csv")
commodities <- csv %>% 
  dplyr::select(date, price_index_2000_base) 
commodities <- commodities[61:444,]  

tr <- fredr(
  series_id = "TOTRESNS",
  observation_start = as.Date("1965-01-01"),
  observation_end = as.Date("1996-12-01")
)

nbr <- fredr(
  series_id = "BOGNONBR",
  observation_start = as.Date("1965-01-01"),
  observation_end = as.Date("1996-12-01")
) 

ff <- fredr(
  series_id = "FEDFUNDS",
  observation_start = as.Date("1965-01-01"),
  observation_end = as.Date("1996-12-01")
) 

br <- fredr(
  series_id = "BORROW",
  observation_start = as.Date("1965-01-01"),
  observation_end = as.Date("1996-12-01")
) 

commodities <- commodities %>%
  mutate(date = as.Date(date)) %>% 
  rename(cp = price_index_2000_base)

df_var <- list(
  industrial_prod %>% dplyr::select(date, y = value),
  cpi %>% dplyr::select(date, p = value),
  commodities,
  tr %>% dplyr::select(date, tr = value),
  nbr %>% dplyr::select(date, nbr = value),
  ff %>% dplyr::select(date, ff = value)
) %>% 
  reduce(full_join, by = "date") %>%
  arrange(date) %>%
  drop_na() 

df_var_final <- df_var %>%
  mutate(across(c(y, p, cp, tr, nbr), ~ log(.) - lag(log(.)), .names = "dlog_{.col}")) %>%
  mutate(d_ff = ff - lag(ff)) %>%
  dplyr::select(-c(y, p, cp, tr, nbr, ff)) %>%
  tidyr::drop_na()

var_data <- ts(df_var_final[,-1], 
               start = c(1965, 1), 
               frequency = 12)

plot(var_data, 
     main = "Log-differentiated Time Series", 
     panel = function(...) {
       lines(...) 
       abline(v = c(1979.75, 1984.083), col = "red", lty = 2, lwd = 1.5)
     })
adf.test(var_data[,4])
adf.test(var_data[,5])
adf.test(var_data[,6]) #all monetary ts stationary

# ------------ VAR MODEL ----------------------

lag_selection <- VARselect(var_data, lag.max = 12, type = "const", season = 12)
print(lag_selection$selection)

# LL2008 CHOSE 13 LAGS. 

var_model13 <- VAR(var_data, p = 13, type = "const", season = 12)
roots(var_model13)

var13.serial <- serial.test(var_model13)
var13.serial

plot(var13.serial, names = "dlog_y")
plot(var13.serial, names = "dlog_p")
plot(var13.serial, names = "dlog_cp")
plot(var13.serial, names = "dlog_tr")
plot(var13.serial, names = "dlog_nbr")
plot(var13.serial, names = "d_ff")

var13.arch <- arch.test(var_model13)
var13.arch

var13.norm <- normality.test(var_model13)
var13.norm


var_data_monetary <- var_data[, c("dlog_tr", "dlog_nbr", "d_ff")]

restr_NBR <- matrix(c(NA, NA, NA,  
                      0, NA, 0, 
                      NA,  NA, NA),
                    nrow = 3, ncol = 3, byrow = T)
restr_NBR_TR <- matrix(c(NA, 0, 0,  
                         NA, NA, 0, 
                         NA,  NA, NA),
                       nrow = 3, ncol = 3, byrow = T)


sb1 <- which(df_var_final$date == as.Date("1979-10-01"))
sb2 <- which(df_var_final$date == as.Date("1984-02-01"))

#------------ BREAK VALIDATION -------------------

var_mon13 <- VAR(var_data_monetary, p = 13, type = "const")

two_breaks13 <- id.cv(var_mon13, SB = sb1, SB2 = sb2)
summary(two_breaks13)

fevd13 <- fevd(two_breaks13)
plot(fevd13)

break_1.2_13 <- id.cv(var_mon13, SB = sb1)
summary(break_1.2_13)

fevd1.2_13 <- fevd(break_1.2_13)
plot(fevd1.2_13)

break_2.3_13 <- id.cv(var_mon13, SB = sb2)
summary(break_2.3_13)

logL_2breaks_13 <- two_breaks13$Lik
logL_break1_13  <- break_1.2_13$Lik
logL_break2_13  <- break_2.3_13$Lik

dof_diff <- 3
#is 1st break necessary? H0: no
lr_stat_A_13 <- 2 * (logL_2breaks_13 - logL_break2_13)
p_val_A_13   <- 1 - pchisq(lr_stat_A_13, df = dof_diff)
p_val_A_13

#is 2nd break necessary? H0: no
lr_stat_B_13 <- 2 * (logL_2breaks_13 - logL_break1_13)
p_val_B_13   <- 1 - pchisq(lr_stat_B_13, df = dof_diff)
p_val_B_13


two_breaks_completo_13 <- id.cv(var_model13, SB = sb1, SB2 = sb2)

fevd_completo_13 <- fevd(two_breaks_completo_13, n.ahead = 24)

fevd_data_long_13 <- lapply(names(fevd_completo_13), function(target_var) {
  as.data.frame(fevd_completo_13[[target_var]]) %>%
    mutate(
      horizon = 1:n(),
      target = target_var
    )
}) %>% 
  bind_rows() %>%
  pivot_longer(
    cols = -c(horizon, target), 
    names_to = "shock", 
    values_to = "percentage"
  )

ggplot(fevd_data_long_13, aes(x = horizon, y = percentage, fill = shock)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  facet_wrap(~ target, ncol = 2) + 
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Structural FEVD - 6 Variable Model (p = 13)",
    x = "Horizon (Months)",
    y = "Percentage Variance Explained",
    fill = "Shock Origin"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

#------------ TESTING MONEY MARKET MODELS -----------------

model_NBR_13 <- id.cv(var_mon13, SB = sb1, SB2 = sb2, restriction_matrix = restr_NBR)
model_NBR_TR_13 <- id.cv(var_mon13, SB = sb1, SB2 = sb2, restriction_matrix = restr_NBR_TR)
summary(model_NBR_13)
summary(model_NBR_TR_13)
summary(two_breaks13)

#---------- IRF analysis ---------------

irf_unres_13 <- irf(two_breaks13, n.ahead = 24)
irf_nbr_13   <- irf(model_NBR_13, n.ahead = 24)
irf_nbr_tr_13 <- irf(model_NBR_TR_13, n.ahead = 24)


plot(irf_unres_13)
plot(irf_nbr_13)
plot(irf_nbr_tr_13)

boot_irf_unres_13 <- wild.boot(
  two_breaks13,
  n.ahead = 24,
  nboot = 500
)

boot_irf_nbr_13 <-  wild.boot(
  model_NBR_13,
  n.ahead = 24,
  nboot = 500
)

boot_irf_nbr_tr_13 <-  wild.boot(
  model_NBR_TR_13,
  n.ahead = 24,
  nboot = 500
)

plot(boot_irf_unres_13)
plot(boot_irf_nbr_13)
plot(boot_irf_nbr_tr_13)

#too many parameters to be estimated without restrictions: 13 LAGS SEEM TO PROVIDE
#INADEQUATE IRFs

#----------------------------------------------------------------

# MORE PARSIMONIOUS MODEL BASED ON AIC IC

var_model2 <- VAR(var_data, p = 2, type = "const", season = 12)
roots(var_model2)

var2.serial <- serial.test(var_model2)
var2.serial
plot(var2.serial, names = "dlog_y")
plot(var2.serial, names = "dlog_p")
plot(var2.serial, names = "dlog_cp")
plot(var2.serial, names = "dlog_tr")
plot(var2.serial, names = "dlog_nbr")
plot(var2.serial, names = "d_ff")

var2.arch <- arch.test(var_model2)
var2.arch

var2.norm <- normality.test(var_model2)
var2.norm

#------------ BREAK VALIDATION -------------------
VARselect(var_data_monetary, type = "const", season = 12)

var_mon2 <- VAR(var_data_monetary, p = 2, type = "const")

two_breaks2 <- id.cv(var_mon2, SB = sb1, SB2 = sb2)
summary(two_breaks2)


break_1.2_2 <- id.cv(var_mon2, SB = sb1)
summary(break_1.2_2)


break_2.3_2 <- id.cv(var_mon2, SB = sb2)
summary(break_2.3_2)

logL_2breaks_2 <- two_breaks2$Lik
logL_break1_2  <- break_1.2_2$Lik
logL_break2_2  <- break_2.3_2$Lik

#is 1st break necessary? H0: no
lr_stat_A_2 <- 2 * (logL_2breaks_2 - logL_break2_2)
p_val_A_2   <- 1 - pchisq(lr_stat_A_2, df = dof_diff)
p_val_A_2

#is 2nd break necessary? H0: no
lr_stat_B_2 <- 2 * (logL_2breaks_2 - logL_break1_2)
p_val_B_2   <- 1 - pchisq(lr_stat_B_2, df = dof_diff)
p_val_B_2

var_completo2 <- VAR(var_data, p = 2, type = "const", season = 12)

two_breaks_completo_2 <- id.cv(var_completo2, SB = sb1, SB2 = sb2)

fevd_completo_2 <- fevd(two_breaks_completo_2, n.ahead = 24)

fevd_data_long_2 <- lapply(names(fevd_completo_2), function(target_var) {
  as.data.frame(fevd_completo_2[[target_var]]) %>%
    mutate(
      horizon = 1:n(),
      target = target_var
    )
}) %>% 
  bind_rows() %>%
  pivot_longer(
    cols = -c(horizon, target), 
    names_to = "shock", 
    values_to = "percentage"
  )

ggplot(fevd_data_long_2, aes(x = horizon, y = percentage, fill = shock)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  facet_wrap(~ target, ncol = 2) + 
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Structural FEVD - 6 Variable Model (p = 2)",
    x = "Orizzonte Temporale (Mesi)",
    y = "Percentuale di Varianza Spiegata",
    fill = "Origine dello Shock"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

#------------ TESTING MONEY MARKET MODELS -----------------

model_NBR_2 <- id.cv(var_mon2, SB = sb1, SB2 = sb2, restriction_matrix = restr_NBR)
model_NBR_TR_2 <- id.cv(var_mon2, SB = sb1, SB2 = sb2, restriction_matrix = restr_NBR_TR)
summary(model_NBR_2)
summary(model_NBR_TR_2)
summary(two_breaks2)

#---------- IRF analysis ---------------

irf_unres_2 <- irf(two_breaks2, n.ahead = 24)
irf_nbr_2   <- irf(model_NBR_2, n.ahead = 24)
irf_nbr_tr_2 <- irf(model_NBR_TR_2, n.ahead = 24)


plot(irf_unres_2)
plot(irf_nbr_2)
plot(irf_nbr_tr_2)

boot_irf_unres_2 <- wild.boot(
  two_breaks2,
  n.ahead = 24,
  nboot = 500
)

boot_irf_nbr_2 <-  wild.boot(
  model_NBR_2,
  n.ahead = 24,
  nboot = 500
)

boot_irf_nbr_tr_2 <-  wild.boot(
  model_NBR_TR_2,
  n.ahead = 24,
  nboot = 500
)

plot(boot_irf_unres_2)
plot(boot_irf_nbr_2)
plot(boot_irf_nbr_tr_2)

