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
  mutate(across(c(y, p, cp, tr, nbr), ~ log(.))) %>% 
  tidyr::drop_na()

var_data <- ts(df_var_final[,-1], 
               start = c(1965, 1), 
               frequency = 12)

plot.ts(var_data, main = "Time Series")
adf.test(var_data[,1])
adf.test(var_data[,2])
adf.test(var_data[,3])

# ------------ VAR MODEL ----------------------

lag_selection <- VARselect(var_data, lag.max = 12, type = "const", season = 12)
print(lag_selection$selection)

# LL2008 CHOSE 13 LAGS. 

var_model13 <- VAR(var_data, p = 13, type = "const", season = 12)
roots(var_model13)
print(var_model13)
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

var_model3 <- VAR(var_data, p = 3, type = "const",season = 12)
roots(var_model3)

# MORE PARSIMONIOUS MODEL BASED ON AIC IC

var3.serial <- serial.test(var_model3)
var3.serial
plot(var3.serial, names = "y")
plot(var3.serial, names = "p")
plot(var3.serial, names = "cp")
plot(var3.serial, names = "tr")
plot(var3.serial, names = "nbr")
plot(var3.serial, names = "ff")

var3.arch <- arch.test(var_model3)
var3.arch



var3.norm <- normality.test(var_model3)
var3.norm


sb1 <- which(df_var_final$date == as.Date("1979-10-01"))
sb2 <- which(df_var_final$date == as.Date("1984-02-01"))


#------------ BREAK VALIDATION -------------------

var_data_monetary <- var_data[, c("tr", "nbr", "ff")]
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

logL_2breaks <- two_breaks13$Lik
logL_break1  <- break_1.2_13$Lik
logL_break2  <- break_2.3_13$Lik

dof_diff <- 3  #2 break model has three more parameters

#is 1st break necessary? H0: no
lr_stat_A <- 2 * (logL_2breaks - logL_break2)
p_val_A   <- 1 - pchisq(lr_stat_A, df = dof_diff)
p_val_A

#is 2nd break necessary? H0: no
lr_stat_B <- 2 * (logL_2breaks - logL_break1)
p_val_B   <- 1 - pchisq(lr_stat_B, df = dof_diff)
p_val_B
#-----------------------------------

var_completo13 <- VAR(var_data, p = 13, type = "const", season = 12)

two_breaks_completo <- id.cv(var_completo13, SB = sb1, SB2 = sb2)

fevd_completo <- fevd(two_breaks_completo, n.ahead = 24)

fevd_data_long <- lapply(names(fevd_completo), function(target_var) {
  as.data.frame(fevd_completo[[target_var]]) %>%
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

ggplot(fevd_data_long, aes(x = horizon, y = percentage, fill = shock)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  facet_wrap(~ target, ncol = 2) + 
  scale_fill_brewer(palette = "Set2") + 
  labs(
    title = "Structural FEVD - 6 Variable Model (p = 13)",
    x = "Orizzonte Temporale (Mesi)",
    y = "Percentuale di Varianza Spiegata",
    fill = "Origine dello Shock"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")


#------------ TESTING MONEY MARKET MODELS -----------------

restr_NBR <- matrix(c(NA, NA, NA,  
                      0, NA, 0, 
                      NA,  NA, NA),
                      nrow = 3, ncol = 3, byrow = T)
restr_NBR_TR <- matrix(c(NA, 0, 0,  
                         NA, NA, 0, 
                         NA,  NA, NA),
                       nrow = 3, ncol = 3, byrow = T)

model_NBR <- id.cv(var_mon13, SB = sb1, SB2 = sb2, restriction_matrix = restr_NBR)
model_NBR_TR <- id.cv(var_mon13, SB = sb1, SB2 = sb2, restriction_matrix = restr_NBR_TR)
summary(model_NBR)
summary(model_NBR_TR)
summary(two_breaks13)

#---------- IRF analysis ---------------

irf_unres <- irf(two_breaks13, n.ahead = 24)
irf_nbr   <- irf(model_NBR,    n.ahead = 24)
irf_nbr_tr<- irf(model_NBR_TR, n.ahead = 24)

plot(irf_unres)
plot(irf_nbr)
plot(irf_nbr_tr)

boot_irf_unres <- wild.boot(
  two_breaks13,
  n.ahead = 24,
  nboot = 500
)

boot_irf_nbr <-  wild.boot(
  model_NBR,
  n.ahead = 24,
  nboot = 500
)
  
boot_irf_nbr_tr <-  wild.boot(
  model_NBR_TR,
  n.ahead = 24,
  nboot = 500
)
  
plot(boot_irf_unres)
plot(boot_irf_nbr)
plot(boot_irf_nbr_tr)



#------------ BREAK VALIDATION -------------------

var_mon3 <- VAR(var_data_monetary, p = 3, type = "const")

two_breaks3 <- id.cv(var_mon3, SB = sb1, SB2 = sb2)
summary(two_breaks3)

break_1.2_3 <- id.cv(var_mon3, SB = sb1)
summary(break_1.2_3)

break_2.3_3 <- id.cv(var_mon3, SB = sb2)
summary(break_2.3_3)

logL_2breaks.3 <- two_breaks3$Lik
logL_break1.3  <- break_1.2_3$Lik
logL_break2.3  <- break_2.3_3$Lik

#is 1st break necessary? H0: no
lr_stat_A3 <- 2 * (logL_2breaks.3 - logL_break2.3)
p_val_A3   <- 1 - pchisq(lr_stat_A3, df = dof_diff)
p_val_A3

#is 2nd break necessary? H0: no
lr_stat_B3 <- 2 * (logL_2breaks.3 - logL_break1.3)
p_val_B3   <- 1 - pchisq(lr_stat_B3, df = dof_diff)
p_val_B3

#------------ TESTING MONEY MARKET MODELS -----------------

model_NBR_3 <- id.cv(var_mon3, SB = sb1, SB2 = sb2, restriction_matrix = restr_NBR)
model_NBR_TR_3 <- id.cv(var_mon3, SB = sb1, SB2 = sb2, restriction_matrix = restr_NBR_TR)
summary(model_NBR_3)
summary(model_NBR_TR_3)
summary(two_breaks3)

#---------- IRF analysis ---------------

irf_unres3 <- irf(two_breaks3, n.ahead = 24)
irf_nbr3   <- irf(model_NBR_3, n.ahead = 24)
irf_nbr_tr3 <- irf(model_NBR_TR_3, n.ahead = 24)

plot(irf_unres3)
plot(irf_nbr3)
plot(irf_nbr_tr3)

boot_irf_unres3 <- wild.boot(
  two_breaks3,
  n.ahead = 24,
  nboot = 500
)

boot_irf_nbr3 <-  wild.boot(
  model_NBR_3,
  n.ahead = 24,
  nboot = 500
)

boot_irf_nbr_tr3 <-  wild.boot(
  model_NBR_TR_3,
  n.ahead = 24,
  nboot = 500
)

plot(boot_irf_unres3)
plot(boot_irf_nbr3)
plot(boot_irf_nbr_tr3)












  
  
  
