if(!require(pacman)) install.packages('pacman')
pacman::p_load(tidyverse, fixest, haven, marginaleffects, broom, useful, readxl)
theme_set(theme_classic())

dir  <- "~/Library/CloudStorage/Dropbox/gdp-temp/"
dir.out <- paste0(dir, "/out/projection/")
code <- "~/Documents/GitHub/gdp_project/"

set.seed(123)
source(paste0(code, "code/0_funcs.R"))

# SPECIFY PARAMETERS ------------------------------------------------------

max_lags <- 10

# GET BASELINES -----------------------------------------------------------

yrs <- get_years()

# Load regression data 
df.reg <- load_historic_data(paste0(dir, "/replication/"), lags=max_lags, 
                             min.year = 1900, tempvar = 'cru_mwtemp') %>% 
  select(ID, time1, time2, y, g, contains(c('temp1', 'temp2'))) %>% 
  na.omit()

df.bal <- df.reg %>% 
  group_by(ID) %>% 
  add_tally() %>% 
    ungroup() %>% 
  filter(n == max(n))

table(df.bal$n)

df.r <- df.bal
for (var in c("temp1", "temp2", "dtemp1", "dtemp2")) {
  for(ll in 0:max_lags){
    fml   <- as.formula(paste0("l", ll, "_", var, " ~ 1 | ID + time1 + ID[time1] + ID[time2]"))  
    model <- feols(fml, data = df.r)
    df.r[[paste0("l", ll, "_", var)]] <- resid(model)
  }
}
for (var in c("y", "g")) {
  fml   <- as.formula(paste0(var, " ~ 1 | ID + time1 + ID[time1] + ID[time2]"))  
  model <- feols(fml, data = df.r)
  df.r[[var]] <- resid(model, na.rm=F)
}

list(
  run_reg(df.bal, "levels", 2, FE='ID + time1 + ID[time1]'), 
  run_reg(df.r, "levels", 2, FE = NULL)
) %>% 
  coefplot()

# CV ----------------------------------------------------------------------
library(tidymodels)
set.seed(2)

df.cv <- df.r %>% select(-time1, -time2)

cv <- group_vfold_cv(df.cv, v = 5, repeats = 100, group='ID')
# cv <- vfold_cv(df.cv, v = 5, repeats = 10)

pdf <- map_dfr(
  cv$splits, 
  function(x){ 
    
    # Extract training and testing data  
    train_data <- analysis(x)
    test_data  <- assessment(x)
    
    # Run regression and predict
    map_dfr(
      0:max_lags, 
      function(max_lag){
        m.l    <- run_reg(train_data, "levels", max_lag, FE=NULL, cluster = NULL)
        pred.l <- predict(m.l, newdata = test_data) 
        
        out.l <- tibble(rmse = rmse_vec(test_data$y, pred.l), 
                        max_lag = max_lag,
                        fold = x$id, type = "levels")
        
        m.g    <- run_reg(train_data, "growth", max_lag, FE=NULL, cluster = NULL)
        pred.g <- predict(m.g, newdata = test_data) 
        
        out.g <- tibble(rmse = rmse_vec(test_data$g, pred.g), 
                        max_lag = max_lag,
                        fold = x$id, type = "growth")
        
        bind_rows(out.l, out.g) 
      }
    ) 
  }
)

pdf %>% 
  ggplot() + 
  geom_point(aes(x = max_lag, y = rmse))+
  facet_wrap(~type, scales='free') + 
  scale_x_continuous(breaks = 0:max_lags) 

pdf %>% 
  group_by(max_lag, type) %>% 
  summarize(rmse = median(rmse)) %>% 
  ggplot() + 
  geom_line(aes(x = max_lag, y = rmse)) + 
  facet_wrap(~type, scales='free') + 
  scale_x_continuous(breaks = 0:max_lags) 
