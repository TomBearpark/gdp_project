if(!require(pacman)) install.packages('pacman')
pacman::p_load(tidyverse, fixest, haven, marginaleffects, broom, useful)
theme_set(theme_classic())

dir  <- "~/Library/CloudStorage/Dropbox/gdp-temp/replication/"
code <- "~/Documents/GitHub/gdp_project/"

set.seed(123)
source(paste0(code, "code/0_funcs.R"))

# SPECIFY PARAMETERS ------------------------------------------------------
type    <- "levels"
lags    <- 0
spec    <- "poly2"
warming <- 4

# GET BASELINES -----------------------------------------------------------

pre.yrs   <- 1990:2019
proj.yrs  <- 2020:2100
yrs       <- c(pre.yrs, proj.yrs) 

TT        <- length(yrs)

pre.ids  <- 1:length(pre.yrs)
post.ids <- (length(pre.yrs)+1):length(yrs)
yr.ids   <- c(pre.ids, post.ids) 

# Filter to period we are going to use in plotting / getting baseline
df <- df.reg %>% 
  filter(year >= min(pre.yrs))  %>% 
  filter(!is.na(g), !is.na(temp1), !is.na(y)) %>% 
  group_by(ID) %>%
  add_tally() %>%
  ungroup() %>%
  filter(n == max(n))

# Get some averages which are used to define counterfactual
df.base <- df %>% 
  summarize(g = mean(g, na.rm=T), 
            temp = mean(temp1, na.rm=T), 
            pop = mean(pop, na.rm=T),
            .by = ID) %>% 
  select(ID, temp, g, pop) %>% 
  mutate(g = if_else(g > 0.03, 0.03, g))

NN <- length(unique(df$ID))
print(NN)

# Create matrices to store data
base <- gen_mats(NN, TT, df.base$ID, yrs)
proj <- gen_mats(NN, TT, df.base$ID, yrs)

## Format pre-projection data -----------------
pre.proj <- format_pre_proj(pre.ids, pre.yrs, df, base, proj)
base <- pre.proj$base
proj <- pre.proj$proj

## Format post-projection baseline ------------------
post.proj <- format_post_proj_baseline(post.ids, base, proj, warming, TT)
base <- post.proj$base
proj <- post.proj$proj

## Sense check the data 
plotdf_mats(toPlot, base, proj)

# PROJECTION --------------------------------------------------------------

opts <- expand_grid(lags=0:10, type=c("levels", "growth"))

pdf <- pmap(
  opts,
  function(lags, type){
    m  <- run_reg(df.reg, type=type, lags=lags, spec=spec)
    # me.cum <- get_me_cum(m, lags, type=type)
    get_damages(m, type, base, 
                proj, post.ids, TT, lags, df.base, yr.ids, yrs, 
                uncertainty=T)
  }
)

pdf %>% 
  mutate(lags = as.factor(lags)) %>% 
  ggplot(aes(x=year, y=damage)) +
  geom_vline(xintercept=2020, linetype="dashed") +
  geom_line() +
  facet_wrap(~type+lags, nrow=2) + 
  scale_x_continuous(breaks=c(2000, 2050, 2100))
  