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
df.reg <- load_historic_data(paste0(dir, "/replication/"), lags=max_lags)

# Filter to period we are going to use in plotting / getting baseline. 
# Should interpolate missing values to get full panel?

df <- df.reg %>% 
  filter(year >= min(yrs$pre)) %>% 
  filter(!is.na(g), !is.na(temp1), !is.na(y))  %>% 
  group_by(ID) %>%
  add_tally() %>%
  ungroup() %>%
  filter(n == max(n))

# Get some averages which are used to define counter-factual
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
base <- gen_mats(NN, yrs, df.base$ID)
proj <- gen_mats(NN, yrs, df.base$ID)

## Format pre-projection data -----------------
pre.proj <- format_pre_proj(yrs, df, base, proj)
base <- pre.proj$base
proj <- pre.proj$proj

## Format post-projection baseline ------------------
warming <- load_warming(dir, "median") %>% filter(ID %in% df$ID)
# warming <- 4
post.proj <- format_post_proj_baseline(yrs, base, proj, warming)
base <- post.proj$base
proj <- post.proj$proj

# PROJECTION --------------------------------------------------------------

opts <- expand_grid(lags=0:max_lags, type=c("levels", "growth"), 
                    FE = c("ID+time1", 
                           "ID+time1+ID[time1]",
                           "ID+time1+ID[time1]+ID[time2]"))
df.pop <- df.base %>% select(ID, pop)

pdf <- pmap(
  opts,
  function(lags, type, FE){
    
    m  <- run_reg(df.reg, type=type, lags=lags, 
                  FE = FE)
    
    print(m)
    dam.df <- get_damages(m=m, 
                          type=type, 
                          base=base, 
                          proj=proj,
                          yrs=yrs, 
                          lags=lags, 
                          df.pop=df.pop, 
                          uncertainty=F)
    dam.df$central %>% 
      mutate(type=type, lags=lags, FE=FE) %>% 
      select(type, lags, everything())
  }
)

bind_rows(pdf) %>% 
  filter(year == 2100) %>% 
  ggplot() + 
  geom_point(aes(x = lags, y = damage, color = type)) +
  facet_wrap(~FE) +
  scale_x_continuous(breaks = 0:max_lags) +
  theme_bw()
