if(!require(pacman)) install.packages('pacman')
pacman::p_load(tidyverse, fixest, haven, marginaleffects, broom, useful, 
               readxl)
theme_set(theme_classic())

dir  <- "~/Library/CloudStorage/Dropbox/gdp-temp/"
dir.out <- paste0(dir, "/out/projection/")
code <- "~/Documents/GitHub/gdp_project/"

set.seed(123)
source(paste0(code, "code/0_funcs.R"))

# SPECIFY PARAMETERS ------------------------------------------------------

warming_list  <- 0:8
max_lags <- 10
lags     <- 10
type     <- "levels"

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


# PROJECTION --------------------------------------------------------------

df.pop <- df.base %>% select(ID, pop)
m  <- run_reg(df.reg, type=type, lags=lags)

# Uniform warming 
pdf <- map(
  set_names(warming_list),
  function(ww){
    df.warming <- tibble(ID = unique(df$ID), warming = ww)
    
    post.proj <- format_post_proj_baseline(yrs, 
                                           base, proj, ww)
    base <- post.proj$base
    proj <- post.proj$proj
    
    get_damages(m=m, 
                type=type, 
                base=base, 
                proj=proj,
                yrs=yrs, 
                lags=lags, 
                df.pop=df.pop, 
                uncertainty=F)
  }
)

map_dfr(seq_along(pdf), 
        function(ww){
          pdf[[ww]] %>% pluck("central") %>% mutate(warming=ww)
        }) %>% 
  ggplot() + 
  geom_line(aes(x = year, y= damage, color = warming, group=warming)) + 
  scale_color_viridis_c() 

map_dfr(seq_along(pdf), 
        function(ww){
          pdf[[ww]] %>% pluck("central") %>% mutate(warming=warming_list[[ww]])
        })  %>% 
  filter(year == 2100) %>% 
  ggplot() + 
  geom_point(aes(x = warming, y = damage))

# USA only 
plot.ids <- c("USA", "RUS", "SDN")
map_dfr(seq_along(pdf), 
        function(ww){
          base <- pdf[[ww]]$base
          proj <- pdf[[ww]]$proj
          plotdf_mats(plot.ids, base, proj, F) %>% 
            mutate(warming = warming_list[[ww]])
        })  %>% 
  filter(year == 2100, scen=='proj', name == "dam") %>% 
  ggplot() + 
  geom_point(aes(x = warming, y = value)) + 
  geom_line(aes(x = warming, y = value)) + 
  facet_wrap(~ ID, scales='free') 



# RCP8.5 warming
pdf.rcp <- map_dfr(
  c("median", "p10", "p90"), 
  function(scen){
    
    df.warming <- load_warming(dir, scen=scen) %>% filter(ID %in% df$ID)
    
    post.proj <- format_post_proj_baseline(post.ids, proj.yrs, base, proj, df.warming, TT)
    base <- post.proj$base
    proj <- post.proj$proj
    get_damages(m, type, base, 
                proj, post.ids, TT, lags, df.base, yr.ids, yrs, 
                uncertainty=F) %>% 
      pluck("central") %>%
      mutate(scen=scen)
  }
)

Dt <- map_dfr(
  c("median", "p10", "p90"), 
  function(scen){
    
    load_warming(dir, scen=scen) %>% 
      filter(ID %in% df$ID, year == 2100) %>% select(ID, warming) %>% 
      left_join(df.base) %>% 
      mutate(scen=scen)
  }
) %>% 
  summarize(warming = weighted.mean(warming, pop), .by=scen)

pdf.rcp %>% 
  ggplot() + 
  geom_line(aes(x = year, y= damage, color = scen)) + 
  scale_color_viridis_d() 

pdf.rcp %>% 
  filter(year == 2100) %>% 
  left_join(Dt) %>% 
  ggplot() + 
  geom_smooth(aes(x = warming, y = damage), method='lm', se=F, 
              color = 'grey', alpha=.5) +
  geom_point(aes(x = warming, y = damage, color=scen)) +
  geom_point(aes(x = warming, y = damage), data = pdf %>% 
               filter(year == 2100))+
  labs(x="Warming (C)", y="Damage (% of GDP)") + 
  theme(legend.position = "bottom")
