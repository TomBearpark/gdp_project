if(!require(pacman)) install.packages('pacman')
pacman::p_load(tidyverse, fixest, haven, marginaleffects, broom, useful)
theme_set(theme_classic())

dir  <- "~/Library/CloudStorage/Dropbox/gdp-temp/replication/"
code <- "~/Documents/GitHub/gdp_project/"

set.seed(123)
source(paste0(code, "code/0_funcs.R"))

# SPECIFY PARAMETERS ------------------------------------------------------

type    <- "levels"
lags    <- 10
spec    <- "poly2"
warming <- 4

toPlot  <- c("RUS", "CHN", "SDN")

# LOAD DATA ---------------------------------------------------------------

df.reg <- load_historic_data(dir)

plt.temps <- df.reg %>% 
  filter(ID %in% toPlot) %>% 
  summarize(temp = mean(temp1), .by=ID) %>% 
  mutate(temp_proj = temp + warming)

df.reg %>% 
  filter(year == 2019) %>% 
  ggplot() + 
  geom_histogram(aes(x = temp1)) + 
  geom_vline(data = plt.temps, aes(xintercept=temp, color=ID)) + 
  geom_vline(data = plt.temps, aes(xintercept=temp_proj, color=ID), linetype=2)

# RUN REGRESSION ----------------------------------------------------------

m  <- run_reg(df.reg, type=type, lags=lags, spec=spec)
get_me_sep(m, lags, type=type)
get_me(m, lags, type=type)

me.cum <- get_me_cum(m, lags, type=type)
me.cum + 
  geom_vline(data = plt.temps, aes(xintercept=temp, color=ID))+
  geom_vline(data = plt.temps, aes(xintercept=temp_proj, color=ID), linetype=2)

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
b0 <- extract_coefs(m, "temp1")
b1 <- extract_coefs(m, "temp2")

projected <- project(m, type, base, proj, post.ids, TT, lags)
base <- projected$base
proj <- projected$proj

# Get global damages 
dam.out <- global_damages(proj, base, df.base, yr.ids, yrs)


# plotting ----------------------------------------------------------------


plotdf_mats(toPlot, base, proj, plt=T)

pdf.all <- plotdf_mats(df.base$ID, base, proj, plt=F)%>% 
  left_join(df.base) %>% 
  filter(name == 'y') %>% 
  pivot_wider(names_from = scen, values_from = value) %>% 
  mutate(pc_dam = proj - base, 
         tot_dam = pc_dam * pop) 

pdf.all %>% 
  pivot_longer(cols = c(pc_dam, tot_dam), 
               names_to = 'type', 
               values_to = 'dam') %>%
  # filter(ID == "IND") %>% 
  ggplot() + 
  geom_line(aes(x = year, y = dam, group = ID), alpha=.4)+
  facet_wrap(~type, scales='free')

pdf.all %>% 
  filter(year == 2100) %>% 
  pivot_longer(cols = c(pc_dam, tot_dam), 
               names_to = 'type', values_to = 'dam') %>% 
  ggplot() + geom_point(aes(x = temp, y = dam, color = base))+
  # geom_label(aes(x = temp, y = dam, label = ID)) + 
  facet_wrap(~type, scales='free') 

pdf.all %>% 
  filter(year == 2100) %>% 
  mutate(rDam = rank(pc_dam), 
         rTemp = rank(temp)) %>%
  ggplot() + 
  geom_point(aes(x= rTemp, y = rDam, color=log(base))) + 
  scale_color_viridis_c()

# CALCULATE DAMAGES  ------------------------------------------------------
dim(proj$y)
length(df.base$pop)


dam.out %>% 
  ggplot() + 
  geom_vline(xintercept=max(pre.yrs)) + geom_line(aes(x = year, y = damage))
