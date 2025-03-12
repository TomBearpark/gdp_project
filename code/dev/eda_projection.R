if(!require(pacman)) install.packages('pacman')
pacman::p_load(tidyverse, fixest, haven, marginaleffects, broom, useful, 
               readxl)
theme_set(theme_classic())

dir  <- "~/Library/CloudStorage/Dropbox/gdp-temp/"
code <- "~/Documents/GitHub/gdp_project/"

set.seed(123)
source(paste0(code, "code/0_funcs.R"))

# SPECIFY PARAMETERS ------------------------------------------------------

type     <- "growth"
lags     <- 0
spec     <- "poly2"
warming  <- 4
max_lags <- 15

toPlot  <- c("RUS", "USA", "SDN")

# LOAD DATA ---------------------------------------------------------------

df.reg <- load_historic_data(paste0(dir, "/replication/"), lags=max_lags) 

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
type <- "levels"
lags <- 10
m  <- run_reg(df.reg, type=type, lags=lags, spec=spec, global = F, 
              FE = "ID+time1+ID[time1]")

coefplot(m, keep="temp")
get_me(m, lags, type=type)

get_me_cumulative(m, lags, type=type) + 
  geom_vline(data = plt.temps, aes(xintercept=temp, color=ID))+
  geom_vline(data = plt.temps, aes(xintercept=temp_proj, color=ID), linetype=2)

# GET BASELINES -----------------------------------------------------------
yrs <- get_years()

# Filter to period we are going to use in plotting / getting baseline
df <- df.reg %>% 
  filter(year >= min(yrs$pre))  %>% 
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

# Load warming data
df.warming <- load_warming(dir) %>% 
  filter(ID %in% df$ID)

# Plot the warming 
total.warming <- df.warming %>% filter(year %in% c(2019, 2100)) 
total.warming %>% filter(year == 2100) %>% ggplot() + geom_histogram(aes(x = warming))

left_join(df %>% filter(year == 2019), total.warming) %>% 
  ggplot() + geom_point(aes(x = temp1, y = temp)) + geom_abline()

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
post.proj <- format_post_proj_baseline(yrs, base, proj, df.warming)
base <- post.proj$base
proj <- post.proj$proj

## Sense check the data 
plotdf_mats("USA", base, proj)

# PROJECTION --------------------------------------------------------------
df.pop <- df.base %>% select(ID, pop)

damages <- get_damages(m=m, 
                       type=type, 
                       base=base, 
                       proj=proj,
                       yrs=yrs, 
                       lags=lags, 
                       df.pop=df.pop, 
                       uncertainty=T)

plotdf_mats("USA", damages$base, damages$proj)
plot_global_damages(damages$central, damages$uncert)

# plotting ----------------------------------------------------------------
base <- damages$base
proj <- damages$proj
plotdf_mats(toPlot, damages$base, damages$proj, plt=T)

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
