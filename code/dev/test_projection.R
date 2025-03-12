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

df.glob <- 
  read_dta(
    file.path(dir, 
              'replication/bk_micc_replication/data/micc_data.dta')
    ) %>% 
  select(year, gtemp1 = gtmp_noaa_aw)  %>% 
  mutate(id = "global", gtemp2 = gtemp1^2, 
         dgtemp1 = lag(gtemp1), 
         dgtemp2 = lag(gtemp2)) %>% 
  add_lags(vars = c("gtemp","dgtemp"),  sort_df = F, max.p = 2, 
           lags=max_lags) %>% 
  select(-id)
  

# LOAD DATA ---------------------------------------------------------------

df.reg <- load_historic_data(paste0(dir, "/replication/"), lags=max_lags) %>% 
  left_join(df.glob)

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
m1  <- run_reg(df.reg, type=type, lags=lags, spec=spec, global = F, 
              FE = "ID+time1+ID[time1]")

m2  <- run_reg(df.reg, type=type, lags=lags, spec=spec, global = T, 
              FE = "ID+ID[time1]", cluster = c("ID", "time1"))
etable(m1, m2)
coefplot(m2, keep="gtemp")

m3 <- feols(g ~ 
              l(l0_dtemp1, 0:10) + 
              l(l0_dtemp2, 0:10) +
              l(l0_dgtemp1, 0:10)
            
            |
              
        ID + ID[time1]+ID[time2], 
        df.reg, 
        panel.id=c("ID", "year"))

tidy(m3) %>% 
  filter(str_detect(term, "gtemp")) %>% 
  mutate(cum = cumsum(estimate), lag=0:10) %>% 
  ggplot() + 
  geom_point(aes(x = lag, y = cum)) + 
  geom_hline(yintercept = 0)

m3
etable(m1, m2, m3)
coefplot(m2, keep="gtemp")



m
coefplot(m, keep="gtemp")

get_me_sep(m, lags, type=type)
get_me_sep(m, lags, type=type, name='gtemp')

get_me(m, lags, type=type)

me.cum <- get_me_cum(m, lags, type=type)
me.cum + 
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

# SCC
w_mean   <- damages$central %>% filter(year == 2100) %>% pull(damage)  
glob_gdp <- damages$central %>% filter(year == 2100) %>% pull(global_gdp)

-(glob_gdp / 1e6) * (1-exp( w_mean /100 )) / 873 * 0.001 * 0.98^40 / (1 - 0.98)


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
