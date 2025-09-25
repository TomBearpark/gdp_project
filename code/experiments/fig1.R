# YOU NEED THIS PACKAGE: install it the first time you use this code
# install.packages("devtools")
# devtools::install_github("TomBearpark/useful")

if(!require(pacman)) install.packages('pacman')
pacman::p_load(MASS, 
               tidyverse, 
               fixest,
               haven,
               marginaleffects,
               broom,
               useful,
               readxl,
               patchwork, 
               rnaturalearth
               )
theme_set(theme_classic())
set.seed(1)
select <- dplyr::select
# get user info and set directory locations
user <- Sys.info()[["user"]]
if(user == "tombearpark"){
  # dir  <- "~/Library/CloudStorage/Dropbox/gdp-temp/"  
  dir <- '~/Library/CloudStorage/Dropbox/damage_uncertainty/'
  code <- "~/Documents/GitHub/gdp_project/"
}else if(user == "jordan"){
  stop("ADD PATH TO CODE AND DROPBOX")
}

dir.out <- paste0(dir, "/out/projection/")

set.seed(123)
source(paste0(code, "code/0_funcs.R"))

shp <- ne_countries(scale = "large", returnclass = "sf") %>% 
  filter(!is.na(iso_a3)) %>% 
  select(ID = iso_a3_eh, name)

shp %>% ggplot() + geom_sf()

# SPECIFY PARAMETERS ------------------------------------------------------

max_lags     <- 10
type     <- "levels"
toPlot   <- c("USA", "CHN", "SDN")
# Number of bootstraps from statistical uncertainty
Ndraws   <- 250

# GET BASELINES -----------------------------------------------------------

yrs <- get_years(base_start_year = 1990, end_data_year = 2019, 
                 start_proj_year = 2020, end_proj_year = 2095)

# Load regression data 
df.reg <- load_historic_data(paste0(dir, "/data/historic/"), lags=max_lags)

# Filter to period we are going to use in plotting / getting baseline. 
# Should interpolate missing values to get full panel?

vars_to_fill <- c("g", "temp1", "pop", "y")

df <- df.reg %>% 
  filter(year >= min(yrs$pre)) %>% 
  group_by(ID) %>% 
    filter(!if_any(all_of(vars_to_fill), ~ all(is.na(.)))) %>% 
    mutate(across(all_of(vars_to_fill), ~ 
                    if_else(is.na(.), mean(., na.rm = TRUE), .))) %>% 
  ungroup()

# Get some averages which are used to define counter-factual
df.base <- df %>% 
  summarize(g = mean(g, na.rm=T), 
            temp = mean(temp1, na.rm=T), 
            pop = mean(pop, na.rm=T),
            .by = ID) %>% 
  select(ID, temp, g, pop) %>% 
  mutate(g = if_else(g > 0.03, 0.03, g))

stopifnot(nrow(df.base)==nrow(na.omit(df.base)))
NN <- length(unique(df.base$ID))
print(NN)

# Create matrices to store data
base <- gen_mats(NN, yrs, df.base$ID)
proj <- gen_mats(NN, yrs, df.base$ID)

## Format pre-projection data -----------------
pre.proj <- format_pre_proj(yrs, df, base, proj)
base <- pre.proj$base
proj <- pre.proj$proj

# Get climate projections -------------------------------------------------

files <- list.files(path = paste0(dir, "data/cmip6_countrylevel/"), 
                    pattern = "*countrytemp.csv", full.names = F)
mods <- str_remove_all(files, "cmip6_|_countrytemp.csv") 

# Only use climate models that go to 2099

df.clim <- 
  map_dfr(files, 
    function(ff){
      print(ff)
      df.out <- read_csv(paste0(dir, "data/cmip6_countrylevel/", ff)) 
      if(! any(str_detect(names(df.out), "2099"))) return(tibble())
      
      df.out %>% 
        # Australia has multiple values
        filter(!(iso3 == "AUS" & country_name != "Australia")) %>% 
        
        select(iso3, all_of(paste0("temp_", yrs$proj))) %>% 
        filter(iso3 %in% df$ID) %>%
        pivot_longer(cols = starts_with("temp_"), 
                     names_to = "year", 
                     values_to = "temp")  %>% 
        mutate(model = str_remove_all(ff, "cmip6_|_countrytemp.csv")) %>%
        rename(ID = iso3) %>% 
        mutate(temp_proj = temp -273.15, 
               year = str_remove(year, "temp_") %>% as.numeric()) 
      }
    )  %>% 
  arrange(model, ID, year) %>% 
  group_by(model, ID) %>% 
    mutate(warming = temp_proj -first(temp_proj)) %>% 
  ungroup()
  
df.clim  %>% 
  filter(ID %in% toPlot) %>% 
  ggplot() + 
  geom_line(aes(x = year, y = warming, color = model, group=model)) + 
  facet_wrap(~ID, scales='free')

# Some checks to make sure its in the right shape  
stopifnot(all(sort(unique(df.base$ID)) %in% sort(unique(df.clim$ID)) ))

df.clim <- df.clim %>% 
  filter(ID %in% unique(df.base$ID)) %>% 
  arrange(model, ID, year) %>% 
  select(ID, model, year, warming)

# PROJECTION --------------------------------------------------------------

lags <- 0
df.pop <- df.base %>% select(ID, pop)
m      <- run_reg(df.reg, type=type, lags=lags)
models <- unique(df.clim$model)

draws  <- MASS::mvrnorm(n = Ndraws, mu = coef(m), Sigma = vcov(m))

# Loop over models, saving outputs
pdf.central <- list()
pdf.uncert  <- list()

for (i in seq_along(models)) {

  ww <- models[i]
  
  df.warming <- df.clim %>% filter(model == ww)
  post.proj <- format_post_proj_baseline(yrs, base, proj, df.warming)
  base <- post.proj$base
  proj <- post.proj$proj
  
  result <- get_damages(
    m = m, 
    type = type, 
    base = base, 
    proj = proj,
    yrs = yrs, 
    lags = lags, 
    df.pop = df.pop, 
    uncertainty = TRUE,
    reduce_uncert = FALSE, 
    draws = draws
  )
  
  pdf.central[[i]] <- result$central %>% mutate(model = ww)
  pdf.uncert[[i]]  <- result$uncert %>% mutate(model = ww)
  
  rm(result)
  gc()
}

pdf.central <- bind_rows(pdf.central)
pdf.uncert  <- bind_rows(pdf.uncert)
gc()

# Plotting ----------------------------------------------------------------

pdf.central %>% 
  ggplot() + 
  geom_line(aes(x = year, y= damage, color = model, group=model)) + 
  scale_color_viridis_d() 

uncert.2095 <- pdf.uncert %>% filter(year == max(yrs$proj)) %>% 
  left_join(df.pop) %>% 
  group_by(model, draw, year) %>% 
  summarise(damage=weighted.mean(damage, pop)) %>% 
  ungroup()

cc <- pdf.uncert %>% 
  left_join(df.pop) %>% 
  group_by(year, draw, model) %>%
  summarize(damage = weighted.mean(damage, pop)) %>% 
  group_by(year) %>% 
  add_q() %>% 
  mutate(year= as.numeric(year))

bpdf <- 
  bind_rows(
    uncert.2095 %>% 
      mutate(year = 2100, Uncertainty="Both") %>% 
      group_by(year, Uncertainty) %>% 
      add_q(), 
    pdf.central %>% 
      filter(year == 2095) %>% 
      mutate(year = 2105, Uncertainty = "Climate") %>% 
      group_by(year, Uncertainty) %>% 
      add_q(), 
    uncert.2095 %>% 
      group_by(draw) %>% 
      summarize(year = unique(year), damage = mean(damage))  %>% 
      mutate(year = 2110, Uncertainty = "Statistical") %>% 
      group_by(year, Uncertainty) %>% 
      add_q()
    )

ts1 <- ggplot() + 
  geom_hline(yintercept=0, linetype='dashed', color='grey') +
  geom_line(data = 
              pdf.central %>% 
              group_by(year) %>% summarize(damage = mean(damage)) %>% 
              filter(year >= 2020), 
            aes(x = year, y = damage)) + 
  geom_ribbon(data = cc, aes(x = year, ymin=q01, ymax = q99), 
              alpha=.1) + 
  geom_ribbon(data = cc, aes(x = year, ymin=q025, ymax = q975), 
              alpha=.2) + 
  geom_ribbon(data = cc, aes(x = year, ymin=q05, ymax = q95), 
              alpha=.3) + 
  geom_ribbon(data = cc, aes(x = year, ymin=q25, ymax = q75),
              alpha=.4) + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1, scale=1)) + 
  xlab("") + 
  ylab("% global gdp-pc") + 
  geom_errorbar(
    data = bpdf, 
    aes(x = year, ymin = q01, ymax = q99, color=Uncertainty), 
    width=4, alpha=.7
  ) + 
  geom_errorbar(
    data = bpdf, 
    aes(x = year, ymin = q025, ymax = q975, color=Uncertainty), 
    width=3, alpha=.8
  ) + 
  geom_errorbar(
    data = bpdf, 
      aes(x = year, ymin = q05, ymax = q95, color=Uncertainty), 
      width=2, alpha=.9
  ) + 
  geom_errorbar(
    data = bpdf, 
    aes(x = year, ymin = q25, ymax = q75, color=Uncertainty), 
    width=1, alpha=1
  )

ts1

# MAPS

errors <- uncert.2095 %>% 
  group_by(model, year, draw) %>% 
  summarize(damage = weighted.mean(damage, pop)) %>% 
  group_by(model)%>% 
  add_q()

total.error <- uncert.2095 %>% 
  group_by(model, draw) %>% 
  summarize(damage = weighted.mean(damage, pop)) %>% 
  ungroup() %>% 
  add_q()

damages.2095 <- pdf.central %>% 
  filter(year == max(yrs$proj)) %>% 
  arrange(desc(damage)) %>%
  mutate(model = fct_reorder(model, damage)) %>%
  ggplot() + 
  geom_vline(xintercept=0, linetype='dashed', color='grey') +
  geom_point(aes(y = model, x = damage)) + 
  geom_errorbar(data = errors, aes(y = model, xmin = q05, xmax = q95), alpha=.3) + 
  geom_errorbar(data = errors, aes(y = model, xmin = q25, xmax = q75))

# variance map ------------------------------------------------------------
map.df <- uncert.2095 %>% 
  # filter(model == "access_cm2") %>% 
  group_by(ID) %>% 
  add_q()

map <- left_join(shp, map.df) %>% 
  # filter(sd < 20) %>% 
  ggplot() + 
  geom_sf(aes(fill = log(sd)), color=NA) + 
  scale_fill_viridis_c()
map
 
df.base %>% 
  left_join(map.df) %>% 
  filter(pop > 1000000) %>% 
  ggplot() + 
  geom_smooth(aes(x = temp, y = q75-q25), color='black', se=F)   +
  geom_point(aes(x = temp, y =q75-q25, size = pop),alpha=.6) 

df.clim %>% 
  filter(year %in% c(2020, 2090)) %>% 
  pivot_wider(names_from = year, values_from = warming)  %>% 
  mutate(warming = `2090` - `2020`)  %>% 
  group_by(ID) %>% 
  summarize(warming = mean(warming)) %>% 
  left_join(map.df) %>% 
  # filter(pop > 1000000) %>% 
  ggplot() + 
  geom_smooth(aes(x = log(warming), y = log(q75-q25)), color='black', se=F)   +
  geom_point(aes(x = log(warming), y = log(q75-q25)),alpha=.6) 



# Some countries only 
plot.ids <- c("USA", "RUS", "SDN")
map_dfr(seq_along(models), 
        function(ww){
          base <- pdf[[ww]]$base
          proj <- pdf[[ww]]$proj
          plotdf_mats(plot.ids, base, proj, F) %>% 
            mutate(model = models[[ww]])
        })  %>% 
  filter(year == max(yrs$proj), scen=='proj', name == "dam") %>% 
  ggplot() + 
  geom_point(aes(x = model, y = value)) + 
  facet_wrap(~ ID, scales='free') 

# ASSEMBLE PLOT -----------------------------------------------------------

# Panel 1
# - Climate projection uncertainty - time series 
# Panel 2
# - Estimated damage function
# Panel 3
# - Global GDP loss time series 
# Panel 4
# - Map across the globe - uncertainty in GDP loss in 2100
# - Uncertainty likely highest for countries near the peak

p1  

p2 <- get_me_cumulative(m, lags, type=type)

(p1 + p2) / (damages.2095 +  map)



# extra plots -------------------------------------------------------------

toPlot <- c("RUS", "CHN", "USA", "SDN")
p1 <- ggplot() + 
  geom_line(data = df.clim %>% filter(ID %in% toPlot), 
            aes(x = year, y = warming, group=model),alpha=.2) + 
  geom_line(
    data = df.clim %>% 
      filter(ID %in% toPlot) %>% 
      summarise(warming = mean(warming), .by=c(year, ID)), 
    aes(x = year, y = warming), linewidth=1
  ) + 
  facet_wrap(~ID, scales='fixed') + 
  geom_hline(yintercept=0, linetype='dashed', color='grey') +
  theme(legend.position='none')+
  labs(title = "Climate projection uncertainty")
