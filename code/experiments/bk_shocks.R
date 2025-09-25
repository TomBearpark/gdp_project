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


# get BK shocks -----------------------------------------------------------
df <- paste0("/Users/tombearpark/Library/CloudStorage/Dropbox/gdp-temp/replication/", 
              'bk_micc_replication/data/micc_data.dta') %>% 
  read_dta() %>%
  filter(year >= 1960 & year <= 2019) %>% 
  rename(y= lnrgdppc_world_pwt, 
         x = gtmp_noaa_aw_dtfe2s, 
         xraw = gtmp_noaa_aw, 
         w = recessiondates) %>% 
  mutate(id = "world", 
         dy = y - lag(y)) %>% 
  arrange(year)

extract_bk_coefs <- function(M){
  map_dfr(seq_along(M), 
          function(ii){
            broom::tidy(M[[ii]], conf.int=T, conf.level = .9) %>%
              mutate(h = ii-1)
          } 
  )
}
horizon <- 10

M.y <- feols(f(y, 0:horizon) - l(y, 1) ~ 
               l(x, 0:2)  + l(dy, 1:2) + l(w, 0:2), 
             df, panel.id = c("id", "year"), vcov = "hetero")

M.t <- feols(f(xraw, 0:horizon) - l(xraw, 1) ~  
               l(x, 0:2) + l(w, 0:2), 
             df, panel.id = c("id", "year"), vcov = "hetero")

pdf.y <-  extract_bk_coefs(M.y) %>% 
  filter(str_detect(term,  "x") & !str_detect(term, "\\(")) %>% 
  mutate(var = "GDP")

pdf.t <-  extract_bk_coefs(M.t) %>% 
  filter(str_detect(term,  "x") & !str_detect(term, "\\(")) %>% 
  mutate(var = "Temp")

t.shock <- pdf.t$estimate

fig.1 <- ggplot() + geom_line(data = pdf.t, aes(x = h, y = estimate)) + 
  geom_hline(yintercept=0) + 
  ggtitle("Fig 1 - Temperature shock persistence")

# SPECIFY PARAMETERS ------------------------------------------------------

lags     <- 10
type     <- "levels"
toPlot   <- c("USA", "CHN", "SDN")
# Number of bootstraps from statistical uncertainty
Ndraws   <- 20
max_lags <- 10

# GET BASELINES -----------------------------------------------------------

yrs <- get_years(base_start_year = 1990, end_data_year = 2019, 
                 start_proj_year = 2020, end_proj_year = 2050)

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

df.clim <- expand_grid(ID = df.base$ID, 
                       year = 2020:2050) %>% 
  left_join(select(df.base, ID, temp))

df.clim <- map_dfr(unique(df.clim$ID), 
        function(ii){
          df.ii <- df.clim %>% filter(ID == ii) 
          df.ii$temp[2:(length(t.shock)+1)] <- df.ii$temp[1:length(t.shock)] + t.shock
          df.ii
        }
)

# Some checks to make sure its in the right shape  
stopifnot(all(sort(unique(df.base$ID)) %in% sort(unique(df.clim$ID)) ))

df.clim <- df.clim %>% 
  filter(ID %in% unique(df.base$ID)) %>% 
  mutate(warming = temp-first(temp), 
         .by = c("ID"))

df.clim %>% 
  filter(ID %in% toPlot) %>% 
  ggplot() + geom_line(aes(x = year, y = warming, group=ID))  + 
  facet_wrap(~ID, scales='free')

# PROJECTION --------------------------------------------------------------
type <- "levels"
lags <- 10
df.pop <- df.base %>% select(ID, pop)
m      <- run_reg(df.reg, type=type, lags=lags)
get_me_cumulative(m, lags,name='temp', type=type)

fig.2 <- get_me(m, name='temp', lags=lags, type=type, plot=F) %>% 
  filter(temp %in% c(10, 20, 30)) %>% 
  mutate(temp = paste0(temp, "C")) %>% 
  ggplot() + 
  geom_point(aes(x = lag, y = estimate, color = temp)) + 
  geom_errorbar(aes(x = lag, ymin = conf.low, ymax = conf.high, color = temp), width=0.2) +
  geom_hline(yintercept=0) +
  ggtitle("Fig 2 - Marginal effects of temperature on growth") + 
  facet_wrap(~temp) + 
  theme(legend.position = "none") 

draws  <- MASS::mvrnorm(n = Ndraws, mu = coef(m), Sigma = vcov(m))

post.proj <- format_post_proj_baseline(yrs, base, proj, df.clim)
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
  uncertainty = FALSE,
  reduce_uncert = FALSE, 
  draws = draws
)

df.base %>% arrange(temp)
toPlot <- c("GBR", "MEX", "MRT")
fig.3 <- plotdf_mats(toPlot, result$base, result$proj, F) %>% 
  filter(year >2019) %>% 
  plot_proj() + 
  ggtitle("Fig 3. Projected damages from Fig 1 temp shock")

fig.1
fig.2
fig.3
plotdf_mats(toPlot, result$base, result$proj, F)  %>% 
  filter(year>2009) %>%
  filter(name=='g') %>% 
  pivot_wider(names_from = scen, values_from = value) %>% 
  mutate(damage = 100*(proj - base)) %>% 
  ggplot() + 
  geom_line(aes(x = year, y = damage)) + 
  facet_wrap(~ID)
