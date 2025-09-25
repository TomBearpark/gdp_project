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

# shp %>% ggplot() + geom_sf()

# SPECIFY PARAMETERS ------------------------------------------------------

max_lags     <- 10
type     <- "levels"
toPlot   <- c("USA", "CHN", "SDN")
# Number of bootstraps from statistical uncertainty
Ndraws   <- 250

# GET BASELINES -----------------------------------------------------------


get.opt <- function(m, var = 'temp'){
  b1 <- coef(m)[paste0(var, '1')] 
  b2 <- coef(m)[paste0(var, '2')]
  tibble(
    opt = -b1/ (2 * b2), 
    type = if_else(b2 < 0, "max", "min")
    )
}

df.in <- load_historic_data(paste0(dir, "/data/historic/"), lags=max_lags)


df.reg <- df.in %>% 
  filter(!is.na(g), !is.na(temp1))

clim <- df.reg %>% summarize(temp = mean(temp1, na.rm=T), 
                             g = mean(g, na.rm=T), 
                             y = mean(y, na.rm=T), 
                             tmin = min(temp1, na.rm=T),
                             tmax = max(temp1, na.rm=T), 
                             .by = ID) %>% 
  mutate(trange = tmax - tmin)

df.reg$decade <- 10*(floor(df.reg$year/10))
m <- feols(
  g ~ temp1 + temp2 | ID +time1 + ID[time1], df.reg, 
  cluster = "ID"
)


# everyone different ------------------------------------------------------
m1 <- feols(
  g ~ i(ID, temp1) + i(ID, temp2) | time1 + ID + ID[time1], 
  df.reg, 
  cluster = 'ID'
)
get.opt(m, var = 'temp')

c.opt <- map_dfr(
  unique(df.reg$ID), function(ID) {
    get.opt(m1, paste0('ID::', ID, ':temp')) %>% 
      mutate(ID = !!ID)
  }
) %>% 
  left_join(clim) 

c.opt %>% 
  ggplot(aes(x = temp, y = opt)) + 
  geom_point() + 
  geom_label(aes(label = ID)) + 
  facet_wrap(~type) + 
  geom_abline(color = 'red')

feols(opt ~ temp + y + g, c.opt)

c.opt %>% 
  mutate(dist = opt -temp) %>%
  ggplot() + 
  geom_point(aes(x = trange, y =dist)) + 
  ylim(-4, 4)


df.reg %>% 
  filter(ID == "STP")

# decades -----------------------------------------------------------------
m2 <- feols(
  g ~ i(decade, temp1) + i(decade, temp2) | ID + time1 + ID[time1], df.reg
)

d.opt <- map_dfr(
  unique(df.reg$decade), function(decade) {
    tibble(decade = !!decade, 
           opt = get.opt(m2, paste0('decade::', decade, ':temp')))    
  }
)  %>% 
  left_join(summarize(df.reg, temp = mean(temp1, na.rm=T), .by = decade))

d.opt %>% ggplot() + geom_point(aes(x = temp, y = opt))


# continent  -----------------------------------------------------------------
m3 <- feols(
  g ~ i(region, temp1) + i(region, temp2) | region+ID + time1 + ID[time1], df.reg
)

r.opt <- map_dfr(
  unique(df.reg$region), function(region) {
    tibble(region = !!region, 
           opt = get.opt(m3, paste0('region::', region, ':temp')))    
  }
)  %>% 
  left_join(summarize(df.reg, temp = mean(temp1, na.rm=T), .by = region))

r.opt %>% ggplot() + geom_point(aes(x = temp, y = opt)) + 
  geom_label(aes(x = temp, y = opt, label = region))




# everyone different ------------------------------------------------------
m1 <- feols(
  g ~ temp1 + temp2 , 
  df.reg, 
  split = ~ID
)
clim
s.opt <- map_dfr(seq_along(m1), 
        function(ii){
          m <- m1[[ii]]
          tibble(ID = m$model_info$sample$value, 
                 opt = get.opt(m))
        }
        ) %>% 
  left_join(clim)

s.opt %>% 
  ggplot() + 
  geom_point(aes(x = temp, y = opt)) + 
  geom_label(aes(x = temp, y = opt, label = ID)) + 
  geom_abline()


m2 <- feols(
  g ~ i(decade, temp1) + i(decade, temp2) , 
  df.reg, 
  split = ~ID
)
clim
s2.opt <- map_dfr(seq_along(m1), 
                 function(ii){
                   m <- m2[[ii]]
                   map_dfr(
                     unique(df.reg$decade), 
                     function(decade){
                       tibble(ID = m$model_info$sample$value, 
                              decade = decade, 
                              opt = get.opt(m, paste0('decade::', decade, ':temp')))    
                     }
                     
                   )
                 }
) %>% 
  left_join(clim)

s2.opt %>% 
  ggplot() + 
  geom_point(aes(x = temp, y = opt)) + 
  geom_label(aes(x = temp, y = opt, label = ID)) + 
  geom_abline()

s2.opt %>% 
  filter(between(opt, 0, 30)) %>% 
  filter(ID %in% c("USA", "CHN", "RUS", "SDN")) %>% 
  ggplot() + 
  geom_point(aes(x = decade, y = opt)) + 
  facet_wrap(~ID, scales='free')



# check non-parametrics ---------------------------------------------------

clim <- clim %>% mutate(trank = rank(temp))
predict_poly_het(m1, df = df.reg, "ID", fix_psd=F, infer.range=T, step.length=.1) %>% 
  left_join(
    select(c.opt, ID, opt, opt.type = type)
  ) %>% 
  rename(tt = temp) %>% 
  left_join(clim) %>% 
  ggplot() + 
  geom_line(aes(x = tt, y = response+trank/20, 
                group=ID, 
                color = opt.type), alpha=.5)

