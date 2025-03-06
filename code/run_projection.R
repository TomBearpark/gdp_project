if(!require(pacman)) install.packages('pacman')
pacman::p_load(tidyverse, fixest, haven, marginaleffects, broom, useful)
theme_set(theme_classic())

dir  <- "~/Library/CloudStorage/Dropbox/gdp-temp/replication/"
code <- "~/Documents/GitHub/gdp_project/"

set.seed(123)
source(paste0(code, "code/0_funcs.R"))

# SPECIFY PARAMETERS ------------------------------------------------------

type    <- "growth"
lags    <- 0
spec    <- "poly2"
warming <- 2

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

b0 <- coef(m)[sort(str_subset(names(coef(m)), "temp1"))] %>% as.matrix()
b1 <- coef(m)[sort(str_subset(names(coef(m)), "temp2"))] %>% as.matrix()

if(type == "levels"){
  for(tt in 2:TT) {
    base$tcons[,tt] <- base$temp[,tt]-base$temp[,tt-1]
    proj$tcons[,tt] <- proj$temp[,tt]-proj$temp[,tt-1]  
  }
}else{
  base$tcons <- base$temp
  proj$tcons <- proj$temp
}

for(tt in post.ids){
  lls <- tt-0:lags
  
  delta <- (proj$tcons[, lls]   - base$tcons[, lls]) %*% b0 + 
           (proj$tcons[, lls]^2 - base$tcons[, lls]^2) %*% b1
    
  # -- Growth
  proj$g[,tt] <- base$g[, tt] + delta
  
  # -- GDP
  proj$y[,tt] <- (1 + proj$g[, tt]) * proj$y[,tt-1]
    
}

plotdf_mats(toPlot, base, proj, plt=T)

pdf.all <- plotdf_mats(df.base$ID, base, proj, plt=F)%>% 
  left_join(df.base) %>% 
  filter(name == 'y') %>% 
  pivot_wider(names_from = scen, values_from = value) %>% 
  mutate(pc_dam = proj - base, 
         tot_dam = pc_dam * pop) 

pdf.all %>% pivot_longer(cols = c(pc_dam, tot_dam), 
                         names_to = 'type', values_to = 'dam') %>%
  ggplot() + 
  geom_line(aes(x = year, y = dam, group = ID), alpha=.4)+
  facet_wrap(~type, scales='free')

pdf.all %>% 
  filter(year == 2100) %>% 
  pivot_longer(cols = c(pc_dam, tot_dam), 
               names_to = 'type', values_to = 'dam') %>% 
  ggplot() + geom_point(aes(x = temp, y = dam, color =base))+
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

dam.out <- map_dfr(
  yr.ids, 
  function(tt){
    pc.diff    <- 100*(proj$y[,tt]-base$y[,tt]) / base$y[,tt]
    tot.damage <- weighted.mean(pc.diff, df.base$pop)
    # glob.gdp   <- sum(proj$y[,tt]*df.base$pop)
    tibble(year = yrs[tt], 
           damage = tot.damage)
  }
)

dam.out %>% 
  ggplot() + 
  geom_vline(xintercept=max(pre.yrs)) + geom_line(aes(x = year, y = damage))
