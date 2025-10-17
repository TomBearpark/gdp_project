if(!require(pacman)) install.packages('pacman')
pacman::p_load(tidyverse, fixest, haven, marginaleffects, broom, useful, readxl, 
               arrow)
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
df.reg <- load_historic_data(dir, lags=max_lags, 
                             max.year = 2019, old=T)

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
warming <- load_warming(paste0(dir, "data/projection/"), "median",
                        file= 'rcp70') %>% filter(ID %in% df$ID)
post.proj <- format_post_proj_baseline(yrs, base, proj, warming)
base <- post.proj$base
proj <- post.proj$proj

# PROJECTION --------------------------------------------------------------

opts <- expand_grid(lags=0:max_lags, type=c("levels", "growth"))
df.pop <- df.base %>% select(ID, pop)

pdf <- pmap(
  opts,
  function(lags, type){
    print(paste0(type, ", lags=", lags))
    # browser()
    m  <- run_reg(df.reg, type=type, lags=lags)
    
    # Save ME info
    get_me_cumulative(m, lags,name='temp', type=type)
    ggsave(paste0(dir.out, type, "_", lags, "_me_cum", ".pdf"), 
           height = 10, width = 10)
    get_me(m, name='temp', lags=lags, type=type)
    ggsave(paste0(dir.out, type, "_", lags, "_me_sep",".pdf"), 
           height = 10, width = 10)
    
    # Project damages
    dam.df <- get_damages(m=m, 
                          type=type, 
                          base=base, 
                          proj=proj,
                          yrs=yrs, 
                          lags=lags, 
                          df.pop=df.pop, 
                          uncertainty=T, 
                          reduce_uncert = T)
    
    plot_global_damages(dam.df$central, dam.df$uncert) + 
      ggtitle(paste0(type, ", ", lags, " lags"))
    ggsave(paste0(dir.out, type, "_", lags, "_damages",".pdf"), 
           height = 10, width = 10)
    dam.df
  }
)

central <- map_dfr(seq_along(pdf), function(ii) pdf[[ii]]$central) %>% 
  mutate(lags = as.factor(lags))

uncert <- map_dfr(seq_along(pdf), function(ii) pdf[[ii]]$uncert) %>% 
  mutate(lags = as.factor(lags))

# PLOTS -------------------------------------------------------------------

plot_global_damages(central %>% filter(lags == 0), 
                    uncert %>% filter(lags == 0)) + 
  facet_wrap(~type+lags, nrow=1, scales='free')

plot_global_damages(central, 
                    uncert) + 
  ylim(c(-100, 100))+
  facet_wrap(~type+lags, nrow=2, scales='fixed')

  

plot_global_damages(central %>% filter(type == "levels", lags %in% c(0, 5, 10)), 
                    uncert%>% filter(type == "levels", lags %in% c(0, 5, 10))) + 
  facet_wrap(~lags, nrow=1) 


plot_global_damages(central %>% filter(type == "growth", lags %in% c(0, 5, 10)), 
                    uncert%>% filter(type == "growth", lags %in% c(0, 5, 10))) + 
  facet_wrap(~lags, nrow=1) 

ggplot() + 
  geom_hline(yintercept = 0, linetype="dashed", color = "black") +
  geom_point(aes(x = lags, y = damage, color = type), 
             data=central %>% filter(year == 2100), 
             position=position_dodge(width = .5)) +
  # geom_errorbar(aes(x = lags, ymin = q05, ymax=q95, color = type),
  #               position=position_dodge(width = .5), width=.5,
  #            data=uncert %>% filter(year == 2100)) +
  ylab("Damages, %GDP 2100") + 
  xlab("No. lags in model")

ggplot() + 
  geom_vline(xintercept = 2020, linetype="dashed", color = "black") +
  geom_line(data = central, aes(x=year, y = damage, color = lags)) + 
  facet_wrap(~type)+ 
  scale_color_viridis_d()+
  ylab("Damages, %GDP 2100") 
ggsave(paste0(dir.out, "1_combined.png"), 
       height = 4, width = 8)
