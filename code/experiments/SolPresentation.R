if(!require(pacman)) install.packages('pacman')
pacman::p_load(tidyverse, fixest, haven, marginaleffects, broom, useful, readxl, 
               arrow)
theme_set(theme_classic())

dir  <- "~/Library/CloudStorage/Dropbox/gdp-temp/"
dir.out <- paste0(dir, "/presentations/norway/")
code <- "~/Documents/GitHub/gdp_project/"

set.seed(123)
source(paste0(code, "code/0_funcs.R"))

# SPECIFY PARAMETERS ------------------------------------------------------

max_lags <- 10
gen_data <- TRUE

if(gen_data){
for(data in c("old", "new")){
    
  if(data == 'old') old <- TRUE else old <- FALSE
  dir.out <- paste0(dir, "/presentations/norway/", data, "_data/")
  dir.create(dir.out, showWarnings = F, recursive = T)
  
  # GET BASELINES -----------------------------------------------------------
  
  yrs <- get_years()
  
  # Load regression data 
  df.reg <- load_historic_data(dir, lags=max_lags, max.year = 2019, old=old)
  
  # Filter to period we are going to use in plotting / getting baseline. 
  # Should interpolate missing values to get full panel?
  
  # Get old data so we have pop
  df <- load_historic_data(dir, lags=max_lags, max.year = 2019, old=T) %>% 
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
  warming   <- load_warming(paste0(dir, "data/projection/"), 
                            file='rcp70', 
                            "median") %>% 
    filter(ID %in% df$ID)
  
  post.proj <- format_post_proj_baseline(yrs, base, proj, warming)
  base <- post.proj$base
  proj <- post.proj$proj
  
  # TESTING ----------------------------------------------------------------
  # m  <- run_reg(df.reg, type='levels', lags=10, name = c("temp", 'prec'))
  # get_me_cumulative(m, 10, name='temp', type='levels')
  # m
  
  # GET results -----------------------------------------------------------------
  opts   <- expand_grid(lags=0:max_lags, type=c("levels", "growth"))
  df.pop <- df.base %>% select(ID, pop)
  
  # Run models
  MM <- pmap(opts, function(lags, type) 
    run_reg(df.reg, type = type, lags = lags, FE='ID+year+ID[year]')) |>
    set_names(pmap_chr(opts, \(lags, type) paste0(type, lags)))
  
  # Get cumulative marginal effects
  cum.ME <- pmap_dfr(
    opts,
    function(lags, type){
      m <- MM[[paste0(type, lags)]]
      get_me_cumulative(m, lags, name='temp', type=type, plot = F) %>%
        mutate(lags = lags, type = type)
    }
  )
  
  # Get separated marginal effects
  sep.ME <- pmap_dfr(
    opts,
    function(lags, type){
      m <- MM[[paste0(type, lags)]]
      get_me(m, name='temp', lags=lags, type=type, plot = F) %>%
        mutate(lags = lags, type = type)
    }
  )
  
  # Get projected damages 
  pdf <- pmap(
    opts,
    function(lags, type){
      m <- MM[[paste0(type, lags)]]
      get_damages(m=m, 
                  type=type, 
                  base=base, 
                  proj=proj,
                  yrs=yrs, 
                  lags=lags, 
                  df.pop=df.pop, 
                  uncertainty = TRUE, 
                  reduce_uncert = TRUE) 
    }
  )
  
  central <- map_dfr(seq_along(pdf), function(ii) pdf[[ii]]$central) %>% 
    mutate(lags_id = paste0(lags, " lag model")) %>% 
    mutate(type = if_else(type == "levels", "Levels", "Growth"))
  
  uncert <- map_dfr(seq_along(pdf), function(ii) pdf[[ii]]$uncert) %>% 
    mutate(lags_id = paste0(lags, " lag model")) %>% 
    mutate(type = if_else(type == "levels", "Levels", "Growth"))
  
  # Save outputs
  write_csv(sep.ME,  paste0(dir.out, "0_separated_marginal_effects.csv"))
  write_csv(cum.ME,  paste0(dir.out, "0_cumulative_marginal_effect.csv"))
  write_csv(central, paste0(dir.out, "0_projected_damages_central.csv"))
  write_csv(uncert,  paste0(dir.out, "0_projected_damages_uncert.csv"))
  }
}
  # PLOTS -------------------------------------------------------------------
for(data in c("old", "new")){
  dir.out <- paste0(dir, "/presentations/norway/", data, "_data/")
  sep.ME <- read_csv(paste0(dir.out, "0_separated_marginal_effects.csv"))
  cum.ME <- read_csv(paste0(dir.out, "0_cumulative_marginal_effect.csv"))
  central <- read_csv(paste0(dir.out, "0_projected_damages_central.csv"))
  uncert  <- read_csv(paste0(dir.out, "0_projected_damages_uncert.csv"))
  
  ## 1. 0, 1 and 5 lag models: projected damages -------
  plot_global_damages(central %>% filter(lags %in% c(0, 1, 5)), 
                      uncert %>% filter(lags %in% c(0, 1, 5))) + 
    facet_wrap(~type+lags_id, nrow=2, scales='fixed') + 
    xlab("Year") + ylab("Damages (%GDP)") 
  
  ggsave(paste0(dir.out, "1a_proj_damages_0_1_5lags.pdf"), 
         height = 6, width = 10)
  
  plot_global_damages(central %>% filter(lags %in% c(0, 1, 5), 
                                          type == "Levels"), 
                      uncert %>% filter(lags %in% c(0, 1, 5), 
                                        type == "Levels")) + 
    facet_wrap(~type+lags_id, nrow=1, scales='fixed') + 
    xlab("Year") + ylab("Damages (%GDP)") 
  
  ggsave(paste0(dir.out, "1b_levels_proj_damages_0_1_5lags.pdf"), 
         height = 4, width = 10)
  
  plot_global_damages(central %>% filter(lags %in% c(0, 1, 5), 
                                         type == "Growth"), 
                      uncert %>% filter(lags %in% c(0, 1, 5), 
                                        type == "Growth")) + 
    facet_wrap(~type+lags_id, nrow=1, scales='fixed') + 
    xlab("Year") + ylab("Damages (%GDP)") 
  
  ggsave(paste0(dir.out, "1c_levels_proj_damages_0_1_5lags.pdf"), 
         height = 4, width = 10)
  
  ## 2. box_plots  ---------------------------------------------
  
  uncert %>%
    left_join(central) %>% 
    filter(year == 2100) %>% 
    ggplot(aes(x = as.factor(lags))) +
    geom_hline(yintercept = 0, color = "grey") +
    geom_boxplot(
      aes(ymin = q05, lower = q25, middle = damage, upper = q75, ymax = q95),
      stat = "identity",
      alpha = 0.6
    ) +
    # geom_point(data = central %>% filter(year == 2100), 
    #            aes(x = as.factor(lags), y = damage), shape = 2)+
    facet_wrap(~type)+
    labs(y = "2100 Damages (%GDP)", 
         x = "Number of lags in model", 
         caption = paste0("Boxes show the interquartile range (25th-75th percentiles).", 
                          "\nWhiskers show the 5th-95th percentile range.", 
                          "\nCentral line shows the mean estimate.") )
  
  ggsave(paste0(dir.out, "2a_2100_damages_boxplots.pdf"), 
         height = 6, width = 10)
  
  uncert %>%
    left_join(central) %>% 
    filter(year == 2100, type == "Growth") %>% 
    ggplot(aes(x = as.factor(lags))) +
    geom_hline(yintercept = 0, color = "grey") +
    geom_boxplot(
      aes(ymin = q05, lower = q25, middle = damage, upper = q75, ymax = q95),
      stat = "identity",
      alpha = 0.6
    ) +
    labs(y = "2100 Damages (%GDP)", 
         x = "Number of lags in model", 
         caption = paste0("Boxes show the interquartile range (25th-75th percentiles).", 
                          "\nWhiskers show the 5th-95th percentile range.", 
                          "\nCentral line shows the median estimate.") )
  
  ggsave(paste0(dir.out, "2b_growth_2100_damages_boxplots.pdf"), 
         height = 6, width = 10)
  
  uncert %>%
    left_join(central) %>% 
    filter(year == 2100, type == "Levels") %>% 
    ggplot(aes(x = as.factor(lags))) +
    geom_hline(yintercept = 0, color = "grey") +
    geom_boxplot(
      aes(ymin = q05, lower = q25, middle = damage, upper = q75, ymax = q95),
      stat = "identity",
      alpha = 0.6
    ) +
    geom_point(data = central %>% filter(year == 2100, type == "Levels"), 
               aes(x = as.factor(lags), y = damage), shape = 2)+
    labs(y = "2100 Damages (%GDP)", 
         x = "Number of lags in model", 
         caption = paste0("Boxes show the interquartile range (25th-75th percentiles).", 
                          "\nWhiskers show the 5th-95th percentile range.", 
                          "\nCentral line shows the median estimate.") )
  
  ggsave(paste0(dir.out, "2c_levels_2100_damages_boxplots.pdf"), 
         height = 6, width = 10)
  
  ## 3. Cumulative ME plot --------------------------------------------
  # Cumulative marginal effects
  for(temp.val in c(5, 15, 25)){
    cum.ME %>% 
      mutate(type = if_else(type == "levels", "Levels", "Growth")) %>% 
      filter(temp == temp.val) %>% 
      ggplot() + 
      geom_hline(yintercept = 0, linetype="dashed", color = "black") +
      geom_line(aes(x = lags, y = estimate)) + 
      geom_ribbon(aes(x = lags, ymin = conf.low, ymax = conf.high), alpha = .2) + 
      facet_wrap(~type) + 
      ggtitle(paste0("Cumulative ME for ", temp.val, "C country")) + 
      xlab("No. lags in model") + 
      scale_x_continuous(breaks = 0:max_lags) +ylab("")
    
    ggsave(paste0(dir.out, "3_cumulative_me", temp.val, ".pdf"), 
           height = 4, width = 8)
  }
  
  # ## 4. Extra stuff --------------------------------------------
  
  plot_global_damages(central %>% filter(lags == 0), 
                      uncert %>% filter(lags == 0)) + 
    facet_wrap(~type+lags_id, nrow=1, scales='free') + 
    ylab("Damages (%GDP)") + xlab("Year")
  
  ggsave(paste0(dir.out, "4_0_lag_proj_damages.pdf"), 
         height = 4, width = 8)
  
  central %>% 
    mutate(lags_id = fct_reorder(lags_id, lags)) %>%
    ggplot() + 
      geom_vline(xintercept = 2020, linetype="dashed", color = "black") +
      geom_line(aes(x=year, y = damage, color = lags_id, group=lags_id)) + 
      facet_wrap(~type)+ 
      scale_color_viridis_d(name = "Model")+
      ylab("Damages, %GDP 2100") + xlab("Year")
  ggsave(paste0(dir.out, "4_1_proj_damages.pdf"), 
         height = 3.5, width = 8)
}
