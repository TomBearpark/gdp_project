# CONTENTS:

## 1. Running regressions:
  # loading historic data
  # running regression
  # plotting marginal effects

## 2. Prepping data for projection:
  # generating matrices to store projection data
  # formatting pre-projection data
  # formatting post-projection baseline
    # including loading warming, population, and GDP projections

## 3. Projecting:
  # extracting coefficients
  # getting global damages from country level projections
  # projecting: core function is `project()`
  # This is called by the `get_damages()` function, which 
    # calls `project()` and `global_damages()`, and bootstraps to get uncertainty
  # plotting global damages

## 4. Outputs:
  # plotting various stages of the analysis

# RUNNING REGRESSIONS -----------------------------------------------------

load_historic_data <- function(dir, lags=10, max.p=2){
  df.in <- read_rds(file.path(dir, "temp_gdp_world_panel.rds")) %>% 
    rename(temp1 = era_mwtemp, pop = SP.POP.TOTL) %>%
    group_by(ID = ISO3) %>% 
    arrange(ID, year) %>% 
    mutate(y = NY.GDP.PCAP.KD, 
           dy = log(y) - lag(log(y)), 
           temp2 = temp1^2, 
           time1 = year - 1960, 
           time2 = time1^2, 
           Tbar = mean(temp1, na.rm=TRUE), 
           dtemp1 = temp1-lag(temp1), 
           dtemp2 = temp2-lag(temp2)
    ) %>% 
    ungroup()  %>% 
    filter(year >= 1960 & year <= 2019)  %>% 
    arrange(ID, year) %>% 
    group_by(ID) %>% 
    useful::add_lags(vars = c("temp", "dtemp"), 
                     lags=lags, max.p=max.p, 
                     sort_df = FALSE) %>% 
    select(ID, year, time1, time2, y, g=dy, pop, contains("temp"))
}

get_var <- function(type, name="temp"){
  if(type=="levels") var <- paste0("d", name)
  else if(type == "growth") var <- name
  else stop("not implemented")
  return(var)
}

run_reg <- function(df, 
                    type="levels", 
                    lags=1, 
                    global=FALSE, 
                    cluster="ID", 
                    spec="poly2", 
                    FE = "ID + time1 + ID[time1]", 
                    return_ff=F){
  
  var <- get_var(type)
  poly_order <- as.numeric(str_extract(spec, "[0-9]"))
  
  
  if(global){
    control <- get_var(type, 'gtemp')
    if(lags==0) control <- paste0("l0_", control)
    
    poly_control <- 1
  }else{
    control <- poly_control <- NULL
  }
  if(lags==0) var <- paste0("l0_", var)
  
  ff <- useful::build_formula_poly(
    yvar = "g", treat = var, poly_treat = poly_order, control=control, 
    poly_control = poly_control, 
    leads = 0, lags=lags , FE=FE
  )
  if(return_ff) return(ff)
  feols(
    ff, 
    data = df, 
    cluster = cluster
  )
}

plot_me <- function(pdf, type, lags=NULL){
  
  p <- ggplot(pdf, 
              aes(temp, estimate)) + 
    geom_hline(yintercept=0, linetype="dashed") +
    geom_line() + 
    geom_ribbon(aes(temp, ymin = conf.low, ymax = conf.high), 
                alpha = .2) + 
    ggtitle(paste0(type, ", ", lags, " lags")) + xlab("Temp") + ylab("ME")
  
  if(!is.null(lags)) p <- p+facet_wrap(~lag) 
  
  p
}

get_me_sep <- function(m, lags, 
                       name="temp", 
                       type="growth", 
                       xrange=seq(0, 30, by = 5), 
                       id="", 
                       plot=T){
  
  var <- get_var(type=type, name=name)
  
  pdf <- map_dfr(
    0:lags, 
    function(lag){
      map_dfr(
        xrange, 
        function(temp){
          hypotheses(m, 
                     paste0(  "l", lag, "_", var, "1 + ", 
                            "2*l", lag, "_", var, "2*",temp," = 0")) %>% 
            broom::tidy() %>% 
            mutate(temp = !!temp, lag=!!lag)
          
        }
      ) %>% 
        mutate(id = !!id, lag=as.character(lag))
    }
  ) %>% 
    mutate(lag = fct_relevel(lag, as.character(0:lags)))
  if(plot){
    return(plot_me(pdf, type, lags))
  }else{
    return(pdf)
  }
}

get_me_cumulative <- function(m, lags, name="temp", 
                              type="growth", xrange=seq(0, 30, by = 5), id="", 
                              plot=T){
  
  var <- get_var(type, name=name)
  
  pdf <-
      map_dfr(
        xrange, 
        function(temp){
          ff <- ""
          for(lag in 0:lags){
            ff <- paste0(ff, 
                         paste0(  "l", lag, "_", var, "1 + ", 
                                "2*l", lag, "_", var, "2*",temp," + "))
          }
          hypotheses(m, 
                     paste0(ff, " = 0")) %>% 
            broom::tidy() %>% 
            mutate(temp = !!temp, lag="Cumulative")
          
        }
      ) %>% 
        mutate(id = !!id)
  if(plot){
    return(plot_me(pdf, type, lags))
  }else{
    pdf
  }
}

get_me <- function(m, lags, name='temp', 
                   type="growth", xrange=seq(0, 30, by = 5), id="", 
                       plot=T){
   
  
  pdf.sep <- get_me_sep(m, lags, name, type, xrange, id, plot=FALSE) 
  pdf.cum <- get_me_cumulative(m, lags, name, type, xrange, id, plot=FALSE) %>% 
    select(-lag)
  
  if(plot){
    return(plot_me(pdf.sep, type, lags) + 
             geom_line(data= pdf.cum, aes(temp, estimate), 
                       color = 'blue'))
  }else{
    return(pdf)
  }
}

# PREPPING DATA -----------------------------------------------------------

get_years <- function(base_start_year = 1990,
                      end_data_year = 2019, 
                      start_proj_year = 2020,
                      end_proj_year = 2100){
  
  pre.yrs   <- base_start_year:end_data_year
  proj.yrs  <- start_proj_year:end_proj_year
  yrs       <- c(pre.yrs, proj.yrs)
  
  TT        <- length(yrs)
  
  pre.ids  <- 1:length(pre.yrs)
  post.ids <- (length(pre.yrs)+1):length(yrs)
  yr.ids   <- c(pre.ids, post.ids) 
  
  return(list(
    pre = pre.yrs, 
    proj = proj.yrs, 
    yrs = yrs, 
    TT = TT, 
    pre.ids = pre.ids, 
    proj.ids = post.ids, 
    ids = yr.ids
  ))
}

gen_mats <- function(NN, yrs, IDs, 
                     var_names = c("y", "g", "temp", "tcons1", "tcons2")) {
  mat <- matrix(NA, nrow=NN, ncol=yrs$TT)
  rownames(mat) <- IDs
  colnames(mat) <- yrs$yrs
  
  setNames(lapply(var_names, function(x) mat), var_names)
}

check_ID_order <- function(mat, df, yr=2000){
  if('year' %in% names(df)) df <- df %>% filter(year == yr)
  
  stopifnot(rownames(mat)==df %>% arrange(ID) %>% pull(ID))
}

format_pre_proj <- function(yrs, df, base, proj){
  
  check_ID_order(base$temp, df, yrs$pre[1])
  
  for(tt in yrs$pre.ids){
    
    yr <- yrs$pre[tt]
    
    # -- Temperature
    base$temp[, tt] <- proj$temp[, tt] <- df %>% 
      filter(year == yr) %>% arrange(ID) %>% pull(temp1)
    
    # --- Growth
    base$g[, tt] <- proj$g[, tt] <- df %>% 
      filter(year == yr) %>% arrange(ID) %>% pull(g)
    
    # --- GDP 
    base$y[, tt] <- proj$y[, tt] <- df %>% 
      filter(year == yr) %>% arrange(ID) %>% pull(y)
  }
  return(list(base=base, proj=proj))
}

load_warming <- function(dir, 
                         scen = "median"
){
  read_excel(paste0(dir, 
                    'cmip6-x0.25_timeseries_tas', 
                    '_timeseries_annual_2015-2100_median,p10,', 
                    'p90_ssp585_ensemble_all_mean.xlsx'), 
             sheet = scen) %>% 
    select(-name) %>% 
    pivot_longer(cols = -c(code), values_to = 'temp') %>% 
    mutate(year = str_remove(name, "-07") %>% as.numeric()) %>% 
    select(-name) %>%
    arrange(code, year) %>% 
    rename(ID=code) %>% 
    filter(year >= 2019) %>% 
    group_by(ID) %>% 
    mutate(warming = temp-first(temp)) %>% 
    ungroup() %>% 
    arrange(ID)
}

load_gdp_ssp <- function(dir, scen){
  stop("not implemented")
}

load_pop_ssp <- function(dir, scen){
  read_excel(paste0(dir, 
                    '/data/projection/', 
                    'pop-x1_timeseries_popcount_timeseries_annual_2010-2100', 
                    '_mean_ssp585_gpw-v4_rev11_mean.xlsx')) %>% 
    select(-name) %>% 
    pivot_longer(cols = -c(code), values_to = 'pop') %>% 
    mutate(year = str_remove(name, "-07") %>% as.numeric()) %>% 
    select(-name) %>%
    arrange(code, year) %>% 
    rename(ID=code) %>% 
    filter(year >= 2019) %>% 
    group_by(ID) %>%
      complete(year = full_seq(year, 1)) %>%  # Create a complete sequence of years
      mutate(pop = approx(year, pop, year)$y) %>%  # Linear interpolation %>% 
    ungroup()
}

# Warming here is either a scalar, for the total warming over the full projection
# or a dataframe, with annual warming for each year in the projection
format_post_proj_baseline <- function(yrs, base, proj, warming){
  
  if(is.data.frame(warming)) {
    check_ID_order(base$temp, warming, yrs$proj[1])
  }else{
    stopifnot(is.numeric(warming))
  }
  
  for(tt in yrs$proj.ids){
    
    yr <- yrs$proj[tt]
    
    # -- Temperature
    base$temp[,tt] <- df.base$temp  
    
    if(is.numeric(warming)){
      proj$temp[,tt] <- df.base$temp + warming / yrs$TT * (tt-min(yrs$proj.ids))  
    }else if(is.data.frame(warming)){
      proj$temp[,tt] <- df.base$temp + warming %>% 
        filter(year == yrs$yrs[tt]) %>% pull(warming)
    }
    
    # -- Growth
    base$g[,tt] <- df.base$g
    
    # -- GDP
    base$y[,tt] <- (1 + base$g[, tt]) * base$y[,tt-1]
    
  }
  return(list(base=base, proj=proj))
}


# PROJECTION --------------------------------------------------------------
extract_coefs <- function(coefs, var){
  
  # Extract coefficients
  matched_names <- str_subset(names(coefs), var)

  # Extract the numeric part using regex and convert to numeric for sorting
  numeric_order <- as.numeric(str_extract(matched_names, "\\d+"))
  
  # Reorder coefficients based on numeric part
  ordered_names <- matched_names[order(numeric_order)]
  
  # Convert to matrix
  coefs[ordered_names] %>% as.matrix()
}

project <- function(b0, b1, type, base, proj, yrs, lags){
  
  if(type == "levels"){
    
    for(tt in 2:max(yrs$proj.ids)) {
      base$tcons1[,tt] <- base$temp[,tt]-base$temp[,tt-1]
      proj$tcons1[,tt] <- proj$temp[,tt]-proj$temp[,tt-1]  
      
      base$tcons2[,tt] <- base$temp[,tt]^2-base$temp[,tt-1]^2
      proj$tcons2[,tt] <- proj$temp[,tt]^2-proj$temp[,tt-1]^2  
    }
  }else if(type == "growth"){
    base$tcons1 <- base$temp
    proj$tcons1 <- proj$temp
    base$tcons2 <- base$temp^2
    proj$tcons2 <- proj$temp^2
  }else{
    stop("not implemented")
  }
  
  for(tt in yrs$proj.ids){
    lls <- tt-0:lags
    
    delta <- (proj$tcons1[, lls] - base$tcons1[, lls]) %*% b0 + 
      (proj$tcons2[, lls] - base$tcons2[, lls]) %*% b1
    
    # -- Growth
    proj$g[, tt] <- base$g[, tt] + delta
    
    # -- GDP
    proj$y[, tt] <- (1 + proj$g[, tt]) * proj$y[, tt-1]
    
  }
  return(list(base=base, proj=proj))
}

global_damages <- function(proj, 
                           base, 
                           df.base, 
                           yrs, 
                           add_global_gdp=FALSE){
  map_dfr(
    1:yrs$TT, 
    function(tt){
      pc.diff    <- 100*(proj$y[,tt]-base$y[,tt]) / base$y[,tt]
      tot.damage <- weighted.mean(pc.diff, df.base$pop)
      out.df <- tibble(year = yrs$yrs[tt], damage = tot.damage)
      if(add_global_gdp){
        out.df <- out.df %>% 
          mutate(global_gdp = sum(base$y[,tt] * df.base$pop))
      }
      out.df
    }
  )  
}


get_damages <- function(m, 
                        type, 
                        base, 
                        proj, 
                        yrs,
                        lags, 
                        df.pop, 
                        uncertainty=T){
  beta <- coef(m)
  vcov <- vcov(m)
  
  # Get central estimate
  b0 <- extract_coefs(beta, "temp1")
  b1 <- extract_coefs(beta, "temp2")
  
  projected <- project(b0=b0, 
                       b1=b1, 
                       type=type, 
                       base=base, 
                       proj=proj, 
                       yrs=yrs, 
                       lags=lags)
  
  base <- projected$base
  proj <- projected$proj
  
  central <- global_damages(proj, base, df.pop, yrs, 
                            add_global_gdp = TRUE) %>% 
    mutate(type = type, lags=lags)
  
  # Draw from statistical uncertainty
  if(uncertainty){
    draws <- MASS::mvrnorm(n = 500, mu = beta, Sigma = vcov)
    uncert <- 
      map_dfr(
        1:dim(draws)[1], function(kk){
          print(kk)
          beta <- draws[kk,]
          b0 <- extract_coefs(beta, "temp1")
          b1 <- extract_coefs(beta, "temp2")
          
          projected <- project(b0=b0, 
                               b1=b1, 
                               type=type, 
                               base=base, 
                               proj=proj, 
                               yrs=yrs, 
                               lags=lags)
          base <- projected$base
          proj <- projected$proj
          
          # Get global damages 
          global_damages(proj, base, df.pop, yrs) %>% 
            mutate(draw = kk)
        }
      ) %>% 
      group_by(year) %>%
      summarize(q025 = quantile(damage, .025), 
                q05  = quantile(damage, .05),
                q25  = quantile(damage, .25),
                q75  = quantile(damage, .75),
                q95  = quantile(damage, .95),
                q975 = quantile(damage, .975), 
                .groups = 'drop') %>% 
      mutate(type = type, lags=lags)
  }else{
    uncert <- tibble()
  }
  
  return(list(central=central, uncert=uncert, base=base, proj=proj))
}


# OUTPUTS -----------------------------------------------------------------

plot_global_damages <- function(central, uncert, vline=2020){
  ggplot(data=uncert) + 
    geom_vline(xintercept=vline, linetype=2)+
    geom_ribbon(aes(x = year, ymin = q025, ymax = q975), alpha=.1) +
    geom_ribbon(aes(x = year, ymin = q05, ymax = q95), alpha=.3) +
    geom_ribbon(aes(x = year, ymin = q25, ymax = q75), alpha=.5) +
    geom_line(data = central, aes(x = year, y = damage), color='red')
}


plot_proj <- function(pdf, vline=2020){
  pdf %>% 
    ggplot() + 
    geom_vline(xintercept=vline, linetype=2)+
    geom_line(aes(x = year, y = value, color = scen)) + 
    facet_wrap(~ID+name, scales='free')
}

plotdf_mats <- function(toPlot, base, proj, plt=T, vline=2020){
  pdf <- map_dfr(
    toPlot, 
    function(ii){
      bind_rows(
        tibble(
          year = as.numeric(colnames(base$y)), 
          ID = ii, 
          y = base$y[ii, ], 
          g = base$g[ii, ], 
          temp = base$temp[ii, ], 
          dam = 0, 
          scen = "base"
        ), 
        tibble(
          year = as.numeric(colnames(proj$y)), 
          ID = ii, 
          y = proj$y[ii, ], 
          g = proj$g[ii, ], 
          temp = proj$temp[ii, ], 
          dam = (proj$y-base$y)[ii, ], 
          scen = "proj"
        )
      )
    }
  ) %>% 
    pivot_longer(cols = c(temp, y, g, dam)) %>% 
    mutate(name = fct_relevel(name, c("temp", "g", "y", "dam")))
  
  if(plt) return(plot_proj(pdf, vline))
  else return(pdf)
}
