# CONTENTS:
## 1. Running regressions:
  # loading historic data
  # running regression
  # plotting marginal effects

## 2. Prepping data for projection

# RUNNING REGRESSIONS -----------------------------------------------------

load_historic_data <- function(dir){
  df.in <- read_rds(file.path(dir,"temp_gdp_world_panel.rds")) %>% 
    rename(temp1 = era_mwtemp, pop=SP.POP.TOTL) %>%
    group_by(ID=ISO3) %>% 
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
    filter(!is.na(y), !is.na(temp1)) %>% 
    group_by(ID) %>% 
    useful::add_lags(vars = c("temp", "dtemp"), lags=10, max.p=2, 
                     sort_df = FALSE) %>% 
    select(ID, year, time1, time2, y, g=dy, pop, contains("temp"))
}

get_var <- function(type){
  if(type=="levels") var <- "dtemp"
  else if(type == "growth") var <- "temp"
  else stop("not implemented")
  return(var)
}

run_reg <- function(df, 
                    type="levels", 
                    lags=1, 
                    spec="poly2", 
                    FE = "ID + time1 + ID[time1]"){
  
  var <- get_var(type)
  poly_order <- as.numeric(str_extract(spec, "[0-9]"))
  if(lags==0) var <- paste0("l0_", var)
  
  ff <- useful::build_formula_poly(
    yvar = "g", treat = var, poly_treat = poly_order, control=NULL, 
    leads = 0, lags=lags , FE=FE
  )
  feols(
    ff, 
    data = df, 
    cluster = "ID"
  )
}

plot_me <- function(pdf, type, lags=NULL){
  
  p <- ggplot(pdf, 
              aes(temp, estimate)) + 
    geom_hline(yintercept=0, linetype="dashed") +
    geom_line() + 
    geom_ribbon(aes(temp, ymin = conf.low, ymax = conf.high), 
                alpha = .2) + 
    ggtitle(type) + xlab("Temp") + ylab("ME")
  
  if(!is.null(lags)) p <- p+facet_wrap(~lag) 
  
  p
}

get_me_sep <- function(m, lags, type="growth", xrange=seq(0, 30, by = 5), id="", 
                   plot=T){
  
  var <- get_var(type)
  
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
  )
  if(plot){
    return(plot_me(pdf, type, lags))
  }else{
    return(pdf)
  }
}

get_me_cum <- function(m, lags, type="growth", xrange=seq(0, 30, by = 5), id="", 
                   plot=T){
  
  var <- get_var(type)
  
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

get_me <- function(m, lags, type="growth", xrange=seq(0, 30, by = 5), id="", 
                       plot=T){
   
  
  pdf.sep <- get_me_sep(m, lags, type, xrange, id, plot=FALSE) 
  pdf.cum <- get_me_cum(m, lags, type, xrange, id, plot=FALSE) %>% 
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

gen_mats <- function(NN, TT, IDs, yrs){
  mat <- matrix(NA, nrow=NN, ncol=TT)
  rownames(mat) <- IDs
  colnames(mat) <- yrs
  y <- g <- temp <- tcons <- mat
  return(list(y=y, g=g, temp=temp, tcons=tcons))
}

check_ID_order <- function(mat, df, yr){
  stopifnot(rownames(mat)==df %>% 
              filter(year == yr) %>% arrange(ID) %>% 
              pull(ID))
}

format_pre_proj <- function(pre.ids, pre.yrs, df, base, proj){
  
  check_ID_order(base$temp, df, pre.yrs[1])
  
  for(tt in pre.ids){
    
    yr <- pre.yrs[tt]
    
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

format_post_proj_baseline <- function(post.ids, base, proj, warming, TT){
  
  for(tt in post.ids){
    
    # -- Temperature
    base$temp[,tt] <- df.base$temp  
    proj$temp[,tt] <- df.base$temp + warming / TT * (tt-min(post.ids))
    
    # -- Growth
    base$g[,tt] <- df.base$g
    
    # -- GDP
    base$y[,tt] <- (1 + base$g[, tt]) * base$y[,tt-1]
    
  }
  return(list(base=base, proj=proj))
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
