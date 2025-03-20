
### interpolation_st.R ###

### A collection of spatio-temporal interpolation functions             ###
### Inverse Distance Weighting and Ordinary Kriging methods are         ###
### available with or without regression. Additionally, a function is   ###
### available to estimate the spatio-temporal anisotropy of a dataset.  ###

### Note that the response variable is renamed to "y" at the beginning  ###
### of each function. Before running, rename any non-response columns   ###
### labelled "y" in your data.                                          ###

### Additionally, spatial coordinates must always be labelled "east"    ###
### and "north", respectively. Time coordinates are labelled "t".       ###
### Again, make sure to avoid duplicated columns for these coordinates. ###
### All spatial coordinates are assumed as Euclidean (e.g. Irish Grid). ### 




# Function for estimating the spatio-temporal anisotropy of a given value

estimate_anisotropy <- function(df, response = "y", cutoff = NULL, 
                                width = NULL, t_cutoff = 6, t_width = 1,
                                progress = T, plot = F){
  
  # df              -dataframe of response variable and coordinates
  # response        -variable of interest for finding anisotropy
  # cutoff          -spatial cutoff of variogram
  # width           -spatial width of bins in variogram
  # t_cutoff        -temporal cutoff of variogram
  # t_width         -temporal width of variogram (recommend to keep at 1)
  # plot            -plot variogram
  
  # Load in required packages
  require(dplyr); require(sp); require(gstat); require(spacetime)

  # Set name of response value to y
  df <- df %>% dplyr::rename(y = all_of(response))
  
  # Get points and dates
  points <- SpatialPoints(cbind(df$east, df$north))
  dates <- as.Date(df$t, origin = "1970-01-01")
  
  # Create spatio-temporal dataframe
  STIDF <- STIDF(points, dates, data = df)
  STSDF <- as(STIDF, "STSDF")
  
  # Heuristic for cutoff and width if none provided
  
  if (missing(cutoff) | missing(width)){
    cutoff <- sqrt((max(df$east) - min(df$east))^2 +
                     (max(df$north) - min(df$north))^2) / 4
    width <- cutoff / 10
  }
  
  # Empirical VGM
  vgmst <- variogramST(y ~  1, data = STSDF, 
                       tlags = seq(0, t_cutoff, by = t_width),
                       cutoff = cutoff, width = width, progress = T)
  
  plot = T
  if (plot == T){
    print(plot(vgmst, wireframe = T, zlim = c(0, max(vgmst$gamma[-1]) * 1.1)))
  }
  
  # Give anisotropy in km / timestep
  C <- estiStAni(vgmst, c(0.1, 1e7)) / 1000
  
  print(paste0("Estimated Anisotropy :", signif(C, 4), "km per timestep"))
  return(C)
}




# Spatio-Temporal Inverse Distance Weighting (No Regression)

IDW_raw_ST <- function(data, newdata, response = "y", idp = 2, nmax = 5,
                       debug.level = 1, C = NULL, 
                       
                       transformation={function (x) x},
                       reverse_transformation={function (x) x}, ...){
  
  # data            -observed values and covariates
  # newdata         -values to be predicted 
  # response        -variable to be interpolated
  # idp             -inverse distance weighting power
  # nmax            -maximum number of neighbours to use for prediction
  # debug.level     -set to -1 to see progress of kriging
  # C               -anisotropy, relates spatial and temporal distance
  #                 units are in km/timestep (i.e., "t")
  
  # transformation          -transformation of response variable
  # reverse_transformation  -reverse transformation before printing results
  
  require(dplyr); require(sp); require(gstat)
  
  # Set name of response value to y
  data <- data %>% dplyr::rename(y = all_of(response)) 
  newdata <- newdata %>% dplyr::rename(y = all_of(response))
  
  # Transform response variable
  data$y <- transformation(data$y)
  
  # If no anisotropy is provided, estimate it from variogram
  # Not recommended as it is quite slow, run the estimate_anisotropy 
  # function on it's own instead
  
  if (is.null(C)){
    C <- estimate_anisotropy(data)
  }
  
  # Convert km to m
  C <- C * 1000
  
  # Scale separation in time by anisotropy
  data$t <- data$t * C; newdata$t <- newdata$t * C
  
  # Create spatio-temporal objects (using time as a third dimension)
  coordinates(data) <- c("east", "north", "t")
  coordinates(newdata) <- c("east", "north", "t")
  
  # Interpolate value at each point in test set
  idw_df <- idw(y ~ 1, data, newdata, idp = idp, nmax = nmax,
                debug.level = debug.level)
  newdata$y <- idw_df$var1.pred
  
  # Reverse transformation of response variable
  newdata$y <- reverse_transformation(newdata$y)
  
  # Reset response column name
  names(newdata)[names(newdata) == "y"] = response
  
  # Reverse anisotropy
  newdata <- as.data.frame(newdata)
  newdata$t <- newdata$t / C
  
  return(as.data.frame(newdata))
}




# Spatio-Temporal Inverse Distance Weighting

IDW_ST <- function(data, newdata, response = "y", idp = 2, nmax = 5, 
                   debug.level = 1, C = NULL,
                   
                   transformation={function (x) x},
                   reverse_transformation={function (x) x},
                   
                   formula = as.formula("y ~ 1"),
                   regression = NULL, plot_regression = F,
                   
                   multicol_fix = F,
                   second_order = F,  second_order_terms = NULL, ...){
  
  # data                    -observed values and covariates
  # newdata                 -values to be predicted 
  # response                -variable to be interpolated
  # idp                     -inverse distance weighting power
  # nmax                    -maximum number of neighbours to use for prediction
  # debug.level             -set to -1 to see progress of kriging
  #                         set to 0 for no printed output
  # C                       -anisotropy, relates spatial and temporal distance
  #                         units are in km/timestep (i.e., "t")
  
  # transformation          -transformation of response variable
  # reverse_transformation  -reverse transformation before printing results
  
  # formula                 -regression formula
  # regression              -regression method used
  # plot_regression         -plot output from regression
  
  # multicol_fix            -remove multi-collinearity in data 
  # second_order            -adds second-order terms for all variables
  # second_order_terms      -only use specific second order terms
  
  require(dplyr); require(sp); require(gstat)
  
  # Set name of response value to y
  data <- data %>% dplyr::rename(y = all_of(response)) 
  newdata <- newdata %>% dplyr::rename(y = all_of(response))
  
  # Transform response variable
  data$y <- transformation(data$y)
  
  # Interpolate without regression if no regression method provided
  if (is.null(regression)){

    newdata <- IDW_raw_ST(data, newdata, idp = idp, nmax = nmax, 
                          debug.level = debug.level, C = C, ...)
    
    # Reverse transformation of response variable
    newdata$y <- reverse_transformation(newdata$y)
    
    # Reset response column name
    names(newdata)[names(newdata) == "y"] = response
    
    return(as.data.frame(newdata))
  }
  
  else{
    # Remove terms with high multi-collinearity 
    if (multicol_fix == T){
      formula <- corr_check(data, formula, ...)
    }
    
    # Add second order terms to formula
    if (second_order == T){
      f <- second_order_formula(formula, data,
                                terms = second_order_terms, ...)
      
      # formula with second order terms
      formula <- f$formula
    }
    
    # Get residuals and model from regression
    r <- regression(data, formula, plot = plot_regression, ...)
    
    data$res <- r$residuals
    model <- r$model
    
    # Need number of components from model if using PCR
    ncomp <- r$ncomp
    
    if (is.null(C)){
      C <- estimate_anisotropy(data %>% dplyr::select(-c("y")),
                               response = "res")
    }
    
    # Convert km to m for anisotropy
    C <- C * 1000

    # Scale timesteps by anisotropy
    data$t <- data$t * C; newdata$t <- newdata$t * C
    
    # Create spatial objects (with time as third dimension)
    coordinates(data) <- c("east", "north", "t")
    coordinates(newdata) <- c("east", "north", "t")
    
    # Interpolate value at each point in test set
    idw_df <- idw(res ~ 1, data, newdata, idp = idp, nmax = nmax,
                  debug.level = debug.level)
    
    # Reset time variable before regression prediction
    newdata <- as.data.frame(newdata)
    newdata$t <- newdata$t / C
    
    # Add interpolated residuals to predicted values from regression
    newdata$y <- as.vector(predict(model, newdata, ncomp = ncomp)) +
      idw_df$var1.pred
    
    # Reverse transformation of response variable
    newdata$y <- reverse_transformation(newdata$y)
    
    # Reset response column name
    names(newdata)[names(newdata) == "y"] = response
    
    return(as.data.frame(newdata))}
}




# Spatio-Temporal Ordinary Kriging (No Regression)

kriging_raw_ST <- function(data, newdata, response = "y", nmax = Inf, 
                           krig_progress = T, start_date = "1970-01-01", 
                           vgm_readout = F,
                           
                           transformation={function (x) x},
                           reverse_transformation={function (x) x}, ...){
  
  # data                  -observed values and covariates
  # newdata               -values to be predicted 
  # response              -variable to be interpolated
  # nmax                  -maximum number of neighbours to use for prediction
  # krig_progress         -add progress bar to kriging interpolation
  # start_date            -starting date for time indices (not too important)
  # vgm_readout           -print summary of variogram
  
  # transformation          -transformation of response variable
  # reverse_transformation  -reverse transformation before printing results
  
  require(dplyr); require(sp); require(gstat); require(spacetime)
  
  # Set name of response value to y
  data <- data %>% dplyr::rename(y = all_of(response)) 
  newdata <- newdata %>% dplyr::rename(y = all_of(response))
  
  # Transform response variable
  data$y <- transformation(data$y)
  
  # Fit spatio-temporal variogram to data
  stvgm <- spatio_temporal_variogram(data, start_date = start_date, 
                                     vgm_readout = vgm_readout, ...)
  fit_vgmst <- stvgm$fit_vgmst
  stAni <- stvgm$stAni
  
  # Create spatial points for all data
  points <- SpatialPoints(cbind(data$east, data$north))
  newpoints <- SpatialPoints(cbind(newdata$east, newdata$north))
  
  # Create dates for all data
  dates <- as.Date(data$t, origin = start_date)
  
  # Create spatio-temporal sparse dataframes
  STIDF <- STIDF(points, dates, data = data)
  ST <- as(STIDF, "STSDF")
  
  # If only interpolating for one moment in time (eg gridding)
  # Need to add one extra date when manipulating data into ST object
  one_time <- F
  time_idx <- unique(newdata$t)
  
  # Add additional day to time values if only one exists
  if(length(time_idx) == 1){ 
    
    one_time <- T
    time_idx <- c(time_idx, time_idx + 1)
    
    # Just double the dataframe for now
    newdata <- bind_rows(newdata, newdata)
    
    # Get date object of data 
    newdates <- as.Date(time_idx, origin = start_date)
    
    # Newdates must be doubled in length to account for extra time added
    newdates <- rep(newdates, each = nrow(newdata) / 2)
  }
  
  # Create points and dates for newdata
  newpoints <- SpatialPoints(cbind(newdata$east, newdata$north))
  
  # Create date object (standard way)
  if (one_time == F){
    newdates <- as.Date(newdata$t, origin = start_date)
  }

  # Create spatio-temporal sparse dataframe
  newSTIDF <- STIDF(newpoints, newdates, data = newdata)
  newST <- as(newSTIDF, "STSDF")
  
  # Kriging missing entries
  krig_df <- krigeST(y ~ 1, ST, newST,
                     fit_vgmst, stAni = stAni, nmax = nmax,
                     progress = krig_progress)
  
  if (one_time == T){
    newST <- newST[, 1]
    krig_df <- krig_df[1:(nrow(krig_df) / 2), ]
  }
  
  newdata <- as.data.frame(newST)
  
  newdata$y <- krig_df$var1.pred
  
  newdata <- newdata[order(newdata$t), ]
  newdata <- newdata[order(newdata$stno), ]
  
  # Reverse transformation of response variable
  newdata$y <- reverse_transformation(newdata$y)
  
  # Reset response column name
  names(newdata)[names(newdata) == "y"] = response
  
  return(as.data.frame(newdata))
}




# Spatio-Temporal Ordinary Kriging

kriging_ST <- function(data, newdata, response = "y", nmax = Inf, 
                       krig_progress = T, start_date = "1970-01-01",
                       vgm_readout = F,
                       
                       transformation={function (x) x},
                       reverse_transformation={function (x) x},
                       
                       formula = as.formula("y ~ 1"),
                       regression = NULL, plot_regression = F,
                       
                       multicol_fix = F, 
                       second_order = F,  second_order_terms = NULL, ...){
  
  # data                  -observed values and covariates
  # newdata               -values to be predicted 
  # response              -variable to be interpolated
  # nmax                  -maximum number of neighbours to use for prediction
  # krig_progress         -add progress bar to kriging interpolation
  # start_date            -starting date for time indices (not too important)
  # vgm_readout           -print summary of variogram
  
  # transformation          -transformation of response variable
  # reverse_transformation  -reverse transformation before printing results
  
  # formula               -regression formula
  # regression            -regression method used
  # plot_regression       -plot output from regression
  
  # multicol_fix          -remove multi-collinearity in data 
  # second_order          -adds second-order terms for all variables
  # second_order_terms    -only use specific second order terms
  
  # Load in required packages
  require(dplyr); require(sp); require(gstat); require(spacetime)

  # Set name of response value to y
  data <- data %>% dplyr::rename(y = all_of(response)) 
  newdata <- newdata %>% dplyr::rename(y = all_of(response))
  
  # Transform response variable
  data$y <- transformation(data$y)
  
  # Interpolate without regression if no regression method provided
  if (is.null(regression)){
    newdata <- kriging_raw_ST(data, newdata, nmax = nmax, 
                              krig_progress = krig_progress,
                              start_date = start_date,
                              vgm_readout = vgm-readout, ...)
    
    # Reverse transformation of response variable
    newdata$y <- reverse_transformation(newdata$y)
    
    # Reset response column name
    names(newdata)[names(newdata) == "y"] = response
    
    return(as.data.frame(newdata))
  }
  
  else{
    # Remove terms with high multi-collinearity 
    if (multicol_fix == T){
      formula <- corr_check(data, formula, ...)
    }
    
    # Add second order terms to formula
    if (second_order == T){
      f <- second_order_formula(formula, data,
                                terms = second_order_terms, ...)
      
      # formula with second order terms
      formula <- f$formula
    }
    
    # Get residuals and model from regression
    r <- regression(data, formula, plot = plot_regression, ...)
    data$res <- r$residuals
    model <- r$model
    
    # Need number of components from model if using PCR
    ncomp <- r$ncomp
    
    # Fit spatio-temporal variogram to data
    stvgm <- spatio_temporal_variogram(data, response = "res",
                                       start_date = start_date, 
                                       vgm_readout = vgm_readout, ...)
    fit_vgmst <- stvgm$fit_vgmst
    stAni <- stvgm$stAni
    
    # Create spatial points for data
    points <- SpatialPoints(cbind(data$east, data$north))
    
    # Create dates for data
    dates <- as.Date(data$t, origin = start_date)

    # Create spatio-temporal sparse dataframes
    STIDF <- STIDF(points, dates, data = data)
    ST <- as(STIDF, "STSDF")

    # If only interpolating for one moment in time (eg gridding)
    # Need to add one extra date when manipulating data into ST object
    one_time <- F
    time_idx <- unique(newdata$t)
    
    # Add additional day to time values if only one exists
    if(length(time_idx) == 1){ 
      
      one_time <- T
      time_idx <- c(time_idx, time_idx + 1)
      
      # Just double the dataframe for now
      newdata <- bind_rows(newdata, newdata)
      
      # Get date object of data 
      newdates <- as.Date(time_idx, origin = start_date)
      
      # Newdates must be doubled in length to account for extra time added
      newdates <- rep(newdates, each = nrow(newdata) / 2)
    }
    
    # Create points and dates for newdata
    newpoints <- SpatialPoints(cbind(newdata$east, newdata$north))
    
    # Create date object (standard way)
    if (one_time == F){
      newdates <- as.Date(newdata$t, origin = start_date)
    }
    
    # Create spatio-temporal sparse dataframe
    newSTIDF <- STIDF(newpoints, newdates, data = newdata)
    newST <- as(newSTIDF, "STSDF")

    # Kriging missing entries
    krig_df <- krigeST(res ~ 1, ST, newST,
                       fit_vgmst, stAni = stAni, nmax = nmax,
                       progress = krig_progress)
    krig_df <- as.data.frame(krig_df)
    
    if (one_time == T){
      newST <- newST[, 1]
      krig_df <- krig_df[1:(nrow(krig_df) / 2), ]
      newdata <- newdata[1:(nrow(newdata) / 2), ]
    }

    newdata <- as.data.frame(newST)
    
    newdata$y <- as.vector(predict(model, newdata, ncomp = ncomp)) +
      krig_df$var1.pred
    
    newdata <- newdata[order(newdata$t), ]
    newdata <- newdata[order(newdata$stno), ]
    
    # Reverse transformation of response variable
    newdata$y <- reverse_transformation(newdata$y)
    
    # Reset response column name
    names(newdata)[names(newdata) == "y"] = response
    
    return(as.data.frame(newdata))
  }
}
