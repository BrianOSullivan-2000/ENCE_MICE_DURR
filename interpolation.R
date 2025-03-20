
### interpolation.R ###

### A collection of spatial interpolation functions                     ###
### Options include: Predict from just a regression model,              ###
### Inverse Distance Weighting (with or without regression),            ###
### Ordinary/Universal Kriging (with or without regression),            ###
### Regression Kriging using Generalised Least Squares (GLS)            ###

### Note that the response variable is renamed to "y" at the beginning  ###
### of each function. Before running, rename any non-response columns   ###
### labelled "y" in your data. Additionally, spatial coordinates are    ###
### labelled "east" and "north", respectively. For all functions,       ###
### spatial coordinates are assumed Euclidean (e.g. Irish Grid TM75)    ###

 


# Basic prediction function that only uses regression model

predict_from_regression <- function(data, newdata, response = "y",
                                    
                                    transformation={function (x) x},
                                    reverse_transformation={function (x) x},
                                    
                                    formula = as.formula("y ~ 1"), 
                                    regression = MLR, plot_regression = F,
                                    
                                    multicol_fix = F, VIF_thresh = 50,
                                    corr_thresh = 0.7,
                                    
                                    second_order = F, 
                                    second_order_terms = NULL,
                                    
                                    ...){
  
  # data                    -observed values and covariates
  # newdata                 -values to be predicted 
  # response                -variable to be interpolated
  
  # transformation          -transformation of response variable
  # reverse_transformation  -reverse transformation before printing results
  
  # formula                 -regression formula
  # regression              -regression method
  # plot_regression         -plot output from regression
  
  # multicol_fix            -remove multicollinearity in data
  # VIF_thresh              -maximum allowed Variance Inflation Factor
  # corr_thresh             -maximum allowed correlation between variables
  
  # second_order            -include second order terms in formula
  # second_order_terms      -can specify second order terms to use
  #                         (if NULL, all possible second order terms included)
  
  
  # Load in packages
  require(dplyr)
  
  # Set name of response value to y
  data <- data %>% dplyr::rename(y = all_of(response)) 
  newdata <- newdata %>% dplyr::rename(y = all_of(response))
    
  # Transform response variable
  data$y <- transformation(data$y)
  
  # Remove terms with high multi-collinearity 
  if (multicol_fix == T){
    formula <- corr_check(data, formula)
  }
  
  # Add second order terms to formula
  if (second_order == T){
    f <- second_order_formula(formula, data,
                              terms = second_order_terms, ...)
    
    # formula with second order terms
    formula <- f$formula
    
    # Create columns for each second order term
    data <- f$df
    newdata <- newdata %>% dplyr::select(names(data))
    }
  
  # Run the regression model
  r <- regression(data, formula, plot = plot_regression, ...)
  model <- r$model
  
  # Need number of components from model if using PCR
  ncomp <- r$ncomp

  # Get predicted values from regression model
  newdata$y <- as.vector(predict(model,
                                 as.data.frame(newdata), ncomp = ncomp))

  # Reverse transformation of response variable
  newdata$y <- reverse_transformation(newdata$y)
  
  # Reset response column name
  names(newdata)[names(newdata) == "y"] = response

  return(as.data.frame(newdata))
}




# Inverse Distance Weighting without any regression

IDW_raw <- function(data, newdata, response = "y", idp = 2, nmax = 5,
                    debug.level = 1,
                    
                    transformation={function (x) x},
                    reverse_transformation={function (x) x}, ...){
  
  # data                    -observed values and covariates
  # newdata                 -values to be predicted 
  # response                -variable to be interpolated
  # idp                     -inverse distance weighting power
  # nmax                    -maximum number of neighbours to use for prediction
  # debug.level             -set to -1 to see progress bar of interpolation
  #                         set to 0 for no printed output
  
  # transformation          -transformation of response variable
  # reverse_transformation  -reverse transformation before printing results
  
  # Load in packages
  require(dplyr); require(sp); require(gstat)
  
  # Set name of response value to y
  data <- data %>% dplyr::rename(y = all_of(response)) 
  newdata <- newdata %>% dplyr::rename(y = all_of(response))

  # Transform response variable
  data$y <- transformation(data$y)
  
  # Create spatial objects
  coordinates(data) <- c("east", "north")
  coordinates(newdata) <- c("east", "north")

  # Interpolate value at each point in test set
  idw_df <- idw(y ~ 1, data, newdata, idp = idp, nmax = nmax, 
                debug.level = debug.level)
  newdata$y <- idw_df$var1.pred
  
  # Reverse transformation of response variable
  newdata$y <- reverse_transformation(newdata$y)
  
  # Reset response column name
  names(newdata)[names(newdata) == "y"] = response
  
  return(as.data.frame(newdata))
}




# Ordinary Kriging without any regression

kriging_raw <- function(data, newdata, response = "y", nmax = Inf, 
                        debug.level = 1, 
                        
                        transformation = {function (x) x},
                        reverse_transformation = {function (x) x}, ...){
  
  # data             -observed values and covariates
  # newdata          -values to be predicted 
  # response         -variable to be interpolated
  # nmax             -maximum number of neighbours to use for prediction
  # debug.level      -set to -1 to see progress of kriging
  #                  set to 0 for no printed output
  
  # transformation          -transformation of response variable
  # reverse_transformation  -reverse transformation before printing results
  
  # Load in packages
  require(dplyr); require(sp); require(gstat)
  
  # Set name of response value to y
  data <- data %>% dplyr::rename(y = all_of(response)) 
  newdata <- newdata %>% dplyr::rename(y = all_of(response))
  
  # Transform response variable
  data$y <- transformation(data$y)
  
  # Fit spatial variogram to values
  fit_vgm <- spatial_variogram(data, ...)
  
  # Create spatial objects
  coordinates(data) <- c("east", "north")
  coordinates(newdata) <- c("east", "north")
  
  
  if (attr(fit_vgm, "singular") == T){
    warning("Singular Fit - Using IDW")
    krig_df <- idw(y ~ 1, data, newdata, idp = 2, nmax = 5, 
                  debug.level = debug.level)
  }
  
  else{
  # Interpolate values with ordinary kriging
  krig_df <- krige(y ~ 1, data, newdata, fit_vgm, nmax = nmax, 
                   debug.level = debug.level)
  }
  
  newdata$y <- krig_df$var1.pred
  
  # Reverse transformation of response variable
  newdata$y <- reverse_transformation(newdata$y)
  
  # Reset response column name
  names(newdata)[names(newdata) == "y"] = response
  
  return(as.data.frame(newdata))
}




# Inverse Distance Weighting Interpolation

IDW <- function(data, newdata, response = "y", idp = 2, nmax = 5, 
                debug.level = 1,
                
                transformation={function (x) x},
                reverse_transformation={function (x) x},
                
                formula = as.formula("y ~ 1"), regression = NULL,
                plot_regression = F,
                
                multicol_fix = F, 
                second_order = F,  second_order_terms = NULL, ...){
  
  # data                    -observed values and covariates
  # newdata                 -values to be predicted 
  # response                -variable to be interpolated
  # idp                     -inverse distance weighting power
  # nmax                    -maximum number of neighbours to use for prediction
  # debug.level             -set to -1 to see progress bar of interpolation
  #                         set to 0 for no printed output
  
  # transformation          -transformation of response variable
  # reverse_transformation  -reverse transformation before printing results
  
  # formula                 -regression formula
  # regression              -regression method
  # plot_regression         -plot output from regression
  
  # multicol_fix            -remove multicollinearity in data
  # second_order            -include second order terms in formula
  # second_order_terms      -can specify second order terms to use
  #                         (if NULL, all possible second order terms included)
  
  # Load in packages
  require(dplyr); require(sp); require(gstat)
  
  # Set name of response value to y
  data <- data %>% dplyr::rename(y = all_of(response)) 
  newdata <- newdata %>% dplyr::rename(y = all_of(response))
  
  # Transform response variable
  data$y <- transformation(data$y)
  
  # Interpolate without regression if no regression method provided
  if (is.null(regression)){
    
    newdata <- IDW_raw(data, newdata, idp = idp,
                       nmax = nmax, debug.level = debug.level, ...)
    
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
  
  # Create spatial objects
  coordinates(data) <- c("east", "north")
  coordinates(newdata) <- c("east", "north")
  
  # Interpolate value at each point in test set
  idw_df <- idw(res ~ 1, data, newdata, idp = idp, nmax = nmax, 
                debug.level = debug.level)

  # Add interpolated residuals to predicted values from regression
  newdata$y <- as.vector(predict(model,
                                 as.data.frame(newdata),
                                 ncomp = ncomp)) + idw_df$var1.pred
  
  # Reverse transformation of response variable
  newdata$y <- reverse_transformation(newdata$y)
  
  # Reset response column name
  names(newdata)[names(newdata) == "y"] = response
  
  return(as.data.frame(newdata))}
}




# Simple regression kriging (only one iteration)

kriging <- function(data, newdata, response = "y", nmax = Inf, 
                    debug.level = 1,
                    
                    transformation={function (x) x},
                    reverse_transformation={function (x) x},
                    
                    formula = as.formula("y ~ 1"),
                    regression = NULL, plot_regression = F,
                    
                    multicol_fix = F, second_order = F,
                    second_order_terms = NULL, ...){
  
  # data                    -observed values and covariates
  # newdata                 -values to be predicted 
  # response                -variable to be interpolated
  # nmax                    -maximum number of neighbours to use for prediction
  # debug.level             -set to -1 to see progress of kriging
  #                         set to 0 for no printed output
  
  # transformation          -transformation of response variable
  # reverse_transformation  -reverse transformation before printing results
  
  # formula                 -regression formula
  # regression              -regression method
  # plot_regression         -plot output from regression
  
  # multicol_fix            -remove multicollinearity in data
  # second_order            -include second order terms in formula
  # second_order_terms      -can specify second order terms to use
  #                         (if NULL, all possible second order terms included)
  
  
  # Load in packages
  require(dplyr); require(sp); require(gstat)
  
  # Set name of response value to y
  data <- data %>% dplyr::rename(y = all_of(response)) 
  newdata <- newdata %>% dplyr::rename(y = all_of(response))
  
  # Transform response variable
  data$y <- transformation(data$y)
  
  # Interpolate without regression if no regression method provided
  if (is.null(regression)){
    
    newdata <- kriging_raw(data, newdata, nmax = nmax, 
                           debug.level = debug.level, ...)
    
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
  
  # Fit spatial variogram to residuals
  fit_vgm <- spatial_variogram(data, response = "res", ...)
  
  # Create spatial objects
  coordinates(data) <- c("east", "north")
  coordinates(newdata) <- c("east", "north")
  
  # Interpolate residuals with ordinary kriging
  krig_df <- krige(res ~ 1, data, newdata, fit_vgm, nmax = nmax)
  
  # Change back to dataframe object before using predict function
  newdata <- as.data.frame(newdata)
  
  # Add interpolated residuals to predicted values from regression
  newdata$y <- as.vector(predict(model, newdata, ncomp = ncomp)) +
    krig_df$var1.pred
  
  # Reverse transformation of response variable
  newdata$y <- reverse_transformation(newdata$y)
  
  # Reset response column name
  names(newdata)[names(newdata) == "y"] = response
  
  return(as.data.frame(newdata))}
}




# Regression Kriging with iterations using GLS

multi_kriging <- function(data, newdata, response = "y", nmax = Inf,
                          debug.level = 1, tol = 1, iter = Inf, 
                       
                          transformation={function (x) x},
                          reverse_transformation={function (x) x},
                          
                          formula = as.formula("y ~ 1"),
                          regression = MLR, plot_regression = F,
                          
                          multicol_fix = F, 
                          second_order = F, second_order_terms = NULL,
                          
                          cutoff = NULL, width = NULL,
                          
                          GLS_vgm_model = "Mat", 
                          GLS_kappa = 1, GLS_fit_kappa = F,
                          GLS_plot_vgm = F,
                          GLS_plot_regression = F, GLS_readout = F, ...){
  
  # data                    -observed values and covariates
  # newdata                 -values to be predicted 
  # response                -variable to be interpolated
  # nmax                    -maximum number of neighbours to use for prediction
  # debug.level             -set to -1 to see progress of kriging
  #                         set to 0 for no printed output
  # tol                     -tolerance, minimum change that continues iteration
  # iter                    -max number of iterations in GLS step
  
  # transformation          -transformation of response variable
  # reverse_transformation  -reverse transformation before printing results
  
  # formula                 -regression formula
  # regression              -regression method
  # plot_regression         -plot output from regression
  
  # multicol_fix            -remove multicollinearity in data
  # second_order            -include second order terms in formula
  # second_order_terms      -can specify second order terms to use
  #                         (if NULL, all possible second order terms included)
  
  # cutoff                  -cutoff of variogram
  # width                   -binwidth of variogram
  
  # GLS_vgm_model           -variogram model for GLS steps
  # GLS_kappa               -specify kappa parameter for Matern model
  # GLS_fit_kappa           -option to fit kappa parameter (usually not great)
  # GLS_plot_vgm            -plot variogram during GLS step
  # GLS_plot_regression     -plot output from GLS regression
  # GLS_readout             -print summary of GLS regression
  
  
  # Load in packages
  require(sp); require(gstat); require(dplyr)
  
  # Set name of response value to y
  data <- data %>% dplyr::rename(y = all_of(response)) 
  newdata <- newdata %>% dplyr::rename(y = all_of(response))
  
  # Transform response variable
  data$y <- transformation(data$y)
  
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
  
  # If regression has removed terms in the formula, update formula
  formula <- as.formula(paste("y", 
                              paste(r$terms, collapse = " + "),
                              sep = " ~ "))
    
  if (missing(cutoff) | missing(width)){

    # Will have to get cutoff and width outside of variogram function
    cutoff <- sqrt((max(data$east) - min(data$east))^2 + 
                     (max(data$north) - min(data$north))^2)
    width = cutoff / 15
    }
  
  # Fit spatial variogram to residuals
  fit_vgm <- spatial_variogram(data, response = "res",
                               cutoff = cutoff, width = width, ...)
  
  # Create spatial objects
  coordinates(data) <- c("east", "north")
  coordinates(newdata) <- c("east", "north")
  
  # Interpolate residuals with ordinary kriging
  krig_df <- krige(res ~ 1, data, newdata, fit_vgm, nmax = nmax, 
                   debug.level = debug.level)
  
  # Change back to dataframe object before using predict function
  newdata <- as.data.frame(newdata)
  
  # Add interpolated residuals to predicted values from regression
  newdata$y <- as.vector(predict(model, newdata, ncomp = ncomp)) +
    krig_df$var1.pred
  
  # Save current values which are used to check for convergence
  old_values <- newdata$y
  
  # Additional counters to track during iterations
  old_diff <- 0; percentage_change <- 1; iter_count <- 0
  
  # Iterate two step process (1. Regression 2. Kriging)
  # Stop loop once values don't change sufficiently (according to tol) or
  # if more iterations have occurred than specified by iter (to a minimum of 3)
  
  while (((percentage_change >= tol) | (iter_count < 3)) & (iter_count < iter)){
    
    data <- as.data.frame(data)
    newdata <- as.data.frame(newdata)
    
    # Linear regression with correlation structure from variogram fit
    GLS <- GLS(data, formula, fit_vgm, plot = GLS_plot_regression,
               readout = GLS_readout)
    model <- GLS$model
    data$res <- GLS$residuals
    
    # Get parameters of variogram from previous fit
    psill <- fit_vgm[2, ]$psill; nugget <- fit_vgm[1, ]$psill
    range <- fit_vgm[2, ]$range
    
    # Fit spatial variogram to residuals
    fit_vgm <- spatial_variogram(data, response = "res", 
                                 vgm_model = GLS_vgm_model, psill = psill,
                                 nugget = nugget, range = range,
                                 cutoff = cutoff, width = width,
                                 flex_vgm = F, flex_fit = F,
                                 plot_vgm = GLS_plot_vgm,
                                 kappa = GLS_kappa, fit.kappa = GLS_fit_kappa)
  
    # Create spatial object
    coordinates(data) <- c("east", "north")
    coordinates(newdata) <- c("east", "north")
    
    # Interpolate residuals with ordinary kriging at missing points
    krig_df <- krige(res ~ 1, data, newdata, fit_vgm, 
                     nmax = nmax, debug.level = debug.level)
  
    # Add interpolated residuals to regression predicted values
    newdata$y <- predict(model, newdata) + krig_df$var1.pred
    
    # Change in predicted values after current iteration
    diff <- rmse(old_values, newdata$y)
    
    old_values <- newdata$y
    
    percentage_change <- abs(old_diff - diff)/diff
    
    # **TEMP**
    if (percentage_change > 1e5){
      percentage_change = 0
      print("Inf error")
      print(paste0("t = ", data$t[1]))
      }
    
    old_diff <- diff; iter_count <- iter_count + 1
    
  }
  
  # Reverse transformation of response variable
  newdata$y <- reverse_transformation(newdata$y)
  
  # Reset response column name
  names(newdata)[names(newdata) == "y"] = response
  
  print(paste0("t = ", data$t[1]))
  return(as.data.frame(newdata))
}
