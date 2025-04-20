
## Elastic-net Chained Equations primary functions ##

# Author: Brian O'Sullivan
# Please refer to https://doi.org/10.1002/joc.8513 for further details 
# where the method is implemented on the Irish rainfall network

## Primary imputation function

## NOTE: The dataframe must follow a "long-table" format. i.e. each row is a 
## measurement with a spatial ID and time ID. All possible ID combinations
## must be included, missing or non-missing, and duplicate 
## combinations are not allowed

ENCE <- function(df, response = "y", 
                 hyp_cycles = 2, 
                 max_cycles = 16, 
                 init_method = c("mean", "idw"),
                 spatial_id = "stno", time_id = "t",
                 tol = 1,
                 transformation = {function(x) x},
                 reverse_transformation = {function(x) x},
                 ...){
  
  # df                      -Data frame of response variable and covariates
  # hyp_cycles              -Number of cycles to fit hyperparameters
  # max_cycles              -Maximum number of cycles to update values
  # init_method             -Method for getting starting imputed values
  # spatial_id, time_id     -Column names for locations and times
  # tol                     -Tolerance value for convergence
  
  # Load in necessary packages
  require(dplyr); require(tidyr)
  
  init_method <- match.arg(init_method)
  
  # Order dataframe 
  df <- df[order(df[[spatial_id]], df[[time_id]]), ]
  
  # Label missing entries
  df <- df %>% mutate(missing = is.na(.data[[response]]))
  original_spatial_order <- unique(df[[spatial_id]])
  original_times <- df[[time_id]]
  
  # Sort the spatial locations from least missing values to most
  # (This isn't necessary, but can improve imputation accuracy)
  df <- df %>%
    group_by(across(all_of(spatial_id))) %>%
    mutate(.n_missing = sum(missing)) %>%
    arrange(.n_missing) %>%
    dplyr::select(-.n_missing)
  
  # Get starting imputed values
  # Mean imputation (mean of each spatial location)
  if(init_method == "mean"){
    df <- df %>% group_by(across(all_of(spatial_id))) %>%
      mutate(across(all_of(response), ~replace_na(., mean(., na.rm=TRUE))))
  }
  # Impute using Inverse Distance Weighting (idw)
  # refer to functions in "idw.R"
  else if(init_method == "idw"){
    df[df$missing, ] <- IDW_ST(df[!df$missing, ], df[df$missing, ], 
                               response, ...)
  }
  else{stop("Error: invalid starting imputation method selected")}
  
  # Transformation of response
  df[response] <- transformation(df[response])
  
  # Get wide table of missing entry locations
  missing_idx <- df %>% dplyr::select(all_of(c(time_id, spatial_id)), missing) %>%
    pivot_wider(names_from = all_of(spatial_id), 
                values_from = missing) %>%
    dplyr::select(-all_of(time_id)) %>% as.matrix()
  
  # Convert dataframe to wide table
  df <- df %>% dplyr::select(all_of(c(time_id, spatial_id, response))) %>%
    pivot_wider(names_from = all_of(spatial_id), 
                values_from = all_of(response)) %>%
    dplyr::select(-all_of(time_id))
  
  # Add character to beginning of IDs in case they are just numbers
  colnames(missing_idx) <- paste0(spatial_id, colnames(missing_idx))
  names(df) <- paste0(spatial_id, names(df))
  
  # List of lambda and alpha values
  ls <- NA; as <- NA
  # Iterator and convergence tracking
  i <- 1; old_rmse <- .Machine$double.xmax; new_rmse <- .Machine$double.xmax
  
  # Loop through imputation cycle multiple times
  while((i <= max_cycles) & ((old_rmse/new_rmse > tol) | (i == 1))){
    
    # Update imputed values and convergence tracking
    past_values <- df[missing_idx]
    old_rmse <- new_rmse
    
    # Compute imputed values and update lambdas/alphas
    imputed_df <- ENCE_impute(df, missing_idx, ls, as, ...)
    df <- imputed_df$df
    ls <- imputed_df$ls
    as <- imputed_df$as
    
    # New RMSE for convergence
    new_rmse <- rmse(past_values, df[missing_idx])
    
    # Reset lambdas and alphas if still updating hyperparameters
    if (i < hyp_cycles){
      ls <- NA; as <- NA
    }
    print(paste0("Cycle ", i, " completed"))
    i <- i + 1
  }
  if(i == max_cycles){
    print("Maximum number of cycles reached")
  }
  else{
    print(paste0("Converged at cycle: ", i))
  }
  
  # Return df to long format
  names(df) <- substr(names(df), nchar(spatial_id)+1, nchar(names(df)))
  df <- pivot_longer(df, cols = everything(), 
                     names_to = spatial_id, values_to = response)
  df <- df %>%
    mutate(across(all_of(spatial_id), 
                  ~factor(.x, levels = original_spatial_order))) %>%
    arrange(across(all_of(spatial_id)))
  
  # Add times back and undo transformation
  df[time_id] <- original_times
  df[response] <- reverse_transformation(df[response])
  
  return(df)
}




# One single imputation cycle for all spatial locations

ENCE_impute <- function(df, missing_idx, ls = NA, as = NA, ...){
  
  # df              -Input data
  # missing_idx     -Index of all originally missing values
  # ls              -Initial lambda parameters for all imputation models
  # as              -Initial alpha parameters for all imputation models
  
  # Loop through and update each column
  for (column in 1:ncol(df)){
    
    # Get target column, other columns, and missing index
    target <- df[, column]; covariates <- df[, -column]
    missing_values <- missing_idx[, column]
    
    # Impute target station
    imputed_column <- column_impute(covariates, target, missing_values, 
                                    ls[column], as[column], ...)
    
    # Update covariates and hyper  parameters with imputation model output
    target <- imputed_column$target; df[, column] <- target
    ls[column] <- imputed_column$lambda; as[column] <- imputed_column$alpha
    
  return(list("df" = df, "ls" = ls, "as" = as))
}




# Impute all missing values for a column using elastic-net (glmnet)

column_impute <- function(covariates, target, missing_values, 
                          lambda = NA, alpha = NA, 
                          lambdas = c(0.05, 0.1, 0.15, 0.2),
                          alphas = c(0.1, 0.35, 0.65, 0.9),
                          cluster = NULL){
  
  # covariates          -Covariate Columns (X)
  # target              -Target Column (Y)
  # missing_values      -Positions of missing values to be imputed
  # lambda              -Tuning Parameter (rate of regularisation)
  # alpha               -Weighted value between lasso and ridge regression
  # lambdas/alphas      -Hyper parameter values to be checked using CV
  # cluster             -Can provide a cluster for parallel computing
  
  require(glmnet); require(glmnetUtils)
  
  # Some stations might have no missing covariates
  if (any(missing_values)){
    
    # Input matrices for glmnet
    x <- model.matrix(~ ., covariates)
    y <- unlist(target)
    ymiss <- unlist(missing_values)
    
    # Lambda/alpha need to be fit if not provided
    if (is.na(lambda)){
      
      set.seed(222)
      # Fit E-Net with CV to fit lambda and alpha
      model <- cva.glmnet(x = x[!ymiss, ], y = y[!ymiss], 
                          lambda = lambdas, alpha = alphas,
                          outerParallel = cluster)
      
      # Get the best cvm (error) and corresponding lambda for each alpha
      cvms <- unlist(lapply(model$modlist, extract_best_cvm))
      cvms <- matrix(cvms, ncol = 2, byrow = T)
      best_ls <- cvms[,1]; cvms <- cvms[,2]
      
      # Get alpha/lambda pair with lowest cvm
      alpha <- alphas[which(cvms == min(cvms))]
      lambda <- best_ls[which(cvms == min(cvms))]
    }
    
    # If lambda and alpha do not need to be fitted
    else{
      model <- glmnet(x = x[!ymiss, ], y = y[!ymiss],
                      lambda = lambda, alpha = alpha)
    }
    
    # Calculate predicted values for missing covariates
    target[missing_values, ] <- predict(model, x[ymiss,,drop=F], 
                                   s = lambda, alpha = alpha)
  }
  
  return(list("target" = target, "lambda" = lambda, "alpha" = alpha))
}




# Get lowest cvm and corresponding lambda from a `cvfit` object

extract_best_cvm <- function(cvfit){
  lambda <- cvfit$lambda.min
  cvm <- cvfit$cvm[which(cvfit$lambda == lambda)]
  return(c(lambda, cvm))
}




# Quick rmse function

rmse <- function(y1, y2){
  return(sqrt(mean((y1 - y2)^2)))
}
