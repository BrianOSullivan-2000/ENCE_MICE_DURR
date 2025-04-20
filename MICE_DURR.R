
## Multiple Imputation by Chained Equations 
## with Driect Use of Regularised Regression

# Author: Brian O'Sullivan
# Please refer to https://doi.org/10.1002/joc.8513 for further details 
# where the method is implemented on the Irish rainfall network

## Primary imputation function

## NOTE: The dataframe must follow a "long-table" format. i.e. each row is a 
## measurement with a spatial ID and time ID. All possible ID combinations
## must be included, missing or non-missing, and duplicate 
## combinations are not allowed

MICE_DURR <- function(df, response = "y", m = 20, 
                      max_impute_cycles = 5, 
                      max_EM_cycles = 5, 
                      init_method = c("mean", "idw"),
                      spatial_id = "stno", time_id = "t",
                      EM_tol = 1,
                      nclusters = 1,
                      outfile = "MICE_DURR_log.txt",
                      transformation = {function(x) x},
                      reverse_transformation = {function(x) x},
                      ...){
  
  # df                      -Data frame of response variable and covariates
  # max_impute_cycles       -Maximum number of cycles when getting parameters
  # max_EM_cycles           -Maximum number of cycles when getting final values
  # init_method             -Method for getting starting imputed values
  # spatial_id, time_id     -Column names for locations and times
  # EM_tol                  -Tolerance for convergence of final values 
  # nclusters               -Number of clusters if doing parallel computing
  # outfile                 -Location to print outputs from clusters
  
  # Load in necessary packages
  require(dplyr); require(tidyr)
  
  init_method <- match.arg(init_method)
  
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
  
  # Replicate wide table into m copies
  df_copies <- rep(list(df), m)
  # each element in df_copies has two elements: a random seed and the data
  first_seed <- 222
  df_copies <- mapply(list, 222:(222+(m-1)), rep(list(df), m), 
                SIMPLIFY = F, USE.NAMES = F)
  
  if(nclusters > 1){
    require(parallel)
    
    # Make clusters for parallel computing
    cl <- makeCluster(nclusters, outfile = outfile)
    
    # Import necessary functions to each cluster
    clusterExport(cl, c("elastic_net_DURR", "extract_best_cvm", "rmse"),
                  envir=environment())
    
    # Get parameters for every imputation model
    params <- parSapply(cl = cl, df_copies, MICE_get_params, 
                        missing_idx = missing_idx,
                        max_cycles = max_impute_cycles, ...)
    stopCluster(cl)
  }
  else{
    # This provides the option to not do parallel computing
    params <- sapply(df_copies, MICE_get_params,
                     missing_idx = missing_idx,
                     max_cycles = max_impute_cycles, ...)
  }
  # Pool the regression parameters together from each replicated data set
  params <- apply(params, 1, {function(x) apply(do.call(cbind, x), 1, mean)})
  colnames(params) <- names(df)
  
  # Iterator and convergence tracking
  i <- 1; old_rmse <- .Machine$double.xmax; new_rmse <- .Machine$double.xmax
  
  # Update missing values over multiple cycles with fixed regression parameters
  while((i <= max_EM_cycles) & ((old_rmse/new_rmse > EM_tol) | (i == 1))){
  
    # Update imputed values and convergence tracking
    past_values <- df[missing_idx]
    old_rmse <- new_rmse

    # Update imputed values
    df <- MICE_impute_fixed(df, missing_idx, params, ...)
    
    # New RMSE for convergence
    new_rmse <- rmse(past_values, df[missing_idx])
    
    print(paste0("Cycle ", i, " completed"))
    i <- i + 1
  }
  if(i == max_EM_cycles){
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
  
  return(as.data.frame(df))
}




# For each replicated dataset, calculate the regression parameters
# for each imputation model

MICE_get_params <- function(df, missing_idx, max_cycles = 16, 
                            hyp_cycles = 2,
                            lambda_options = c(0.05, 0.1, 0.15, 0.2),
                            alpha_options = c(0.1, 0.35, 0.65, 0.9),
                            impute_tol = 0, ...){
  
  # df              -Input data
  # missing_idx     -Index of all originally missing values
  # max_cycles      -Max number of cycles to update values
  # hyp_cycles      -Number of cycles to fit hyperparameters
  # lambda_options  -Lambda values to consider for elastic-net
  # alpha_options   -Alpha values to consider for elastic-net
  # impute_tol      -Tolerance for convergence between cycles

  require(glmnet)
  
  # Extract the random seed and data
  seed <- df[[1]]
  df <- df[[2]]
  
  # Set random seed
  set.seed(seed)
  
  # Hyperparameters for each model are saved as elements in an array
  lambdas <- NA
  alphas <- NA
  
  # Iterator and convergence tracking
  i <- 1; old_rmse <- .Machine$double.xmax; new_rmse <- .Machine$double.xmax
  
  # Update the imputed values over multiple cycles
  while((i <= max_cycles) & ((old_rmse/new_rmse > impute_tol) | (i == 1))){
    
    # Update imputed values and convergence tracking
    past_values <- df[missing_idx]
    old_rmse <- new_rmse
    
    # Reset lambdas and alphas if still updating hyperparameters
    if (i < hyp_cycles){
      lambdas <- NA; alphas <- NA
    }
    
    # Loop through and update each column
    for (column in 1:ncol(df)){
      
      # Get target column, other columns, and missing index
      target <- df[, column]; covariates <- df[, -column]
      missing_values <- missing_idx[, column]
      
      # Only care if there are missing values to be updated in column
      if(sum(missing_values) > 0){
        
        # Get lambda and alpha values
        if(is.na(lambdas[column])){
          lambda <- lambda_options
          alpha <- alpha_options
        }
        else{
          # If lambdas and alphas are fixed, get previous values
          lambda <- lambdas[column]
          alpha <- alphas[column]
        }
        
        # Fit imputation model for column using DURR, sample imputed 
        # values from predictive distribution
        imputation_model <- elastic_net_DURR(y = unlist(target),
                                            ry = unlist(missing_values),
                                            x = as.matrix(covariates),
                                            lambda = lambda,
                                            alpha = alpha, 
                                            seed = seed,
                                            ...) 
        
        # Update imputed values and hyperparameters
        df[missing_values, column] <- imputation_model$y
        lambdas[column] <- imputation_model$lambda
        alphas[column] <- imputation_model$alpha
      }
    }
    # Check for convergence
    new_rmse <- rmse(past_values, df[missing_idx])
    
    print(paste0("Cycle ", i, " completed"))
    i <- i + 1
  }
  if(i == max_cycles){
    print("Maximum number of cycles reached")
  }
  else{
    print(paste0("Converged at cycle: ", i-1))
  }
  
  # Each set of regression parameters is saved as an element in a list
  params <- list()
      
  # Get regression parameters for each column
  for (column in 1:ncol(df)){
        
    # Get target column, other columns, missing index, and hyperparameters
    target <- df[, column]; covariates <- df[, -column]
    missing_values <- missing_idx[, column]
    lambda <- lambdas[column]; alpha <- alphas[column]
    
    # Only care if there are missing values to be updated in column
    if(sum(missing_values) > 0){
  
      # Input matrices for glmnet
      x <- as.matrix(covariates)
      y <- unlist(target)
      ymiss <- unlist(missing_values)
      
      # Fit model using elastic-net and save regression parameters
      model <-  glmnet(x = x, y = y,
                       lambda = lambda, alpha = alpha)
      params[[column]] <- as.matrix(coef(model))
      
    }
    else{
      params[[column]] <- rep(0, ncol(covariates)+1)
    }
  }
  return(params)
}




# Fit a regression model using elastic-net and DURR then update missing
# values by sampling from a predictive distribution
# A lot of this function follows the Mice.impute.lasso.norm function 
# from the mice library on CRAN: 
# https://cran.r-project.org/web/packages/mice/index.html

elastic_net_DURR <- function (y, ry, x, nfolds = 10, 
                              lambda = c(0.05, 0.1, 0.15, 0.2),
                              alpha = c(0.1, 0.35, 0.65, 0.9), 
                              seed = 222, ...){
  
  # y               -Target column
  # ry              -Indices of missing values to be updated
  # x               -Auxiliary columns (containing covariates in model)
  # lambda/alpha    -Candidate lambda and alpha values to check
  # seed            -Random seed for cv functions
  
  require(glmnet); require(glmnetUtils)
  
  # Bootstrapping of observed values
  n_obs <- sum(!ry)
  set.seed(seed)
  s <- sample(n_obs, n_obs, replace = TRUE)
  
  # Prepare bootstrapped matrices for glmnet
  x_glmnet <- cbind(1, x)
  x_star <- x_glmnet[!ry, , drop = FALSE][s, , drop = FALSE]
  y_star <- y[!ry][s]
  
  if ((length(lambda) > 1) | (length(alpha) > 1)){
    # If more than one lambda/alpha given, check each lambda/alpha combination
    
    # Fit model for each lambda.alpha combination
    set.seed(seed)
    enet <- cva.glmnet(x = x_star, y = y_star, 
                       alpha = alpha, lambda = lambda) 
    
    # Get lambda/alpha combination with best mean cross-validated error
    cvms <- unlist(lapply(enet$modlist, extract_best_cvm))
    cvms <- matrix(cvms, ncol = 2, byrow = T)
    best_lambdas <- cvms[,1]; cvms <- cvms[,2]
    
    # Get alpha and lambda pair that gives lowest cvm
    alpha <- alpha[which(cvms == min(cvms))]
    lambda <- best_lambdas[which(cvms == min(cvms))]
  }
  else{
    # Fit model using elastic-net
    enet <- glmnet(x = x_star, y = y_star, alpha = alpha, lambda = lambda)     
  }
  
  # We need the regression error for sampling
  s2hat <- mean((predict(enet, x_star, s = lambda, alpha = alpha) - y_star)^2)
  
  # Sample imputed values using DURR
  set.seed(seed)
  y_imp <- predict(enet, x_glmnet[ry, ], s = lambda, alpha = alpha) + 
    rnorm(sum(ry), 0, sqrt(s2hat))
  
  return(list("y" = as.vector(y_imp), "lambda" = lambda, "alpha" = alpha))
}




# Quick function to get the lambda/alpha combination with the lowest
# mean cross-validated error for an elastic-net fitted model
extract_best_cvm <- function(model){
  lambda <- model$lambda.min
  cvm <- model$cvm[which(model$lambda == lambda)]
  return(c(lambda, cvm))
}




# Update missing values with fixed regression parameters

MICE_impute_fixed <- function(df, missing_idx, params, ...){
  
  # Impute each column
  for (column in 1:ncol(df)){
    
    # Get covariates (other columns) and positions of missing values
    covariates <- as.matrix(df[, -column])
    missing_values <- missing_idx[, column]
    
    # Only care if there are missing values to be updated in column
    if(sum(missing_values) > 0){
    
      # If a column has no missing values just skip it
      if (sum(missing_values) > 0){
        
        # Get fixed regression model parameters for current column
        coefs <- params[, column]
        
        # Linear model -> Imputed values calculated with X %*% coefs
        if(sum(missing_values) == 1){
          # Edge case if only one missing value in column
          X <- c(1, covariates[missing_values, ])
          df[missing_values, column] <- as.numeric(X) %*% coefs 
        }
        else{
          df[missing_values, column] <- 
            cbind(rep(1, nrow(df[missing_values, ])),
                  as.matrix(covariates)[missing_values, ]) %*% coefs
        }
      }
    }
  }
  return(df)
}




# Quick rmse function

rmse <- function(y1, y2){
  return(sqrt(mean((y1 - y2)^2)))
}
