
### regression.R ###

### Various wrappers for different regression methods.              ### 
### Examples: Multiple Linear Regression (MLR), Elastic-Net         ###
### Regularisation (E-Net), Principal Component Regression (PCR)    ###
### Generalised Least Squares Regression (GLS)                      ###
### Each regression method returns the regression model (model),    ###
### the model residuals (residuals), and model coefficients (coef). ###

### Additional regression tools are included, see below.            ###




# Simple Step-wise Multiple Linear Regression

MLR <- function(df, formula, step_method = "forward", plot = F, readout = F, ...){

  # df                      -data frame of response variable and covariates
  # formula                 -regression formula
  # step_method             -method to carry out stepwise regression
  #                         (forward, backward, both, none)
  #                         if none, no step-wise regression will be calculated
  # plot                    -plot the output of the regression model
  # readout                 -print a summary of the model

  # MASS package required for step-wise
  require(MASS)

  # Standard Regression  
  model <- lm(formula = formula, data = df)

  if (tolower(step_method) != "none"){
    
    # Stepwise Regression
    model <- stepAIC(model, method = step_method , trace = F)
  }
  
  # Plot and Summary
  if (plot == T){ plot(model) }
  if (readout == T){ print(summary(model)) }

  return(list("model" = model, "residuals" = model$residuals,
              "terms" = names(coef(model))[-1], "coef" = model$coefficients))
}




# Elastic-Net Regression (Lasso and Ridge)

Elastic_Net <- function(df, formula, alpha = 0.1, lambda = NA,
                        plot = F, readout = F, ...){
  
  # df                      -data frame of response variable and covariates
  # formula                 -regression formula
  # alpha                   -Elastic Net penalty
  # lambda                  -scaling parameter of model 
  #                         Manually choosing lambda not recommended
  #                         instead find lambda with cv.glmnet in function
  
  # NOTE: alpha is a hyperparameter describing how much of the regression
  # is attributed to lasso (alpha = 1) and ridge (alpha = 0)
  # it's optimal value must be determined using cross validation
  
  # plot                    -plot the output of the regression model
  # readout                 -print a summary of the model
  
  # Elastic-Net Regularisation is done through the glmnet package
  require(glmnet); require(glmnetUtils)

  # E-Net model (CV selects best scaling parameter, lambda)
  
  # Get lambda from CV if unspecified
  if (is.na(lambda)){
    
    # Fix folds for CV to ensure consistency
    set.seed(222)
    foldid <- sample(rep(1:10, ceiling(nrow(df)/10)))
    foldid <- head(foldid, nrow(df))
    
    cvfit <-
      glmnetUtils::cv.glmnet(formula, df, alpha = alpha,
                             use.model.frame = T, 
                             lambda = seq(0.05, 0.25, 0.05),
                             foldid = foldid)
    lambda <- cvfit$lambda.min
  }
  
  else{
    cvfit <- glmnetUtils::glmnet(formula, df, alpha = alpha, lambda = lambda,
                                 use.model.frame = T)
  }
  
  # Response vector is determined from LHS of formula
  y <- eval(parse(text = paste0(obj_name(df),"$",as.character(formula)[2])))
  
  # Calculate residuals
  res <- y - as.vector(predict(cvfit, df, s = lambda))
  
  # Plotting
  if (plot == T){
    
    # Fitted values vs. residuals
    plot(predict(cvfit, df), res)
    
    # QQ-plot
    qqnorm(res); qqline(res)
  }
  
  # Summary
  if (readout == T){
    # Readout just prints model coefficients in this case
    print(coef(glmnet(formula, df, alpha = alpha,
                      lambda = lambda, use.model.frame = T)))
    print(paste0("Lambda ", lambda))
  }
  
  return(list("model" = cvfit, "residuals" = res, "coef" = coef(cvfit),
              "terms" = rownames(coef(cvfit))[as.vector(coef(cvfit) !=
                                                          0)][-1],
              "lambda" = lambda))
}




# Generalised Least Squares Regression
# Uses spatial correlation structure obtained from a variogram (fit_vgm)

GLS <- function(df, formula, fit_vgm, plot = F, readout = F, ...){
  
  # df              -dataframe (response variable and explanatory variables)
  # formula         -regression formula 
  # fit_vgm         -spatial variogram (covariance structure of residuals)
  # plot            -plot the output of the regression model
  # readout         -print a summary of the model
  
  
  # need gls function from nlme package
  require(nlme)
  
  # Create spatial object
  sdf <- df
  coordinates(sdf) <- c("east", "north")
  
  # Get distances between as points
  ds <- spDists(sdf)
  
  # Get covariances between all points using variogram
  V <- variogramLine(fit_vgm, dist_vector = ds, covariance = T)
  
  # Convert covariance matrix to a correlation matrix
  CM <- cov2cor(V)
  
  # Turn correlation matrix into a "corStruct" object 
  # (this is used in gls function)
  C <- corSymm(CM[lower.tri(CM)], fixed = T)
  
  # GLS regression 
  model <- gls(formula, df, correlation = C)

  # Plot
  if (plot == T){
    
    # Fitted values vs. residuals
    print(plot(model))
    
    # QQ-plot
    qqnorm(model$res); qqline(model$res)
  }
  
  # Summary
  if (readout == T){
    print(model$coefficients)
  }
  
  return(list("model" = model, "residuals" = model$res,
              "coef" = model$coefficients))
}




# Principal Component Regression
# Take the covariates (east, north, elev) and combine the ones that have high
# covariance into principal components. These components are used for regression

PCR <- function(df, formula, center = T, scale = T, 
                plot = F, readout = F, elbow_plot = F, validation = "CV",
                ncomp = NA, plsr = F, ...){
  
  # df              -dataframe (response variable and explanatory variables)
  # formula         -regression formula
  # center          -center the data around a mean before making components
  # scale           -scale the data (standardise)
  #                 highly recommended to set to True
  
  # plot            -plot the output of the regression model
  # readout         -print a summary of the model
  # elbow_plot      -create an elbow plot of the data
  # validation      -include validation method for selecting ncomp
  # ncomp           -number of principal components
  # plsr            -can opt to use plsr instead of pcr (partial least squares)
  
  # Done with pls package
  require(pls)
  
  # Create PLSR regression model
  if (plsr == T){
    model <- plsr(formula = formula, data = df,
                  center = center, scale = scale)
  }
  
  # Or create PCR regression model
  else{
    model <- pcr(formula = formula, data = df, 
                 center = center, scale = scale, validation = validation)
  }
  
  # Find optimal number of components if not specified
  # Option to output an elbow plot
  if (is.na(ncomp)){
    ncomp <- selectNcomp(model, plot = elbow_plot)
    }
  
  # Get the residuals which correspond to the best number of components
  res <- as.numeric(model$residuals[,,ncomp])
  
  # Plotting
  if (plot == T){
    
    # Fitted values vs. residuals
    plot(predict(model, df, ncomp = ncomp), res)
    
    # QQ-plot
    qqnorm(res); qqline(res)
  }
  
  # Summary of model
  if (readout == T){
    print(summary(model))
  }

  return(list("model" = model, "residuals" = res, 
              "terms" = rownames(coef(model)),
              "coef" = model$coefficients, "ncomp" = ncomp))
}




# This function removes multicollinearity in a linear regression model
# (i.e. covariates with high correlation)
# NOTE: this function will remove second-order terms and interactions

corr_check <- function(df, formula, 
                       keep_vars = c("east", "north"),
                       VIF_thresh = 50, corr_thresh = 0.7, 
                       verbose = F, ...){
  
  # df             -dataframe (response variable and explanatory variables)
  # formula        -regression formula
  # step_method    -method for stepwise regression (forward, backward, both)
  # keep_vars      -these variables will always be kept in the formula
  #                For example, you will commonly always need easting/northing
  
  # VIF_thresh     -threshold Variance Inflation Factor
  # corr_thresh    -threshold correlation between variables
  # verbose        -print the change in the formula
  
  # Load in packages
  require(regclass); require(MASS)
  
  if (verbose == T){
    print("Starting Formula"); print(formula)
    }
  
  # Get the origingal regression model (original formula)
  model <- lm(formula = formula, data = df)
  
  # If there are any terms in the formula that don't exist in the dataframe,
  # remove them (second order terms, interactions)
  vars <- unlist(strsplit(as.character(formula)[3],  " + ", fixed = T))
  vars <- vars[vars %in% names(df)]
  vars <- vars[!(vars %in% keep_vars)]

  # Change formula to only include vars
  formula <- as.formula(
    paste(as.character(formula)[2],
          paste(vars, collapse = " + "),
          sep = " ~ ")
  )
  
  # These are used to track how the model changes in the while loop
  old_vars <- vars; no_change <- F
  
  # If there are no variables with a VIF above the threshold, or if
  # the model doesn't change from one iteration to the next, 
  # terminate and return the current model
  
  while (any(VIF(model) > VIF_thresh) & no_change == F){
    
    # Find variables with VIF over the threshold
    multicol_vars <- VIF(model)[VIF(model) > VIF_thresh]
    multicol_vars <- multicol_vars[names(multicol_vars) %in% old_vars]
    
    # Select variable with largest VIF
    if (length(multicol_vars) > 0){
      big_var <- names(multicol_vars)[multicol_vars == max(multicol_vars)]
    }
    else{big_var <- NULL} # This will only be necessary on final loop
    
    # Get the correlation between big_var and all other variables in model
    cors <- cor(df[big_var], df[old_vars])

    # Identify the variables with a correlation higher than the threshold
    other_vars <- colnames(cors)[abs(cors) >= corr_thresh]
    
    # Only concerned with variables still in model
    other_vars <- other_vars[other_vars %in% old_vars]
    
    # Names of variables to check for highest adjR^2
    correlated_vars <- c(big_var, other_vars)
    
    # For each variable, run a regression to see which is best
    adjust_formula <- function(i){
      
      # Name of correlated variable of interest
      correlated_var <- correlated_vars[i]
      
      # Remove all correlated variables except the current one
      new_vars <- unlist(strsplit(as.character(formula)[3],  " + ", fixed = T))
      new_vars <- new_vars[!(new_vars %in% 
                               correlated_vars[correlated_vars != 
                                                 correlated_var])]
      new_vars <- c(keep_vars, new_vars)
      
      # Make a new regression formula
      new_formula <- as.formula(
        paste(as.character(formula)[2],
              paste(new_vars, collapse = " + "),
              sep = " ~ ")
      )

      # Get adjR^2 of new formula
      model <- lm(formula = new_formula, data = df)
      
      return(summary(model)$adj.r.squared)
    }
    
    # Get adjusted R squared values for each correlated variable
    adjrsqs <- unlist(lapply(1:length(correlated_vars), adjust_formula))
    
    # Get the best performing correlated variable 
    best_var <- correlated_vars[adjrsqs == max(adjrsqs)]
    
    # Sometimes the adjR^2s are equal, can just pick one
    if (length(best_var) > 1){best_var = best_var[1]}
    
    # Remove all correlated variables except the best one
    new_vars <- unlist(strsplit(as.character(formula)[3],  " + ", fixed = T))
    new_vars <- new_vars[!(new_vars %in% 
                             correlated_vars[correlated_vars != best_var])]
    
    # The kept variables are removed from the next loop 
    # to avoid double counting (added back to formula later)
    new_vars <- new_vars[!(new_vars %in% keep_vars)]
    
    # Recreate the formula
    if (!(is.null(keep_vars))){
      
    # Always include keep_variables
    formula <- as.formula(
      paste(as.character(formula)[2],
            paste(paste(keep_vars, collapse = " + "),
                  paste(new_vars, collapse = " + "),
                  sep = " + "),
            sep = " ~ ")
      )
    }
    
    else{
      
      formula <- as.formula(
        paste(as.character(formula)[2],
              paste(paste(new_vars, collapse = " + "),
                    sep = " + "),
              sep = " ~ ")
        )
      }
  
    # If model hasn't changed, end the loop
    if (all(old_vars %in% new_vars)){
      no_change = T
    }
    
    # Update variables to be checked
    old_vars <- new_vars
  }
  
  if (verbose == T){
    print("Final Formula"); print(formula)
  }
  
  return(formula)
}




# This function adds each term in a regression formula to 
# the dataframe as columns (design variables) and changes the format of
# the formula to accomodate this (the formula doesn't change)

set_up_regression <- function(df, formula, vars = NA,
                              keep_vars = c("east", "north", "stno"), 
                              ...){
  
  # df              -dataframe (response variable and explanatory variables)
  # formula         -regression formula
  # vars            -only variables to be included (all included if NA)
  # keep_vars       -variables to always keep even if not in formula
  
  # LHS of formula (always kept)
  response <- as.character(formula)[2]
  df$y <- unlist(df[response])
  
  # If specific variables are mentioned, only those variables are included
  if (any(!is.na(vars))){
    df <- df %>% dplyr::select(c(response, vars, 
                                        keep_vars[!(keep_vars %in% vars)]))
  }
  
  # If the responses are NA (e.g. in test set) need to change them to 
  # non-NA value for model.matrix function
  nas_in_response <- F
  if (all(is.na(df$y))){
    nas_in_response = T
    df$y <- 1
  }
  
  if (any(!is.na(keep_vars))){
    # Get values for each term in the dataframe, and add as columns
    df <- cbind(df$y, df[keep_vars],
                as.data.frame(model.matrix(formula, df))[-1])
  }

  # Quick fix of first column name
  names(df)[1] <- c(response)
  
  # Remove placeholder response values
  if (nas_in_response == T){
    df[response] <- NA
  }
  
  return(list("df" = df, "formula" = formula))
}




# Creates formula of second order terms from formula of first order terms

second_order_formula <- function(formula, df = NA, terms = NULL, 
                                 non_formula_vars = "stno", 
                                 keep_vars = c("east", "north", "t"), ...){

  # formula             -regression formula  
  # df                  -dataframe (response and covariates)
  # terms               -specific terms to make second order 
  #                     (or interaction terms)
  # non_formula_vars    -some columns are not used as covariates (e.g. station
  #                     labels) they are treated separately
  # keep_vars           -some columns need to always be kept even if they
  #                     are not present in formula
  
  # If keep_vars not in dataframe originally, ignore them
  if (any(!is.na(df))){
    keep_vars <- keep_vars[keep_vars %in% names(df)]
    }
  
  # Get listed variables from formula
  vars <- unlist(strsplit(as.character(formula)[3], " + ", fixed = T))
  
  # If specific terms are listed, only want second-order terms specified
  if (any(!is.null(terms))){
  
    # Get variables that already exist in formula (get made second-order)
    second_order_vars <- terms[terms %in% vars]
    
    # Terms not in formula are considered as interactions
    interactions <- terms[!(terms %in% vars)]
    
    # If there are interactions, an additional plus sign needs
    # to be added to the formula
    plus_sign <- ifelse(any(!is.na(interactions)), " + ", "")

    # Create a formula with first order terms,
    # second order terms, and interactions
    
    if (any(!is.na(second_order_vars))){
    formula <- as.formula(paste(as.character(formula)[2],
                                paste(paste(paste(vars, collapse = " + "),
                                            paste("I(",
                                                  paste(second_order_vars,
                                                        collapse = "^2) + I("),
                                                  "^2)", sep = ""),
                                            sep = " + "),
                                      plus_sign,
                                      paste(interactions,
                                            collapse = " + "),
                                      sep = ""),
                                sep = " ~ ")
    )}
    
    # Small change if only including interactions (edge-case)
    else{
    formula <- as.formula(paste(as.character(formula)[2],
                                paste(paste(paste(vars, collapse = " + "),
                                            sep = " + "),
                                      plus_sign,
                                      paste(interactions,
                                            collapse = " + "),
                                      sep = ""),
                                sep = " ~ ") 
    )}
  }
  
  else{
    
    # This makes the first term in formula (i.e. "(. - stno)*(. - stno)" )
    first_term <- paste(non_formula_vars, collapse = " - ")
    first_term <- paste("(. -", first_term, ")*(. - ", first_term, ") +")
    
    # Make second order formula (with terms for each variable)
    formula <- as.formula(
      paste(as.character(formula)[2],
            paste(first_term,
                  "I(",
                  paste(vars, collapse = "^2) + I("),
                  "^2)"),
            sep = "~")
    )
  }
  
  # Option to edit data frame
  if (any(!is.na(df))){
    # We only need columns for these variables
    df <- df %>% dplyr::select(c(all_of(keep_vars[!(keep_vars %in% vars)]),
                                 as.character(formula)[2], 
                                 all_of(non_formula_vars), all_of(vars)))
  }
  
  return(list("formula" = formula, "df" = df))
}