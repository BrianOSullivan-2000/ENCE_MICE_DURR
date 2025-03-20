
### tools.R ###
### Miscellaneous tools ###




# Import many useful packages
import_everything <- function(){
  pkgs <- list("Rcpp", "plyr", "ggplot2", "Metrics", "sp", "sf", "gstat", 
               "MASS", "dplyr", "tidyr", "glmnet", "pls", "glmnetUtils", 
               "pbapply", "ggnewscale", "palleteer", "nlme", "regclass", 
               "spacetime",  "mice", "R.utils", "fitdistrplus", "mgcv", "VIM", 
               "parallel", "abind", "paletteer", "lattice", "plotly", "fields", 
               "proxy") 
  suppressMessages(suppressWarnings(invisible(
    lapply(pkgs, require, character.only = TRUE))))
}




# Simple checkpoint for debugging

cp <- function(i = 1){print(paste("Checkpoint", i))}




# Returns object name as string

obj_name <- function(v1) {
  deparse(substitute(v1))
}




# Returns specific column of dataframe

df_col <- function(df_name, col = NA, just_df = F){
  
  # df_name   -name of dataframe (string)
  # col       -name of desired column (string)
  # just_df   -just returns the entire dataframe if true
  
  # Just return the dataframe
  if (just_df == T){
    return(eval(parse(text = df_name)))
  }
  
  # return specific column (i.e. "df_name$col")
  else{
    return(eval(parse(text = paste0(df_name,"$",col))))
  }
}




# Put dataframe in order according to specified values (columns of dataframe)

df_order <- function(df, values){
  
  # df          -dataframe
  # values      -vector of values to sort by (e.g. c("stno", "year"))
  
  # Will order by each value in sequence
  for (value in values){
    df <- df[order(eval(parse(text = paste0(obj_name(df),"$",value)))), ]
  }
  
  return(df)
}




# Convert coordinates in longitude and latitude to easting and northing
# NOTE: long/lat coordinates are assumed to have WGS84 epsg:4326 projection
# projection for easting northing is Irish Grid TM 75 epsg:29903

lonlat_to_eastnorth <- function(longitudes, latitudes, reverse = F){
  
  require(sf); require(sp)
  
  # Get coordinates into one object
  coords <- data.frame(lon = longitudes, lat = latitudes)
  
  # If reverse is true, convert east/north to long/lat instead
  if (reverse == T){
    crs1 <- CRS("+init=epsg:29903")
    crs2 <- CRS("+init=epsg:4326")
  }
  
  # Otherwise, convert long/lat to east/north
  else{
    
    # long/lat projection
    crs1 <- CRS("+init=epsg:4326")
    
    # east/north projection
    crs2 <- CRS("+init=epsg:29903")
  }
  
  # Create spatial object with first projection, transform to second projection
  sp <- SpatialPoints(coords = coords, proj4string = crs1)
  sp <- spTransform(sp, CRSobj = crs2)
  
  # Organise transformed data
  coords <- data.frame(east = sp$lon, north = sp$lat)
  return(coords)
}




# Splits into training and testing data

tt_split <- function(df, test_index){
  
  # df                -dataframe to split
  # test_index        -indices of test data
  
  # Train and testing dataframes
  train_df <- df[-test_index, ]
  test_df <- df[test_index, ]
  
  # Remove measured values from test set
  test_values <- test_df$y
  
  if ("t" %in% colnames(test_df)){
    test_values <- test_df[order(test_df$stno, test_df$t), ]$y
  }
  test_df$y <- NA
  
  return(list("train_df" = train_df, "test_df" = test_df,
              "test_values" = test_values))
}




# Predict values at given indices using specified function

predict_each_value <- function(test_index, df, func, response, 
                               extra_data = data.frame(), split_var = NULL, 
                               ...){
  
  # test_index      -rows in dataframe to be predicted
  # df              -dataframe (response variable and explanatory variables)
  # func            -function to predict response variables
  # response        -variable to be predicted (rain, temperature, etc.)
  # extra_data      -additional data to add to training set
  # split_var        -remove rows from extra_data if split_var matches test set
  # ...             -additional arguments for prediction function
  
  # Split data into training and test sets
  s <- tt_split(df, test_index)
  
  train_df <- s$train_df
  test_df <- s$test_df
  
  # Add extra data to training data (eg data from other days)
  if(!missing(extra_data)){
    # Set name of response value to y
    extra_data <- extra_data %>% dplyr::rename(y = all_of(response))
    
    if(!is.null(split_var)){
      extra_data <- extra_data[!(unlist(extra_data[split_var]) %in% 
                                   unlist(test_df[split_var])), ]
    }
    train_df <- bind_rows(train_df, extra_data)
  }
  
  # Mark test data (for functions where all data gets puts together)
  train_df$test <- 0
  test_df$test <- 1
  
  # Output copy of test set with predicted values
  pred_df <- func(train_df, test_df, ...)
  
  if(any(pred_df$test == 0)){
    pred_df <- pred_df[pred_df$test == 1, ]
  }
  
  if ("t" %in% colnames(pred_df)){
    pred_df <- pred_df[order(pred_df$stno, pred_df$t), ]
  }
  
  # Of course you can't trust an RMSE that still has NA, but might
  # be useful for debugging
  if (is.na(rmse(s$test_values, pred_df$y))){
    print("Warning: there are NAs among predicted values")
    print("RMSE")
    print(rmse(s[!is.na(pred_df$y), ]$test_values, 
               pred_df[!is.na(pred_df$y), ]$y))
  }
  
  else{
    print("RMSE")
    print(rmse(s$test_values, pred_df$y))
  }
  
  return(as.data.frame(pred_df))
}




# Generic Cross Validation Function (Leave One Out or k-fold)

CV <- function(df, func, response, jobID = NULL,
               LOOCV = T, kfolds = 10, split_var = NULL, 
               check_by_year = F, check_by_year_var = NA,
               check_by_years = F, nyears = 10, nstations = 20, 
               check_by_final_year = F,
               
               progress = T, scale = F, scale_vars = NA,
               
               plot = F, plot_title = "", plot_dst = NA, 
               get_results = F, get_first = F, ...){
  
  # df                      -dataframe (response and explanatory variables)
  # func                    -function to predict response variables
  # response                -variable to be predicted (rain, temperature, etc.)
  
  # LOOCV                   -if true, use Leave One Out Cross Validation
  # kfolds                  -if LOOCV is false, number of folds for k-fold CV
  # split_var               -can split into training and testing sets 
  #                         according to specific variables (e.g. stno)
  # check_by_year           -If true, an entire year(s) will be removed for a 
  #                         single station to be checked during CV
  # check_by_year_var       -Remove years for randomly selected values of var
  # check_by_years          -Run cross validation on sequence of years for
  #                         randomly selected stations
  # nyears                  -number of years to run check_by_years
  # nstations               -number of stations to run check_by_years
  # check_by_final_year     -check_by_years, but only the last year considered
  
  # progress                -include progress bar when running
  # scale                   -scale explanatory variables
  # scale_vars              -can choose to only scale certain variables
  
  # plot                    -plot observed values against predicted
  # plot_title              -title for results plot
  # plot_dst                -destination to save plot file
  # get_results             -return results (RMSE, R2, predicted dataframe)
  # get_first               -can just check RMSE of first fold
  # ...                     -additional arguments for prediction function
  
  # Load in packages
  require(dplyr); require(Metrics); require(pbapply)
  
  # Print name of job (for HPC)
  if (!is.null(jobID)){
    print(paste0("Running job: ", jobID))
  }
    
  # Start a timer for CV
  st_CV = Sys.time()
  
  # Scale values
  if (scale == T){
    
    # Scale variables specified by scale_var
    if (!all(is.na(scale_vars))){
      df[scale_vars] <- scale(df[scale_vars])
    }
    
    # Scale all variables if scale_var is empty
    # (except response variable and station number)
    else{
      df[!(names(df) %in% c(response, "stno"))] <- 
        scale(df[!(names(df) %in% c(response, "stno"))])
    }
  }
  
  # Set name of response value to y
  df <- df %>% dplyr::rename(y = all_of(response))
  missing_indices <- as.vector(is.na(df[,"y"]))
  
  # When a column is specified in split_var, all instances of each unique
  # value in that column will be removed when splitting into testing and 
  # training sets (for example, all of a stations measurements will be added
  # to the testing set at the same time)
  
  if (!is.null(split_var)){
    
    # Get each unique split variable (e.g. each unique station)
    split_vals <- unlist(df[split_var])
    vals <- unique(split_vals)
    
    # Split unique values chunks if LOOCV is not selected
    if (LOOCV == F){
      rows <- sample(vals)
      vals <- split(rows,
                    cut(seq_along(rows),
                        kfolds,
                        labels = FALSE))
    }
    
    # Test indices are grouped together according to unique values
    test_indices <- lapply(vals, {function(val, x) which(x %in% val)},
                           x = split_vals)
  }
  
  else if (check_by_year){
    
    # Get each possible combination of year and station
    years <- unique(unlist(df["year"]))
    split_vals <- unique(unlist(df[check_by_year_var]))
    vals_by_year <- expand.grid(years, split_vals)
    names(vals_by_year) <- c("year", check_by_year_var)
    
    # Split year/var combinations randomly into k folds
    if (LOOCV == F){
      rows <- sample(nrow(vals_by_year))
      vals <- split(rows,
                    cut(seq_along(rows),
                        kfolds,
                        labels = FALSE))
      vals <- lapply(seq(kfolds), {function(x) vals_by_year[vals[[x]], ]})
    }
    
    # Get indices corresponding to each year/var combination in data
    js <- lapply(vals,
                 {function(x) join.keys(df, x, by = c("year", 
                                                      check_by_year_var))})
    
    test_indices <- lapply(js, {function(j) which(j$x %in% j$y)})
    
    # Don't check values with NA response
    test_indices <- lapply(test_indices, {function(t_i) 
      t_i[t_i %in% which(!is.na(df["y"]))]})
  }
  
  else if (check_by_years){
    
    # Randomly select stations to check
    random_stations <- sample(x = unique(df$stno), size = nstations)
    
    # Get sequences of years starting at random points for each station
    years <- seq(min(unique(df$year)), max(unique(df$year)))
    start_years <- sample(seq(min(df$year), max(df$year)  - (nyears - 1)), 
                          size = nstations, replace = T)
    years_to_remove <- lapply(start_years, 
                              {function(start_year) 
                                seq(start_year, start_year + (nyears - 1))})
    
    # Get corresponding indices to run tests on for each station/years group
    test_indices <- 
      lapply(seq(nstations), 
             {function(i) which((df$stno %in% random_stations[i]) & 
                                  (df$year %in% years_to_remove[[i]]))})
    
    # Only want to validate against non-missing values
    test_indices <- lapply(test_indices, {function(t_i)
      t_i[t_i %in% which(!is.na(df["y"]))]})
  }
  
  
  else if (check_by_final_year){
    
    # Years to consider (sequence from end of time period)
    years <- seq(max(unique(df$year)) - nyears + 1, max(unique(df$year)))
    
    # Stations with complete data during target period
    candidate_stations <- (df[df$year %in% years, ] %>% 
                             group_by(stno) %>% filter(sum(is.na(y)) == 0))$stno
    
    # Randomly select stations to test against
    random_stations <- sample(x = candidate_stations, size = nstations)
    
    # Get corresponding indices to run tests on for each station/years group
    test_indices <- 
      lapply(seq(nstations), 
             {function(i) which((df$stno %in% random_stations[i]) & 
                                  (df$year %in% years))})
  }
  
  
  # Primary method for getting test indices (i.e. split_var = NULL)
  else{
    # Indices are each row of dataframe if LOOCV is True
    if (LOOCV == T){
      test_indices <- 1:nrow(df)
      
      # Don't include any missing entries in CV
      if (any(is.na(df[,"y"]))){
        test_indices <- test_indices[!is.na(df[,"y"])] 
      }
    }
    
    else{
      # Randomly select k equally sized chunks of dataframe for k-fold CV 
      rows <- sample(1:nrow(df))
      
      if (any(is.na(df[,"y"]))){
        rows <- sample((1:nrow(df))[!is.na(df[,"y"])])
      }
      
      test_indices <<- split(rows,
                            cut(seq_along(rows),
                                kfolds,
                                labels = FALSE))
    }
  }
  
  # Get RMSE of only one fold (ie get_first = 5 for fifth fold)
  if(get_first){
    test_indices <- test_indices[get_first]
  }
  
  # With progress bar
  if (progress == T){
    
    # For each set of test indices, predict response variable using 
    # the rest of the data
    pred_df <- pbapply::pblapply(test_indices, predict_each_value,
                                 df = df, func = func, response = response, 
                                 split_var = split_var, ...)
  }
  
  # No progress bar
  else{
    
    # For each set of test indices, predict response variable using 
    # the rest of the data  
    pred_df <- lapply(test_indices, predict_each_value,
                      df = df, func = func, response = response, 
                      split_var = split_var, ...)
  }
  
  # Group all predicted rows together
  pred_df <- bind_rows(pred_df)
  
  # Put predicted dataframe back in order of station
  if (LOOCV == F){
    pred_df <- pred_df[order(pred_df$stno), ]
    
    # Order by station AND time if using spatio-temporal data
    if ("t" %in% colnames(pred_df)){
      pred_df <- pred_df[order(pred_df$stno, pred_df$t), ]
    }
  }
  
  # Now there are predicted values, remove any rows that were initially missing
  # if (length(missing_indices) > 0){
  #   df <- df[!missing_indices, ]
  # }
  
  df <- df[unlist(test_indices), ]
  df <- df[order(df$stno, df$t), ]
  
  # Root Mean Squared Error between observed and predicted values
  RMSE <- rmse(df$y, pred_df$y)
  
  # Normalised Root Mean Squared Error
  RMSEr <- RMSE / sd(df$y)
  
  # Root Relative Squared Error between observed and predicted values
  RRSE <- rrse(df$y, pred_df$y)
  
  # Coefficient of determination
  R2 <- cor(df$y, pred_df$y) ^ 2
  
  # Plot observed vs. predicted values
  if (plot == T){
    
    # Save plot is a filename is specified
    if(!is.na(plot_dst)){
      jpeg(file = plot_dst)
    }
    
    plot(df$y, pred_df$y, xlab = "Observed", ylab = "Predicted",
         main = plot_title)
    abline(a=0, b=1)
    
    if(!is.na(plot_dst)){
      dev.off()
    }
  }
  
  print(noquote(paste("Root Mean Squared Error = ", round(RMSE, 4))))
  print(noquote(paste("Relative Root Mean Squared Error = ",
                      round(RMSEr, 4))))
  print(noquote(paste("Pearson Correlation Coefficient (R-squared) = ",
                      round(R2, 4))))
  et_CV <- Sys.time()
  print(noquote(paste("Total time for CV: ", et_CV - st_CV, 
                      " ", units(et_CV - st_CV))))
  
  # Return results if get_results is true
  if (get_results == T){
    
    # Dataframe of predicted values
    pred_df <- as.data.frame(pred_df)
    names(pred_df)[names(pred_df) == "y"] <- paste0("predicted_", response)
    pred_df <- cbind(pred_df, df$y)
    names(pred_df)[ncol(pred_df)] <- paste0("observed_", response)
    df <- pred_df
    
    return(list("RMSE" = RMSE, "RMSEr" = RMSEr, 
                "RRSE" = RRSE, "R2" = R2, "df" = df))
  }
}




# Simple function to plot gridded data over Ireland

plot_grid <- function(grid, response = "y", crs = CRS("+init=epsg:29903"),
                      
                      plot_ireland = T, inside_border = F,
                      ireland_dsn = "datasets/counties.json",
                      
                      added_components = NULL, xlim = NULL, ylim = NULL,
                      color_scale = "grDevices::Blues", reverse = F,
                      
                      plot_dsn = NULL, height = 15, width = 12, 
                      dpi = 300, display = T){
  
  # grid              -dataframe containing response variable and coordinates
  # response          -variable to be plotted
  # crs               -projection of grid (recommend default Irish Grid)
  
  # plot_ireland      -add Irish Border to plot
  # inside_border     -only include points within border
  # ireland_dsn       -location of Ireland shapefile
  
  # added_components  -additional components can be added to ggplot
  #                   (Argument must be passed as list of each component)
  # xlim, ylim        -plot limits (in longitude and latitude)
  # color_scale       -color scheme to use for plot 
  #                   (see paletteer_d_names for options)
  # reverse           -Change change color direction of color scheme
  
  # plot_dsn          -destination to save plotted image
  # height            -plot height
  # width             -plot width
  # dpi               -plot resolution
  # display           -option to print plot
  
  # Load in packages
  require(dplyr); require(ggplot2); library(paletteer); require(sf)
  require(ggnewscale)
  
  # Add Irish Border if specified
  if (plot_ireland == T){
    ireland <- read_sf(dsn = ireland_dsn)
    ireland <- st_transform(x = ireland, crs = crs)
  }
  
  # Confine grid within Ireland's Border
  if (inside_border == T){
    
    if(missing(ireland)){
      ireland <- read_sf(dsn = ireland_dsn)
      ireland <- st_transform(x = ireland, crs = crs)
    }
    
    # Turn grid into spatial object
    coordinates(grid) <- c("east", "north")
    grid <- st_as_sf(grid)
    st_crs(grid) <- crs
    
    # Only keep points within border
    s <- st_contains(ireland, grid)
    grid <- grid[unlist(s), ]
    
    # Turn back into standard dataframe
    coords <- st_coordinates(grid)
    grid <- as.data.frame(grid) %>% dplyr::select(-c("geometry"))
    grid$east <- coords[, 1]; grid$north <- coords[, 2]
  }
  
  # Make ggplot component for Ireland's borders
  if (plot_ireland == T){
    ireland <- geom_sf(data = ireland, col = "black", fill = NA)
  }
  
  # Warning, must include both x and y limits or neither
  if (xor(!missing(xlim), !missing(ylim))){
    stop("If specifying plot limits, limits are required for longitude
         AND latitude")
  }
  
  # Convert longitude and latitude to easting and northing
  if(!missing(xlim) & !missing(ylim)){
    limits <- lonlat_to_eastnorth(xlim, ylim)
    xlim <- limits$east; ylim <- limits$north
  }
  
  # Set direction of color scheme
  d <- ifelse(reverse, -1, 1)
  
  # Truth value indicating if additional components are included
  addcomp <- ifelse(missing(added_components), F, T)
  
  # Save response as an expression (which is then evaluated in ggplot)
  response_var <- enquo(response)
  
  # Create the plot
  grid_plot <- ggplot() +
    
    # Plot the grid
    geom_tile(data = grid, aes(x = east, y = north, fill = {{response}})) +
    scale_fill_paletteer_c(color_scale, direction = d) +
    
    # Simple theme settings
    theme(plot.title = element_text(size = 20, hjust = 0.5),
          panel.background = element_rect(fill = "grey80"),
          panel.grid = element_blank(), 
          axis.title = element_text(size = 12), 
          axis.text = element_blank()) +
    
    # Additional components (irish border, plot limits, etc.)
    {if(plot_ireland) ireland} +
    {if(addcomp) added_components} +
    {if(!missing(xlim)) xlim(xlim)} +
    {if(!missing(ylim)) ylim(ylim)}
  
  # Print and save 
  if(display){print(grid_plot)}
  if(!missing(plot_dsn)){ggsave(plot_dsn, 
                                height = height, width = width, dpi = dpi)}
}




# Visualise a matrix, designed to identify any patterns 
# or matrix structure, you may need to round the matrix before plotting
# (i.e. round(M, 10) should be sufficient)

matrix_map <- function(M, no_unique = F){
  
  # M         - Matrix to visualise
  # no_unique - print the number of unique values in the matrix
  
  # Load in packages
  require(plot.matrix)
  
  # 702 letter pairs (You shouldn't need more that 702)
  LETTERS702 <- c(LETTERS, sapply(LETTERS, function(x) paste0(x, LETTERS)))
  
  # Group equal values using factors
  f <- factor(M, levels = unique(as.numeric(M))) 
  levels(f) <- LETTERS702[1:length(levels(f))]
  
  if(no_unique){print(length(levels(f)))}
  
  # Make factor matrix and plot
  f <- matrix(f, nrow = nrow(M))
  plot(f, col = rainbow, key = NULL, digits = 2)
}




# Visualise two matrices side by side for comparison 
# You may need to round the matrix before plotting
# (i.e. round(M, 10) should be sufficient)

matrix_map_compare <- function(A, B){
  
  # A, B  - Two matrices to be plotted (must have same dimensions)
  
  # Load in packages
  require(plot.matrix)
  
  # 702 letter pairs (You shouldn't need more that 702)
  LETTERS702 <- c(LETTERS, sapply(LETTERS, function(x) paste0(x, LETTERS)))
  
  # Apeend matrices together with blank column in between
  M <- cbind(A, rep(NA, nrow(A)), B)
  
  # Group equal values together using factors
  f <- factor(M, levels = unique(as.numeric(M))) 
  levels(f) <- LETTERS702[1:length(levels(f))]
  
  # Plot factor matrix
  f <- matrix(f, nrow = nrow(M))
  plot(f, col = rainbow, key = NULL, digits = 2)
}

