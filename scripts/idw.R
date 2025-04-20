## Inverse distance weighting functionality

## This function is mostly a wrapper of the idw function from the 
## gstat package

## A spatio-temporal distance metric is considered, where the anisotropy
## between spatial distances and temporal distances can be estimated

IDW_ST <- function(data, newdata, response, 
                   x_name = "lon", y_name = "lat", t_name = "t",
                   idp = 2, nmax = 8, C = NULL, 
                   transformation={function (x) x},
                   reverse_transformation={function (x) x}, ...){
  
  # data            -dataframe of observed values and covariates
  # newdata         -dataframe for values to be predicted 
  # response        -variable to be interpolated
  # x_name, y_name, t_name        -labels for each coordinate
  # idp             -inverse distance weighting power
  # nmax            -maximum number of neighbours to use for prediction
  # C               -anisotropy, relates spatial and temporal distance
  
  # transformation          -transformation of response variable
  # reverse_transformation  -reverse transformation before final results
  
  require(dplyr); require(sp); require(gstat)
  
  # Transform response variable
  data[response] <- transformation(data[response])
  
  # If no anisotropy is provided, estimate it from a variogram
  # Not recommended as it is quite slow
  
  if (is.null(C)){
    # Get points and dates
    points <- SpatialPoints(cbind(data[[x_name]], data[[y_name]]))
    dates <- as.Date(as.vector(data[[t_name]]), origin = "1970-01-01")
    C <- estimate_anisotropy(df, points, dates, response, ...)
  }
  
  # Scale separation in time by anisotropy
  data[t_name] <- data[t_name] * C
  newdata[t_name] <- newdata[t_name] * C
  
  # Create spatio-temporal objects (using time as a third dimension)
  coordinates(data) <- c(x_name, y_name, t_name)
  coordinates(newdata) <- c(x_name, y_name, t_name)
  
  # Interpolate value at each point in test set
  idw_df <- idw(as.formula(paste0(response, " ~ 1")), 
                data, newdata, 
                idp = idp, nmax = nmax)
  newdata <- as.data.frame(newdata)
  newdata[response] <- idw_df$var1.pred
  
  # Reverse transformation of response variable
  newdata[response] <- reverse_transformation(newdata[response])
  
  # Reverse anisotropy
  newdata[t_name] <- newdata[t_name] / C
  
  return(newdata)
}




# Function for estimating the spatio-temporal anisotropy of a given value

estimate_anisotropy <- function(df, points, dates, response, 
                                cutoff = NULL, width = NULL, 
                                t_cutoff = 6, t_width = 1, ...){
  
  # df              -dataframe of response variable and coordinates
  # response        -variable of interest for finding anisotropy
  # points          -spatial locations
  # dates           -temporal locations
  # cutoff          -spatial cutoff of variogram
  # width           -spatial width of bins in variogram
  # t_cutoff        -temporal cutoff of variogram
  # t_width         -temporal width of variogram (recommend to keep at 1)
  
  # Load in required packages
  require(sp); require(gstat); require(spacetime)
  
  # Create spatio-temporal dataframe
  STIDF <- STIDF(points, dates, data = df)
  STSDF <- as(STIDF, "STSDF")
  
  # Heuristic for cutoff and width if none provided
  if (missing(cutoff) | missing(width)){
    cutoff <- 
      sqrt((max(points@coords[, 1]) - 
              min(points@coords[, 1]))^2 +
             (max(points@coords[, 2]) - min(points@coords[, 2]))^2) / 4
    width <- cutoff / 10
  }
  
  # Empirical VGM
  vgmst <- variogramST(as.formula(paste0(response, " ~ 1")), 
                       data = STSDF, 
                       tlags = seq(0, t_cutoff, by = t_width),
                       cutoff = cutoff, width = width, progress = T)
  
  # Calculate anisotropy
  C <- estiStAni(vgmst, c(0.1, 1e7))

  return(C)
}