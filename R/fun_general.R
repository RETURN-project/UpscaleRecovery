#' Load RData file and returns it. This function is a substitute for the load function, allowing to assign a user defined variable name when loading a RData file.
#'
#' @param fileName the path to the Rdata file that needs to be loaded
#'
#' @return R object
#' @export
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

#' Convert time series to annual frequency: The toAnnualTS function converts a time series with n observations per year to an annual time series (one observation per year). The main concept is to select observations per year closest to a given day of year that have no missing value (NA). Here, the day of year for which the seasonality is maximum is being used.
#'
#' @param tsseas vector of observations (time series) representing the seasonal component of the time series to be converted
#' @param tsi vector of observations (time series) that needs to be converted to an annual time series
#' @param obspyr number of observations per year of the time series to be converted
#'
#' @return vector of observations (time series) with annual observation frequency
#' @export
 toAnnualTS <- function(tsseas, tsi, obspyr){
    seasi <- rowMeans(matrix(tsseas, nrow = obspyr),na.rm=T)# average seasonality
    smax <- which(seasi == max(seasi, na.rm=T))# yearly observation number with max seas
    tsmi <- matrix(tsi, nrow = obspyr)
    dst <- abs(matrix(rep(1:obspyr,times = (length(tsi)/obspyr)), nrow = obspyr)-smax)# distance of observations to seasonal max
    dst[is.na(tsmi)] <- NA #set distance of NA observations equal to NA, these can not be selected
    rsel <- as.matrix(apply(dst, 2, which.min))# row numbers of observations to be selected, i.e. those closest to seasonal max
    toNA <- unlist(lapply(rsel, identical, integer(0)))# years without observation: assign temporary the first observation of the year
    rsel[toNA] <- 1
    rsel <- unlist(rsel) # rows to be selected
    csel <- 1:dim(tsmi)[2] # columns to be selected
    tsyr <- tsmi[rsel + nrow(tsmi) * (csel - 1)]# get values of yearly time series
    tsyr[toNA] <- NA# years without observation: set value to NA
    tsyr
}

#' Create regular time series
#'
#' @param tsi vector of observations
#' @param dts dates associated with the observations. This should be a Date object.
#' @param fun function used to aggregate observations to monthly observations. This should be 'max' or 'mean'.
#' @param resol desired temporal resolution of the output. This could be 'monthly' or 'daily'
#'
#' @return a vector with a regular time series object
#' @export
#' @import zoo
#' @import bfast
toRegularTS <- function(tsi, dts, fun, resol){
  tsi <- as.numeric(tsi)
  if(resol == 'monthly'){
    z <- zoo(tsi, dts) ## create a zoo (time) series
    if(fun == 'max'){
      mz <- as.ts(aggregate(z, as.yearmon, mmax)) ## max
    }
    if(fun == 'mean'){
      mz <- as.ts(aggregate(z, as.yearmon, mean)) ## mean
    }
  }else if (resol == 'daily'){
    mz <- bfastts(tsi, dts, type = "irregular")
  }else if (resol == 'quart'){
    z <- zoo(tsi, dts) ## create a zoo (time) series
    if(fun == 'max'){
      mz <- as.ts(aggregate(z, as.yearqtr, mmax)) ## max
    }
    if(fun == 'mean'){
      mz <- as.ts(aggregate(z, as.yearqtr, mean)) ## mean
    }
  }
   return(mz)
 }


#' Convert a raster stack with irregular time observations to regular time series
#'
#' @param x stack of observations, the first raster is a mask (with values 0 and 1). The pixels of the mask that have a value 1 are processed, the pixels with a 0 value are not processed.
#' @param dts dates of the observations
#' @param fun unction used to aggregate observations to monthly observations. This should be 'max' or 'mean'.
#' @param resol desired temporal resolution of the output. This could be 'monthly' or 'daily'
#'
#' @return stack with regular observations
#' @export
toRegularTSStack <- function(x, dts, fun, resol)
{
   msk <- x[,1]
   msk[is.na(msk)] <- 0
   x <- x[,-1]
   mask <- apply(x, 1, FUN = function(x) { sum(is.na(x)) / length(x) } )
   i <- ((mask < 1) & (msk == 1))
  len <- length(toRegularTS(dts, dts, fun=fun, resol = resol))
  res <- matrix(NA, length(i), len)
  if(sum(i) == 1) {
    res[i,] <- toRegularTS(x[i,], dts, fun=fun, resol = resol)
  } else if(sum(i) > 1) {
    res[i,] <- t(apply(x[i,], 1, toRegularTS, dts, fun, resol))
  }
  res
}

#' Helper function for the toRegularTSStack function
#'
#' @param x vector of observations
#'
#' @return the maximum value of the vector
#' @export
mmax <- function(x) {
  if(length(which.max(x)) == 0) {
    out <- NA
  } else {
    out <- as.numeric(x[which.max(x)])
  }
  return(out)
}


#' Convert CCI fire stack to a stack with a predefined temporal resolution and containing the value 0 when no fire is present and 1 if a fire is present
#'
#' @param x stack of CCI fire imagery, the first layer contains a mask (here the value 0 denotes that the pixel should not be included, the value 1 denotes that the pixel should be included), followed by n fire confidence layers and n fire doy layers
#' @param dts dates associated with the stack
#' @param resol the desired temporal resolution of the output data
#' @param thres threshold on the fire confidence layer, only fires having a confidence higher than the threshold are included
#'
#' @return a stack with with a predefined temporal resolution and containing the value 0 when no fire is present and 1 if a fire is present
#' @export
createFireStack <- function(x, dts, resol, thres)
{
  msk <- x[,1]
  x <- x[,-1]
  i <- (msk == 1)
  strtyr <- format(dts[1], "%Y")
  endyr <- format(dts[length(dts)], "%Y")
  if(resol == 'monthly'){
    len <- length(dts)}else if(resol == 'daily'){
      len <- length(seq(as.Date(paste0(strtyr,'-01-01')), as.Date(paste0(endyr,'-12-31')), by = "1 day"))
    }else if (resol == 'quart'){
      len <- length(seq(as.Date(paste0(strtyr,'-01-01')), as.Date(paste0(endyr,'-12-31')), by = "3 months"))
    }
  res <- matrix(NA, length(i), len)
  if(sum(i) == 1) {
    res[i,] <- toFireTS(x[i,], dts, resol = resol, thres = thres)
  } else if(sum(i) > 1) {
    res[i,] <- t(apply(x[i,], 1, toFireTS, dts, resol, thres))
  }
  res
}

#' Create a time series with the desired temporal resolution that has the value 0 if no fire occurs and the value 1 if a fire occurs
#'
#' @param x a time series, the first n values contain the fire confidence values, the subsequent n values the doy of the fires
#' @param dts date values associated with the time series observations (Date object)
#' @param resol desired temporal resolution of the output
#' @param thres threshold on the fire confidence (value between 0 and 100). The higher the confidence value, the more likely the fire really occured. Only fires with a value above the threshold are taken into account.
#'
#' @return vector with 0 if no fire occurs and 1 if a fire occurs
#' @export
toFireTS <- function(x, dts, resol, thres = 95){
  x <- as.numeric(x)
  len <- length(x)
  cl <- x[1:(len/2)] # confidence values
  jd <- x[(1+(len/2)):len] # doy of fire
  if(resol == 'monthly'){
    out <- rep(0,length(cl)) # initialise a vector of zeros
    out[cl>thres & jd>0] <- 1 # dates with high fire confidence and doy of fire > 0 are set to 1
  }else if (resol == 'daily'){
    # create daily time series and associated dates
    strtyr <- format(dts[1], "%Y") # start year of the fire dataset
    endyr <- format(dts[length(dts)], "%Y") # end year of the fire dataset
    tsdts <- seq(as.Date(paste0(strtyr,'-01-01')), as.Date(paste0(endyr,'-12-31')), by = "1 day")  # all potential fire dates
    out <- rep(0,length(tsdts))# initialise output vector with zeros

    # get timing of fires
    ind <- which(cl>thres & jd>0) # find observations with high fire confidence and doy of fire > 0
    fireyr <- format(dts[ind], "%Y") # year of the observed fires
    firedoy <- jd[ind] # doy of the observed fires
    firedate <- as.Date(paste0(fireyr,'-',firedoy),'%Y-%j')# create fire observation dates
    out[tsdts %in% firedate] = 1# set fire observation dates to 1
  }else if (resol == 'quart'){
    outm <- rep(0,length(cl)) # initialise a vector of zeros
    outm[cl>thres & jd>0] <- 1 # dates with high fire confidence and doy of fire > 0 are set to 1
    out <- toRegularTS(outm, dts, 'max', 'quart')
  }
  out
}

#' Calculate the Normalized Burn ratio
#'
#' @param NIR Near Infrored reflectance
#' @param SWIR Shortwave Infrared reflectance
#'
#' @return Normalized Burn Ratio
#' @export
calcNBR <- function(NIR,SWIR){
  nbr <- (NIR-SWIR)/(NIR+SWIR)
  return(nbr)
}

#' Mask an image using its quality layer, quality parameters set the quality preferences for pixels that should NOT be masked.
#'
#' @param im the image that should be masked
#' @param qaim the quality layer
#' @param valid quality parameter - valid data, values can be 0 (valid) or 1 (no data)
#' @param cloud quality parameter - cloud state, values can be 0 (clear), 1 (less confident cloud, i.e. buffered cloud 300m), 2 (confident, opaque cloud), or 3 (cirrus)
#' @param shadow quality parameter - cloud shadow flag, values can be 0 (no shadow), or 1 (shadow)
#' @param snow quality parameter - snow flag, values can be 0 (no snow), or 1 (snow)
#' @param water quality parameter - water flag, values can be 0 (no water), or 1 (water)
#' @param aero quality parameter - aerosol state, values can be 0 (estimated, best quality), 1 (interpolated, mid quality), 2 (high, aerosol optical depth > 0.6, use with caution)
#' @param subzero quality parameter - subzero flag, values can be 0 (no subzero values), or 1 (subzero values)
#' @param sat quality parameter - saturation flag, values can be 0 (no saturation), or 1 (saturation)
#' @param sunZen quality parameter - high sun zenith flag, values can be 0 (no high sun zenith), or 1 (sun elevation < 15 degrees, use with caution)
#' @param illum quality parameter - illumination state, values can be 0 (good, incidence angle <55°, best quality for topographic correction), 1 (medium, incidence angle 55°-80°, low quality for top. correction), 2 (poor, incidence angle > 80°, low quality for top. correction), or 3 (shadow, incidence angle > 90°, no topographic correction applied)
#' @param slope quality parameter - slope flag, values can be 0 (cosine correction applied) or 1 (enchanced C-correction applied)
#' @param wvp quality parameter - water vapor flag, values can be 0 (measured, best quality, only for Sentinel-2) or 1 (fill, scene average, only for Sentinel-2)
#'
#' @return a masked image
#' @export
#'
mskQA <- function(im, qaim, valid = 0, cloud = 0, shadow = 0, snow = 0, water = 0, aero = 0, subzero = 0, sat = 0, sunZen = 0, illum = 0, slope = 0, wvp = 0){
  # QA values of pixels that shouldn't be masked
  toKeep <- getQAvals(valid, cloud, shadow, snow, water, aero, subzero, sat, sunZen, illum, slope, wvp)
  # make raster with 0 values for pixels that should be masked and 1 for pixels that should be kept
  msk <- qaim
  msk[] <- 0
  for(i in 1:length(toKeep)){
    msk[qaim == toKeep[i]] <- 1
  }
  # mask values
  im[msk == 0] <- NA
  return(im)
}

#' Generate QA values for all possible combinations of quality parameter values. Here, multiple values for one parameter can be provided.
#'
#' @param valid quality parameter - valid data, values can be 0 (valid) or 1 (no data)
#' @param cloud quality parameter - cloud state, values can be 0 (clear), 1 (less confident cloud, i.e. buffered cloud 300m), 2 (confident, opaque cloud), or 3 (cirrus)
#' @param shadow quality parameter - cloud shadow flag, values can be 0 (no shadow), or 1 (shadow)
#' @param snow quality parameter - snow flag, values can be 0 (no snow), or 1 (snow)
#' @param water quality parameter - water flag, values can be 0 (no water), or 1 (water)
#' @param aero quality parameter - aerosol state, values can be 0 (estimated, best quality), 1 (interpolated, mid quality), 2 (high, aerosol optical depth > 0.6, use with caution)
#' @param subzero quality parameter - subzero flag, values can be 0 (no subzero values), or 1 (subzero values)
#' @param sat quality parameter - saturation flag, values can be 0 (no saturation), or 1 (saturation)
#' @param sunZen quality parameter - high sun zenith flag, values can be 0 (no high sun zenith), or 1 (sun elevation < 15 degrees, use with caution)
#' @param illum quality parameter - illumination state, values can be 0 (good, incidence angle <55°, best quality for topographic correction), 1 (medium, incidence angle 55°-80°, low quality for top. correction), 2 (poor, incidence angle > 80°, low quality for top. correction), or 3 (shadow, incidence angle > 90°, no topographic correction applied)
#' @param slope quality parameter - slope flag, values can be 0 (cosine correction applied) or 1 (enchanced C-correction applied)
#' @param wvp quality parameter - water vapor flag, values can be 0 (measured, best quality, only for Sentinel-2) or 1 (fill, scene average, only for Sentinel-2)
#'
#' @return vector of QA values for the combinations of the given parameter values
#' @export
#'
getQAvals <- function(valid = 0, cloud = 0, shadow = 0, snow = 0, water = 0, aero = 0, subzero = 0, sat = 0, sunZen = 0, illum = 0, slope = 0, wvp = 0){
  # make all possible combinations of settings
  cmbs <- expand.grid(valid, cloud, shadow, snow, water, aero, subzero, sat, sunZen,illum, slope, wvp)
  names(cmbs) = c('valid', 'cloud', 'shadow', 'snow', 'water', 'aero', 'subzero', 'sat', 'sunZen', 'illum', 'slope', 'wvp')
  # convert to dataframe
  mcmbs <- as.data.frame(cmbs)
  names(mcmbs) = c('valid', 'cloud', 'shadow', 'snow', 'water', 'aero', 'subzero', 'sat', 'sunZen', 'illum', 'slope', 'wvp')
  # get QA values of all inputs
  out <- mapply(getQA,mcmbs$valid, mcmbs$cloud, mcmbs$shadow, mcmbs$snow, mcmbs$water, mcmbs$aero, mcmbs$subzero, mcmbs$sat, mcmbs$sunZen, mcmbs$illum, mcmbs$slope, mcmbs$wvp)
  return(out)
}

#' Generate the QA value associated with specific quality settings
#'
#' @param valid quality parameter - valid data, values can be 0 (valid) or 1 (no data)
#' @param cloud quality parameter - cloud state, values can be 0 (clear), 1 (less confident cloud, i.e. buffered cloud 300m), 2 (confident, opaque cloud), or 3 (cirrus)
#' @param shadow quality parameter - cloud shadow flag, values can be 0 (no shadow), or 1 (shadow)
#' @param snow quality parameter - snow flag, values can be 0 (no snow), or 1 (snow)
#' @param water quality parameter - water flag, values can be 0 (no water), or 1 (water)
#' @param aero quality parameter - aerosol state, values can be 0 (estimated, best quality), 1 (interpolated, mid quality), 2 (high, aerosol optical depth > 0.6, use with caution)
#' @param subzero quality parameter - subzero flag, values can be 0 (no subzero values), or 1 (subzero values)
#' @param sat quality parameter - saturation flag, values can be 0 (no saturation), or 1 (saturation)
#' @param sunZen quality parameter - high sun zenith flag, values can be 0 (no high sun zenith), or 1 (sun elevation < 15 degrees, use with caution)
#' @param illum quality parameter - illumination state, values can be 0 (good, incidence angle <55°, best quality for topographic correction), 1 (medium, incidence angle 55°-80°, low quality for top. correction), 2 (poor, incidence angle > 80°, low quality for top. correction), or 3 (shadow, incidence angle > 90°, no topographic correction applied)
#' @param slope quality parameter - slope flag, values can be 0 (cosine correction applied) or 1 (enchanced C-correction applied)
#' @param wvp quality parameter - water vapor flag, values can be 0 (measured, best quality, only for Sentinel-2) or 1 (fill, scene average, only for Sentinel-2)
#'
#' @return QA value for the given parameter value setting
#' @export
#'
getQA <- function(valid = 0, cloud = 0, shadow = 0, snow = 0, water = 0, aero = 0, subzero = 0, sat = 0, sunZen = 0, illum = 0, slope = 0, wvp = 0){
  # convert the two bit numbers
  bcloud <- paste(as.numeric(intToBits(cloud)[1:2]), collapse = '')
  baero <- paste(as.numeric(intToBits(aero)[1:2]), collapse = '')
  billum <- paste(as.numeric(intToBits(illum)[1:2]), collapse = '')
  # combine bits to bitword
  bitword <- paste0(0,wvp,slope,billum,sunZen,sat,subzero,baero,water,snow,shadow,bcloud, valid)
  # convert bitword to integer
  out <- BinToDec(bitword)
  return(out)
}

#' Convert a binary number to integer
#'
#' @param x the binary number
#'
#' @return an integer
#' @export
#'
#' @examples
BinToDec <- function(x){
  sum(2^(which(rev(unlist(strsplit(as.character(x), "")) == 1))-1))
}

#' Fix the time span of a multi-temporal raster stack to a user defined time window
#'
#' @param br raster stack
#' @param starttime start date of the desired time span (given as a numeric vector containing the year, month and day)
#' @param endtime end date of the desired time span (given as a numeric vector containing the year, month and day)
#' @param tempRes temporal resolution of the raster stack
#' @param dtsbr vector of dates (Date object) associated with the raster stack
#'
#' @return raster stack with adjusted time span
#' @import lubridate
#' @export
#'
setPeriod <- function(br, starttime, endtime,tempRes, dtsbr){
  tres <- c('1 month', '3 months', '1 day', '1 year')
  names(tres) <- c('monthly', 'quarterly', 'daily', 'yearly')
  # make sure that image stack covers time period of interest
  rstNA <- br[[1]]
  rstNA[] <- NA # empty image
  startyr <- as.Date(paste0(starttime[1],'-',starttime[2],'-',starttime[3])) # create date object from start date
  endyr <- as.Date(paste0(endtime[1],'-',endtime[2],'-',endtime[3]))# create date object from end date
  dtstot <- seq(startyr, endyr, by = tres[tempRes])

  out <- br
  out <- out[[which((dtsbr>=min(dtstot)) & (dtsbr<=max(dtstot)))]]# remove observations that fall outside the study period
  npre <- sum(dtstot<min(dtsbr))# missing dates at the start of the study period
  npost <- sum(dtstot>max(dtsbr))# missing dates at the end of the study period
  if(npre>0){for(i in 1:npre){out <- addLayer(rstNA, out)}}
  # else{# add observations at the beginning if needed
    # remove observations before the start of the study period if needed
  #   ind <- which(dtsbr>=min(dtstot))
  #   dtsbr <- dtsbr[ind]
  #   out <- out[[ind]]
  # }
  if(npost>0){for(i in 1:npost){out <- addLayer(out,rstNA)}}
  # else{# add observations at the end if needed
  #   # remove observations after the end of the study period if needed
  #   ind <- which(dtsbr<=max(dtstot))
  #   out <- out[[ind]]
  # }
  names(out) <- dtstot
  return(out)
}
