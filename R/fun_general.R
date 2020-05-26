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

#' The TScompress function compresses a time series vector. The purpose is to avoid storing many NA values.
#'
#' @param ts a vector that will be compressed
#'
#' @return a vector containing, in order, the length of the input vector, the number of observations without NA, the observations without NA, and the values of the observations that are no NA
#' @export
TScompress <- function(ts){
    c(length(ts), length(which(!is.na(ts))), which(!is.na(ts)), ts[!is.na(ts)])
}

#' The TSdecompress recovers a vector that has been compressed by the TScompress function in its original format.
#'
#' @param ts a vector compressed be the TScompress function
#'
#' @return the vector restored in its original format
#' @export
TSdecompress <- function(ts){
    vec <- rep(NA,ts[1])
    if(ts[2] > 0){
        vec[ts[3:(ts[2]+2)]] <- ts[(ts[2]+3):((2*ts[2])+2)]
    }
    vec
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
  # len <- length(tsi)
  # tdist <- tsi[1:(len/2)]
  # tsi <- tsi[(1+(len/2)):len]
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
  if(sum(i, na.rm = T)>0){
    if(sum(i) == 1) {
      res[i,] <- toFireTS(x[i,], dts, resol = resol, thres = thres)
    } else if(sum(i) > 1) {
      res[i,] <- t(apply(x[i,], 1, toFireTS, dts, resol, thres))
    }
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

#' Mask low quality observations for Landsat 4, 5, or 7 time series.
#'
#' @param x a vector, where the first n observations represent Landsat observations,
#' the next n observations represent the associated SR_CLOUD_QA values,
#' the next n observations represent the QA_PIXEL values,
#' and the final n observations represent the QA_RADSAT values
#'
#' @return vector of Landsat observations for which low quality observations are masked
#' @export
#'
maskL47 <- function(x){
  len <- length(x)/4
  vals <- x[1:len]
  cld <- x[(len+1):(2*len)]
  px <- x[((2*len)+1):(3*len)]
  sat <- x[((3*len)+1):(4*len)]
  # SR_CLOUD_QA mask: cloud (bit 1), cloud shadow (bit 2), adjacent to cloud (bit 3), snow (bit 4), water (bit 5)
  m <- sapply(cld,function(x){ as.integer(intToBits(x))})# unpack bit words
  vals[(m[2,] == 1) | (m[3,] == 1) | (m[4,] == 1) | (m[5,] == 1) | (m[6,] == 1)] <- NA
  # QA_PIXEL mask: fill (bit 0), dilated cloud (bit 1), cloud (bit 3), cloud shadow (bit 4)
  m <- sapply(px,function(x){ as.integer(intToBits(x))})
  vals[(m[1,] == 1) | (m[2,] == 0) | (m[3,] == 1) | (m[4,] == 1) | (m[5,] == 1)| (m[6,] == 1)] <- NA
  # QA_RADSAT radionetric saturation quality (shows which sensor bands were saturated during data capture, yielding unusable data): saturation band 4 (bit 3) and band 7 (bit 6)
  m <- sapply(sat,function(x){ as.integer(intToBits(x))})
  vals[(m[5,] == 1) | (m[8,] == 1)] <- NA

  return(vals)
}

#'  Mask low quality observations for a stack of Landsat 4, 5, or 7 images.
#'
#' @param x a raster stack, where the first image is a mask (only pixels with a value equal to 1 are processed),
#' the subsequent n images represent Landsat observations,
#' the next n images represent the QA_PIXEL values,
#' and the final n images represent the QA_RADSAT values
#'
#' @return stack of images for which low quality observations are masked
#' @export
#'
maskL47Stack <- function(x)
{
  msk <- x[,1]
  x <- x[,-1]
  i <- (msk == 1)
  len <- dim(x)[2]/4
  res <- matrix(NA, length(i), len)
  if(sum(i) == 1) {
    res[i,] <- maskL47(x[i,])
  } else if(sum(i) > 1) {
    res[i,] <- t(apply(x[i,],1, maskL47))
  }
  res
}


#' Mask low quality observations for Landsat 8 time series.
#'
#' @param x a vector, where the first n observations represent Landsat observations,
#' the next n observations represent the QA_PIXEL values,
#' and the final n observations represent the QA_RADSAT values
#'
#' @return vector of Landsat observations for which low quality observations are masked
#' @export
#'
maskL8 <- function(x){
  len <- length(x)/3
  vals <- x[1:len]
  px <- x[(len+1):(2*len)]
  sat <- x[((2*len)+1):(3*len)]
  # QA_PIXEL mask: fill (bit 0), dilated cloud (bit 1), cloud (bit 3), cloud shadow (bit 4)
  m <- sapply(px,function(x){ as.integer(intToBits(x))})
  vals[(m[1,] == 1) | (m[2,] == 0) | (m[3,] == 1) | (m[4,] == 1) | (m[5,] == 1)| (m[6,] == 1)] <- NA
  # QA_RADSAT radionetric saturation quality (shows which sensor bands were saturated during data capture, yielding unusable data): saturation band 5 (bit 5) and band 7 (bit 7)
  m <- sapply(sat,function(x){ as.integer(intToBits(x))})
  vals[(m[6,] == 1) | (m[8,] == 1)] <- NA

  return(vals)
}

#'  Mask low quality observations for a stack of Landsat 8 images.
#'
#' @param x a raster stack, where the first image is a mask (only pixels with a value equal to 1 are processed),
#' the subsequent n images represent Landsat observations,
#' the next n images represent the associated SR_CLOUD_QA values,
#' the next n images represent the QA_PIXEL values,
#' and the final n images represent the QA_RADSAT values
#'
#' @return stack of images for which low quality observations are masked
#' @export
#'
maskL8Stack <- function(x)
{
  msk <- x[,1]
  x <- x[,-1]
  i <- (msk == 1)
  len <- dim(x)[2]/3
  res <- matrix(NA, length(i), len)
  if(sum(i) == 1) {
    res[i,] <- maskL8(x[i,])
  } else if(sum(i) > 1) {
    res[i,] <- t(apply(x[i,],1, maskL8))
  }
  res
}

#' Calculate Normalized Burn Ratio
#'
#' @param NIR Reflectance NIR band (band 4 for landsat 4-7 and band 5 for landsat 8)
#' @param SWIR Reflectance SWIR band (band 7 of landsat)
#'
#' @return
#' @export
#'
calcNBR <- function(NIR, SWIR){
  return ((NIR - SWIR) / (NIR + SWIR))
}
