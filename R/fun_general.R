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

#' Set and generate folder structure to store data for the FORCE processing workflow
#'
#' @param forcefolder the main folder where all data needs to be stored (full path)
#'
#' @return generates a folder structure
#' @export
#'
setFolders <- function(forcefolder){
  if (!dir.exists(forcefolder)){
    stop('directory does not exist')
  }
  tmpfolder <- file.path(forcefolder, 'temp')
  l1folder <- file.path(forcefolder, 'level1')
  l2folder <- file.path(forcefolder, 'level2')
  queuefolder <- file.path(forcefolder, 'level1')
  queuefile <- 'queue.txt'
  demfolder <- file.path(forcefolder, 'misc','dem')
  wvpfolder <- file.path(forcefolder, 'misc','wvp')
  logfolder <- file.path(forcefolder, 'log')
  paramfolder <- file.path(forcefolder, 'param')
  paramfile <- 'l2param.prm'
  lcfolder <- file.path(forcefolder, 'misc','lc')# raw land cover data
  tcfolder <- file.path(forcefolder, 'misc','tc')# raw tree cover data
  firefolder <- file.path(forcefolder, 'misc','fire')# raw fire data
  S2auxfolder <- file.path(forcefolder, 'misc', 'S2')# auxiliary S2 data (eg tile grid)

  demlogfile <- file.path(logfolder,'DEM.txt')
  wvplogfile <- file.path(logfolder,'WVP.txt')
  landsatlogfile <- file.path(logfolder, 'Landsat.txt')
  lclogfile <- file.path(logfolder, 'LC.txt')
  firelogfile <- file.path(logfolder,'fire.txt')
  tclogfile <- file.path(logfolder, 'tc.txt')
  Sskiplogfile <- file.path(logfolder, 'Sskip.txt')
  Ssuccesslogfile <- file.path(logfolder, 'Ssuccess.txt')
  Smissionlogfile <- file.path(logfolder, 'Smission.txt')
  Sotherlogfile <- file.path(logfolder, 'Sother.txt')

  if(!dir.exists(forcefolder)){dir.create(forcefolder)}
  if(!dir.exists(tmpfolder)){dir.create(tmpfolder)}
  if(!dir.exists(l1folder)){dir.create(l1folder)}
  if(!dir.exists(file.path(l1folder,'landsat'))){dir.create(file.path(l1folder,'landsat'))}
  if(!dir.exists(file.path(l1folder,'sentinel'))){dir.create(file.path(l1folder,'sentinel'))}
  if(!dir.exists(l2folder)){dir.create(l2folder)}
  if(!dir.exists(queuefolder)){dir.create(queuefolder)}
  if(!dir.exists(demfolder)){dir.create(demfolder, recursive = TRUE)}
  if(!dir.exists(wvpfolder)){dir.create(wvpfolder, recursive = TRUE)}
  if(!dir.exists(logfolder)){dir.create(logfolder)}
  if(!dir.exists(S2auxfolder)){dir.create(S2auxfolder)}

  if(!file.exists(demlogfile)){file.create(demlogfile)}# logfile for DEM
  if(!file.exists(wvplogfile)){file.create(wvplogfile)}# logfile for WVP
  if(!file.exists(landsatlogfile)){file.create(landsatlogfile)}# logfile for DEM
  if(!file.exists(lclogfile)){file.create(lclogfile)}# logfile for DEM
  if(!file.exists(firelogfile)){file.create(firelogfile)}# logfile for DEM
  if(!file.exists(tclogfile)){file.create(tclogfile)}# logfile for DEMS
  if(!file.exists(Sskiplogfile)){file.create(Sskiplogfile)}# logfile for skipped scenes
  if(!file.exists(Ssuccesslogfile)){file.create(Ssuccesslogfile)}# logfile for successful scenes
  if(!file.exists(Smissionlogfile)){file.create(Smissionlogfile)}# logfile for scenes with an unknown mission
  if(!file.exists(Sotherlogfile)){file.create(Sotherlogfile)}# logfile for scenes with an unrecoginized processing status
  if(!file.exists(file.path(queuefolder,queuefile))){file.create(file.path(queuefolder,queuefile))}# generate a queue file
  if(!dir.exists(paramfolder)){dir.create(paramfolder)}
  if(!dir.exists(lcfolder)){dir.create(lcfolder)}
  if(!dir.exists(tcfolder)){dir.create(tcfolder)}
  if(!dir.exists(firefolder)){dir.create(firefolder)}

  out <- c(tmpfolder, l1folder, l2folder, queuefolder, queuefile, demfolder, wvpfolder, logfolder, paramfolder, paramfile,
           lcfolder, tcfolder, firefolder, S2auxfolder, demlogfile, wvplogfile, landsatlogfile, lclogfile, firelogfile, tclogfile, Sskiplogfile, Ssuccesslogfile, Smissionlogfile, Sotherlogfile)
  names(out) <- c('tmpfolder', 'l1folder', 'l2folder', 'queuefolder', 'queuefile', 'demfolder', 'wvpfolder', 'logfolder', 'paramfolder', 'paramfile',
                  'lcfolder', 'tcfolder', 'firefolder', 'S2auxfolder', 'demlogfile', 'wvplogfile', 'landsatlogfile', 'lclogfile', 'firelogfile', 'tclogfile','Sskiplogfile', 'Ssuccesslogfile', 'Smissionlogfile', 'Sotherlogfile')
  return(out)

}

#' Extract the extent of each grid tile that covers an area of interest
#'
#' @param l2folder directory of the level2 data cube
#' @param ext extent of the area of interest
#'
#' @return list of extents
#' @import raster
#' @import sp
#' @export
#'
getGrid <- function(cubefolder, ext){
  system(paste0("force-tabulate-grid ", cubefolder, " ", ext[3]," ", ext[4]," ", ext[1]," ", ext[2], " shp"), intern = TRUE, ignore.stderr = TRUE)
  # load shapefile
  p <- shapefile(file.path(cubefolder, 'shp',"grid.shp"))
  # transform crs to crs of interest
  p_wgs <- spTransform(p, CRS("+proj=longlat +datum=WGS84"))
  # extent of each polygon/tile
  elist <- lapply(1:length(p_wgs), function(i) extent(p_wgs[i,]))
  names(elist) <- p_wgs$Tile_ID
  unlink(file.path(cubefolder, 'shp'), recursive = T)
  return(elist)
}

#' Convert time series to annual frequency: The toAnnualTS function converts a time series with n observations per year to an annual time series (one observation per year). The main concept is to select observations per year closest to a given day of year that have no missing value (NA). Here, the day of year for which the seasonality is maximum is being used.
#'
#' @param tsseas vector of observations (time series) representing the seasonal component of the time series to be converted
#' @param tsi vector of observations (time series) that needs to be converted to an annual time series
#' @param obspyr number of observations per year of the time series to be converted
#' @param dtmax maximum time (expressed in year, so 1/12 equals one month) between selected observation and the seasonal maximum
#'
#' @return vector of observations (time series) with annual observation frequency
#' @export
toAnnualTS <- function(tsseas, tsi, obspyr, dtmax = 1/12){
  seasi <- rowMeans(matrix(tsseas, nrow = obspyr),na.rm=T)# average seasonality
  smax <- which(seasi == max(seasi, na.rm=T))# yearly observation number with max seas
  tsmi <- matrix(tsi, nrow = obspyr)
  dst <- abs(matrix(rep(1:obspyr,times = (length(tsi)/obspyr)), nrow = obspyr)-smax[1])# distance of observations to seasonal max
  dst[is.na(tsmi)] <- NA #set distance of NA observations equal to NA, these can not be selected
  dst[dst > floor(obspyr*dtmax)] <- NA
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
    } else if(fun == 'mean'){
      mz <- as.ts(aggregate(z, as.yearmon, mmean)) ## mean
    }}else if(resol == 'daily'){
    mz <- bfastts(tsi, dts, type = "irregular")
  }else if(resol == 'quarterly'){
    z <- zoo(c(NA,tsi), c(as.Date(paste0(format(min(dts),'%Y'),'-01-01')),dts)) ## create a zoo (time) series
    if(fun == 'max'){
      mz <- as.ts(aggregate(z, as.yearqtr, mmax)) ## max
    } else if(fun == 'mean'){
      mz <- as.ts(aggregate(z, as.yearqtr, mmean)) ## mean
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
   msk <- x[1]# mask value for the time series
   msk[is.na(msk)] <- 0
   x <- x[-1]# time series
   missval <- sum(is.na(x)) / length(x) # fraction of missing values
   len <- length(toRegularTS(dts, dts, fun=fun, resol = resol))# length of output time series
   if((missval < 1) & (msk == 1)){
     res <- toRegularTS(x, dts, fun=fun, resol = resol)
     if(length(res) != len){res <- rep(NA,len)}
   } else{
     res <- rep(NA,len)
   }
  return(res)
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

#' Helper function for the toRegularTSStack function
#'
#' @param x vector of observations
#'
#' @return the mean value of the vector
#' @export
mmean <- function(x) {
  if(length(which.max(x)) == 0) {
    out <- NA
  } else {
    out <- as.numeric(mean(x, na.rm = T))
  }
  return(out)
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
#' @import terra
#'
mskQA <- function(im, qaim, valid = 0, cloud = 0, shadow = 0, snow = 0, water = 0, aero = 0, subzero = 0, sat = 0, sunZen = 0, illum = 0, slope = 0, wvp = 0){
  # QA values of pixels that shouldn't be masked
  toKeep <- getQAvals(valid, cloud, shadow, snow, water, aero, subzero, sat, sunZen, illum, slope, wvp)
  # make raster with 0 values for pixels that should be masked and 1 for pixels that should be kept
  msk <- qaim
  values(msk) <- rep(0,size(qaim))
  # msk[] <- 0
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
  bcloud <- paste(as.numeric(rev(intToBits(cloud)[1:2])), collapse = '')
  baero <- paste(as.numeric(rev(intToBits(aero)[1:2])), collapse = '')
  billum <- paste(as.numeric(rev(intToBits(illum)[1:2])), collapse = '')
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
#' @param br terra raster stack
#' @param starttime start date of the desired time span (Date object)
#' @param endtime end date of the desired time span (Date object)
#' @param tempRes temporal resolution of the raster stack
#' @param dtsbr vector of dates (Date object) associated with the raster stack
#'
#' @return terra raster stack with adjusted time span
#' @import lubridate
#' @export
#'
setPeriod <- function(br, startyr, endyr,tempRes, dtsbr){
  tres <- c('1 month', '3 months', '1 day', '1 year')
  names(tres) <- c('monthly', 'quarterly', 'daily', 'yearly')
  # make sure that image stack covers time period of interest
  rstNA <- br[[1]]
  values(rstNA) <- rep(NA, ncell(rstNA)) # empty image
  dtstot <- seq(startyr, endyr, by = tres[tempRes])

  out <- br
  out <- out[[which((dtsbr>=min(dtstot)) & (dtsbr<=max(dtstot)))]]# remove observations that fall outside the study period
  npre <- sum(dtstot<min(dtsbr))# missing dates at the start of the study period
  npost <- sum(dtstot>max(dtsbr))# missing dates at the end of the study period
  if(npre>0){for(i in 1:npre){out <- c(rstNA, out)}}
  if(npost>0){for(i in 1:npost){out <- c(out,rstNA)}}

  # if(npre>0){for(i in 1:npre){out <- c(rstNA, out)}}
  # else{# add observations at the beginning if needed
    # remove observations before the start of the study period if needed
  #   ind <- which(dtsbr>=min(dtstot))
  #   dtsbr <- dtsbr[ind]
  #   out <- out[[ind]]
  # }
  # if(npost>0){for(i in 1:npost){out <- addLayer(out,rstNA)}}
  # else{# add observations at the end if needed
  #   # remove observations after the end of the study period if needed
  #   ind <- which(dtsbr<=max(dtstot))
  #   out <- out[[ind]]
  # }
  names(out) <- dtstot
  return(out)
}

#' Set terra raster values outside area of interest to 0
#'
#' @param fmask a terra raster layer
#' @param ext extent of area of interest (vector with xmin, xmax, ymin, ymax), expressed as longitude and latitude
#'
#' @return terra raster layer
#' @export
#' @import terra
#'
maskAOI <- function(fmask, ext){
  # first generate a raster with value 0 outside the area of interest and 1 inside the area of interest
  z <- cbind(object=1, part=1, rbind(c(ext[1],ext[4]),
                                     c(ext[2],ext[4]),
                                     c(ext[2],ext[3]),
                                     c(ext[1],ext[3])), hole=0)
  colnames(z)[3:4] <- c('x','y')
  z <- data.frame(z)
  aoi <- terra::vect(z, type="polygons", crs = CRS("+init=epsg:4326"))
  aoi <- terra::project(aoi,showP4(crs(fmask)))
  aoi <- rasterize(aoi, fmask, field=1, background=0)
  fmask[aoi == 0] <- 0# set all values of the mask layer outside the area of interest to 0
  return(fmask)
}

#' Prepare NBR stack
#' - crop temporal extent to study period
#' - mask low quality observations using QAI
#' - calculate NBR
#' - aggregate Sentinel-2 to Landsat resolution
#' - aggregate observations to temporal resolution of interest
#'
#' @param img_list list of terra rasters with spectral information
#' @param qai_list list of terra rasters with quality assurance information
#' @param fmask terra raster with a processing mask (0 = not processed, 1 = to be processed)
#' @param lsens the sensor types associated with the images ()
#' @param ldts dates (Date object) associated with the images
#' @param tempFun function used to temporally aggregate time series (mean or max)
#' @param tempRes temporal resolution of interest ('daily', 'quarterly', 'monthly' or 'yearly')
#' @param starttime start date (Date object) of study period
#' @param endtime end date (Date object) of study period
#'
#' @return
#' @export
#' @import lubridate
#' @import terra
#'
#' @examples
prepareNBRstack <- function(img_list, qai_list, fmask, lsens, ldts, tempFun, tempRes, startdt, enddt){
  img_list <- img_list[(ldts >= startdt) & (ldts <= enddt)]
  qai_list <- qai_list[(ldts >= startdt) & (ldts <= enddt)]
  lsens <- lsens[(ldts >= startdt) & (ldts <= enddt)]
  ldts <- ldts[(ldts >= startdt) & (ldts <= enddt)]

  for(i in 1:length(img_list)){# iterate over the BOA layers (each layer is associated with a particular date)
    # open image
    im <- img_list[[i]]#rast(file.path(tilefolder,img[i]))# load the image with spectral data
    qaim <- qai_list[[i]]#rast(file.path(tilefolder,gsub("BOA", "QAI", img[i])))# load the Quality Assurance layer
    # calculate NBR
    nbri <- switch(1+ ((lsens[i] == 'SEN2A') | (lsens[i] == 'SEN2B')), calcNBR(im[[4]], im[[6]]), calcNBR(im[[8]], im[[10]]))
    # mask low quality data
    nbrf <- mskQA(nbri, qaim, slope = c(0,1), aero = c(0,1))

    # resample S2 to the spatial resolution of Landsat
    if ((lsens[i] == 'SEN2A') | (lsens[i] == 'SEN2B')){
      nbrf <- terra::aggregate(nbrf, fact = 3, fun=mmean)
    }
    # mask data using forest fire mask
    nbrf[fmask!=1] <- NA
    # add image to stack
    if(i == 1){
      st <- nbrf
    }else{
      st <- c(st, nbrf)
    }
    rm(nbri, nbrf, im, qaim)
  }
  # assign dates to stack
  names(st) <- ldts

  # generate a regular image stack
  if(tempRes == 'monthly'){
    ldts <- rollback(ldts, roll_to_first = TRUE, preserve_hms = TRUE)
  }
  tsbr <- app(c(fmask,st), fun = function(x){toRegularTSStack(x, ldts, fun = tempFun, resol = tempRes)})# tsbr <- calc(stack(fmask,st), function(x){toRegularTSStack(x, ldts, fun = tempFun, resol = tempRes)})#
  tres <- c('1 month', '3 months', '1 day', '1 year')
  names(tres) <- c('monthly', 'quarterly', 'daily', 'yearly')

  dtsbr <- switch(((tempRes == 'quarterly')+1),seq(min(ldts), max(ldts), by = tres[tempRes]),
    seq(as.Date(paste0(format(min(ldts),'%Y'),'-01-01')), max(ldts), by = tres[tempRes]))

  # make sure that image stack covers time period of interest
  if(tempRes == 'monthly' || tempRes == 'quarterly'){
    if (tempRes == 'monthly'){
      startdt <- rollback(startdt, roll_to_first = TRUE, preserve_hms = TRUE)
    }else if (tempRes == 'quarterly'){
      startdt <- rollback_quart(startdt)
      }
    endtdt <- rollback(enddt, roll_to_first = TRUE, preserve_hms = TRUE)
  }
  tsNBR <- setPeriod(tsbr, startdt, enddt,tempRes, dtsbr)

}

#' Rollback dates to first day of quarterly series
#'
#' @param dts dates to be converted (date object)
#'
#' @return date object
#' @export
#' @import plyr
#'
rollback_quart <- function(dts){
  mnths <- revalue(format(dts, '%m'), c("01" = "01", "02" = "01", "03" = "01",
                                        "04" = "04", "05" = "04", "06" = "04",
                                        "06" = "06", "08" = "07", "09" = "07",
                                        "10" = "10", "11" = "10", "12" = "10"))
  out <- as.Date(paste0(format(dts,'%Y'),'-',mnths,'-01'))
  return(out)
}
