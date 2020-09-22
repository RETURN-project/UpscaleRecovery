context("Data preparation")

test_that("Calculate Normalized Burn Ratio", {

  nir <- 0.2
  swir <- 0.8
  nbr <- (nir - swir)/(nir+swir)
  expect_equal(calcNBR(nir,swir), nbr)
})

test_that("Calculate QA", {
  expect_equal(getQA(illum = 1, slope = 1, wvp = 1), 28672)
})

test_that("convert binary to integer", {
  expect_equal(BinToDec(c(1,0)), 2)
  expect_equal(BinToDec(c(1,0,1)), 5)
})

test_that("change time span of raster stack", {
  r1 <- raster(vals = 1:12, nrows = 3, ncols = 4)
  r2 <- raster(vals = 13:24, nrows = 3, ncols = 4)
  r3 <- raster(vals = 25:36, nrows = 3, ncols = 4)
  r4 <- raster(vals = 37:48, nrows = 3, ncols = 4)
  st <- stack(r1,r2,r3,r4)
  tempRes <- 'monthly'
  dts <- seq(as.Date('2000-01-01'), as.Date('2000-04-01'), by = '1 month')

  # case 1 - append start and end
  starttime <- c(1999,8,1)
  endtime <- c(2001,1,1)
  st1 <- setPeriod(st, starttime, endtime,tempRes, dts)
  exp <- c(NA, NA, NA, NA, NA,  1, 13, 25, 37, NA, NA, NA, NA, NA, NA, NA, NA, NA)
  expect_equal(as.numeric(st1[1]), exp)

  # case 2 - shorten start and extend end
  starttime <- c(2000,2,1)
  endtime <- c(2000,3,1)
  st2 <- setPeriod(st, starttime, endtime,tempRes, dts)
  exp <- c(13, 25)
  expect_equal(as.numeric(st2[1]), exp)

})

test_that("Get combination of QA values", {
  vals <- getQAvals(valid = c(0,1), aero = c(0,3))
  exp <- c(getQA(), getQA(valid = 1), getQA(aero = 3), getQA(valid = 1, aero = 3))
  expect_equal(vals, exp)
})

test_that("Mask using QA values",{
  im <- raster(vals = 1:12, nrows = 3, ncols = 4)
  qaim <- raster(vals = c(0, 1, 4, 9, 11, 12, 128, 132, 966, 8192, 8324, 28795), nrows = 3, ncols = 4)
  mskd <- mskQA(im, qaim, valid = 0, cloud = c(0,1), shadow = 0, snow = 0, water = 0, aero = c(0,1), subzero = 0, sat = 0, sunZen = 0, illum = 0, slope = c(0,1), wvp = 0)
  exp <- c(1,NA,3,NA,NA,NA,7,8,NA,10,11,NA)
  expect_equal(mskd[], exp)
})

test_that("dense to annual", {
  tsseas <- rep(c(1,4,5,2),2)
  tsi <- c(2,5,7,4,5,6,7,3,4,5,6,7,2,1,3,4,5,3,6,9)
  obspyr <- 4

  tsa <- toAnnualTS(tsseas, tsi, obspyr)
  annual <- tsi[seq(3,20,by=4)]

  diff <- sum(tsa - annual)

  expect_equal(diff, 0, tolerance = 1e-4)
})

test_that("dense to annual with missing values", {
  tsseas <- rep(c(1,4,5,2),2)
  tsi <- c(2,5,NA,4,5,6,7,3,4,5,6,7,2,1,3,4,5,3,6,9)
  obspyr <- 4

  tsa <- toAnnualTS(tsseas, tsi, obspyr, dtmax = 1/12)
  annual <- tsi[seq(3,20,by=4)]

  expect_equal(tsa, annual, tolerance = 1e-4)
})

test_that("dense to annual with varying intra-annual selection period", {
  tsseas <- rep(c(1,4,5,2),2)
  tsi <- c(2,5,NA,4,5,6,7,3,4,5,6,7,2,1,3,4,5,NA,NA,9)
  obspyr <- 4

  tsa1 <- toAnnualTS(tsseas, tsi, obspyr, dtmax = 1/12)
  annual1 <- tsi[seq(3,20,by=4)]
  tsa2 <- toAnnualTS(tsseas, tsi, obspyr, dtmax = 1/2)
  annual2 <- c(5, 7, 6, 3, 9)

  expect_equal(tsa1, annual1, tolerance = 1e-4)
  expect_equal(tsa2, annual2, tolerance = 1e-4)
})

test_that("Temporal aggregation", {
  library(raster)
  library(lubridate)
  # mask
  m <- raster(ncol=3, nrow=3, vals=c(1,1,0,1,1,1,1,1,1))
  # rasters of time series observations
  r1<- raster(ncol=3, nrow=3, vals=rep(1, 9))
  r2<- raster(ncol=3, nrow=3, vals=rep(2, 9))
  r3<- raster(ncol=3, nrow=3, vals=rep(1, 9))
  r4<- raster(ncol=3, nrow=3, vals=rep(20, 9))
  r5<- raster(ncol=3, nrow=3, vals=rep(1, 9))
  r6<- raster(ncol=3, nrow=3, vals=rep(30, 9))
  r7<- raster(ncol=3, nrow=3, vals=c(-12,1,1,1,1,1,1,1,1))
  r8<- raster(ncol=3, nrow=3, vals=c(-12,1,1,1,1,1,1,1,1)+10)
  r9<- raster(ncol=3, nrow=3, vals=c(-11,1,1,-5,-3,-3,-3,-3,-3))
  r10<- raster(ncol=3, nrow=3, vals=c(-11,1,1,-5,-3,-3,-3,-3,-3)-10)
  r11<- raster(ncol=3, nrow=3, vals=c(-10,1,1,-4,-2,-2,-2,-2,-2))
  r12<- raster(ncol=3, nrow=3, vals=c(-10,1,1,-4,-2,-2,-2,-2,-2)-20)
  r13<- raster(ncol=3, nrow=3, vals=c(-9,1,1,-3,-1,-1,-1,-1,-1))
  r14<- raster(ncol=3, nrow=3, vals=c(-9,1,1,-3,-1,-1,-1,-1,-1)-30)
  r15<- raster(ncol=3, nrow=3, vals=c(-8,1,1,-2,0,0,0,0,0))
  r16<- raster(ncol=3, nrow=3, vals=c(-8,1,1,-2,0,0,0,0,0)-40)
  r17<- raster(ncol=3, nrow=3, vals=c(-7,1,1,-1,1,1,1,1,1))
  r18<- raster(ncol=3, nrow=3, vals=c(-7,1,1,-1,1,1,1,1,1)-50)
  r19<- raster(ncol=3, nrow=3, vals=c(-6,1,1,0,1,1,1,1,1))
  r20<- raster(ncol=3, nrow=3, vals=c(-6,1,1,0,1,1,1,1,1)-60)
  #create stack
  st <- stack(m,r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,
              r11,r12,r13,r14,r15,r16,r17,r18,r19,r20)
  # dates associated with observations
  dts <- as.Date(c('2001-01-02','2001-01-03','2001-02-02','2001-02-04','2001-03-02','2001-03-04','2001-04-02','2001-04-05','2001-05-02','2001-05-12',
                   '2001-06-02','2001-06-03','2001-07-02','2001-07-22','2001-08-02','2001-08-12','2001-09-02','2001-09-22','2001-12-02','2001-10-02'))
  # names(st) <- dts
  dts <- as.Date(dts, format = "X%Y.%m.%d") ## needed as input in the helper function of get_m_agg

  # create time series of monthly max
  brmomax <- calc(st, function(x){toRegularTSStack(x, dts, fun = 'max', resol = 'monthly')})
  names(brmomax) <- as.Date(toRegularTS(dts, dts, fun = 'max', resol = 'monthly'))
  # create time series of monthly mean values
  brmomean <- calc(st, function(x){toRegularTSStack(x, dts, fun = 'mean', resol = 'monthly')})
  names(brmomean) <- as.Date(toRegularTS(dts, dts, fun = 'mean', resol = 'monthly'))
  # create daily time series
  brday <- calc(st, function(x){toRegularTSStack(x, dts, fun = 'max', resol = 'daily')})
  names(brday) <- date_decimal(as.numeric(time(bfastts(rep(1,length(dts)), dts, type = "irregular"))))
  # create quarterly time series
  brquart <- calc(st, function(x){toRegularTSStack(x, dts, fun = 'max', resol = 'quart')})
  names(brquart) <- as.Date(toRegularTS(dts, dts, fun = 'mean', resol = 'quart'))


  # convert to matrix'
  mmomax <- raster::as.matrix(brmomax)
  mmomean <- raster::as.matrix(brmomean)
  mday <- raster::as.matrix(brday)
  mquart <-  raster::as.matrix(brquart)

  # case 1 - monthly max
  expect_equal(as.numeric(mmomax[1,]), c(2,20,30,-2,-11,-10,-9,-8,-7,-66,NA,-6), tolerance = 1e-4)
  # case 2 - monthly mean
  expect_equal(as.numeric(mmomean[1,]), c(1.5,10.5,15.5,-7.0,-16.0,-20.0,-24.0,-28.0,-32.0,-66.0,NA,-6.0), tolerance = 1e-4)
  # case 3 - daily
  expect_equal(as.numeric(mday[1,c(1,2,32,34,60,62,91,94,121,131,152,153,182,202,213,223,244,264,274,335)]),
               c(1,2,1,20,1,30,-12,-2,-11,-21,-10,-30,-9,-39,-8,-48,-7,-57,-66,-6), tolerance = 1e-4)
  # case 4 - quarterly
  expect_equal(as.numeric(mquart[1,]),c(30,-2, -7,-6))
})


test_that("Prepare fire time series", {
  library(raster)
  #mask
  m <- raster(ncol=3, nrow=3, vals=c(1,1,0,1,1,1,1,1,1))
  # confidence
  cl1<- raster(ncol=3, nrow=3, vals=rep(0, 9))
  cl2<- raster(ncol=3, nrow=3, vals=c(96,0,0,0,0,0,15,0,0))
  cl3<- raster(ncol=3, nrow=3, vals=c(0,97,0,0,0,0,0,0,0))
  cl4<- raster(ncol=3, nrow=3, vals=c(0,0,97,0,0,0,0,0,0))
  cl5<- raster(ncol=3, nrow=3, vals=c(0,0,0,99,0,0,0,0,0))
  cl6<- raster(ncol=3, nrow=3, vals=c(0,0,0,0,98,0,0,0,0))
  cl7<- raster(ncol=3, nrow=3, vals=c(0,0,0,0,0,96,0,0,0))
  cl8<- raster(ncol=3, nrow=3, vals=c(99,0,0,0,0,0,99,0,0))
  cl9<- raster(ncol=3, nrow=3, vals=c(0,0,0,0,0,0,0,99,97))
  cl10<- raster(ncol=3, nrow=3, vals=c(0,0,0,0,0,0,0,0,0))
  cl11<- raster(ncol=3, nrow=3, vals=c(0,0,0,0,0,0,0,0,0))
  cl12<- raster(ncol=3, nrow=3, vals=c(0,0,0,0,0,0,0,0,0))

  # day of observation
  jd1<- raster(ncol=3, nrow=3, vals=c(0,0,0,0,0,0,0,0,0))
  jd2<- raster(ncol=3, nrow=3, vals=c(33,0,0,0,0,0,0,0,0))# 2 feb
  jd3<- raster(ncol=3, nrow=3, vals=c(0,62,0,0,0,0,0,0,0))# 3 mar
  jd4<- raster(ncol=3, nrow=3, vals=c(0,0,94,0,0,0,0,0,0))# 4 apr
  jd5<- raster(ncol=3, nrow=3, vals=c(0,0,0,125,0,0,0,0,0))#5 may
  jd6<- raster(ncol=3, nrow=3, vals=c(0,0,0,0,157,0,0,0,0))#6jun
  jd7<- raster(ncol=3, nrow=3, vals=c(0,0,0,0,0,188,0,0,0))#7jul
  jd8<- raster(ncol=3, nrow=3, vals=c(220,0,0,0,0,0,232,0,0))#8aug & 20 aug
  jd9<- raster(ncol=3, nrow=3, vals=c(0,0,0,0,0,0,0,252,254))#9sep 11sep
  jd10<- raster(ncol=3, nrow=3, vals=c(0,0,0,0,0,0,0,0,0))
  jd11<- raster(ncol=3, nrow=3, vals=c(0,0,0,0,0,0,0,0,0))
  jd12<- raster(ncol=3, nrow=3, vals=c(0,0,0,0,0,0,0,0,0))

  #create stack
  st <- stack(m,cl1,cl2,cl3,cl4,cl5,cl6,cl7,cl8,cl9,cl10,cl11,cl12,
              jd1,jd2,jd3,jd4,jd5,jd6,jd7,jd8,jd9,jd10,jd11,jd12)
  dts <- seq(as.Date(paste0(2001,'-01-01')), as.Date(paste0(2001,'-12-31')), by = "1 month")

  firemo <- calc(st, function(x){createFireStack(x, dts, resol = 'monthly', thres = 95)})
  fireday <- calc(st, function(x){createFireStack(x, dts, resol = 'daily', thres = 95)})

  mfmo <- raster::as.matrix(firemo)
  mfday <- raster::as.matrix(fireday)

  d1 <- rep(0,365)
  d1[c(33,220)] <- 1
  d7 <- rep(0,365)
  d7[232] <- 1

  # case 1 - monthly
  expect_equal(as.numeric(mfmo[1,]), c(0,1,0,0,0,0,0,1,0,0,0,0), tolerance = 1e-4)
  expect_equal(as.numeric(mfmo[7,]), c(0,0,0,0,0,0,0,1,0,0,0,0), tolerance = 1e-4)
  expect_equal(sum(is.na(mfmo[3,])), 12, tolerance = 1e-4)
  # case 2 - daily
  expect_equal(as.numeric(mfday[1,]), d1, tolerance = 1e-4)
  expect_equal(as.numeric(mfday[7,]), d7, tolerance = 1e-4)

})


