context("Calculate recovery")

test_that("Frazier - annual - too short time series", {

  tsio <- c(rep(0,1), seq(-1, 0), rep(0,1))
  tdist <- 2
  obspyr <- 1
  shortDenseTS <- FALSE
  nPre <- 1
  nDist <- 1
  nPostMin <- 1
  nPostMax <- 1

  metrics <- calcFrazier(tsio, tdist, obspyr, shortDenseTS, nPre, nDist, nPostMin, nPostMax)

  expect_equal(metrics$RRI, NA)
  expect_equal(metrics$R80P, NA)
  expect_equal(metrics$YrYr, NA)
})

test_that("Frazier - annual", {

  tsio <- c(rep(1,2), seq(-5, -1), rep(-2,1))
  tdist <- 3
  obspyr <- 1
  shortDenseTS <- FALSE
  nPre <- 1
  nDist <- 1
  nPostMin <- 1
  nPostMax <- 1

  metrics <- calcFrazier(tsio, tdist, obspyr, shortDenseTS, nPre, nDist, nPostMin, nPostMax)
  pre <- 1
  dnbr <- 6
  ari <- 4

  rri <- ari/dnbr
  r80p <- -1/(0.8*pre)
  yryr <- (-2 + 5)/5

  expect_equal(metrics$RRI, rri, tolerance = 1e-4)
  expect_equal(metrics$R80P, r80p, tolerance = 1e-4)
  expect_equal(metrics$YrYr, yryr, tolerance = 1e-4)
})

test_that("Frazier - dense", {

  tsio <- c(rep(1,24), seq(-5, -1, length.out=60), rep(-2,12))
  tdist <- 25
  obspyr <- 12
  shortDenseTS <- TRUE
  nPre <- 2
  nDist <- 12
  nPostMin <- 4
  nPostMax <- 5

  metrics <- calcFrazier(tsio, tdist, obspyr, shortDenseTS, nPre, nDist, nPostMin, nPostMax)
  pre <- 1
  dist <- mean(tsio[25:36])
  post <- mean(tsio[73:85])
  dnbr <- pre-dist
  ari <- post-dist

  rri <- ari/dnbr
  r80p <- post/(0.8*pre)
  yryr <- (post - dist)/4.5


  expect_equal(metrics$RRI, rri, tolerance = 1e-4)
  expect_equal(metrics$R80P, r80p, tolerance = 1e-4)
  expect_equal(metrics$YrYr, yryr, tolerance = 1e-4)
})

test_that("Frazier - segmented", {

  tsio <- c(rep(1,24), seq(-5, -1, length.out=60), rep(-2,12))
  tdist <- rep(0,96)
  tdist[25] <- 1
  obspyr <- 12
  shortDenseTS <- TRUE
  nPre <- 2
  nDist <- 12
  nPostMin <- 4
  nPostMax <- 5
  h <- 0.1
  timeThres <- 2
  slpThres <- 2

  metrics <- calcSegRec(tsio, tdist, maxBreak=T, obspyr, h, shortDenseTS, nPre, nDist, nPostMin, nPostMax, timeThres, slpThres)
  pre <- 1
  dist <- mean(tsio[25:36])
  post <- mean(tsio[73:85])
  dnbr <- pre-dist
  ari <- post-dist

  rri <- ari/dnbr
  r80p <- post/(0.8*pre)
  yryr <- (post - dist)/4.5
  sl <- 4/60

  expect_equal(metrics$RRI, rri, tolerance = 1e-4)
  expect_equal(metrics$R80P, r80p, tolerance = 1e-4)
  expect_equal(metrics$YrYr, yryr, tolerance = 1e-4)
  expect_equal(metrics$Sl, sl, tolerance = 1e-2)
})

test_that("Calc recovery indicators from stack using yearly, raw observations", {
  #mask
  library('raster')
  m1 <- raster(ncol=3, nrow=3, vals = c(1,0,0,1,1,1,1,1,1))
  # fire: 1 equals fire
  f1 <- raster(ncol=3, nrow=3, vals = c(0,0,0,0,0,0,0,0,0))
  f2 <- raster(ncol=3, nrow=3, vals = c(0,0,0,0,0,0,0,0,0))
  f3 <- raster(ncol=3, nrow=3, vals = c(0,0,0,0,0,0,0,0,0))
  f4 <- raster(ncol=3, nrow=3, vals = c(1,0,0,0,0,0,0,0,0))
  f5 <- raster(ncol=3, nrow=3, vals = c(0,0,0,1,1,1,1,1,1))
  f6 <- raster(ncol=3, nrow=3, vals = c(0,0,0,0,0,0,0,0,0))
  f7 <- raster(ncol=3, nrow=3, vals = c(1,0,0,0,0,0,0,0,0))
  f8 <- raster(ncol=3, nrow=3, vals = c(0,0,0,0,0,0,0,0,0))
  f9 <- raster(ncol=3, nrow=3, vals = c(0,0,0,0,0,0,0,0,0))
  f10 <- raster(ncol=3, nrow=3, vals = c(0,0,0,0,0,0,0,0,0))
  # raster observations
  r1<- raster(ncol=3, nrow=3, vals=rep(1, 9))
  r2<- raster(ncol=3, nrow=3, vals=rep(1, 9))
  r3<- raster(ncol=3, nrow=3, vals=rep(1, 9))
  r4<- raster(ncol=3, nrow=3, vals=c(-12,1,1,1,1,1,1,1,1))
  r5<- raster(ncol=3, nrow=3, vals=c(-11,1,1,-5,-3,-3,-3,-3,-3))
  r6<- raster(ncol=3, nrow=3, vals=c(-10,1,1,-4,-2,-2,-2,-2,-2))
  r7<- raster(ncol=3, nrow=3, vals=c(-9,1,1,-3,-1,-1,-1,-1,-1))
  r8<- raster(ncol=3, nrow=3, vals=c(-8,1,1,-2,0,0,0,0,0))
  r9<- raster(ncol=3, nrow=3, vals=c(-7,1,1,-1,1,1,1,1,1))
  r10<- raster(ncol=3, nrow=3, vals=c(-6,1,1,0,1,1,1,1,1))
  #create stack
  st <- stack(m1,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,r1,r2,r3,r4,r5,r6,r7,r8,r9,r10)

  # as.matrix(st)

  # Calculate stability indicators (RRI, R80P, YrYr, Sl)
  out <- calc(st, function(x){calcRecoveryStack(x, maxBreak=T, obspyr=1, inp = 'raw', shortDenseTS = FALSE,
                                                nPre = 2, nDist = 12, nPostMin = 4, nPostMax = 6, h = 0.15, timeThres = 2, slpThres = 2)})
  names(out) <- c('RRI', 'R80p', 'YrYr', 'Slope', 'missingVal', 'loglik', 'AIC')
  mout <- raster::as.matrix(out)
  # observations that were masked
  msked <- sum(is.na(mout[2:3,]))

  # case 1
  pre <- 1 # value 2 years prior
  dnbr <- 13 # drop due to disturbance
  ari <- -7+12 # max 4-5 years post disturbance - value at disturbance
  # case 1 - rri
  c1rri <- ari/dnbr
  # case 1 - r80p
  c1r80p <- -7/(0.8*pre) # max 4-5 years post disturbance / 0.8 * pre
  # case 1 - YrYr
  c1yryr <- (-7+12)/5 # 5 years post - disturbance value / 5

  # case 2
  pre <- 1 # value 2 years prior
  dnbr <- 6 # drop due to disturbance
  ari <- 5 # max 4-5 years post disturbance - value at disturbance
  # case 2 - rri
  c2rri <- ari/dnbr
  # case 2 - r80p
  c2r80p <- 0/(0.8*pre) # max 4-5 years post disturbance / 0.8 * pre
  # case 2 - YrYr
  c2yryr <- (0+5)/5 # 5 years post - disturbance value / 5

  # case 3
  pre <- 1 # value 2 years prior
  dnbr <- 4 # drop due to disturbance
  ari <- 4 # max 4-5 years post disturbance - value at disturbance
  # case 2 - rri
  c3rri <- ari/dnbr
  # case 2 - r80p
  c3r80p <- 1/(0.8*pre) # max 4-5 years post disturbance / 0.8 * pre
  # case 2 - YrYr
  c3yryr <- (1+3)/5 # 5 years post - disturbance value / 5

  # masked time series should have NA value for the recovery indicators
  expect_equal(msked, 14, tolerance = 1e-4)
  # case 1 - multiple disturbance dates
  expect_equal(as.numeric(mout[1,1]), c1rri, tolerance = 1e-4)
  expect_equal(as.numeric(mout[1,2]), c1r80p, tolerance = 1e-4)
  expect_equal(as.numeric(mout[1,3]), c1yryr, tolerance = 1e-4)
  # case 2
  expect_equal(as.numeric(mout[4,1]), c2rri, tolerance = 1e-4)
  expect_equal(as.numeric(mout[4,2]), c2r80p, tolerance = 1e-4)
  expect_equal(as.numeric(mout[4,3]), c2yryr, tolerance = 1e-4)
  # case 3
  expect_equal(as.numeric(mout[5,1]), c3rri, tolerance = 1e-4)
  expect_equal(as.numeric(mout[5,2]), c3r80p, tolerance = 1e-4)
  expect_equal(as.numeric(mout[5,3]), c3yryr, tolerance = 1e-4)

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
