context("Calculate recovery")

test_that("Frazier - annual - too short time series", {

  tsio <- c(rep(0,1), seq(-1, 0), rep(0,1))
  tdist <- 2
  obspyr <- 1
  shortDenseTS <- FALSE
  nPre <- 2
  nDist <- 1
  nPostMin <- 4
  nPostMax <- 5

  metrics <- calcFrazier(tsio, tdist, obspyr, shortDenseTS, nPre, nDist, nPostMin, nPostMax)

  expect_equal(metrics$RRI, NA)
  expect_equal(metrics$R80P, NA)
  expect_equal(metrics$YrYr, NA)
})

test_that("Frazier - annual", {

  tsio <- c(rep(1,2), seq(-5, -1), rep(-2,1),1)
  tdist <- 3
  obspyr <- 1
  shortDenseTS <- FALSE
  nPre <- 2
  nDist <- 0
  nPostMin <- 4
  nPostMax <- 5

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
  nDist <- 1
  nPostMin <- 4
  nPostMax <- 5

  metrics <- calcFrazier(tsio, tdist, obspyr, shortDenseTS, nPre, nDist, nPostMin, nPostMax)
  pre <- 1
  dist <- mean(tsio[25:36])
  post <- max(tsio[73:84])
  dnbr <- pre-dist
  ari <- post-dist

  rri <- ari/dnbr
  r80p <- post/(0.8*pre)
  yryr <- (mean(tsio[73:84]) - dist)/(4*12)


  expect_equal(metrics$RRI, rri, tolerance = 1e-4)
  expect_equal(metrics$R80P, r80p, tolerance = 1e-4)
  expect_equal(metrics$YrYr, yryr, tolerance = 1e-4)
})

test_that("Frazier - segmented", {

  tsio <- c(rep(1,24), seq(-5, -1, length.out=60), rep(-2,12))
  tdist <- 25
  tdist[25] <- 1
  obspyr <- 12
  shortDenseTS <- TRUE
  nPre <- 2
  nDist <- 1
  nPostMin <- 4
  nPostMax <- 5
  h <- 0.1
  timeThres <- 2

  metrics <- calcSegRec(tsio, tdist, maxBreak=T, obspyr, h, shortDenseTS, nPre, nDist, nPostMin, nPostMax, timeThres)
  pre <- 1
  dist <- mean(tsio[25:36])
  post <- mean(tsio[73:84])
  dnbr <- pre-dist
  ari <- post-dist

  rri <- ari/dnbr
  r80p <- post/(0.8*pre)
  yryr <- (mean(tsio[73:84]) - dist)/(4*12)

  expect_equal(metrics$RRI, rri, tolerance = 1e-4)
  expect_equal(metrics$R80P, r80p, tolerance = 1e-4)
  expect_equal(metrics$YrYr, yryr, tolerance = 1e-4)
})

test_that("Frazier - segmented annual - long", {

  tsio <- c(rep(1,8), seq(-5, 0, by = 0.5), rep(0,8))
  tdist <- 9
  obspyr <- 1
  shortDenseTS <- FALSE
  nPre <- 2
  nDist <- 0
  nPostMin <- 4
  nPostMax <- 5
  h <- 0.2
  seas <- F

  metrics <- calcSegRec(tsio, tdist, maxBreak = T, obspyr, h, shortDenseTS, nPre, nDist, nPostMin, nPostMax, timeThres = 2,seas = F)
  pre <- 1
  dnbr <- 6
  ari <- 2.5

  rrim <- ari/dnbr
  r80pm <- tsio[14]/(0.8*pre)
  yryrm <- (tsio[14] - tsio[9])/5

  expect_equal(metrics$RRI, rrim, tolerance = 1e-4)
  expect_equal(metrics$R80P, r80pm, tolerance = 1e-4)
  expect_equal(metrics$YrYr, yryrm, tolerance = 1e-4)
})

test_that("Frazier - segmented annual - testing the preconditions", {

  tsio <- c(rep(1,8), seq(-5, 0, by = 0.5), rep(0, 8))
  tdist <- 9
  obspyr <- 1
  shortDenseTS <- FALSE
  nPre <- 2
  nDist <- 0
  nPostMin <- 4
  nPostMax <- 5
  h <- 0.2
  seas <- F
  # max time span between break and disturbance
  metrics <- calcSegRec(tsio, 2, maxBreak = F, obspyr, h, shortDenseTS, nPre, nDist, nPostMin, nPostMax, timeThres = 2,seas = F)
  expect_equal(metrics$RRI, NA, tolerance = 1e-4)
  expect_equal(metrics$R80P, NA, tolerance = 1e-4)
  expect_equal(metrics$YrYr, NA, tolerance = 1e-4)
  # pos break
  tsio <- c(rep(1,8), seq(13, 0, by = -0.5), rep(0, 8))
  metrics <- calcSegRec(tsio, 2, maxBreak = F, obspyr, h, shortDenseTS, nPre, nDist, nPostMin, nPostMax, timeThres = 2,seas = F)
  expect_equal(metrics$RRI, NA, tolerance = 1e-4)
  expect_equal(metrics$R80P, NA, tolerance = 1e-4)
  expect_equal(metrics$YrYr, NA, tolerance = 1e-4)

  # neg recovery
  tsio <- c(rep(1,8), seq(-5, -10, by = -0.5), rep(0, 8))
  metrics <- calcSegRec(tsio, 2, maxBreak = F, obspyr, h, shortDenseTS, nPre, nDist, nPostMin, nPostMax, timeThres = 2,seas = F)
  expect_equal(metrics$RRI, NA, tolerance = 1e-4)
  expect_equal(metrics$R80P, NA, tolerance = 1e-4)
  expect_equal(metrics$YrYr, NA, tolerance = 1e-4)

  # second break within recovery period
  tsio <- c(rep(1,8), seq(-5, -2, by = 0.5),seq(-5, -3, by = 0.5), rep(0, 8))
  metrics <- calcSegRec(tsio, 2, maxBreak = F, obspyr, h, shortDenseTS, nPre, nDist, nPostMin, nPostMax, timeThres = 2,seas = F)
  expect_equal(metrics$RRI, NA, tolerance = 1e-4)
  expect_equal(metrics$R80P, NA, tolerance = 1e-4)
  expect_equal(metrics$YrYr, NA, tolerance = 1e-4)

})

test_that("Frazier - segmented annual - short", {

  tsio <- c(rep(1,8), seq(-5, 0, by = 0.5), rep(0,8))
  tdist <- 9
  obspyr <- 1
  shortDenseTS <- FALSE
  nPre <- 2
  nDist <- 0
  nPostMin <- 1
  nPostMax <- 1
  h <- 0.2
  seas <- F

  metrics <- calcSegRec(tsio, tdist, maxBreak = T, obspyr, h, shortDenseTS, nPre, nDist, nPostMin, nPostMax, timeThres = 2,seas = F)
  pre <- 1
  dnbr <- 6
  ari <- 0.5

  rrim <- ari/dnbr
  r80pm <- tsio[10]/(0.8*pre)
  yryrm <- (tsio[10] - tsio[9])/1

  expect_equal(metrics$RRI, rrim, tolerance = 1e-4)
  expect_equal(metrics$R80P, r80pm, tolerance = 1e-4)
  expect_equal(metrics$YrYr, yryrm, tolerance = 1e-4)
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
                                                nPre = 2, nDist = 0, nPostMin = 4, nPostMax = 5, h = 0.15, timeThres = 2, seas=F)})
  names(out) <- c('RRI', 'R80p', 'YrYr', 'missingVal', 'loglik', 'AIC')
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
  expect_equal(msked, 12, tolerance = 1e-4)
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

