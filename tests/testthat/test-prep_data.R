context("Data preparation")

test_that("Calculate Normalized Burn Ratio", {

  nir <- 0.2
  swir <- 0.8
  nbr <- (nir - swir)/(nir+swir)
  expect_equal(calcNBR(nir,swir), nbr)
})

test_that("Calculate QA", {
  expect_equal(getQA(illum = 2, slope = 1, wvp = 1), 28672)
})

test_that("convert binary to integer", {
  expect_equal(BinToDec(c(1,0)), 2)
  expect_equal(BinToDec(c(1,0,1)), 5)
})

test_that("change time span of raster stack", {
  empty_rast <- rast(nrows = 3, ncols = 4)

  r1 <- empty_rast; values(r1) <- 1:12
  r2 <- empty_rast; values(r2) <-  13:24
  r3 <- empty_rast; values(r3) <-  25:36
  r4 <- empty_rast; values(r4) <-  37:48
  st <- c(r1,r2,r3,r4)
  tempRes <- 'monthly'
  dts <- seq(as.Date('2000-01-01'), as.Date('2000-04-01'), by = '1 month')

  # case 1 - append start and end
  starttime <- as.Date('1999-8-1')
  endtime <- as.Date('2001-1-1')
  st1 <- setPeriod(st, starttime, endtime,tempRes, dts)
  exp <- c(NA, NA, NA, NA, NA,  1, 13, 25, 37, NA, NA, NA, NA, NA, NA, NA, NA, NA)
  expect_equal(as.numeric(st1[1]), exp)

  # case 2 - shorten start and extend end
  starttime <- as.Date('2000-2-1')
  endtime <- as.Date('2000-3-1')
  st2 <- setPeriod(st, starttime, endtime,tempRes, dts)
  exp <- c(13, 25)
  expect_equal(as.numeric(st2[1]), exp)

})

test_that("Get combination of QA values", {
  vals <- getQAvals(valid = c(0,1), aero = c(0,3))
  exp <- c(getQA(), getQA(valid = 1), getQA(aero = 3), getQA(valid = 1, aero = 3))
  expect_equal(vals, exp)
  vals <- getQAvals(illum = 2, slope = 1, wvp = 1)
  expect_equal(vals, 28672)
})

test_that("Mask using QA values",{
  library(terra)
  empty_rast <- rast(nrows = 3, ncols = 4)
  im <- empty_rast; values(im) <- 1:12
  qaim <- empty_rast; values(qaim) <- c(0, 1, 4, 9, 11, 12, 128, 132, 966, 8192, 8324, 28795)
  mskd <- mskQA(im, qaim, valid = 0, cloud = c(0,1), shadow = 0, snow = 0, water = 0, aero = c(0,1), subzero = 0, sat = 0, sunZen = 0, illum = 0, slope = c(0,1), wvp = 0)
  exp <- c(1,NA,3,NA,NA,NA,7,8,NA,10,11,NA)
  expect_equal(as.numeric(mskd[,]), exp)
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
  library(terra)
  library(lubridate)
  # mask
  empty_rast <- rast(ncol=3, nrow=3)
  m <- empty_rast; values(m) <- c(1,1,0,1,1,1,1,1,1)
  # rasters of time series observations
  r1<- empty_rast; values(r1) <- rep(1, 9)
  r2<- empty_rast; values(r2) <- rep(2, 9)
  r3<- empty_rast; values(r3) <- rep(1, 9)
  r4<- empty_rast; values(r4) <- rep(20, 9)
  r5<- empty_rast; values(r5) <- rep(1, 9)
  r6<- empty_rast; values(r6) <- rep(30, 9)
  r7<- empty_rast; values(r7) <- c(-12,1,1,1,1,1,1,1,1)
  r8<- empty_rast; values(r8) <- c(-12,1,1,1,1,1,1,1,1)+10
  r9<- empty_rast; values(r9) <- c(-11,1,1,-5,-3,-3,-3,-3,-3)
  r10<- empty_rast; values(r10) <- c(-11,1,1,-5,-3,-3,-3,-3,-3)-10
  r11<- empty_rast; values(r11) <- c(-10,1,1,-4,-2,-2,-2,-2,-2)
  r12<- empty_rast; values(r12) <- c(-10,1,1,-4,-2,-2,-2,-2,-2)-20
  r13<- empty_rast; values(r13) <- c(-9,1,1,-3,-1,-1,-1,-1,-1)
  r14<- empty_rast; values(r14) <- c(-9,1,1,-3,-1,-1,-1,-1,-1)-30
  r15<- empty_rast; values(r15) <- c(-8,1,1,-2,0,0,0,0,0)
  r16<- empty_rast; values(r16) <- c(-8,1,1,-2,0,0,0,0,0)-40
  r17<- empty_rast; values(r17) <- c(-7,1,1,-1,1,1,1,1,1)
  r18<- empty_rast; values(r18) <- c(-7,1,1,-1,1,1,1,1,1)-50
  r19<- empty_rast; values(r19) <- c(-6,1,1,0,1,1,1,1,1)
  r20<- empty_rast; values(r20) <- c(-6,1,1,0,1,1,1,1,1)-60
  #create stack
  st <- c(m,r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,
              r11,r12,r13,r14,r15,r16,r17,r18,r19,r20)
  # dates associated with observations
  dts <- as.Date(c('2001-01-02','2001-01-03','2001-02-02','2001-02-04','2001-03-02','2001-03-04','2001-04-02','2001-04-05','2001-05-02','2001-05-12',
                   '2001-06-02','2001-06-03','2001-07-02','2001-07-22','2001-08-02','2001-08-12','2001-09-02','2001-09-22','2001-12-02','2001-10-02'))
  # names(st) <- dts
  dts <- as.Date(dts, format = "X%Y.%m.%d") ## needed as input in the helper function of get_m_agg

  # create time series of monthly max
  brmomax <- app(st, function(x){toRegularTSStack(x, dts, fun = 'max', resol = 'monthly')})
  names(brmomax) <- as.Date(toRegularTS(dts, dts, fun = 'max', resol = 'monthly'))
  # create time series of monthly mean values
  brmomean <- app(st, function(x){toRegularTSStack(x, dts, fun = 'mean', resol = 'monthly')})
  names(brmomean) <- as.Date(toRegularTS(dts, dts, fun = 'mean', resol = 'monthly'))
  # create daily time series
  brday <- app(st, function(x){toRegularTSStack(x, dts, fun = 'max', resol = 'daily')})
  names(brday) <- date_decimal(as.numeric(time(bfastts(rep(1,length(dts)), dts, type = "irregular"))))
  # create quarterly time series
  brquart <- app(st, function(x){toRegularTSStack(x, dts, fun = 'max', resol = 'quart')})
  names(brquart) <- as.Date(toRegularTS(dts, dts, fun = 'mean', resol = 'quart'))

  # case 1 - monthly max
  expect_equal(as.numeric(brmomax[1,][1,]), c(2,20,30,-2,-11,-10,-9,-8,-7,-66,NA,-6), tolerance = 1e-4)
  # case 2 - monthly mean
  expect_equal(as.numeric(brmomean[1,][1,]), c(1.5,10.5,15.5,-7.0,-16.0,-20.0,-24.0,-28.0,-32.0,-66.0,NA,-6.0), tolerance = 1e-4)
  # case 3 - daily
  expect_equal(as.numeric(brday[1,][1,c(1,2,32,34,60,62,91,94,121,131,152,153,182,202,213,223,244,264,274,335)]),
               c(1,2,1,20,1,30,-12,-2,-11,-21,-10,-30,-9,-39,-8,-48,-7,-57,-66,-6), tolerance = 1e-4)
  # case 4 - quarterly
  expect_equal(as.numeric(brquart[1,][1,]),c(30,-2, -7,-6))
})


test_that("Create NBR stack", {
  library(terra)
  empty_rast <- rast(nrows = 3, ncols = 3)
  empty_rast2 <- rast(nrows = 9, ncols = 9)

  lnd5_1 <- empty_rast; values(lnd5_1) <- rep(1,9)
  lnd5_2 <- empty_rast; values(lnd5_2) <- rep(100,9)
  lnd5_3 <- empty_rast; values(lnd5_3) <- rep(55,9)
  lnd5_4 <- empty_rast; values(lnd5_4) <- 1:9
  lnd5_5 <- empty_rast; values(lnd5_5) <- rep(44,9)
  lnd5_6 <- empty_rast; values(lnd5_6) <- 11:19
  lnd5 <- c(lnd5_1, lnd5_2, lnd5_3, lnd5_4, lnd5_5, lnd5_6)

  sen2_1 <- empty_rast2; values(sen2_1) <- rep(5, times = 81)
  sen2_2 <- empty_rast2; values(sen2_2) <- rep(5, times = 81)
  sen2_3 <- empty_rast2; values(sen2_3) <- rep(5, times = 81)
  sen2_4 <- empty_rast2; values(sen2_4) <- rep(5, times = 81)
  sen2_5 <- empty_rast2; values(sen2_5) <- rep(5, times = 81)
  sen2_6 <- empty_rast2; values(sen2_6) <- rep(5, times = 81)
  sen2_7 <- empty_rast2; values(sen2_7) <- rep(5, times = 81)
  sen2_8 <- empty_rast2; values(sen2_8) <- c(rep(rep(41:43, each =3), times = 3),
                                             rep(rep(44:46, each =3), times = 3),
                                             rep(rep(47:49, each =3), times = 3))
  sen2_9 <- empty_rast2; values(sen2_9) <- rep(5, times = 81)
  sen2_10 <- empty_rast2; values(sen2_10) <- c(rep(rep(1:3, each =3), times = 3),
                                             rep(rep(4:6, each =3), times = 3),
                                             c(rep(7:9, each =3),
                                                 c(7,100,7,8,8,8,9,9,9),
                                                 rep(7:9, each =3)))
  sen2 <- c(sen2_1, sen2_2, sen2_3, sen2_4, sen2_5, sen2_6, sen2_7, sen2_8, sen2_9, sen2_10)

  lnd8_1 <- empty_rast; values(lnd8_1) <- rep(1,9)
  lnd8_2 <- empty_rast; values(lnd8_2) <- rep(100,9)
  lnd8_3 <- empty_rast; values(lnd8_3) <- rep(55,9)
  lnd8_4 <- empty_rast; values(lnd8_4) <- 81:89
  lnd8_5 <- empty_rast; values(lnd8_5) <- rep(44,9)
  lnd8_6 <- empty_rast; values(lnd8_6) <- 11:19
  lnd8 <- c(lnd8_1, lnd8_2, lnd8_3, lnd8_4, lnd8_5, lnd8_6)

  lnd7_1 <- empty_rast; values(lnd7_1) <- rep(1,9)
  lnd7_2 <- empty_rast; values(lnd7_2) <- rep(100,9)
  lnd7_3 <- empty_rast; values(lnd7_3) <- rep(55,9)
  lnd7_4 <- empty_rast; values(lnd7_4) <- 61:69
  lnd7_5 <- empty_rast; values(lnd7_5) <- rep(44,9)
  lnd7_6 <- empty_rast; values(lnd7_6) <- 11:19
  lnd7 <- c(lnd7_1, lnd7_2, lnd7_3, lnd7_4, lnd7_5, lnd7_6)

  lnd5_qai <- empty_rast; values(lnd5_qai) <- c(0,0,63, 64,0,8292,0,8256,2)
  sen2_qai <- empty_rast2; values(sen2_qai) <- c(0,0,0,0,0,0,63,63,63,
                                                0,0,0,0,0,0,63,63,63,
                                                0,0,0,0,0,0,63,63,63,
                                                64,64,64,0,0,0,8292,8292,8292,
                                                64,2,64,0,0,0,8292,8292,8292,
                                                64,64,64,0,0,0,8292,8292,8292,
                                                0,0,0,8256,8256,8256,2,2,2,
                                                0,2,0,8256,8256,8256,2,2,2,
                                                0,0,0,8256,8256,8256,2,2,2)
  lnd8_qai <- empty_rast; values(lnd8_qai) <- c(0,0,63, 64,0,8292,0,8256,0)
  lnd7_qai <- empty_rast; values(lnd7_qai) <- c(0,0,64, 64,0,8292,0,8256,2)

  img_list <- list(lnd5,sen2,lnd8,lnd7)
  qai_list <- list(lnd5_qai,sen2_qai,lnd8_qai,lnd7_qai)
  fmask <- empty_rast; values(fmask) <- c(1,0,1,1,NA,NaN,1,1,1)
  lsens <- c('LND05', 'SEN2A', 'LND08', 'LND07')
  ldts <- as.Date(c('2000-2-4','2000-3-15','2000-3-20','2001-8-2'))
  tempFun <- 'mean'
  tempRes <- 'monthly'
  starttime <- as.Date('2000-1-1')
  endtime <- as.Date('2001-12-31')

  st <- prepareNBRstack(img_list, qai_list, fmask, lsens, ldts, tempFun, tempRes, starttime, endtime)

  expect_equal(as.numeric(st[,][1,]),c(NA,(1-11)/(1+11),mean(c((41-1)/(42),(81-11)/(81+11))),rep(NA,16),(61-11)/(61+11),rep(NA,4)))
  expect_equal(sum(is.na(st[,][2,])),24)
  expect_equal(as.numeric(st[,][3,]),c(rep(NA,19),(63-13)/(63+13),rep(NA,4)))
  expect_equal(as.numeric(st[,][4,]),c(NA,(4-14)/(4+14),mean(c((44-4)/(48),(84-14)/(84+14))),rep(NA,16),(64-14)/(64+14),rep(NA,4)))
  expect_equal(sum(is.na(st[,][5,])),24)
  expect_equal(sum(is.na(st[,][6,])),24)
  expect_equal(as.numeric(st[,][7,]),c(NA,(7-17)/(7+17),mean(c((47-7)/(54),(87-17)/(87+17))),rep(NA,16),(67-17)/(67+17),rep(NA,4)))
  expect_equal(as.numeric(st[,][8,]),c(NA,(8-18)/(8+18),mean(c((48-8)/(56),(88-18)/(88+18))),rep(NA,16),(68-18)/(68+18),rep(NA,4)))
  expect_equal(as.numeric(st[,][9,]),c(NA,NA,(89-19)/(89+19),rep(NA,21)))

  tempRes <- 'quarterly'
  ldts <- as.Date(c('2000-2-4','2000-3-15','2000-5-20','2001-8-2'))
  st <- prepareNBRstack(img_list, qai_list, fmask, lsens, ldts, tempFun, tempRes, starttime, endtime)
  expect_equal(as.numeric(st[,][1,]),c(mean(c((1-11)/(1+11),(41-1)/42)),(81-11)/(81+11),NA,NA,NA,NA,(61-11)/(61+11),NA))
  expect_equal(sum(is.na(st[,][2,])),8)
  expect_equal(as.numeric(st[,][3,]),c(rep(NA,6),(63-13)/(63+13),NA))
  expect_equal(as.numeric(st[,][4,]),c(mean(c((4-14)/(4+14),(44-4)/48)),(84-14)/(84+14),rep(NA,4),(64-14)/(64+14),NA))
  expect_equal(sum(is.na(st[,][5,])),8)
  expect_equal(sum(is.na(st[,][6,])),8)
  expect_equal(as.numeric(st[,][7,]),c(mean(c((7-17)/(7+17),(47-7)/54)),(87-17)/(87+17),rep(NA,4),(67-17)/(67+17),NA))
  expect_equal(as.numeric(st[,][8,]),c(mean(c((8-18)/(8+18),(48-8)/56)),(88-18)/(88+18),rep(NA,4),(68-18)/(68+18),NA))
  expect_equal(as.numeric(st[,][9,]),c(NA,(89-19)/(89+19),rep(NA,6)))


  starttime <- as.Date('2000-1-1')
  endtime <- as.Date('2001-6-29')
  st <- prepareNBRstack(img_list, qai_list, fmask, lsens, ldts, tempFun, tempRes, starttime, endtime)
  expect_equal(as.numeric(st[,][1,]),c(mean(c((1-11)/(1+11),(41-1)/42)),(81-11)/(81+11),NA,NA,NA,NA))
  expect_equal(sum(is.na(st[,][2,])),6)
  expect_equal(sum(is.na(st[,][3,])),6)
  expect_equal(as.numeric(st[,][4,]),c(mean(c((4-14)/(4+14),(44-4)/48)),(84-14)/(84+14),rep(NA,4)))
  expect_equal(sum(is.na(st[,][5,])),6)
  expect_equal(sum(is.na(st[,][6,])),6)
  expect_equal(as.numeric(st[,][7,]),c(mean(c((7-17)/(7+17),(47-7)/54)),(87-17)/(87+17),rep(NA,4)))
  expect_equal(as.numeric(st[,][8,]),c(mean(c((8-18)/(8+18),(48-8)/56)),(88-18)/(88+18),rep(NA,4)))
  expect_equal(as.numeric(st[,][9,]),c(NA,(89-19)/(89+19),rep(NA,4)))

})


