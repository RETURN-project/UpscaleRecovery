context("Simulation of benchmarking time series")

test_that("Exponential decay", {

  ts <- seq(0, 25, by = 0.1)
  pert <- 1.5
  tpert <- 0
  thalf <- 1

  ys <- exponential(ts, pert = pert, tpert = tpert, thalf = thalf)

  expect_equal(abs(min(ys)), 0, tolerance = 1e-4)
  expect_equal(abs(max(ys)), abs(pert), tolerance = 1e-4)

})

test_that("Linear decay", {

  ts <- seq(0, 5, by = 0.01)
  pert <- 1.5
  ys <- piecewise(ts, pert = pert)

  expect_equal(abs(min(ys)), 0, tolerance = 1e-4)
  expect_equal(abs(max(ys)), abs(pert), tolerance = 1e-4)
})

test_that("Realistic decay (deterministic)", {

  ts <- seq(0, 25, by = 0.1)
  pert <- 1.5
  tpert <- 0
  thalf <- 1
  offset <- 2

  ys_exp <- exponential(ts, offset = offset, pert = pert, tpert = tpert, thalf = thalf)
  ys_rea <- realistic(ts, offset = offset, pert = pert, tpert = tpert, thalf = thalf)

  # In the absence of noise, realistic and exponential decays should be identical
  expect_equal(ys_rea, ys_exp, tolerance = 1e-4)
})

test_that("Disturbance simulation",{
  # source('../R/fun_simulate.R')

  # disturbance and recovery before end of time series
  distT1 <- 5
  distRec1 <- 12/2
  distMag1 <- -10
  nobs1 <- 20

  ts1 <- piecewise(1:nobs1, pert = distMag1, tpert = distT1, thalf = distRec1)

  # recovery after end of time series
  distT2 <- 5
  distRec2 <- 20/2
  distMag2 <- -10
  nobs2 <- 20

  ts2 <- piecewise(1:nobs2, pert = distMag2, tpert = distT2, thalf = distRec2)

  # disturbance after end of time series
  distT3 <- 25
  distRec3 <- 16/2
  distMag3 <- -10
  nobs3 <- 20

  ts3 <- piecewise(1:nobs3, pert = distMag3, tpert = distT3, thalf = distRec3)

  # tests
  expect_equal(distMag1,min(ts1))# magnitude disturbance
  expect_equal(distT1,which(ts1 == min(ts1)))# timing disturbance
  expect_equal(distRec1,length(which(ts1 < 0)), tolerance = 1)# recovery period

  expect_equal(distMag2,min(ts2))# magnitude disturbance
  expect_equal(distT2,which(ts2 == min(ts2)))# timing disturbance
  expect_equal(distMag2/((2*distRec2)), (ts2[distT2] - ts2[distT2 + 1]))# recovery slope

  expect_equal(nobs3,length(which(ts3 == 0)))# no disturbance and recovery
})



