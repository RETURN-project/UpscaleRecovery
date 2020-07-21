#' Calculate recovery metrics from a time series with known disturbance date. The calcFrazier function derives the RRI, R80P and YrYr recovery indicators, defined by Frazier et al. (2018). The indicators are originally developped for annual long-term time series of optical vegetation indices (the indicators are shown in the figures below). Yet, in order to be able to derive the indicators as well for dense and/or short time series, a modified version is suggested. Here, the user can define the time period before, during and after the disturbance that is used to derive the indicators. To reduce the interference of the seasonal pattern of dense time series, the chosen time period should cover blocks of n years. Moreover, given the potentially high noise levels of dense time series, the mean value instead of the maximum value was used in the formulas. (Frazier, R. J., Coops, N. C., Wulder, M. A., Hermosilla, T., & White, J. C. (2018). Analyzing spatial and temporal variability in short-term rates of post-fire vegetation return from Landsat time series. Remote Sensing of Environment, 205, 32-45.)
#'
#' @param tsio vector of observations (time series with a fixed observation frequency)
#' @param tdist observation number of disturbance, indicating the timing of the disturbance
#' @param obspyr number of observations per year
#' @param shortDenseTS TRUE or FALSE. In case TRUE, the metrics are adjusted to be compatible with short, dense time series
#' @param nPre If shortDenseTS is TRUE, number of years prior to the disturbance used to calculate the pre-disturbance value
#' @param nDist If shortDenseTS is TRUE, number of years used to quantify the time series value during the disturbance
#' @param nPostMin If shortDenseTS is TRUE,  the post-disturbance condition is quantified starting from nPostMin years after the disturbance
#' @param nPostMax If shortDenseTS is TRUE, max number of years after the disturbance used to quantify the post-disturbance condition
#'
#' @return a list containing the RRI recovery indicator, R80p recovery indicator and YrYr recovery indicator
#' @export
#'
calcFrazier <- function(tsio, tdist, obspyr, shortDenseTS, nPre, nDist, nPostMin, nPostMax){
  # check if there are enough observations before and after the disturbance to calculate the metrics
  if((tdist>((nPre*obspyr))) & (tdist < (length(tsio)-(nPostMax*obspyr)+1)) & (sum(!is.na(tsio))>2)){
    # translate parameters to those needed for the recovery functions
    ys <- tsio # response
    ts <- seq(1,length(tsio))# time

    # disturbance dates
    if (obspyr == 1 | nDist == 0){
      tpert <- seq(tdist,tdist+(nDist*obspyr))
    }else{
      tpert <- seq(tdist,tdist+(nDist*obspyr)-1)
    }
    # pre-disturbance period
    ts_pre <- seq(tdist-(nPre*obspyr),tdist-1)
    # post-disturbance period
    if (obspyr == 1 | nPostMin == nPostMax){
      ts_post <-  seq(tdist +(nPostMin*obspyr), tdist +(nPostMax*obspyr))
    }else{
      ts_post <-  seq(tdist +(nPostMin*obspyr), tdist +(nPostMax*obspyr)-1)
    }
    # timing between disturbance and post-disturbance state
    deltat <- switch(shortDenseTS + 1, (nPostMax*obspyr), ts_post-tdist)

    # derive recovery indicators
    RRI <- rri(ts,ys,tpert,ts_pre, ts_post)
    R80P <- r80p(ts,ys,r = 0.8,ts_pre, ts_post)
    YrYr <- yryr(ts,ys,tpert, deltat)
    # make list of recovery indicators as output of the function
    lst <- list(RRI, R80P, YrYr)
    names(lst) <- c('RRI', 'R80P', 'YrYr')
    # give NA as output if not able to calculate the recovery indicators
  }else{
    lst <- list(NA, NA, NA)
    names(lst) <- c('RRI', 'R80P', 'YrYr')
  }
  lst
}

#' Post-disturbance slope and recovery metrics derived from BFAST0n trend segments. The calcBFASTrec function derives a set of recovery indicators after fitting a segmented trend in the time series. Using the breakpoints function of the strucchange package, a segmented trend is fitted (hereafter called BFAST0n trend segments). The detected break showing the largest change (in absolute values) is assumed to represent the disturbance. Using the segmented trend and detected disturbance date, the RRI, R80p, YrYr and the slope of the post-disturbance trend segment are derived as recovery indicators.
#'
#' @param tsio vector of observations (time series)
#' @param obspyr number of observations in one year
#' @param h This parameter defines the minimal segment size either given as fraction relative to the sample size or as an integer giving the minimal number of observations in each segment.
#' @param shortDenseTS TRUE or FALSE. In case FALSE, the metrics follow closely the definitions given by Frazier et al
#' @param nPre number of years prior to the disturbance used to calculate the pre-disturbance value
#' @param nDist number of years used to quantify the time series value during the disturbance
#' @param nPostMin min number of years after the disturbance used to quantify the recovery
#' @param nPostMax max number of years after the disturbance used to quantify the recovery
#' @param tdist the timing of the disturbance [observation number]
#' @param maxBreak only for recovery indicators derived from piecewise regression: if maxbreak is true, the maximum break in the segmented series is used as disturbance date to calculate the recovery indicators. If maxbreak is false, the break closest to the provided disturbance timing is used to calculate recovery.
#' @param seas should a seasonal comonent be used in the piecewise regression?
#' @param timeThres threshold on the duration between the disturbance date and date of the detected break [years]
#'
#' @return a list containing  the RRI, R80p, YrYr recovery indicator derived from the BFAST0n trend segments and slope of the trend segment after the disturbance (sl).
#' @export
#' @import strucchange
#' @import stats
#' @import bfast
calcSegRec <- function(tsio, tdist, maxBreak, obspyr, h, shortDenseTS, nPre, nDist, nPostMin, nPostMax, timeThres, seas = F){
  # Create time series object, needed as input for BFAST
  tsi <- ts(tsio, frequency = obspyr)
  # Convert the time series object into a dataframe, needed for the breakpoints function
  if(obspyr>1){
    datapp <- bfastpp(tsi, order = 1, lag = NULL, slag = NULL,
                      na.action = na.omit, stl = 'none')
  }else if(!seas){
    datapp <- data.frame(response = tsio, trend = seq(1:length(tsio)))
  }else{stop('No seasonal term allowed for time series with one observation per year or less.')}

  nreg <- switch(seas+1, 2, 5)
  # Test if enough observations are available to fit piecewise model
  if(floor(length(tsio[is.na(tsio)==F]) * h) > nreg){

      # Apply BFAST0n on time series: find breaks in the regression
    if (seas){
      bp <- breakpoints(response ~ trend + harmon, data = datapp, h = h)#, breaks = breaks
    } else{
      bp <- breakpoints(response ~ trend, data = datapp, h = h)##, breaks = breaks
    }
      # Check if BFAST0n found breakpoints
      if(is.na(bp$breakpoints[1])){# no breakpoint found
        # tr <- fitted(bp, 0)
        # sl <- (tr[2] - tr[1])
        frz <- list(NA, NA, NA, NA, NA)
        names(frz) <- c('RRI', 'R80P', 'YrYr', 'loglik', 'AIC')
      }else{# at least one breakpoint found
        # Extract BFAST trend component and breaks
        cf <- coef(bp)
        # Extract BFAST trend component and breaks
        tbp <- bp$breakpoints #observation number of break
        #tr <- rep(NA,length(tsi))
        indna <- which(is.na(tsi)==F)
        tbp <- indna[tbp]   # correct observation number for missing values
        totbp <- tbp

        #Derive trend component without missing values
        bpf <- c(0, tbp, length(tsi))
        trf <- rep(NA,length(tsi))
        for(ti in 1:(length(bpf)-1)){
          trf[(bpf[ti]+1):bpf[ti+1]] <- cf[ti,1] + ((cf[ti,2]*((bpf[ti]+1):bpf[ti+1])))
        }

        # Get information criteria
        bp_loglik <- logLik(bp)
        bp_aic <- AIC(bp)[length(tbp) + 1]

        if(maxBreak){
          # Find the major break
          dbr <- trf[tbp+1]-trf[tbp]
          tbp <- tbp[which(abs(dbr) == max(abs(dbr)))]
        }else{
          # Use the break closest to the disturbance date
          dbr <- tbp-tdist
          tbp <- tbp[which(abs(dbr) == min(abs(dbr)))]
        }

        # check the time period between the break and the fire
        timeChck <- ((min(abs(dbr))/obspyr) < timeThres)

        # check the typology of the segments:
        # constant pre-disturbance period
        # preChck <- (abs(trf[tbp] - trf[tbp-1]) < slpThres)
        # positive post-disturbance slope
        postChck <- ((trf[tbp+3] - trf[tbp+2]) > 0)
        # negative break
        distChck <- ((trf[tbp+1] - trf[tbp]) < 0)
        # no negative break in recovery period
        brkthres <- 1+(nPostMax*obspyr) #post-disturbance period used to assess recovery
        if(any((totbp>tbp) & (totbp<(brkthres+tbp)))){
          postbr <- totbp[(totbp>tbp) & (totbp<(brkthres+tbp))]
          postdbr <- trf[postbr+1]-trf[postbr]
          # postsl <- trf[postbr+3]-trf[postbr+2]
          brkChck <- !any((postdbr<0))# & (postsl>slpThres)
        }else{
          brkChck <- TRUE
        }

        if(timeChck & postChck & distChck & brkChck){
          # Calculate Frazier recovery metrics on BFAST trend component
          frz <- calcFrazier(as.numeric(trf), (tbp+1), floor(obspyr), shortDenseTS, nPre, nDist, nPostMin, nPostMax)
          # Calculate the post-disturbance slope of the BFAST trend component (first segment after break)
          sl <- (trf[tbp+3] - trf[tbp+2])# Calculate Frazier recovery metrics on BFAST trend component
          frz <- calcFrazier(as.numeric(trf), (tbp+1), floor(obspyr), shortDenseTS, nPre, nDist, nPostMin, nPostMax)
          # Calculate the post-disturbance slope of the BFAST trend component (first segment after break)
          # sl <- (trf[tbp+3] - trf[tbp+2])
          frz <- c(frz, bp_loglik, bp_aic)
          names(frz) <- c('RRI', 'R80P', 'YrYr', 'loglik', 'AIC')

        }else{
          frz <- list(NA, NA, NA, NA, NA)
          names(frz) <- c('RRI', 'R80P', 'YrYr', 'loglik', 'AIC')
        }
      }
    }else{
      frz <- list(NA, NA, NA, NA, NA, NA)
      names(frz) <- c('RRI', 'R80P', 'YrYr', 'Sl', 'loglik', 'AIC')
    }

    frz
}

#' Calculate recovery for a single time series
#'
#' @param tsi vector: the first n values contain the timing of the disturbances and the next n values the observations for which the recovery indicators should be computed
#' @param maxBreak (only for recovery indicators derived after piecewise regression): if maxbreak is true, the maximum break in the segmented series is used as disturbance date to calculate the recovery indicators. If maxbreak is false, the break closest to the provided disturbance timing is used to calculate recovery.
#' @param obspyr the number of observations per year
#' @param inp the preprocessing applied to the time series before computing the recovery indicators: segmented (for piecewise regression), smooth (time series smoothing using loess), or raw (no preprocessing)
#' @param shortDenseTS In case FALSE, the metrics follow closely the definitions given by Frazier et al
#' @param nPre number of years prior to the disturbance that are used to derive the pre-disturbance condition
#' @param nDist number of observations used to derive the disturbance state
#' @param nPostMin start of the post-disturbance period: number of years after the disturbance
#' @param nPostMax end of the post-disturbance period: number of years after the disturbance
#' @param h h parameter of the breakpoints function in the strucchange package
#' @param seas should a seasonal comonent be used in the piecewise regression?
#' @param timeThres only relevant for piecewise regression: threshold on the duration between the disturbance date and date of the detected break [years]
#'
#' @return the RRI, R80P, YrYr and the slope of the pos-disturbance segment
#' @export
#' @import stats
calcRecoveryTS <- function(tsi, maxBreak, obspyr, inp = 'segmented', shortDenseTS = TRUE,
                         nPre = 2, nDist = 12, nPostMin = 4, nPostMax = 6, h = 0.15, timeThres, seas){
  # if (frq == 'annual'){
  #   #convert time series to annual values by selecting date closest to seasonal max
  #   tsi <- toAnnualTS(tsseas, tsi, obspyr)
  #   tdist <- ceiling(tdist/obspyr)
  #   obspyr <- 1
  # }
  len <- length(tsi)
  tdist <- tsi[1:(len/2)]# data about the timing of the disturbance (0 equals no disturbance happened, 1 equals a disturbance took place)
  tsi <- tsi[(1+(len/2)):len]# data about the vegetation response

  # Find the major break given by the disturbance dataset
  tdist <- which(tdist == 1)# get all dates where a disturbance took place
  if(length(tdist) > 0){
    dbr <- rep(NA,length(tdist))
    for(i  in 1:length(tdist)){
      if((tdist[i] > (2*obspyr)) & (tdist[i] < (length(tsi)-obspyr+1))){
        # calculate magnitude break (difference between 2 years pre and 1 year post disturbance) for each disturbance date
        dbr[i] <- mean(tsi[(tdist[i]+1):(tdist[i]+(1*obspyr))], na.rm = T)-mean(tsi[(tdist[i]-2*obspyr):tdist[i]-1], na.rm = T)
      }else{
        dbr[i] <- 0
      }
    }
    tdist <- tdist[which(abs(dbr) == max(abs(dbr)))]# find largest break
    mval <- sum(is.na(tsi))/length(tsi) #fraction of missing values

    if (inp == 'smoothed'){
      # smooth the time series using a loess filter
      df <- data.frame(dat = tsi,tm = 1:length(tsi))
      ls <- loess(dat  ~ tm, df, span = 0.2)#, span = 0.5
      tmps <- predict(ls, data.frame(tm = 1:length(tsi)), se = TRUE)
      tsi <- tmps$fit
    }

    if((inp == 'smoothed') | (inp == 'raw')){
      # calculate recovery indicators
      tmp <- calcFrazier(tsi, tdist, obspyr, shortDenseTS, nPre, nDist, nPostMin, nPostMax)
      outp <- c(tmp$RRI,tmp$R80P, tmp$YrYr, mval,NA,NA)
    }
    if(inp == 'segmented'){
      # calculate recovery indicators after piecewise regression
      tmp <- calcSegRec(tsi, tdist, maxBreak, obspyr, h, shortDenseTS, nPre, nDist, nPostMin, nPostMax, timeThres, seas)
      outp <- c(tmp$RRI,tmp$R80P, tmp$YrYr, mval, tmp$loglik, tmp$AIC)
    }
  }else{outp <- c(NA,NA,NA,NA,NA,NA)}
  outp
}

#' Calculate recovery indicators from a time series stack
#'
#' @param st raster stack, the first raster is a mask (pixels with value 0 are not considered, pixels with value 1 are included), the next n rasters represent disturbances (0 if no disturbance occours, 1 if a disturbance occurs) and the last n rasters contain the time series observations
#' @param maxBreak (only for recovery indicators derived after piecewise regression): if maxbreak is true, the maximum break in the segmented series is used as disturbance date to calculate the recovery indicators. If maxbreak is false, the break closest to the provided disturbance timing is used to calculate recovery.
#' @param obspyr the number of observations per year
#' @param inp the preprocessing applied to the time series before computing the recovery indicators: segmented (for piecewise regression), smooth (time series smoothing using loess), or raw (no preprocessing)
#' @param shortDenseTS In case FALSE, the metrics follow closely the definitions given by Frazier et al
#' @param nPre number of years prior to the disturbance that are used to derive the pre-disturbance condition
#' @param nDist number of observations used to derive the disturbance state
#' @param nPostMin start of the post-disturbance period: number of years after the disturbance
#' @param nPostMax end of the post-disturbance period: number of years after the disturbance
#' @param h h parameter of the breakpoints function in the strucchange package
#' @param timeThres only relevant for piecewise regression: threshold on the duration between the disturbance date and date of the detected break [years]
#' @param seas only relevant for piecewise regression: should a seasonality term be used?
#'
#' @return a raster with the RRI, R80P, YrYr and the slope of the pos-disturbance segment
#' @export
#'
calcRecoveryStack <- function(st, maxBreak, obspyr, inp = 'segmented', shortDenseTS = TRUE,
                              nPre = 2, nDist = 12, nPostMin = 4, nPostMax = 6, h = 0.15, timeThres, seas) {

  # only consider time series with a mask value of 1
  msk <- (st[,1] == 1)
  st <- st[,-1]
  out <- matrix(NA, length(msk), 6)
  if(sum(msk)>1){
    out[msk,] <- t(apply(st[msk,], 1, calcRecoveryTS, maxBreak=maxBreak, obspyr=obspyr, inp=inp, shortDenseTS=shortDenseTS,
                         nPre=nPre, nDist=nDist, nPostMin=nPostMin, nPostMax=nPostMax, h=h, timeThres, seas))
  }
  if(sum(msk)==1){
    out[msk,] <- calcRecoveryTS(st[msk,], maxBreak, obspyr, inp, shortDenseTS,
                              nPre, nDist, nPostMin, nPostMax, h, timeThres, seas)
  }
  out
}

# toAnnualStack <- function(){
#
# }
