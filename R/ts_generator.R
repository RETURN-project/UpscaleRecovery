#' Piecewise linear decay
#'
#' @param t Times to simulate
#' @param offset Offset
#' @param pert Perturbation intensity (signed)
#' @param tpert Perturbation timing
#' @param thalf Perturbation half-life
#' @param noise Strength of the additive white noise (standard deviation)
#'
#' @return The time series
#' @export
piecewise <- function(t, offset = 0, pert = 0, tpert = 0, thalf = 1, noise = 0) {
  m <- -pert / (2 * thalf) # Slope of the transitory regime
  ttrans <- 2*thalf # Duration of the transitory regime
  y <- offset                             * (t < tpert) +
    (offset + pert + m *(t - tpert))    * (t >= tpert) * (t <= tpert + ttrans) + # Transitory regime
    offset                             * (t > tpert + ttrans)

  y <- y + rnorm(length(t), sd = noise) # Add the noise

  return(y)
}


#' Exponential decay
#'
#' @param t Times to simulate
#' @param offset Offset
#' @param pert Perturbation intensity (signed)
#' @param tpert Perturbation timing
#' @param thalf Perturbation half-life
#' @param noise Strength of the additive white noise (standard deviation)
#'
#' @return The time series
#' @export
exponential <- function(t, offset = 0, pert = 0, tpert = 0, thalf = 1, noise = 0) {
  r <- log(2)/thalf # Translate the half-life to a multiplicative constant
  y <- offset + pert * exp(-r*(t-tpert)) * (t >= tpert)

  y <- y + rnorm(length(t), sd = noise) # Add the noise

  return(y)
}


#' Realistic time series simulation
#'
#' Simulates a return to equilibrium after a perturbation under the influence of a stochastic differential equation where:
#'
#' 1. The deterministic dynamics are given by an exponential decay with the given half life
#' 2. The stochastic dynamics have the given infinitesimal standard deviation
#' 3. Gaussian noise is added to the final time series to simulate measurement errors
#'
#' @param t Times to simulate
#' @param offset Offset
#' @param pert Perturbation intensity (signed)
#' @param tpert Perturbation timing
#' @param thalf Perturbation half-life
#' @param noise Strength of the stochastic term in the differential equation (standard deviation of the integrated time series)
#'
#' @return The time series
#' @export
realistic <- function(t, offset = 0, pert = 0, tpert = 0, thalf = 1, noise = 0) {

  # Translate parameters to the language of differential equations
  y0 <- pert # The perturbation represents the initial condition
  r <- log(2)/thalf # Translate the half-life to a multiplicative constant
  sigma <- noise * sqrt(2 * log(2) / thalf) # Infinitesimal standard deviation
  # The infinitesimal standard deviation `sigma` yields a standard deviation of magnitude `noise` after integration
  # Reference: https://math.stackexchange.com/questions/2558659/expectation-and-variance-of-stochastic-differential-equations

  # Pose the differential equation dy = f(y,t) dt + g(y,t) dW
  #
  # With
  # f(y, t) = -r * y (exponential decay)
  # and
  # g(y, t) = s (white noise)
  #
  # Unfortunately, the package sde uses a very obscure syntax. Instead of functions it expects
  # expressions depending on x and t as an input. When those object contain, additionally,
  # parameters, we need to pass them via the substitute command.
  #
  # e.g: substitute(a + x, list(a = 2)) returns 2 + x
  # 2 + x is an object of class call, that must be converted to expression
  f <- as.expression(
    substitute(-r * x,
               list(r = r))
  )
  g <- as.expression(
    substitute(s,
               list(s = sigma))
  )


  # Solve
  sol <- sde::sde.sim(X0 = y0,
                      T = max(t),
                      N = length(t)-1,
                      drift = f,
                      sigma = g,
                      sigma.x = 0.0,
                      method = 'euler')


  # Shift the time series (only if tpert is not zero)
  if(tpert != 0) {
    # Shift the original time series
    ts <- time(sol) # Store the original times
    i <- min(which(ts >= tpert)) # Find the index corresponding to tpert
    sol <- lag(sol, -i + 1) # Displace the time series, so it begins at tpert

    # Create the time series before tpert (only dynamic noise around equilibrium)
    tfill <- seq(t[1], t[i-1], by = 1 / frequency(sol))
    fill <- ts(realistic(tfill, noise = noise), start = t[1], end = t[i-1], frequency = frequency(sol))

    # Paste both time series together
    sol <- ts(c(fill, sol), start = start(fill), frequency = frequency(fill))

    # Trim the tail, so the displaced time series is equal in size to the original
    sol <- head(sol, length(sol) - length(fill))
  }

  # Extract the state only (so the output has the same structure as in piecewise and exponential)
  sol <- as.data.frame(sol)
  y <- as.numeric(sol$x)

  # Don't forget to add the offset
  y <- y + offset

  return(y)
}
