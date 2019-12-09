
#' Parzen window smoothing of power spectral densities
#'
#' @param freq frequency
#' @param u,bandwidth param of window.
#'   these two are related with each other as
#'   'u <- 280/bandwidth/151'
#' @param power power spectral density (PSD)
#' @param dat1,dat2 time-series data
#' @param dt time interval
#' @param smooth if FALSE, it returns raw periodgram

#' @export
parzen_window <- function(freq, u) {

    W <- 0.75*u*(sin(0.5*pi*u*freq)/(0.5*pi*u*freq))^4
    W[freq == 0] <- 0.75*u

    return(W)
}

#' @export
parzen_window_smoothing <- function(freq, power, bandwidth) {

    df <- freq[2] - freq[1]
    power_smth <- rep(0, length(power))

    # bandwidth to u
    u <- 280/bandwidth/151

    # convolution of power and window
    power_smth <- convolve(power, parzen_window(freq, u))*df

    return(list(freq = freq, power = power_smth))
}

#' @export
psd_parzen <- function (dat1, dat2, dt, bandwidth = 0.2, smooth = TRUE) {

    # fourier transform
    N <- length(dat1)

    ft1 <- fft(dat1)/N
    ft2 <- fft(dat2)/N
    frq <- (1:(N/2))/N/dt

    # periodgram
    Sf <- Conj(ft1)*ft2*N*dt
    Sf <- 2*Sf[2:(N/2+1)]  # two-sided to one-sided

    # smoothing
    if (smooth == TRUE) {
        Sf_smooth <- parzen_window_smoothing(frq, Sf, bandwidth)$power
    } else {
        Sf_smooth <- Sf
    }

    return(list(freq = frq, power = Sf_smooth))
}

