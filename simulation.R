xlim <- 1000 ## 1km
ylim <- 1000 ## 1km
A <- 2.5
TRIALS <- 10000
U <- 0.05 ## 
noise <- `^`(10, -5)/U

source("coverage_base.R")

decibelToWatt = function(D) {
    `^`(10, D/10)
}

## Verified - function to find euclidean distance between 2 points
distance = function(x1, y1, x2, y2) {
    dist(rbind(c(x1, y1), c(x2, y2)))[1]
}

## SINR calculation for the ith tower, assuming interference from all other towers
sinr = function(i, x, y, g, h, N, R) {
    num <- h[i] * `^`(R[i], -1 * A)
    denom <- noise
    for(j in 1:N){
        if(j != i){
            denom <- denom + g[j] * `^`(R[j], -1 * A)
        }
    }
    return (num/denom)
}

coverageForT = function(L, T) {
    ## Number of towers is a poisson distribution
    Ns <- rpois(TRIALS, L)
    num <- 0
    denom <- TRIALS

    for(N in Ns){
        if(N != 0) {
            ## Location of all towers is uniformly random
            x <- runif(N, 0, xlim)
            y <- runif(N, 0, ylim)

            ## Transmit power is exponential
            h <- rexp(N, U)
            ## Rayleigh inteference - fading effect is exponential
            g <- rexp(N, U)

            R <- mapply(distance, x, y, MoreArgs = list(x2 = xlim/2, y2 = ylim/2))

            ## Find SINR at center for each tower
            I <- c(1:N)
            sinrs <- unlist(lapply(I, sinr, x = x, y = y, g = g, h = h, N = N, R = R))

            m <- max(sinrs)

            if(m > T)
                num <- num + 1
        }
        
    }
    return(num/denom)
}

plotCoverageActual = function(L, d){
    decibel = seq(-10, 20, d)
    watts = sapply(decibel, decibelToWatt)
    y1 = sapply(watts, coverageForT, L = L)
    y2 = sapply(watts, coverage, L =  L, SNR = (1/(U*noise)), alpha = A)
    header = paste("Coverage Probability (alpha, Noise, U, Lambda, Trials) = (", A, ",", noise, ",", U, ", ", L, ", ", TRIALS, ")", sep = " ")
    plot(decibel, y1, type = "l", xlim = c(-10, 20), ylim = c(0.1, 1), xlab = "SINR Threshold (dB)", ylab = "Probability of Coverage", main = header, col = "blue")
    lines(decibel, y2, col = "red")
}
