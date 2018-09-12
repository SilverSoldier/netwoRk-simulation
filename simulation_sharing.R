xlim  <- 1000
ylim  <- 1000
TRIALS  <- 10000
A  <- 2.5
U  <- 0.05
noise  <- `^`(10, -5)/U

source("coverage_sharing.R")

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

coverageActualAS = function(T, L){
    N1s  <-  rpois(TRIALS, L[1])
    N2s  <-  rpois(TRIALS, L[2])
    N3s  <-  rpois(TRIALS, L[3])

    num <- 0
    denom <- TRIALS

    for(i in (1:TRIALS)) {
        N1 <- N1s[i]
        N2 <- N2s[i]
        N3 <- N3s[i]

        if(N1 == 0 || N2 == 0 || N3 == 0)
            next;

        N  <- N1 + N2 + N3
        ## distribute each operator's towers uniformly
        x  <- runif(N, 0, xlim)
        y  <- runif(N, 0, ylim)
        
        ## Transmit power is exponential
        h  <- rexp(N, U)

        ## Fading is exponential
        g  <- rexp(N, U)

        R <- mapply(distance, x, y, MoreArgs = list(x2 = xlim/2, y2 = ylim/2))

        ## find max SINR over all towers
        I <- c(1:N)
        sinrs <- unlist(lapply(I, sinr, x = x, y = y, g = g, h = h, N = N, R = R))

        m <- max(sinrs)

        if(m > T)
            num <- num + 1
    }
    return (num/denom)
}

plotCoverageActualAS = function(L, d){
    decibel = seq(-10, 20, d)
    watts = sapply(decibel, decibelToWatt)

    y1 = sapply(watts, coverageActualAS, L = L)
    y2 = sapply(watts, coverageAO, lambdas = L, owner = 1, SNR = (1/(U*noise)), A = A)

    header = paste("Coverage Probability (alpha, Noise, U, Lambda, Trials, SNR) = (", A, ",", noise, ",", U, ", ", "(", paste(L[1], L[2], L[3], sep = ","), ")", ", ", TRIALS, ", ", (1/(U*noise)), ")", sep = " ")
    plot(decibel, y1, type = "l", xlim = c(-10, 20), ylim = c(0.1, 1), xlab = "SINR Threshold (dB)", ylab = "Probability of Coverage", main = header, col = "blue")
    lines(decibel, y2, col = "red")
    legend("top", legend=c("Actual", "Theoretical"), col=c("blue", "red"), lty=c(1, 1))
}

coverageActualBAO = function(T, L){
    N1s  <-  rpois(TRIALS, L[1])
    N2s  <-  rpois(TRIALS, L[2])
    N3s  <-  rpois(TRIALS, L[3])

    num <- 0
    denom <- TRIALS

    for(i in (1:TRIALS)) {
        N1 <- N1s[i]
        N2 <- N2s[i]
        N3 <- N3s[i]

        if(N1 == 0 || N2 == 0 || N3 == 0)
            next;

        N  <- N1 + N2 + N3
        ## distribute each operator's towers uniformly
        x  <- runif(N, 0, xlim)
        y  <- runif(N, 0, ylim)
        
        ## Transmit power is exponential
        h  <- rexp(N, U)

        ## Fading is exponential
        g  <- rexp(N, U)

        R <- mapply(distance, x, y, MoreArgs = list(x2 = xlim/2, y2 = ylim/2))

        ## find SINR over all towers totally
        I <- c(1:N)
        sinrs <- unlist(lapply(I, sinr, x = x, y = y, g = g, h = h, N = N, R = R))

        ## find max own operator SINR
        m1 <- max(sinrs[1:N1])
        m <- max(sinrs)

		## Minimum connection to own operator tower and max total SINR > T
        if(m1 > decibelToWatt(-20) && m > T)
            num <- num + 1
    }
    return (num/denom)
}

plotCoverageActualBAO = function(d){
    L  <- c(5, 10, 7.5)

    decibel = seq(-10, 20, d)
    watts = sapply(decibel, decibelToWatt)

    y1 = sapply(watts, coverageActualBAO, L = L)
    y2 = sapply(watts, coverageBAO, L = L, owner = 1, SNR = (1/(U*noise)), A = A)

    header = paste("Coverage Probability (alpha, Noise, U, Lambda, Trials, SNR) = (", A, ",", noise, ",", U, ", ", "(", paste(L[1], L[2], L[3], sep = ","), ")", ", ", TRIALS, ", ", (1/(U*noise)), ")", sep = " ")
    plot(decibel, y1, type = "l", xlim = c(-10, 20), ylim = c(0.1, 1), xlab = "SINR Threshold (dB)", ylab = "Probability of Coverage", main = header, col = "blue")
    lines(decibel, y2, col = "red")
    legend("top", legend=c("Actual", "Theoretical"), col=c("blue", "red"), lty=c(1, 1))
}

plotCoverageBAOAS = function(d){
    L  <- c(5, 10, 7.5)

    decibel = seq(-10, 20, d)
    watts = sapply(decibel, decibelToWatt)

    y1 = sapply(watts, coverageActualBAO, L = L)
    y2 = sapply(watts, coverageActualAS, L = L)
    y3 = sapply(watts, coverageBAO, L = L, owner = 1, SNR = (1/(U*noise)), A = A)
    y4 = sapply(watts, coverageAO, lambdas = L, owner = 1, SNR = (1/(U*noise)), A = A)

    header = paste("Coverage Probability (alpha, Noise, U, Lambda, Trials, SNR) = (", A, ",", noise, ",", U, ", ", "(", paste(L[1], L[2], L[3], sep = ","), ")", ", ", TRIALS, ", ", (1/(U*noise)), ")", sep = " ")
    plot(decibel, y1, type = "l", xlim = c(-10, 20), ylim = c(0.1, 1), xlab = "SINR Threshold (dB)", ylab = "Probability of Coverage", main = header, col = "blue")
    lines(decibel, y2, col = "red")
	lines(decibel, y3, col = "blue", lty = 2)
	lines(decibel, y4, col = "red", lty = 2)
    legend("top", legend=c("Base Assisted Offloading Simulation", "Autonomous Selection Simulation"), col=c("blue", "red", "blue", "red"), lty=c(1, 1, 2, 2))
}
