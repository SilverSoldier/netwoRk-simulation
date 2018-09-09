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

        ## find SINR over all towers totally
        I <- c(1:N)
        sinrs <- unlist(lapply(I, sinr, x = x, y = y, g = g, h = h, N = N, R = R))

        m <- max(sinrs)

        if(m > T)
            num <- num + 1
    }
    return (num/denom)
}

plotCoverageAS = function(d){
    L  <- c(5, 10, 7.5)

    decibel = seq(-10, 20, d)
    watts = sapply(decibel, decibelToWatt)

    y1 = sapply(watts, coverageActualAS, L = L)
    y2 = sapply(watts, coverageAO, lambdas = L, owner = 1, SNR = (1/(U*noise)), A = A)

    header = paste("Coverage Probability (alpha, Noise, U, Lambda, Trials, SNR) = (", A, ",", noise, ",", U, ", ", "(", paste(L[1], L[2], L[3], sep = ","), ")", ", ", TRIALS, ", ", (1/(U*noise)), ")", sep = " ")
    plot(decibel, y1, type = "l", xlim = c(-10, 20), ylim = c(0.1, 1), xlab = "SINR Threshold (dB)", ylab = "Probability of Coverage", main = header, col = "blue")
    lines(decibel, y2, col = "red")
}

coverageBAO = function(){
    N1s  <-  rpois(TRIALS, L[1])
    N2s  <-  rpois(TRIALS, L[2])
    N3s  <-  rpois(TRIALS, L[3])

    num  <- 0
    denom  <- TRIALS

    for(i in (1:TRIALS)) {
        N1  <- N1s[i]
        N2  <- N2s[i]
        N3  <- N3s[i]

        ## distribute each operator's towers uniformly
        x1  <- runif(N1, 0, xlim)
        x2  <- runif(N2, 0, xlim)
        x3  <- runif(N3, 0, xlim)

        y1  <- runif(N1, 0, ylim)
        y2  <- runif(N2, 0, ylim)
        y3  <- runif(N3, 0, ylim)

        ## Transmit power is exponential
        h1  <- rexp(N1, U)
        h2  <- rexp(N2, U)
        h3  <- rexp(N3, U)

        ## Fading is exponential
        g1  <- repx(N1, U)
        g2  <- rexp(N2, U)
        g3  <- rexp(N3, U)

        R1  <- mapply(distance, x1, y1, x2 = xlim/2, y2 = ylim/2)
        R2  <- mapply(distance, x2, y2, x2 = xlim/2, y2 = ylim/2)
        R3  <- mapply(distance, x3, y3, x2 = xlim/2, y2 = ylim/2)

        ## find max SINR for owner operator tower (assuming first operator)
        sinrs <- unlist(lapply(I, sinr, x = x1, y = y1, g = g1, h = h1, N = N1, R = R1))

        m  <- max(sinrs)
        if(m > T){
            ## find max SINR from all operators
        }
    }
}

