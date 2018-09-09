source("coverage_base.R")

coverageBAO = function(T, n, lambdas, owner, SNR){
    ## For non-owner operator towers, p(O_i) * cov(O_i) * cov(O_owner), for owner, p(O_owner) * cov(0_owner)
    cov <- 0
    denom <- Reduce("+", lambdas)
    ownerCov <- coverage(T, lambdas[owner], SNR)
    for (i in 1:n) {
        if (i != owner) {
            cov <- cov + (lambdas[i]/denom) * coverage(T, lambdas[i], SNR) * ownerCov
        }
        else {
            cov <- cov + (lambdas[owner]/denom) * ownerCov
        }
    }
    return(cov)
}

coverageAO = function(T, lambdas, owner, SNR, A){
    ## For all operator towers, p(O_i) * cov(O_i)
    cov <- 0
    denom <- Reduce("+", lambdas)
    n  <- length(lambdas)
    for(i in 1:n){
        cov <- cov + (lambdas[i]/denom)*coverage(T, lambdas[i], A, SNR)
    }
    return(cov)
}

plotCoverageAO = function () {
    decibel = seq(-10, 20, 0.1)
    watts = sapply(decibel, decibelToWatt)
    L <- c(5, 10)
    y1 = sapply(watts, coverageAO, lambdas = L, owner = 1, A = 2.5, SNR = `^`(10, 5))
    y2 = sapply(watts, coverage, L =  5, SNR = `^`(10, 5), alpha = 2.5)
    ## y3 = sapply(watts, coverage, L =  0.45, SNR = 10, alpha = 2.5)
    plot(decibel, y1, type = "l", xlim = c(-10, 20), ylim = c(0.1, 1), xlab = "SINR Threshold (dB)", ylab = "Probability of Coverage", main = "Coverage Probability (alpha = 4, No noise, SNR = 10)", col = "red")
    ## lines(decibel, y2, type = "l", col = "blue")
    ## lines(decibel, y3, type = "l", col = "green")
}
    
plotCoverageBAO = function () {
    decibel = seq(-10, 20, 0.1)
    watts = sapply(decibel, decibelToWatt)
    L <- c(0.25, 0.45)
    y1 = sapply(watts, coverageBAO, lambdas = L, n = length(L), owner = 1, SNR = 10)
    y2 = sapply(watts, coverage, L =  0.25, SNR = 10)
    y3 = sapply(watts, coverage, L =  0.45, SNR = 10)
    plot(decibel, y1, type = "l", xlim = c(-10, 20), ylim = c(0.1, 1), xlab = "SINR Threshold (dB)", ylab = "Probability of Coverage", main = "Coverage Probability (alpha = 4, No noise, SNR = 10)", col = "red")
    lines(decibel, y2, type = "l", col = "blue")
    lines(decibel, y3, type = "l", col = "green")
    
}
