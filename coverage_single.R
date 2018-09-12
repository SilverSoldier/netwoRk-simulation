innerFuncRho = function(u, alpha) {
    1/(1 + `^`(u, alpha/2))
}

rho = function(T, alpha) {
    res <- `^`(T, 2/alpha)
    lowLim <- 1/res
    res <- res * integrate(innerFuncRho, lowLim, Inf, alpha = alpha)$value
    return(res)
}

innerFuncCoverage = function(v, T, L, alpha, SNR){
    exp(-1 * pi * L * v * (1 + rho(T, alpha)) - T*(`^`(v, alpha/2))/SNR)
}

coverage = function(T, L, alpha, SNR) {
    pi * L * integrate(innerFuncCoverage, lower = 0, upper = Inf, T = T, L = L, alpha = alpha, SNR = SNR)$value
}
    
decibelToWatt = function(D) {
    `^`(10, D/10)
}

plotCoverage = function (A, L, SNR) {
    decibel = seq(-10, 20, 0.1)
    watts = sapply(decibel, decibelToWatt)
    y1 = sapply(watts, coverage, L =  4, SNR = SNR, alpha = A)
    header = paste("Coverage Probability (alpha, SNR, Lambda) = (", A, ",", SNR, ", ", L, ")", sep = " ")
    y2 = sapply(watts, coverage, L =  7, SNR = SNR, alpha = A)
    plot(decibel, y1, type = "l", xlim = c(-10, 20), ylim = c(0.1, 1), xlab = "SINR Threshold (dB)", ylab = "Probability of Coverage", main = header, col = "blue")
    lines(decibel, y2, col = "red")
}
