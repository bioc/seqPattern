.scan.sequence.with.pwm <- function(pwm, seq, minScore){
    pwm.match <- matchPWM(pwm = pwm, subject = seq, min.score = minScore)
    return(start(pwm.match))
}

.get.scanning.score <- function(pwm, seq){
    PWMscoreStartingAt(pwm, subject = seq,
        starting.at = c(1:(length(seq) - ncol(pwm) + 1)))
}

