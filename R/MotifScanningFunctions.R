.scan.sequence.with.pwm <- function(pwm, seq, minScore){
    
    nc <- nchar(minScore)
    if (substr(minScore, nc, nc) == "%"){
        perc.threshold <- substr(minScore, 1L, nc-1L)
        min.score <- minScore(pwm)
        max.score <- maxScore(pwm)
        score.threshold = min.score + as.double(perc.threshold)/100 *
            (max.score-min.score)
    }else{
        score.threshold <- minScore
    }

    pwm.match <- matchPWM(pwm = pwm, subject = seq, min.score = score.threshold)
    return(start(pwm.match))

}

.get.scanning.score <- function(pwm, seq){
    PWMscoreStartingAt(pwm, subject = seq,
        starting.at = c(1:(length(seq) - ncol(pwm) + 1)))
}

