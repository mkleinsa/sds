
resamp = function(x, n) {
    sq = 1:nrow(x)
    sq = sq / max(sq)
    sqz = rep(TRUE, length(sq))
    for(i in 1:nrow(x)) {
        sqz[i] = rbernoulli(1, p = sq[i])
    }
    y = x[sqz, ]
    y[["time"]] = (max(y[["time"]]) - nrow(y) + 1):(max(y[["time"]]))
    # for(i in 1:nrow(y)) {
    #     nrp = sample(0:5, size = 1, replace = FALSE)
    #     if(nrp != 0) {
    #         y = rbind(y, y[i,][rep(1, nrp), ])
    #     }
    # }
    # y = y %>% arrange(day)
    return(list(y, sqz))
}