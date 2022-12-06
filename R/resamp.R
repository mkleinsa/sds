
resamp = function(x, n) {
    sq = 1:nrow(x)
    sq = rep(0.5, nrow(x)) #sq / max(sq)
    sqz = rep(FALSE, length(sq))
    while(sum(sqz) < n) {
        pos = floor(runif(1, 0, 1)*nrow(x)) + 1
        sqz[pos] = rbernoulli(1, p = sq[pos])
        if(sum(sqz) == n) {
            break;
        }
    }
    # for(i in 1:nrow(x)) {
    #     pos = floor(runif(1, 0, 1)*nrow(x)) + 1
    #     sqz[pos] = rbernoulli(1, p = sq[pos])
    #     if(sum(sqz) == n) {
    #         break;
    #     }
    # }
    y = x[sqz, ]
    # y[["time"]] = (max(y[["time"]]) - nrow(y) + 1):(max(y[["time"]]))
    y2 = y
    # for(i in 1:nrow(y)) {
    #     nrp = sample(0:5, size = 1, replace = FALSE)
    #     if(nrp != 0) {
    #         y2 = rbind(y2, y[i,][rep(1, nrp), ])
    #     } else {
    #         y2 = y[-i, ]
    #     }
    # }
    y2 = y2 %>% arrange(time)
    return(list(y = y2, sqz = sqz))
}