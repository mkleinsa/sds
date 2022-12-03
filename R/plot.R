
#' Plot a seir model
#' 
#'
#' 
#'
#' 
#'
#' 
#' @importFrom ggplot2 ggplot geom_col geom_line theme_classic theme aes element_text
#' @importFrom purrr rbernoulli
#' @export 
plot.seirMod = function(fit, CI = FALSE) {

    N = fit[["N"]]
    data = fit[["data"]]
    model = fit[["model"]]
    logli = fit[["logli"]]
    start = fit[["start"]]

    # plot result
    time_points = seq(min(data$date), max(data$date), le = 100)
    I0 = data[["I"]][1] # initial number of infected
    R0 = data[["R"]][1] # initial number of removed
    param_hat = fit[["coef"]]
    # model's best predictions:
    best_predictions = sir_1(beta = param_hat["beta"], gamma = param_hat["gamma"], I0 = I0, R0 = R0, 
            times = 1:nrow(data), N = N)$I

    p = ggplot(data = data.frame(obs = data[["I"]], pred = best_predictions, time = data[["time"]])) + 
        geom_col(aes(x = time, y = obs), fill = "#00274C") + 
        geom_line(aes(x = time, y = pred), color = "#ff0505", linewidth = 1.4) + 
        theme_classic() +
        theme(
            axis.title.x = element_text(size = 18),
            axis.title.y = element_text(size = 18),
            axis.text = element_text(size = 18),
            legend.text = element_text(size = 18))

    if(CI == TRUE) {
        resamp_est = matrix(nrow = 250, ncol = 2)
        resamp_dat = matrix(nrow = 250, ncol = 2)
        resamp_day = matrix(nrow = 250, ncol = 1)
        resamp_stretch = rexp(250, rate = 1)
        resamp_stretch = (resamp_stretch - min(resamp_stretch)) / (max(resamp_stretch) - min(resamp_stretch))
        resamp_stretch = resamp_stretch + 1
        for(i in 1:nrow(resamp_est)) {
            dat = resamp(data, 12)
            dat$I = round(dat$I * (resamp_stretch[i]), digits = 0)
            resamp_day[i, ] = dat[["time"]][1]
            resamp_dat[i, ] = as.numeric(dat[1, 2:3])
            estimates = mle2(minuslogl = logli, start = start, method = "Nelder-Mead", 
                data = list(dat = dat, N = N))
            resamp_est[i,] = exp(coef(estimates))
        }

        best_preds = vector(mode = "list", length = 250)
        for(i in 1:250) {
            best_preds[[i]] = sir_1(beta = resamp_est[i, 1], gamma = resamp_est[i, 2], I0 = resamp_dat[i, 1], R0 = resamp_dat[i, 2], 
                times = 1:(floor(1.2*nrow(data))), N = N)$I
        }
        for(i in 1:250) {
            prds = best_preds[[i]]
            dat = cbind.data.frame(time = resamp_day[i,1]:(resamp_day[i,1] + length(prds) - 1), pred = prds)
            p = p + geom_line(data = dat, aes(x = time, y = pred), color = "#FFCB05", linewidth = 0.4, alpha = 0.08)
        }
    }

    return(p)
}