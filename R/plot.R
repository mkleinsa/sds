
#' Plot a SEIR model
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
plot.seirMod = function(fit, CI = FALSE, forecast = 0) {

    N = fit[["N"]]
    mle2_data = fit[["data"]]
    data = fit[["data"]][["dat"]]
    model = fit[["model"]]
    logli = fit[["logli"]]
    start = fit[["start"]]
    fixed_data = fit[["fixed_data"]]

    # plot result
    time_points = seq(min(data$date), max(data$date), le = 100)
    I0 = data[["I"]][1] # initial number of infected
    R0 = data[["R"]][1] # initial number of removed
    param_hat = fit[["coef"]]
    # model's best predictions: 

    plot_times = 1:(nrow(data) + forecast)

    if(model == "SIR") {
        best_predictions = sir_1(beta = param_hat["beta"], gamma = param_hat["gamma"], I0 = I0, R0 = R0, 
            times = plot_times, N = N)$I
    } else {
        best_predictions = seir_1(beta = param_hat["beta"], gamma = param_hat["gamma"], I0 = I0, R0 = R0, 
            times = plot_times, N = N, De = fixed_data$De)$I
    }

    p = ggplot(data = data.frame(obs = c(data[["I"]], rep(0, forecast)), pred = best_predictions, time = plot_times)) + 
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
        # resamp_stretch = rexp(250, rate = 1)
        resamp_stretch = rnorm(250, mean = 1, sd = 0.25)
        resamp_jitter = rnorm(250, mean = 1, sd = 0.25)
        #resamp_stretch = (resamp_stretch - min(resamp_stretch)) / (max(resamp_stretch) - min(resamp_stretch))
        resamp_stretch = (resamp_stretch - 0.25) / (1.75 - 0.25)
        resamp_stretch = resamp_stretch - mean(resamp_stretch) + 1.0
        #resamp_jitter = (resamp_jitter - 0.5) / (1.5 - 0.5)
        resamp_jitter = resamp_jitter - mean(resamp_jitter)
        #print(mean(resamp_stretch))
        # resamp_stretch = resamp_stretch + 1
        for(i in 1:nrow(resamp_est)) {
            #dat = resamp(data, 15)
            dat$y = data
            dat[["y"]][["I"]] = round(((dat[["y"]][["I"]] + resamp_jitter[i])*resamp_stretch[i]), digits = 0) #* resamp_stretch[i]
            resamp_day[i, ] = dat[["y"]][["time"]][1]
            resamp_dat[i, ] = as.numeric(dat[["y"]][1, 2:3])
            estimates = mle2(minuslogl = logli, start = start, method = "Nelder-Mead", 
                data = mle2_data)
            resamp_est[i,] = exp(coef(estimates))
        }

        best_preds = vector(mode = "list", length = 250)
        for(i in 1:250) {
            if(model == "SIR") {
                best_preds[[i]] = sir_1(beta = resamp_est[i, 1], gamma = resamp_est[i, 2], I0 = resamp_dat[i, 1], R0 = resamp_dat[i, 2], 
                times = plot_times, N = N)$I
            } else if (model == "SEIR"){
                # best_preds[[i]] = seir_1(beta = resamp_est[i, 1], gamma = resamp_est[i, 2], I0 = resamp_dat[i, 1], R0 = resamp_dat[i, 2], 
                #     times = 1:(floor(1.2*nrow(data))), N = N, De = fixed_data$De)$I
                best_preds[[i]] = seir_1(beta = resamp_est[i, 1], gamma = resamp_est[i, 2], I0 = resamp_dat[i, 1], R0 = resamp_dat[i, 2], 
                    times = plot_times, N = N, De = fixed_data$De)$I
            }
        }
        for(i in 1:250) {
            prds = best_preds[[i]]
            print(i)
            print(any(is.na(prds)))
            pdat = cbind.data.frame(time = plot_times, pred = prds)
            # pdat = cbind.data.frame(time = dat[["y"]][["time"]], pred = prds)
            p = p + geom_line(data = pdat, aes(x = time, y = pred), color = "#FFCB05", linewidth = 0.4, alpha = 0.08)
        }
    }

    return(p)
}