
#' Plot a S(E)IR model
#' 
#' Plot observed data, best fit, forecasts, and confidence intervals via simulation.
#' 
#' @param fit seirMod model fit object. Result of calling seir() from sds package.
#' @param CI TRUE/FALSE simulate confidence bands by resampling procedure?
#' @param forecast number of days to forecast (i.e. make predictions). Defaults to 0, or no prediction.
#' @param nsim number of simulations to run if CI = TRUE. Default = 250.
#' @param mse_interval TRUE/FALSE, use smape (estimated model error) (if CI = TRUE) as s.d. of normal distribution in resampling simulation?
#' @param mse user specified error of simulation (s.d. of normal distribution)
#' 
#'
#' @return a ggplot object. Object can be stored and further updated just like a ggplot. See examples from seir() function.
#'
#' 
#' @importFrom ggplot2 ggplot geom_col geom_line theme_classic theme aes element_text scale_y_continuous ylab xlab
#' @importFrom scales number_format cut_short_scale
#' @importFrom purrr rbernoulli
#' @export 
plot.seirMod = function(fit, CI = FALSE, forecast = 0, nsim = 250, mse_interval = TRUE, mse = NULL) {

    # CI = TRUE
    # forecast = 5
    # nsim = 250
    # mse_interval = TRUE 
    # mse = NULL

    N          = fit[["N"]]
    mle2_data  = fit[["data"]]
    data       = fit[["data"]][["dat"]]
    model      = fit[["model"]]
    logli      = fit[["logli"]]
    start      = fit[["start"]]
    fixed_data = fit[["fixed_data"]]

    # plot result
    time_points = seq(min(data$date), max(data$date), le = 100)
    I0 = data[["I"]][1] # initial number of infected
    R0 = data[["R"]][1] # initial number of removed
    param_hat = fit[["coef"]]
    # param_hat[["beta"]] = log(param_hat[["beta"]])
    # param_hat[["gamma"]] = log(param_hat[["gamma"]])
    # model's best predictions: 

    plot_times = 1:(nrow(data) + forecast)

    if(model == "SIR") {
        best_predictions = sir_1(beta = param_hat["beta"], gamma = param_hat["gamma"], I0 = I0, R0 = R0, 
            times = plot_times, N = N)[["I"]]
    } else if(model == "SEIR") {
        best_predictions = seir_1(beta = param_hat["beta"], gamma = param_hat["gamma"], I0 = I0, R0 = R0, 
            times = plot_times, N = N, De = fixed_data[["De"]])[["I"]]
    } else if(model == "SEIR-De") {
        best_predictions = seir_1(beta = param_hat["beta"], gamma = param_hat["gamma"], I0 = I0, R0 = R0, 
            times = plot_times, N = N, De = param_hat["De"])[["I"]]
    } else if(model == "SIR-Demography") {
        best_predictions = sir_1_demography(beta = param_hat["beta"], gamma = param_hat["gamma"], I0 = I0, R0 = R0, 
            times = plot_times, N = N, lambda = fixed_data[["lambda"]], mu = fixed_data[["mu"]])[["I"]]
    } else if(model == "SEIR-Demography") {
        best_predictions = seir_1_demography(beta = param_hat["beta"], gamma = param_hat["gamma"], I0 = I0, R0 = R0, 
            times = plot_times, N = N, De = fixed_data[["De"]], lambda = fixed_data[["lambda"]], mu = fixed_data[["mu"]])[["I"]]
    } else if(model == "SEIR-De-Demography") {
        best_predictions = seir_1_de_demography(beta = param_hat["beta"], gamma = param_hat["gamma"], I0 = I0, R0 = R0, 
            times = plot_times, N = N, De = param_hat["De"], lambda = fixed_data[["lambda"]], mu = fixed_data[["mu"]])[["I"]]
    }

    resids_fit = best_predictions[1:nrow(data)] - data[["I"]]
    resids_fit = resids_fit - mean(resids_fit)
    # (resids_fit - mean(resids_fit)) / sd(resids_fit)
    # hist(resids_fit)
    # hist((resids_fit - mean(resids_fit)) / sd(resids_fit))
    model_rmse = sqrt(sum(( (best_predictions[1:nrow(data)] - data[["I"]]) /  1.0 )^2) / nrow(data)) #(data[["I"]])
    # model_rmse = (1 / nrow(data)) * ( sum( abs(best_predictions[1:nrow(data)] - data[["I"]]) / ((abs(best_predictions[1:nrow(data)]) + abs(data[["I"]])) / 2)  ) )
    # model_error = smape(pred = best_predictions[1:nrow(data)], obs = data[["I"]], n = nrow(data))

    p = ggplot(data = data.frame(obs = c(data[["I"]], rep(0, forecast)), pred = best_predictions, time = plot_times)) + 
        geom_col(aes(x = time, y = obs), fill = "#00274C") + 
        geom_line(aes(x = time, y = pred), color = "#ff0505", linewidth = 1.4) + 
        theme_classic() +
        theme(
            axis.title.x = element_text(size = 18),
            axis.title.y = element_text(size = 18),
            axis.text = element_text(size = 18),
            legend.text = element_text(size = 18)) + 
        scale_y_continuous(labels = scales::number_format(scale = 1, trim = TRUE, big.mark = ",", scale_cut = scales::cut_short_scale()))

    if(CI == TRUE) {
        if(model %in% c("SIR", "SEIR", "SIR-Demography", "SEIR-Demography")) {
            resamp_est = matrix(nrow = nsim, ncol = 2)
        } else {
            resamp_est = matrix(nrow = nsim, ncol = 3)
        }
        resamp_dat = matrix(nrow = nsim, ncol = 2)
        resamp_day = matrix(nrow = nsim, ncol = 1)
        # resamp_stretch = rexp(250, rate = 1)
        # resamp_stretch = rnorm(nsim, mean = 1, sd = model_error)
        resamp_jitter = rnorm(nsim, mean = 1, sd = median(model_rmse / (data[["I"]])))
        #resamp_stretch = (resamp_stretch - min(resamp_stretch)) / (max(resamp_stretch) - min(resamp_stretch))
        #resamp_stretch = (resamp_stretch - 0.25) / (1.75 - 0.25)
        # resamp_stretch = resamp_stretch - mean(resamp_stretch) + 1.0
        #resamp_jitter = (resamp_jitter - 0.5) / (1.5 - 0.5)
        resamp_jitter = resamp_jitter - mean(resamp_jitter) + 1.0
        #print(mean(resamp_stretch))
        # resamp_stretch = resamp_stretch + 1
        for(i in 1:nrow(resamp_est)) {
            resid_samp_pos = sample(1:nrow(data), size = 1)
            resid_samp_prop = resids_fit[resid_samp_pos] / data[["I"]][resid_samp_pos]
            dat = data
            # dat[["I"]] = round(((dat[["I"]] + (resamp_jitter[i]/(dat[["I"]])))), digits = 0) # *resamp_stretch[i] + resamp_jitter[i]
            
            # recent:
            #dat[["I"]] = round(((dat[["I"]] * resamp_jitter[i])), digits = 0)
            dat[["I"]] = round(dat[["I"]] + resid_samp_prop*dat[["I"]], digits = 0)
#sample(c(-1, 1), size = 1)*

            #dat[["I"]] = round(((dat[["I"]] + resamp_jitter[i])), digits = 0)
            resamp_day[i, ] = dat[["time"]][1]
            resamp_dat[i, ] = as.numeric(dat[1, 2:3])
            # estimates = mle2(minuslogl = logli, start = start, method = "Nelder-Mead", data = mle2_data)
            estimates = seir(data = dat, N = N, model = model, start = start, fixed_data = fixed_data)
            ### data = data???? really???
            resamp_est[i, ] = estimates[["coef"]]
        }

        best_preds = vector(mode = "list", length = nsim)
        for(i in 1:nsim) {
            if(model == "SIR") {
                best_preds[[i]] = sir_1(beta = resamp_est[i, 1], gamma = resamp_est[i, 2], I0 = resamp_dat[i, 1], R0 = resamp_dat[i, 2], 
                    times = plot_times, N = N)[["I"]]
            } else if(model == "SEIR") {
                best_preds[[i]] = seir_1(beta = resamp_est[i, 1], gamma = resamp_est[i, 2], I0 = resamp_dat[i, 1], R0 = resamp_dat[i, 2], 
                    times = plot_times, N = N, De = fixed_data[["de"]])[["I"]]
            } else if(model == "SEIR-De") {
                best_preds[[i]] = seir_1_de(beta = resamp_est[i, 1], gamma = resamp_est[i, 2], I0 = resamp_dat[i, 1], R0 = resamp_dat[i, 2], 
                    times = plot_times, N = N, De = esamp_est[i, 3])[["I"]]
            } else if(model == "SIR-Demography") {
                best_preds[[i]] = sir_1_demography(beta = resamp_est[i, 1], gamma = resamp_est[i, 2], I0 = resamp_dat[i, 1], R0 = resamp_dat[i, 2], 
                    times = plot_times, N = N, lambda = fixed_data[["lambda"]], mu = fixed_data[["mu"]])[["I"]]
            } else if(model == "SEIR-Demography") {
                best_preds[[i]] = seir_1_demography(beta = resamp_est[i, 1], gamma = resamp_est[i, 2], I0 = resamp_dat[i, 1], R0 = resamp_dat[i, 2], 
                    times = plot_times, N = N, De = fixed_data[["De"]], lambda = fixed_data[["lambda"]], mu = fixed_data[["mu"]])[["I"]]
            } else if(model == "SEIR-De-Demography") {
                best_preds[[i]] = seir_1_de_demography(beta = resamp_est[i, 1], gamma = resamp_est[i, 2], I0 = resamp_dat[i, 1], R0 = resamp_dat[i, 2], 
                    times = plot_times, N = N, De = resamp_est[i, 3], lambda = fixed_data[["lambda"]], mu = fixed_data[["mu"]])[["I"]]
            }
        }

        q_preds = vector(mode = "numeric", length = nsim)
        for(i in 1:nsim) {
            b_preds = best_preds[[i]]
            lg = length(b_preds)
            q_preds[i] = b_preds[lg]
        }

        q_tiles = quantile(q_preds, p = c(0.025, 0.975))
        print(hist(q_preds))
        lower = which.min(abs(q_tiles[1] - q_preds))
        upper = which.min(abs(q_tiles[2] - q_preds))

        for(i in 1:nsim) {
            prds = best_preds[[i]]
            pdat = cbind.data.frame(time = plot_times, pred = prds)
            p = p + geom_line(data = pdat, aes(x = time, y = pred), color = "#FFCB05", linewidth = 0.5, alpha = 0.08)
        }

        # Re-plot the primary fit above the uncertainty interval layers
        p = p + 
            geom_line(data = data.frame(obs = c(data[["I"]], rep(0, forecast)), pred = best_predictions, time = plot_times), 
                      aes(x = time, y = pred), color = "#ff0505", linewidth = 1.4) + 
            geom_line(data = data.frame(time = plot_times, pred = best_preds[[lower]]), aes(x = time, y = pred), 
                      color = "#ff0505", linewidth = 1.4, linetype = 2) + 
            geom_line(data = data.frame(time = plot_times, pred = best_preds[[upper]]), aes(x = time, y = pred), 
                      color = "#ff0505", linewidth = 1.4, linetype = 2)
    }

    p = p + 
        ylab("Active Cases") + 
        xlab("Time (days)")

    return(p)
}