#' Susceptible-(Exposed)-Infectious-Recovered
#' 
#' Fit a S(E)IR (SIR or SEIR) model or give a custom ODE model.
#'
#' @examples
#' 
#' # Data set
#' data(indiaAprilJune)
#' N = 1366e6
#' 
#' # Basic SIR
#' fit = seir(data = indiaAprilJune, N = N, model = "SIR", start = list(beta = log(1e-2), gamma = log(1e-5)))
#' fit$coef
#' coefs = fit$coef
#' R0 = coefs[["beta"]] / coefs[["gamma"]]
#' R0
#' plot.seirMod(fit)
#' 
#' # Basic SIR with demography (fixed birth and death rate)
#' lambda = mu = 1 / (69.416 * 365)
#' fit = seir(data = indiaAprilJune, N = N, model = "SIR-Demography", start = list(beta = log(1e-2), gamma = log(1e-5)),
#'            fixed_data = list(lambda = lambda, mu = mu))
#' fit$coef
#' coefs = fit$coef
#' R0 = coefs[["beta"]] / coefs[["gamma"]]
#' R0
#' plot.seirMod(fit)
#' 
#' # Basic SEIR
#' fit = seir(data = indiaAprilJune, N = N, model = "SEIR", start = list(beta = log(1e-2), gamma = log(1e-5)),
#'            fixed_data = list(De = 5.2))
#' fit$coef
#' # R0 = ???
#' # R0
#' plot.seirMod(fit)
#' 
#' # Basic SEIR with demography
#' lambda = mu = 1 / (69.416 * 365)
#' fit = seir(data = indiaAprilJune, N = N, model = "SEIR-Demography", start = list(beta = log(1e-2), gamma = log(1e-5)),
#'            fixed_data = list(De = 5.2, lambda = lambda, mu = mu))
#' fit$coef
#' # R0 = ???
#' # R0
#' plot.seirMod(fit)
#' 
#' # SEIR with inverse duration modeled
#' # Estimate beta, gamma, and De (inverse duration of infection in SEIR)
#' fit = seir(data = indiaAprilJune, N = N, model = "SEIR-De", start = list(beta = log(1e-2), gamma = log(1e-5), De = 5.2))
#' 
#' plot.seirMod(fit)
#' 
#' # SEIR with inverse duration modeled and demography (birts and deaths)
#' # Estimate beta, gamma, and De (inverse duration of infection in SEIR)
#' data(indiaAprilJuly)
#' fit = seir(data = indiaAprilJuly, N = N, model = "SEIR-De-Demography", start = list(beta = log(1e-2), gamma = log(1e-5), De = 5.2),
#'            fixed_data = list(lambda = lambda, mu = mu))
#' 
#' plot.seirMod(fit)
#' 
#' # Plots
#' plot.seirMod(fit)
#' p = plot.seirMod(fit, CI = TRUE)
#' library(ggplot2)
#' p + xlim(0, 63) + ylim(0, 1.5e5)
#' p2 = plot.seirMod(fit, CI = TRUE, forecast = 20)
#' p2 + xlim(0, 75) + ylim(0, 4e5)
#' 
#' # Estimate beta, gamma, and De (inverse duration of infection in SEIR)
#' fit = seir(data = indiaAprilJune, N = N, model = "SEIR-De")
#' 
#' 
#' @importFrom deSolve ode
#' @importFrom bbmle mle2 coef
#' 
seir = function(data, model = c("SIR", "SEIR", "custom", "SEIR-De", "SIR-Demography", "SEIR-Demography", "SEIR-De-Demography"), 
        N, start = NULL, fixed_data = NULL) {
    if(model == "SIR") {
        if(is.null(start)) {
            stop("In SIR, starting values must be provided for beta and gamma, e.g. start = list(beta = log(1e-2), gamma = log(1e-5)).")
        } else if(!all(c("beta", "gamma") %in% names(start))) {
            stop("Start requires a named list of starting values for beta and gamma, e.g. start = list(beta = log(1e-2), gamma = log(1e-5)).")
        }
        logli = logli_sir_1
        # start_default = list(beta = log(1e-2), gamma = log(1e-5))
    } else if(model == "SEIR") {
        # start_default = list(beta = log(1e-2), gamma = log(1e-5))
        if(is.null(fixed_data)) {
            stop("In SEIR with fixed inverse duration, De, you must specify the argument fixed_data = list(De = X), where X is an 
                  appropriate De for the problem. To treat De as a parameter, use model = 'SEIR-De' and give De starting values.")
            # fixed_data = list(De = 5.2)
        }
        if(is.null(start)) {
            stop("In SEIR, starting values must be provided for beta and gamma, e.g. start = list(beta = log(1e-2), gamma = log(1e-5)).")
        } else if(!all(c("beta", "gamma") %in% names(start))) {
            stop("Start requires a named list of starting values for beta and gamma, e.g. start = list(beta = log(1e-2), gamma = log(1e-5)).")
        }
        logli = logli_seir_1
    } else if(model == "SEIR-De") {
        if(is.null(start)) {
            stop("In SEIR-De, starting values must be provided for beta, gamma, and De, e.g. start = list(beta = log(1e-2), gamma = log(1e-5), De = 5.2).")
        } else if(!all(c("beta", "gamma", "De") %in% names(start))) {
            stop("Start requires a named list of starting values for beta, gamma, and De, e.g. start = list(beta = log(1e-2), gamma = log(1e-5), De = 5.2).")
        }
        logli = logli_seir_1_de
        # start_default = list(beta = log(1e-2), gamma = log(1e-5), De = 5.2)
    } else if(model == "SIR-Demography") {
        if(is.null(start)) {
            stop("In SIR, starting values must be provided for beta and gamma, e.g. start = list(beta = log(1e-2), gamma = log(1e-5)).")
        } else if(!all(c("beta", "gamma") %in% names(start))) {
            stop("Start requires a named list of starting values for beta and gamma, e.g. start = list(beta = log(1e-2), gamma = log(1e-5)).")
        }
        if(is.null(fixed_data)) {
            stop("In SIR-Demography, fixed_data for lambda, mu must be specified, e.g. ")
            # fixed_data = list(De = 5.2)
        } else if(!all(c("lambda", "mu") %in% names(fixed_data))) {
            stop("In SIR-Demography, fixed_data for both lambda, mu must be specified, e.g. ")
        }
        logli = logli_sir_1_demography
        # start_default = list(beta = log(1e-2), gamma = log(1e-5))
    } else if(model == "SEIR-Demography") {
        if(is.null(fixed_data)) {
            stop("In SEIR-Demography with fixed inverse duration, De, you must specify the argument fixed_data = list(De = X. lambda = , mu = ).
                  To treat De as a parameter, use model = 'SEIR-De-Demography' and give De starting values.")
            # fixed_data = list(De = 5.2)
        }
        if(is.null(start)) {
            stop("In SEIR-Demography, starting values must be provided for beta and gamma, e.g. start = list(beta = log(1e-2), gamma = log(1e-5)).")
        } else if(!all(c("beta", "gamma") %in% names(start))) {
            stop("Start requires a named list of starting values for beta and gamma, e.g. start = list(beta = log(1e-2), gamma = log(1e-5)).")
        }
        logli = logli_seir_1_demography
        # start_default = list(beta = log(1e-2), gamma = log(1e-5))
    } else if(model == "SEIR-De-Demography") {
        logli = logli_seir_1_de_demography
        # start_default = list(beta = log(1e-2), gamma = log(1e-5), De = 5.2)
    } else if(model == "custom") {
        sirmod = custom_s(model_details)
    } else {
        stop("Not a model. Use one of SIR, SEIR, custom for a custom model.")
    }

    data[["time"]] = 1:nrow(data)

    if(model == "SIR") {
        data = list(dat = data, N = N)
    } else if(model == "SEIR") {
        data = list(dat = data, N = N, De = fixed_data[["De"]])
    } else if(model == "SEIR-De") {
        data = list(dat = data, N = N)
    } else if(model == "SIR-Demography") {
        data = list(dat = data, N = N, lambda = fixed_data[["lambda"]], mu = fixed_data[["mu"]])
    } else if(model == "SEIR-Demography") {
        data = list(dat = data, N = N, De = fixed_data[["De"]], lambda = fixed_data[["lambda"]], mu = fixed_data[["mu"]])
    } else if(model == "SEIR-De-Demography") {
        data = list(dat = data, N = N, lambda = fixed_data[["lambda"]], mu = fixed_data[["mu"]])
    }
    
    estimates = mle2(minuslogl = logli, start = start, method = "Nelder-Mead", data = data)

    if(model == "SIR") {
        coef = exp(bbmle::coef(estimates))
    } else if(model == "SEIR") {
        coef = exp(bbmle::coef(estimates))
    } else if(model == "SEIR-De") {
        coef = bbmle::coef(estimates)
        coef = c(exp(coef[c("beta", "gamma")]), coef["De"])
    } else if(model == "SIR-Demography") {
        coef = exp(bbmle::coef(estimates))
    } else if(model == "SEIR-Demography") {
        coef = exp(bbmle::coef(estimates))
    } else if(model == "SEIR-De-Demography") {
        coef = bbmle::coef(estimates)
        coef = c(exp(coef[c("beta", "gamma")]), coef["De"])
    }

    result = list(data = data, coef = coef, N = N, model = model, logli = logli, start = start, fixed_data = fixed_data)
    class(result) = "seirMod"
    return(result)
}