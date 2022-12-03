#' Susceptible-(Exposed)-Infectious-Recovered
#' 
#' Fit a S(E)IR model or give a ODE model.
#' 
#'
#' 
#'
#' 
#'
#' 
#'
#' 
#'
#' @examples
#' data(indiaAprilJune)
#' N = 1366e6
#' fit = seir(data = indiaAprilJune, N = N, model = "SIR")
#' fit$coef
#' plot.seirMod(fit)
#' plot.seirMod(fit, CI = TRUE)
#' 
#' @importFrom deSolve ode
#' @importFrom bbmle mle2 coef
#' 
seir = function(data, model = c("SIR", "SEIR", "custom"), N, start = NULL) {
    if(model == "SIR") {
        logli = logli_sir_1
        start_default = list(beta = log(1e-2), gamma = log(1e-5))
    } else if(model == "SEIR") {
        sirmod = seir_1
    } else if(model == "custom") {
        sirmod = custom_s(model_details)
    } else {
        stop("Not a model. Use one of SIR, SEIR, custom for a custom model.")
    }

    data[["time"]] = 1:nrow(data)

    if(is.null(start)) {
        start = start_default
    }

    estimates = mle2(minuslogl = logli, start = start, method = "Nelder-Mead", 
        data = list(dat = data, N = N))

    result = list(data = data, coef = exp(bbmle::coef(estimates)), N = N, model = model,
        logli = logli, start = start)
    class(result) = "seirMod"

    return(result)
}