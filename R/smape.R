#' Calculate SMAPE (symmetric? mean prediction error)
#' 
#' Used in the plot function to simulate models for confidence intervals
#' 
#'
#' 
#'
#' 
#'
#' 
#'
#' 
smape = function(pred, obs, n) {
    model_smape = (1 / n) * (sum(abs(pred - obs) / ((abs(pred) + abs(obs)) / 2)))
    return(model_smape)
}