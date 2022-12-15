#' Visualize S(E)IR Models for Disease Dynamics
#' 
#' Visualize S(E)IR Models for Disease Dynamics
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
#' 
#'
#'
#'
#' 
#' @examples 
#' vis_seir(model = "SIR", beta = 2, gamma = 1, I0 = 5, R0 = 2, time = 1:70, N = 50000)
#'
#' @importFrom ggplot2 scale_colour_manual
vis_seir = function(model = "SIR", beta, gamma, I0, R0, times, N) {

    # beta = 2
    # gamma = 0.05
    # I0 = 5
    # R0 = 2
    # N = 5e4
    # times = 1:nrow(india)

    predictions = sir_1(beta = beta, gamma = gamma, I0 = I0, R0 = R0, times = times, N = N)

    p = ggplot(predictions) + 
        aes(x = time, y = S, color = "S") + 
        geom_line(size = 1.2) + 
        geom_line(aes(y = I, color = "I"), size = 1.2) + 
        geom_line(aes(y = R, color = "R"), size = 1.2) + 
        theme_classic() + 
        scale_colour_manual("", breaks = c("S", "I", "R"),
            values = c("S"="black", "I"="red", "R"="green")) + 
        xlab("Time (days)") + 
        ylab("Count") + 
        theme(
            axis.title.x = element_text(size = 18),
            axis.title.y = element_text(size = 18),
            axis.text = element_text(size = 18),
            legend.text = element_text(size = 18)
        )

    return(p)

}