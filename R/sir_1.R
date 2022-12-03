# SIR model via mle
sir_1 = function(beta, gamma, I0, R0, times, N) {
  sir_equations = function(time, variables, parameters) {
    with(as.list(c(variables, parameters)), {
      dS = -beta * I * S/N
      dI = beta * I * S/N - gamma * I
      dR = gamma * I
      return(list(c(dS, dI, dR)))
    })
  }
  parameters_values = c(beta = beta, gamma = gamma)
  S0 = N - I0 - R0
  initial_values = c(S = S0, I = I0, R = R0)
  out = ode(initial_values, times, sir_equations, parameters_values, method = "rk4")
  as.data.frame(out)
}

logli_sir_1 = function(beta, gamma, N, dat) {
  I0 = dat$I[1]
  R0 = dat$R[1]
  times = dat[["time"]]
  beta = exp(beta)
  gamma = exp(gamma)
  predictions = sir_1(beta = beta, gamma = gamma, I0 = I0, R0 = R0, 
        times = times, N = N)
  ## negative of log likelihood
  -sum(dpois(x = dat$I, lambda = predictions$I, log = TRUE))
}