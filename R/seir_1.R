# SEIR model via mle
seir_1 = function(beta, gamma, I0, R0, times, N, De) {
  sir_equations = function(time, variables, parameters) {
    with(as.list(c(variables, parameters)), {
      dS = -(beta * I * (S/N))
      dE =  (beta * I * (S/N)) - ((1/De) * E)
      dI =  ((1/De) * E) - (gamma * I)
      dR =  (gamma * I)
      return(list(c(dS, dE, dI, dR)))
    })
  }
  parameters_values = c(beta = beta, gamma = gamma)
  S0 = N - I0 - R0
  E0 = 3 * I0
  initial_values = c(S = S0, E = E0, I = I0, R = R0)
  out = ode(initial_values, times, sir_equations, parameters_values, method = "rk4")
  as.data.frame(out)
}

logli_seir_1 = function(beta, gamma, N, dat, De) {
  I0 = dat$I[1]
  R0 = dat$R[1]
  times = dat[["time"]]
  beta = exp(beta)
  gamma = exp(gamma)
  predictions = seir_1(beta = beta, gamma = gamma, I0 = I0, R0 = R0, 
        times = times, N = N, De = De)
  ## negative of log likelihood
  -sum(dpois(x = dat$I, lambda = predictions$I, log = TRUE))
}

seir_1_demography = function(beta, gamma, I0, R0, times, N, lambda, mu, De) {
  sir_equations = function(time, variables, parameters) {
    with(as.list(c(variables, parameters)), {
      dS = -(beta * I * (S/N)) + (lambda * N) - (mu * S)
      dE =  (beta * I * (S/N)) - ((1/De) * E) - (mu * E)
      dI =  ((1/De) * E) - (gamma * I) - (mu * I)
      dR =  (gamma * I) - (mu * R)
      return(list(c(dS, dE, dI, dR)))
    })
  }
  parameters_values = c(beta = beta, gamma = gamma)
  S0 = N - I0 - R0
  E0 = 3 * I0
  initial_values = c(S = S0, E = E0, I = I0, R = R0)
  out = ode(initial_values, times, sir_equations, parameters_values, method = "rk4")
  as.data.frame(out)
}

logli_seir_1_demography = function(beta, gamma, N, dat, lambda, mu, De) {
  I0 = dat$I[1]
  R0 = dat$R[1]
  times = dat[["time"]]
  beta = exp(beta)
  gamma = exp(gamma)
  predictions = seir_1_demography(beta = beta, gamma = gamma, I0 = I0, R0 = R0, 
        times = times, N = N, lambda = lambda, mu = mu, De = De)
  ## negative of log likelihood
  -sum(dpois(x = dat$I, lambda = predictions$I, log = TRUE))
}

# seir_1_de_demography = function(beta, gamma, I0, R0, times, N, lambda, mu, De) {
#   sir_equations = function(time, variables, parameters) {
#     with(as.list(c(variables, parameters)), {
#       dS = -(beta * I * (S/N)) + (lambda * N) - (mu * S)
#       dE =  (beta * I * (S/N)) - ((1/De) * E) - (mu * E)
#       dI =  ((1/De) * E) - (gamma * I) - (mu * I)
#       dR =  (gamma * I) - (mu * R)
#       return(list(c(dS, dE, dI, dR)))
#     })
#   }
#   parameters_values = c(beta = beta, gamma = gamma, De = De)
#   S0 = N - I0 - R0
#   E0 = 3 * I0
#   initial_values = c(S = S0, E = E0, I = I0, R = R0)
#   out = ode(initial_values, times, sir_equations, parameters_values, method = "rk4")
#   as.data.frame(out)
# }

# logli_seir_1_de_demography = function(beta, gamma, N, dat, lambda, mu, De) {
#   I0 = dat$I[1]
#   R0 = dat$R[1]
#   times = dat[["time"]]
#   beta = exp(beta)
#   gamma = exp(gamma)
#   predictions = seir_1_de_demography(beta = beta, gamma = gamma, I0 = I0, R0 = R0, 
#         times = times, N = N, lambda = lambda, mu = mu, De = De)
#   ## negative of log likelihood
#   -sum(dpois(x = dat$I, lambda = predictions$I, log = TRUE))
# }
