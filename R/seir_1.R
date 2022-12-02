seir_1 = function(beta, gamma, I0, R0, times, N, lambda, mu, De) {
  sir_equations = function(time, variables, parameters) {
    with(as.list(c(variables, parameters)), {
      dS = -beta * I * S/N + lambda * N - mu * S
      dE =  beta * I * S/N - 1/De * E - mu* E
      dI =  1/De * E - gamma * I - mu * I
      dR =  gamma * I - mu*R
      return(list(c(dS, dE, dI, dR)))
    })
  }
  parameters_values = c(beta  = beta, gamma = gamma)
  S0= N-I0-R0
  E0= 3*I0
  initial_values = c(S = S0, E=E0, I = I0, R = R0)
  out = ode(initial_values, times, sir_equations, parameters_values, method = "rk4")
  as.data.frame(out)
}