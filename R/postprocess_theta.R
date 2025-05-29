postprocess_theta <- function(theta_par, NXS, NXO, NE, NV, XS, XO, outcomeS, outcomeC) {
  # Etapa 1: Construir nomes dos parâmetros
  param_names <- c(colnames(XS), colnames(XO))

  # Parte sigma
  if (NE == 1) {
    param_names <- c(param_names, "sigma")
  } else {
    param_names <- c(param_names, "interceptS", colnames(outcomeS))
  }

  # Parte rho
  if (NV == 1) {
    param_names <- c(param_names, "correlation")
  } else {
    param_names <- c(param_names, "interceptC", colnames(outcomeC))
  }

  names(theta_par) <- param_names

  # Etapa 2: Reparametrização
  idx_sigma_start <- NXS + NXO + 1
  idx_rho_start <- idx_sigma_start + ifelse(NE == 1, 1, NE + 1)

  theta_final <- c(
    theta_par[1:(NXS + NXO)],

    # sigma
    if (NE == 1) {
      exp(theta_par[idx_sigma_start])
    } else {
      theta_par[idx_sigma_start:(idx_sigma_start + NE)]
    },

    # rho
    if (NV == 1) {
      tanh(theta_par[idx_rho_start])
    } else {
      theta_par[idx_rho_start:(idx_rho_start + NV)]
    }
  )

  return(theta_final)
}
