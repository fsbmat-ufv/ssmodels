#' Post-process Parameter Vector for Generalized Heckman Models
#'
#' @description
#' Internal helper function that assigns meaningful names to a vector of estimated
#' parameters and applies transformations to \code{sigma} and \code{rho} parameters
#' to obtain interpretable estimates.
#'
#' @details
#' The parameter vector \code{theta_par} is structured as follows:
#' \itemize{
#'   \item The first \code{NXS} elements are coefficients for the selection equation.
#'   \item The next \code{NXO} elements are coefficients for the outcome equation.
#'   \item The next \code{NE} elements are coefficients (or log-sigma if \code{NE == 1}) for the scale model.
#'   \item The next \code{NV} elements are coefficients (or atanh(rho) if \code{NV == 1}) for the correlation model.
#' }
#' For interpretation, the function applies:
#' \itemize{
#'   \item \code{exp()} transformation for \code{sigma} if \code{NE == 1}.
#'   \item \code{tanh()} transformation for \code{rho} if \code{NV == 1}.
#' }
#'
#' @param theta_par A numeric vector containing the estimated parameters.
#' @param NXS Integer. Number of covariates in the selection equation.
#' @param NXO Integer. Number of covariates in the outcome equation.
#' @param NE Integer. Number of covariates (or 1 for intercept-only) in the scale model.
#' @param NV Integer. Number of covariates (or 1 for intercept-only) in the correlation model.
#' @param XS Design matrix for the selection equation (used for naming).
#' @param XO Design matrix for the outcome equation (used for naming).
#' @param outcomeS Design matrix or variable for the scale (variance) model.
#' @param outcomeC Design matrix or variable for the correlation model.
#'
#' @return
#' A named numeric vector with:
#' \itemize{
#'   \item Transformed \code{sigma} and \code{rho} values (if needed),
#'   \item Meaningful names assigned to all parameters.
#' }
#'
#' @keywords internal
#' @export
postprocess_theta <- function(theta_par, NXS, NXO, NE, NV, XS, XO, outcomeS, outcomeC) {

  # Step 1: Construct parameter names using column names from input matrices
  param_names <- c(colnames(XS), colnames(XO))

  # Append names for scale parameters (sigma)
  if (NE == 1) {
    param_names <- c(param_names, "sigma")
  } else {
    param_names <- c(param_names, "interceptS", colnames(outcomeS))
  }

  # Append names for correlation parameters (rho)
  if (NV == 1) {
    param_names <- c(param_names, "correlation")
  } else {
    param_names <- c(param_names, "interceptC", colnames(outcomeC))
  }

  # Assign names to the full parameter vector
  names(theta_par) <- param_names

  # Step 2: Define indices for sigma and rho blocks
  idx_sigma_start <- NXS + NXO + 1
  idx_rho_start <- idx_sigma_start + ifelse(NE == 1, 1, NE + 1)

  # Step 3: Reconstruct the transformed parameter vector
  theta_final <- c(
    theta_par[1:(NXS + NXO)],

    # Apply transformation to sigma if needed
    if (NE == 1) {
      exp(theta_par[idx_sigma_start])
    } else {
      theta_par[idx_sigma_start:(idx_sigma_start + NE)]
    },

    # Apply transformation to rho if needed
    if (NV == 1) {
      tanh(theta_par[idx_rho_start])
    } else {
      theta_par[idx_rho_start:(idx_rho_start + NV)]
    }
  )

  return(theta_final)
}
