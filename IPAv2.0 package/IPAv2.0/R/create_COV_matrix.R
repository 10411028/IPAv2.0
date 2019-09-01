#' @title Creating the Covariance Matrix
#'
#' @description
#' This function creates the covariance matrix based on the given data.
#'
#' @param sd_LogRT The standard deviation of the residuals of the model.
#' @param cov_data A data frame containing the covariance data.
#' The data frame should have the following columns: Name, Predicted.RT, measured.RT, measured.mz, theoretical.mz, adduct.
#'
#' @return The covariance matrix.
#'
#' @author Francesco Del Carratore \email{francescodc87@@gmail.com}
#'
#' @seealso trainModel create_COV_matrix build_model
#'
#' @export


"create_COV_matrix" <- function(sd_LogRT, cov_data){
  ppm <- ((cov_data$theoretical.mz - cov_data$measured.mz) / cov_data$theoretical.mz) * 10^6
  sd_ppm <- sd(ppm)
  DeltaLogRT <- log(cov_data$measured.RT) - log(cov_data$Predicted.RT)
  sd_DeltaLogRT <- sd(DeltaLogRT)
  covariance <- cov(DeltaLogRT, ppm) # Covariance between mass and RT
  cov_matrix <- matrix(c(sd_ppm^2, covariance, covariance, sd_DeltaLogRT^2), ncol = 2)
  colnames(cov_matrix) <- c("ppm", "LogRT")
  rownames(cov_matrix) <- c("ppm", "LogRT")
  return(cov_matrix)
}
