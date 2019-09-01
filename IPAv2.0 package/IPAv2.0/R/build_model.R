#' @title Traing the prediction model, calculating covariance matrix and updatinf the database.
#'
#' @description
#' This function creates the covariance matrix based on the given data.
#'
#' @param train_data A matrix contains information about the predictors for the standard. The data frame should have the following columns: RT, logD, pos, Phosphates, HBD_MW, neg, Rotarable_Bonds.
#' @param cov_data A data frame containing the covariance data.
#' The data frame should have the following columns: Name, Predicted.RT, measured.RT, measured.mz, theoretical.mz, adduct.
#' @param DB The user's metabolite database.
#'
#' @return The covariance matrix, the standard deviation of the models residuals and the updated database.
#'
#' @author Francesco Del Carratore \email{francescodc87@@gmail.com}
#'
#' @seealso trainModel create_COV_matrix build_model
#'
#' @export


"build_model" <- function(train_data, DB, cov_data){
  model <- trainModel(train_data)
  sd_LogRT <- sd(model$residual)
  cov_matrix <- create_COV_matrix(sd_LogRT, cov_data)
  new_DB <- updateDB(model, DB)
  out <- list("cov_matrix" = cov_matrix, "updated_DB" = new_DB, "sd_LogRT" = sd_LogRT)

  return(out)
}
