#' @title Update the database with the new predicted LogRT
#'
#' @description
#' This function updates a user defined database with the new prediction for LogRT.
#'
#' @param data A matrix contains information about the predictors for the standard. The data frame should have the following columns: RT, logD, pos, Phosphates, HBD_MW, neg, Rotarable_Bonds.
#' @param model A trained prediction model
#'
#' @return The updated database.
#'
#' @author Francesco Del Carratore \email{francescodc87@@gmail.com}
#'
#' @seealso trainModel create_COV_matrix build_model
#'
#' @export


"updateDB" <- function(model, data){
  LogRT <- predict(model, data)
  LogRT <- data.frame(LogRT)
  final_data <- cbind(data, LogRT)

  return(final_data)
}
