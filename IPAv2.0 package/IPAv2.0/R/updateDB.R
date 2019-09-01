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


"updateDB" <- function(model, DB){
  data <- data.frame(DB[,15:21])
  data[] <- lapply(data, function(x) as.numeric(as.character(x)))
  
  LogRT <- predict(model, data)
  LogRT <- data.frame(LogRT)
  
  new_DB <- data.frame(DB[,1:14])
  names <- colnames(new_DB)
  names[1] <- "KEGG.id"
  new_DB[10] <- LogRT
  colnames(new_DB) <- names
  new_DB <- cbind(new_DB, data.frame(DB[,15:21]))
  new_DB <- as.matrix(new_DB)

  return(new_DB)
}
