#' @title Training Linear Regression Model
#'
#' @description
#' This function trains a Linear Regression based on the standards that the user provides.
#'
#' @param data A matrix contains information about the predictors for the standard. The data frame should have the following columns: RT, logD, pos, Phosphates, HBD_MW, neg, Rotarable_Bonds.
#'
#' @return A trained Linear Regression Model.
#'
#' @author Francesco Del Carratore \email{francescodc87@@gmail.com}
#'
#' @seealso trainModel create_COV_matrix build_model
#'
#' @export

"trainModel" <- function (data){
  RT <- log(data[,1]) #calculate the log(RT)
  log_data <- cbind(RT, data[,c("logD", "pos", "Phosphates", "HBD_MW", "neg", "Rotarable_Bonds")]) #combine log(RT) and the rest of the data

  model <- lm(RT ~ logD  + pos  + Phosphates + HBD_MW + neg + Rotarable_Bonds, data=log_data)

  return(model)
}
