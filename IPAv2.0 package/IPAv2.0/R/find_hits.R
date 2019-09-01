#' @title Finding putative annotations
#'
#' @description
#' This function finds all the possible annotations given an user-defined ppm and LogRT window.
#'
#' @param adducts.matrix A matrix contains information about all possible adducts.
#' Columns are: KEGG.id (e.g. 'C00002'), adduct (e.g. 'M+H'), RT (previously known
#' retention time for the compound, NA if unknown), formula (e.g. C10H17N5O13P3),
#' theoretical mz (e.g. 508.003570124), charge (e.g. -1) mono ('mono'
#' if monoisotopical form, 'iso' otherwise)
#' @param dataset A matrix containing the measured data, organized in 3 colums: mz, RT and Int
#' @param ppm.thr A numerical value indicating the maximum accuracy value to be considered
#' @param RT.thr A numerical value indicating the maximum difference allowed
#' between measured LogRT and predicted LogRT
#' @param isotopes A matrix containing infomation about isotopes
#' @param iso.threshold A numerical value indicating the probability below which
#' isotope peaks can be omitted
#' @param corr.matrix A matrix containing the correlation values between peaks
#' @param corr.thr A numerical value expressing the treshold used to consider the correlations significant
#' @param relation.id A vector containg class labels of the previously grouped peaks
#' @param v A logical value indicating if the progress will be shown (default TRUE)
#' @param IT A number inticating after how many iteration an update should be shown (default 120)
#'
#' @return A list containing the matrix of the posterior probabilities, the id.masses vector, the all.formulas dataframe and
#' the allsampcomp matrix containing all the assignments for each iteration of the Gibbs sampler
#'
#' @author Francesco Del Carratore \email{francescodc87@@gmail.com}
#'
#' @seealso find.Hits IPAposteriors
#'
#' @import Matrix
#' @import enviPat
#' @import stringr
#' @importFrom utils combn
#' @export

"find.hits" <- function(adducts.matrix, dataset, ppm.thr, RT.thr, isotopes, iso.threshold = 1, corr.matrix = NULL, corr.thr = 0.75, relation.id = NULL, v = T, IT = 500) {
  library(dplyr)
  cat("Finding hits... \n")
  all.adducts.matrix.hits <- matrix(NA, 0, 9)
  colnames(all.adducts.matrix.hits) <- c("KEGG.id", "adduct", "LogRT", "Formula", "theoretical mz", "charge", "mono", "abundance", "mz.id")
  temp_data <- matrix(data = NA, ncol = 2, nrow=2000)
  # Seperate NAs
  copy_adducts.matrix <- adducts.matrix
  temp <- data.frame("ind"= 1:length(adducts.matrix[,1]))
  copy_adducts.matrix <- cbind(copy_adducts.matrix, temp)
  #put NAs to different table in order to check them only with mass
  NA_adductsRT <- copy_adducts.matrix[is.na(copy_adducts.matrix[,3]),]
  copy_adducts.matrix <- copy_adducts.matrix[!is.na(copy_adducts.matrix[,3]),]
  copy_adducts.matrix[,3] <- as.numeric(as.character(copy_adducts.matrix[,3]))
  # For the rest of the data we will check both m/z and LogRT.
  mz_LogRT_adducts <- copy_adducts.matrix
  mz_LogRT_adducts <- as.matrix(mz_LogRT_adducts)
  NA_adductsRT <- as.matrix(NA_adductsRT)
  for (k in 1:nrow(dataset)) {
    mz <- as.numeric(dataset[k, 1])
    RT <- as.numeric(dataset[k, 2])
    deltaM <- ppm.thr * mz * 1e-06
    deltaRT <- RT.thr

    ind <- NULL
    if (length(NA_adductsRT[,1]) > 0){
      temp1 <- NA_adductsRT[which(abs(as.numeric(NA_adductsRT[, 5]) - mz) <= deltaM),8]
      ind <- union_all(temp1, ind)
    }
    if (length(mz_LogRT_adducts[,1]) > 0){
      temp2 <- mz_LogRT_adducts[which((abs(as.numeric(mz_LogRT_adducts[, 5]) - mz) <= deltaM) & (abs(as.numeric(mz_LogRT_adducts[,3]) - log(RT)) <= deltaRT)),8]
      ind <- union_all(temp2, ind)
    }
    ind <- as.numeric(ind)
    if (length(ind) > 0) {
      for (h in 1:length(ind)) {
        ### for each hit I will compute the isotopes and look for hits!
        tmp <- t(c(adducts.matrix[ind[h], ], abundance = 100, mz.id = k))
        isotable <- isopattern(isotopes = isotopes, adducts.matrix[ind[h], 4], charge = as.numeric(adducts.matrix[ind[h], 6]), verbose = FALSE, threshold = iso.threshold)[[1]]
        isotable <- formula_isopattern(isotable, isotopes)
        if (nrow(isotable) > 1) {
          isotable <- isotable[order(isotable[, 2], decreasing = T), ]
        }
        ### find hits
        flag.iso = TRUE
        count = 1
        while (flag.iso & nrow(isotable) > 1 & count < nrow(isotable)) {
          count = count + 1
          deltaM.iso <- ppm.thr * (isotable[count, 1]) * 1e-06
          if (is.null(corr.matrix) & is.null(relation.id)) {
            ind.iso <- which(abs(as.numeric(dataset[, 1]) - isotable[count, 1]) <= deltaM.iso & abs(log(as.numeric(dataset[, 2])) - log(RT)) <= deltaRT)
          } else if (is.null(relation.id)) {
            ind.iso <- which(abs(as.numeric(dataset[, 1]) - isotable[count, 1]) <= deltaM.iso & as.numeric(corr.matrix[k, ]) >= corr.thr)
          } else {
            ind.iso <- which(abs(as.numeric(dataset[, 1]) - isotable[count, 1]) <= deltaM.iso & relation.id == relation.id[k])
          }
          if (length(ind.iso) == 0) {
            flag.iso <- FALSE
          } else {
            for (i in 1:length(ind.iso)) {
              tmp.iso <- c(adducts.matrix[ind[h], 1:3], rownames(isotable)[count], isotable[count, 1], adducts.matrix[ind[h], 6], "iso", isotable[count, 2], mz.id = ind.iso[i])
              tmp <- rbind(tmp, tmp.iso)
            }
          }
        }
        all.adducts.matrix.hits <- rbind(all.adducts.matrix.hits, tmp)
      }
    }
    if (v) {
      if (k%%IT == 0) {
        # Print on the screen some message
        cat(paste0(round(k/nrow(dataset), 1) * 100, "%", "\n"))
      }
    }
  }
  rownames(all.adducts.matrix.hits) <- NULL
  cat("Mapping results... \n")
  id.masses <- as.numeric(unique(all.adducts.matrix.hits[, 9]))
  all.formulas <- unique(all.adducts.matrix.hits[, 1:8], MARGIN = 1)
  hit.matrix <- Matrix(0, nrow = length(id.masses), ncol = nrow(all.formulas))
  for (k in 1:length(id.masses)) {
    ind <- which(all.adducts.matrix.hits[, 9] == id.masses[k])
    tmp <- unique(all.adducts.matrix.hits[ind, 1:8])
    ind <- which(duplicated(rbind(all.formulas, tmp), fromLast = T))
    hit.matrix[k, ind] <- 1
    if (v) {
      if (k%%IT == 0) {
        # Print on the screen some message
        cat(paste0(round(k/length(id.masses), 1) * 100, "%", "\n"))
      }
    }

  }
  out <- list(id.masses = id.masses, all.formulas = all.formulas, hit.matrix = hit.matrix)
  return(out)
}
