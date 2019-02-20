#' @title Generic calibration data available to psifluo package
#' @description This calibration file corresponds to Doug Campbell's PSI_F3500_Calibration.csv
#' @param ./data-raw/PSI_F3500_Calibration.csv path to the calibration file
#' @return Returns a calibration table that can be used by the psifluo functions in default mode.
#' @keywords internal
generic_cali<-read.csv('./data-raw/PSI_F3500_Calibration.csv')
devtools::use_data(generic_cali)
