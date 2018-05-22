#make generic calibration data available to psifluo package

generic_cali<-read.csv('./data-raw/PSI_F3500_Calibration.csv')
devtools::use_data(generic_cali)
