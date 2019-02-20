#' @title Wrapper function to import data from PSI FRRf and fit induction/re-opening models to the dataset
#' @description This function will import a PSI Fluorwin output file and will attempt to fit FRRf induction/re-opening models to the dataset. It can process files generated with Doug Campbell's double FRRf protocol and with the PSI RLC wizard (not implemented yet).
#' @param file_name Path to the file to be imported
#' @param out_name Name the file to be exported with the results. Only works if the user selects out_files=TRUE.
#' @param calib_file Path to the calibration file to convert lamp voltage settings to photons. A^-2^ . s^-1^. If not provided the default is a calibration function from a Doug Campbell PSI FL3500 instrument. If you are not using his machine the parameters sigmaPSII, cumulative excitation and actinic PAR values associated to an individual machine will be wrongly calibrated.
#' @param sti_model Parameters to fit to the induction curve. Default is "all". Options are: "all", "sigma_rho" and "sigma". "All" - fo, fm, sigma and rho are fitted; "sigma_rho" -  sigma and rho are fitted, fo and fm are fixed at user-defined levels, defaults are fo =  1st point of the fluorescence data vector y, fm = average of the maximum 5 points of the fluorescence data vector y;  "sigma" - fo, fm and rho are fixed, only sigma is fitted.
#' @param str_model Type of re-opening model. Options: "tau1", "tau2", "tau3" and "tau1_subset" (the user selects a subset of the data starting from the first point of the re-opening curve). Default is tau1, fitting re-opening with a single phase exponential decay with lifetime of tau1.
#' @param subset_time selects a subset of the re-opening dataset based on time (micro seconds) from the first point. To be used only with the tau1_subset model.
#' @param protocol Type of protocol used to produce the PSI file. Currently only "campbell" is implemented.
#' @param turn_plots_off Logical, plots off. Default = FALSE
#' @param plot_matrix Logical, should plots be plotted to a plot matrix. Default = TRUE
#' @param out_files Logical, should results (fitted parameters) be exported to a csv file. Default = FALSE
#' @param out_pdf Logical, should plot be exported to a pdf file. Default = FALSE
#' @param dec decimal separator, default= "."
#' @return Returns a list with 3 elements. The first element (data) contains five data frames and one numeric vector: 1 - all.data  = raw data input; 2 - meta.data = acquisiton settings; flashlet_energy = the calculated flashlet energy from the calibration file; 3 - data.sti  = induction only; 4 - data.str  = re-opening only; 5 - measuring_steps = voltage and equivalent actinic PAR at each step. The second element (sti_parameters) contains a data frame with fitted parameters for the induction step and respective SE, p-values.The third element (str_parameters) contains a data frame with fitted parameters for the re-opening step and respective SE, p-values.
#' @keywords external
#' @export

process_psi<-function(file_name,
                         out_name = "results",
                         sti_model = "all",
                         str_model = "tau2",
                         calib_file = NA,
                         protocol = "campbell",
                         turn_plots_off = FALSE,
                         plot_matrix = FALSE,
                         out_files = FALSE,
                         out_pdf = FALSE,
                         dec = ".",
                         subset_time = 600){

################################################################################
#read psi data file and stores input data  inside an object "a"

a<-read_psi_file_campbell(filename=file_name,calib_file=calib_file,dec=dec)
#print(a)
#checking if the user wants plots or not
if(turn_plots_off == 'FALSE'){
  plots = 'TRUE'
}else{plots = 'FALSE'}


#storing the number of datasets
number_datasets<-length(a$measuring_steps$step)

################################################################################
#fit induction models

#creating a data.frame that will store the fitted values for each sti data series
#it gets the dimensions from the number of light levels and the number of flashlets
#fitted_values_sti<-data.frame(matrix(ncol = 2*a$meta.data$Light_curves_repeats, nrow = length(a$data.sti[,1])))

fitted_values_sti<-data.frame(matrix(ncol = 2*number_datasets, nrow = length(a$data.sti[,1])))


#creating a data.frame to store the fitted parameters
parameters_sti<-setNames(data.frame(matrix(ncol = 16, nrow = 0)), c('fo_sti','fm_sti','sigma_sti','rho_sti','fo_se_sti','fm_se_sti','sigma_se_sti','rho_se_sti','fo_tvalue_sti','fm_tvalue_sti','sigma_tvalue_sti','rho_tvalue_sti','fo_p_sti','fm_p_sti','sigma_p_sti','rho_p_sti'))


###############################################################################
#####OLD CODE
#counter<-1
#for (i in 2:(2*a$meta.data$Light_curves_repeats+1))
#{
#  a1<-fit_sti(x = a$flashlet_energy,y = a$data.sti[,i],fit_model = sti_model,plots = plots)
#  legend('topleft', legend=names(a$data.sti[i]),bty='n',cex=1)
#  parameters_sti[counter,]<-a1$fitted_parameters
#  fitted_values_sti[,counter]<-a1$fitted_values
#  counter<-counter+1
#}
###############################################################################
#print(a$data.sti)
#print(number_datasets)
#print(number_datasets+1)

#storing the number of datasets light and dark measurements
number_datasets_LD<-number_datasets*2

counter<-1

for (i in 2:(number_datasets_LD+1))
{
  a1<-fit_sti(x = a$flashlet_energy,y = a$data.sti[,i],fit_model = sti_model,plots = plots)
  legend('topleft', legend=names(a$data.sti[i]),bty='n',cex=1)
  parameters_sti[counter,]<-a1$fitted_parameters
  fitted_values_sti[,counter]<-a1$fitted_values
  counter<-counter+1
}
###############################################################################

#print(parameters_sti)

#Adding a column to the sti_parameters with the measuring step (i.e. light vs dark)
names_steps<-names(a$data.sti)

names_steps<-names_steps[-c(1,length(names_steps))]

#print(names_steps)

parameters_sti$steps<-names_steps

#print(parameters_sti)
################################################################################
#fit relaxation models

#creating a data.frame that will store the fitted values for each str data series
#it gets the dimensions from the number of light levels and the number of flashlets

if (str_model == 'tau1_subset'){
  x_subset<-a$data.str$time[a$data.str[,1]<=subset_time]
  fitted_values_str<-data.frame(matrix(ncol = 2*number_datasets, nrow = length(x_subset)))}else{
  fitted_values_str<-data.frame(matrix(ncol = 2*number_datasets, nrow = length(a$data.str[,1])))
  }




parameters_str<-setNames(data.frame(matrix(ncol = 30, nrow = 0)), c('fo_str','fm_str','tau1_str','alpha1_str','tau2_str','alpha2_str','tau3_str','fo_se_str','fm_se_str','tau1_se_str','alpha1_se_str','tau2_se_str','alpha2_se_str','tau3_se_str','fo_tvalue_str','fm_tvalue_str','tau1_tvalue_str','alpha1_tvalue_str','tau2_tvalue_str','alpha2_tvalue_str','tau3_tvalue_str','fo_p_str','fm_p_str','tau1_p_str','alpha1_p_str','tau2_p_str','alpha2_p_str','tau3_p_str','alpha3_str','alpha3_se_str'))



counter<-1
for (i in 2:(number_datasets_LD+1))
{
  a1<-fit_str(a$data.str[,1],a$data.str[,i],
              tau_model=str_model,plots=plots,subset_time = subset_time)
  legend('bottomleft', legend=names(a$data.str[i]),bty='n',cex=1)
  parameters_str[counter,]<-a1$fitted_parameters
  fitted_values_str[,counter]<-a1$fitted_values
  counter<-counter+1
}

#Adding a column to the str_parameters with the measuring step (i.e. light vs dark)
parameters_str$steps<-names_steps



################################################################################
#preparing outputs
output<-list(a,parameters_sti,parameters_str)
names(output)<-c('data','sti_parameters','str_parameters')


################################################################################
#prepare output files if the user chose to have them
if (out_files==TRUE){
  #write.table(as.data.frame(a[1]),file="data_psi.txt",col.names = TRUE,row.names = FALSE,sep=";")
  write.table(parameters_sti,file=paste(out_name,"_parameters_sti.txt",sep=""),col.names = TRUE,row.names = FALSE,sep=";")
  write.table(parameters_str,file=paste(out_name,"_parameters_str.txt",sep=""),col.names = TRUE,row.names = FALSE,sep=";")
   }

#option of turning all plotting off
if (turn_plots_off=='FALSE'){
################################################################################
#plot plot mosaic if plot_matrix is TRUE
################################################################################


if (plot_matrix==TRUE){

  #number of plots
  num_plots=2*a$meta.data$Light_curves_repeats

  #closing ALL open devices
    graphics.off()
################################
#plotting STI data
################################

dev.set()
dev.new()
nf<-layout(matrix(c(1:(ceiling(num_plots/4)*4)),nrow=ceiling(num_plots/4),ncol=4,byrow=TRUE))
#nf<-layout(matrix(c(1:(ceiling(no_light_steps/4)*4)),nrow=ceiling(no_light_steps/4),ncol=4,byrow=TRUE))

par(oma=c(4,4,3,4),mar=c(0,0,0,0), xpd=NA,tcl=-0.3,bg="white",cex=0.8,cex.axis=0.9,cex.lab=0.9,bty="o",las=1,mgp=c(3,0.5,0),adj=0.5)

for (i in 1:num_plots){
    plot(a$data.sti[,i+1],xlab="", ylab="", bty="l",xaxt='n',yaxt='n',pch=21,col="blue")
  lines(fitted_values_sti[,i],type="l",col="red",lwd=1.5)
  legend('topleft', legend=names(a$data.sti[i+1]),bty='n',cex=0.7)

  #this is for the legend at the bottom right
  # Light: Print optimized parameter values on the plot ####
  leg.L.model<-paste("Model: ",sti_model,sep="");
  #leg.L.main<-"Parameter \u00B1 SE" # \u00B1 is the unicode character for the plus/minus symbol
  leg.L.Fo<-tryCatch({paste("Fo =", round(parameters_sti$fo_sti[i], digits=3), "\u00B1", round(parameters_sti$fo_se_sti[i], digits=3), sep=" ")},error=function(Fit.Fail){return("Fit Fail")})
  leg.L.Fm<-tryCatch({paste("Fm =", round(parameters_sti$fm_sti[i], digits=3), "\u00B1", round(parameters_sti$fm_se_sti[i], digits=3), sep=" ")},error=function(Fit.Fail){return("Fit Fail")})
  leg.L.sigma<-tryCatch({paste("sigma =", round(parameters_sti$sigma_sti[i], digits=0), "\u00B1", round(parameters_sti$sigma_se_sti[i], digits=0), sep=" ")},error=function(Fit.Fail){return("Fit Fail")})
  leg.L.rho<-tryCatch({paste("rho =", round(parameters_sti$rho_sti[i], digits=3), "\u00B1", round(parameters_sti$rho_se_sti[i], digits=3), sep=" ")},error=function(Fit.Fail){return("Fit Fail")})

  legend("bottomright",legend=c(leg.L.model,leg.L.Fo,leg.L.Fm,leg.L.sigma,leg.L.rho),bty="n",text.font=c(2,1,1,1,1),adj=c(0, 0.5),cex = 0.7)


  }


################################
#plotting STR data
################################

dev.new()
nf<-layout(matrix(c(1:(ceiling(num_plots/4)*4)),nrow=ceiling(num_plots/4),ncol=4,byrow=TRUE))

par(oma=c(4,4,3,4),mar=c(0,0,0,0), xpd=NA,tcl=-0.3,bg="white",cex=0.8,cex.axis=0.9,cex.lab=0.9,bty="o",las=1,mgp=c(3,0.5,0),adj=0.5)

for (i in 1:num_plots){
  plot(a$data.str[,i+1],xlab="", ylab="", bty="l",xaxt='n',yaxt='n',cex=1,pch=21,col="blue")
  lines(fitted_values_str[,i],type="l",col="red",lwd=1.5)
  legend('bottomleft', legend=names(a$data.str[i+1]),bty='n',cex=0.7)

  # Light: Print optimized parameter values on the plot ####
  leg.L.model<-paste("Model: ",str_model,sep="");
  #leg.L.main<-"Parameter \u00B1 SE" # \u00B1 is the unicode character for the plus/minus symbol
  leg.L.Fo<-tryCatch({paste("Fo =", round(parameters_str$fo_str[i], digits=3), "\u00B1", round(parameters_str$fo_se_str[i], digits=3), sep=" ")},error=function(Fit.Fail){return("Fit Fail")})
  leg.L.Fm<-tryCatch({paste("Fm =", round(parameters_str$fm_str[i], digits=3), "\u00B1", round(parameters_str$fm_se_str[i], digits=3), sep=" ")},error=function(Fit.Fail){return("Fit Fail")})
  leg.L.tau1<-tryCatch({paste("tau1 =", round(parameters_str$tau1_str[i], digits=5), "\u00B1", round(parameters_str$tau1_se_str[i], digits=5), sep=" ")},error=function(Fit.Fail){return("Fit Fail")})
  leg.L.alpha1<-tryCatch({paste("alpha1 =", round(parameters_str$alpha1_str[i], digits=2), "\u00B1", round(parameters_str$alpha1_se_str[i], digits=2), sep=" ")},error=function(Fit.Fail){return("Fit Fail")})
  leg.L.tau2<-tryCatch({paste("tau2 =", round(parameters_str$tau2_str[i], digits=5), "\u00B1", round(parameters_str$tau2_se_str[i], digits=5), sep=" ")},error=function(Fit.Fail){return("Fit Fail")})
  leg.L.alpha2<-tryCatch({paste("alpha2 =", round(parameters_str$alpha2_str[i], digits=2), "\u00B1", round(parameters_str$alpha2_se_str[i], digits=2), sep=" ")},error=function(Fit.Fail){return("Fit Fail")})
  leg.L.tau3<-tryCatch({paste("tau3 =", round(parameters_str$tau3_str[i], digits=5), "\u00B1", round(parameters_str$tau3_se_str[i], digits=5), sep=" ")},error=function(Fit.Fail){return("Fit Fail")})

  if(is.na(parameters_str$alpha1_str[i]) && is.na(parameters_str$alpha2_str[i])){
    legend("topright",legend=c(leg.L.model,leg.L.Fo,leg.L.Fm,leg.L.tau1),bty="n",text.font=c(2,1,1,1),adj=c(0, 0.5),cex=0.7)
  }
  if(!is.na(parameters_str$alpha1_str[i])&is.na(parameters_str$tau3_str[i])){

    legend("topright",legend=c(leg.L.model,leg.L.Fo,leg.L.Fm,leg.L.tau1,leg.L.alpha1,leg.L.tau2),bty="n",text.font=c(2,1,1,1,1,1),adj=c(0, 0.5),cex=0.7)
  }
  if(!is.na(parameters_str$alpha2_str[i])){

    legend("topright",legend=c(leg.L.model,leg.L.Fo,leg.L.Fm,leg.L.tau1,leg.L.alpha1,leg.L.tau2,leg.L.alpha2,leg.L.tau3),bty="n",text.font=c(2,1,1,1,1,1,1,1),adj=c(0, 0.5),cex=0.7)
  }



}


}

#end of turn_plots_off
}
################################################################################
#output pdf with the plots if the user decided


#open pdf device if the user wants to store the plots in a single pdf file
if (out_pdf==TRUE){
  pdf("output_plots.pdf")

  #number of plots
  num_plots=2*a$meta.data$Light_curves_repeats

  for (i in 1:num_plots){
    plot(a$data.sti[,i+1],xlab="Flashlets", ylab="Fluorescence", bty="l",cex=1,pch=21,col="blue",las=1)
    lines(fitted_values_sti[,i],type="l",col="red",lwd=1.5)
    legend('topleft', legend=names(a$data.sti[i+1]),bty='n')
    #this is for the legend at the bottom right
    # Light: Print optimized parameter values on the plot ####
    leg.L.model<-paste("Model: ",sti_model,sep="");
    #leg.L.main<-"Parameter \u00B1 SE" # \u00B1 is the unicode character for the plus/minus symbol
    leg.L.Fo<-tryCatch({paste("Fo =", round(parameters_sti$fo_sti[i], digits=3), "\u00B1", round(parameters_sti$fo_se_sti[i], digits=3), sep=" ")},error=function(Fit.Fail){return("Fit Fail")})
    leg.L.Fm<-tryCatch({paste("Fm =", round(parameters_sti$fm_sti[i], digits=3), "\u00B1", round(parameters_sti$fm_se_sti[i], digits=3), sep=" ")},error=function(Fit.Fail){return("Fit Fail")})
    leg.L.sigma<-tryCatch({paste("sigma =", round(parameters_sti$sigma_sti[i], digits=0), "\u00B1", round(parameters_sti$sigma_se_sti[i], digits=0), sep=" ")},error=function(Fit.Fail){return("Fit Fail")})
    leg.L.rho<-tryCatch({paste("rho =", round(parameters_sti$rho_sti[i], digits=3), "\u00B1", round(parameters_sti$rho_se_sti[i], digits=3), sep=" ")},error=function(Fit.Fail){return("Fit Fail")})

    legend("bottomright",legend=c(leg.L.model,leg.L.Fo,leg.L.Fm,leg.L.sigma,leg.L.rho),bty="n",text.font=c(2,1,1,1,1),adj=c(0, 0.5),cex = 1)

  }

  for (i in 1:num_plots){
    plot(a$data.str[,i+1],xlab="Time (log(seconds))", ylab="Fluorescence", bty="l",cex=1,pch=21,col="blue",las=1)
    lines(fitted_values_str[,i],type="l",col="red",lwd=1.5)
    legend('bottomleft', legend=names(a$data.str[i+1]),bty='n')


    # Light: Print optimized parameter values on the plot ####
    leg.L.model<-paste("Model: ",str_model,sep="");
    #leg.L.main<-"Parameter \u00B1 SE" # \u00B1 is the unicode character for the plus/minus symbol
    leg.L.Fo<-tryCatch({paste("Fo =", round(parameters_str$fo_str[i], digits=3), "\u00B1", round(parameters_str$fo_se_str[i], digits=3), sep=" ")},error=function(Fit.Fail){return("Fit Fail")})
    leg.L.Fm<-tryCatch({paste("Fm =", round(parameters_str$fm_str[i], digits=3), "\u00B1", round(parameters_str$fm_se_str[i], digits=3), sep=" ")},error=function(Fit.Fail){return("Fit Fail")})
    leg.L.tau1<-tryCatch({paste("tau1 =", round(parameters_str$tau1_str[i], digits=5), "\u00B1", round(parameters_str$tau1_se_str[i], digits=5), sep=" ")},error=function(Fit.Fail){return("Fit Fail")})
    leg.L.alpha1<-tryCatch({paste("alpha1 =", round(parameters_str$alpha1_str[i], digits=2), "\u00B1", round(parameters_str$alpha1_se_str[i], digits=2), sep=" ")},error=function(Fit.Fail){return("Fit Fail")})
    leg.L.tau2<-tryCatch({paste("tau2 =", round(parameters_str$tau2_str[i], digits=5), "\u00B1", round(parameters_str$tau2_se_str[i], digits=5), sep=" ")},error=function(Fit.Fail){return("Fit Fail")})
    leg.L.alpha2<-tryCatch({paste("alpha2 =", round(parameters_str$alpha2_str[i], digits=2), "\u00B1", round(parameters_str$alpha2_se_str[i], digits=2), sep=" ")},error=function(Fit.Fail){return("Fit Fail")})
    leg.L.tau3<-tryCatch({paste("tau3 =", round(parameters_str$tau3_str[i], digits=5), "\u00B1", round(parameters_str$tau3_se_str[i], digits=5), sep=" ")},error=function(Fit.Fail){return("Fit Fail")})

    if(is.na(parameters_str$alpha1_str[i]) && is.na(parameters_str$alpha2_str[i])){
      legend("topright",legend=c(leg.L.model,leg.L.Fo,leg.L.Fm,leg.L.tau1),bty="n",text.font=c(2,1,1,1),adj=c(0, 0.5),cex=1)
    }
    if(!is.na(parameters_str$alpha1_str[i])&is.na(parameters_str$tau3_str[i])){

      legend("topright",legend=c(leg.L.model,leg.L.Fo,leg.L.Fm,leg.L.tau1,leg.L.alpha1,leg.L.tau2),bty="n",text.font=c(2,1,1,1,1,1),adj=c(0, 0.5),cex=1)
    }
    if(!is.na(parameters_str$alpha2_str[i])){

      legend("topright",legend=c(leg.L.model,leg.L.Fo,leg.L.Fm,leg.L.tau1,leg.L.alpha1,leg.L.tau2,leg.L.alpha2,leg.L.tau3),bty="n",text.font=c(2,1,1,1,1,1,1,1),adj=c(0, 0.5),cex=1)
    }



  }


  #end of pdf output
  dev.off()
}







return(output)

}



