#' @title Function to import and fit model
#' @description  Main processing scrip for PSI fluorometer generated files. Script cycles through all FluorWin output files in a specified folder and calculates FRR paramaters.
#' @param sti the user decides if fit should be done in the induction points
#' @param str the user decides if fit should be done in the relaxation points
#' @param sti_model decide which sti model to fit c("fofmsig","fofmsigp","sigp","sig"). "fofmsigp", which fits fo fm sigmaPSII and rho; "fofmsig", which fits fo fm 
#' and sigmaPSII and sets rho to 0; "sigp", which fits sigmaPSII and p, and calculates fo and fm as the mean of points specified in fo_pts and fm_pts; "sig", which fits sigmaPSII and sets p to 0, calcualtes fo and fm as the mean of points specified in fo_pts and fm_pts
#' @param str_model decide which str model to fit
#' @param protocol the current supported protocols are Campbell and RLC
#' @return list all the objects the function returns plus the output files
#' @export 
fit_model_fofm<-function(file_name,
                         out_name="results",
                         samples=1,
                         num_steps=11,
                         sti=TRUE,
                         str=TRUE,
                         sti_model="fofmsig",
                         str_model="tau2",
                         calib_file="Bruno_Calibration.csv",
                         protocol="Campbell"){

# List of settings that will be used throughout the script
settings=list() 
settings$run_sti<-sti
settings$sti_model<-sti_model  
settings$PSI.Calibration.File <- calib_file
# set the maximum value for Sigma PSII
settings$sigma.PSII.max.lim <- 10000 #why do we need this?
settings$run_str<-str
settings$str_model<-str_model


###########################################  
##PARSING FUNCTION
#Added a parameter to select the input file

if (protocol=='RLC'){
  
  read_psi(file_name,out_name,samples,num_steps,calib_file)
  
}else{
  psiworx.read(filename=file_name)
}

#these are objects that are needed for the plotting function
#and that are produced by the parsing function
no_light_steps<-local(no_light_steps,env = e2)

rep_names<-samples

nam<-out_name

###############so we should now have the following data frame available for subsequent analyses
#all.data - raw data
#meta.data - meta information
#cor.data - corrected data with identifier columns
#data.sti - corrected sti data named well
#data.str - corrected str data named well
#time.data - unformatted time data from raw Flurowin output.


#format the time data
#the time format is determined by the computer...
#the first thing is to check if there's an AM or PM string the time.data[,3] object
#it will fail if the time is in third type format

a<-grepl("AM|PM", time.data[,3])

if (a[1]==FALSE){
  formatted.time <- strptime(time.data[,3], format="%d/%m/%Y %H:%M:%S")  
}else{
  formatted.time <- strptime(time.data[,3], format="%m/%d/%Y %I:%M:%S %p")  
}

#print(formatted.time)

#Load PSI light level calibration .csv file - we need to make this!! ##################******************$$############!!!!!!!!!!!!!!!!!!!!!
#BJ: why is this file imported twice? here and in settings$PSI.Calibration.File 
cali <- read.csv(calib_file,header=T)
cali<-assign("cali",cali, envir = e2)

#Extract the power level from the calibration table based on meta data flashlet voltage.
#it comes from the function pw_level()
   
power<-pw_level(protocol)

#print(power)    
############################################ Fit STI

# Power per flash is determined from the calibration .csv and is specific to instument used
incident<-rep(power, times=nrow(data.sti))

#print(incident)

# initialize data frame
sti.fit.names<-c('fo_sti','fm_sti','sigma_sti','rho_sti','fo_se_sti','fm_se_sti','sigma_se_sti','rho_se_sti','fo_tvalue_sti','fm_tvalue_sti','sigma_tvalue_sti','rho_tvalue_sti','fo_p_sti','fm_p_sti','sigma_p_sti','rho_p_sti','model_sti')
sti.fit <- as.data.frame(matrix(nrow=(ncol(data.sti)-1),ncol=length(sti.fit.names)))
names(sti.fit)=sti.fit.names[order(sti.fit.names)]
rm(sti.fit.names)

# open blank .pdf file where all sti plots will be saved, each focus file whill have an individual .pdf output - this uses the output name previously input by the user at the beginning of the script.
#BJ: include a parameter to activate/desactivate the pdf output (true by default) 
pdf(paste(nam,"_sti_fit.pdf",sep=""))

# loop through each column of fluorescence data and apply the STI model####

for(c in 2:ncol(data.sti)){ #note we run from 2nd column this avoiding first column which is time.
    x<-incident
    y<-data.sti[,c]
    wm<-settings$sti_model
    sigma_max<-settings$sigma.PSII.max.lim
    titl<-names(data.sti)[c]

###################################
####FITTING FUNCTION 
fit_sti(x,y,fo.opt,fm.opt,sigma.opt,rho.opt,sigma_max,wm,sti.fit,c)        
    
#Output object needed for the plotting function
#it is stored in the Environment e2
res<-local(res, env = e2)

res$model_sti=wm;
names(sti.fit)<-names(res)
#add res into sti.fit data frame in appropriate row
sti.fit[c-1,]<-res

###################################
####PLOTING FUNCTION 
psi_plots(y, fit1,fit1u,fit1l,titl,res)  
  

}


#turn off pdf catcher for STI curves
dev.off() # puts plots in .pdf created, adobe cannot be open


################################################################################
#screen version
#temporary fix to plot the sti on the screen, this is far from being ideal 
#we are repeating the same calculations and fits that were previously done for the pdf
#maybe it would be better to first create an object with all the data needed for the
#plottings and then do one pdf and one screen version. The non-fitted would be
#marked as NA


#create the layout using the number of light steps
dev.new()
nf<-layout(matrix(c(1:(ceiling(no_light_steps/4)*4)),nrow=ceiling(no_light_steps/4),ncol=4,byrow=TRUE))
par(oma=c(4,4,3,4),mar=c(0,0,0,0), xpd=NA,tcl=-0.3,bg="white",cex=0.8,cex.axis=0.9,cex.lab=0.9,bty="o",las=1,mgp=c(3,0.5,0),adj=0.5)
#layout.show(nf)
#cycle through the light steps

for(c in 2:ncol(data.sti)){ #note we run from 2nd column this avoiding first column which is time.
  x<-incident
  y<-data.sti[,c]
  wm<-settings$sti_model
  sigma_max<-settings$sigma.PSII.max.lim
  titl<-names(data.sti)[c]
  ###################################
  ####FITTING FUNCTION 
  fit_sti(x,y,fo.opt,fm.opt,sigma.opt,rho.opt,sigma_max,wm,sti.fit,c)        
  #Output object needed for the plotting function
  #it is stored in the Environment e2
  res<-local(res, env = e2)
  res$model_sti=wm;
  names(sti.fit)<-names(res)
  #add res into sti.fit data frame in appropriate row
  sti.fit[c-1,]<-res
  ###################################
  ####PLOTING FUNCTION 
  psi_plots_screen(y, fit1,fit1u,fit1l,titl,res)  
}
#dev.off()
#par(def.par)




#checking if the user also wants to plot the relaxation step
if (settings$run_str == TRUE) {
  
################################## FIT the STR model ####################################

# open blank .pdf file where all plots will be saved, each focus file whill have an individual .pdf output

pdf(file=paste(nam,"_str_fit.pdf",sep=""))

    # initialize data frame
    str.fit.names<-c('fo_str','fm_str','tau1_str','alpha1_str','tau2_str','alpha2_str','tau3_str','fo_se_str','fm_se_str','tau1_se_str','alpha1_se_str','tau2_se_str','alpha2_se_str','tau3_se_str','fo_tvalue_str','fm_tvalue_str','tau1_tvalue_str','alpha1_tvalue_str','tau2_tvalue_str','alpha2_tvalue_str','tau3_tvalue_str','fo_p_str','fm_p_str','tau1_p_str','alpha1_p_str','tau2_p_str','alpha2_p_str','tau3_p_str','model_str')
    str.fit <- as.data.frame(matrix(nrow=(ncol(data.str)-1),ncol=length(str.fit.names)))
    names(str.fit)=str.fit.names[order(str.fit.names)]
    rm(str.fit.names)

   
#loop through columns (starting at second to avoid the time column) of data.str i.e. the relaxation curves
for(d in 2:ncol(data.str)){
        #define our new x, y etc for str model inputs
        x<-data.str$time
        y<-data.str[,d]
        wm<-settings$str_model

# change time to seconds
x<-x/1e6;

# First guesses and limits for relaxation fit
  fg=list();
  fg$fo.opt=min(y);
  fg$fm.opt=max(y);
  fg$tau1.opt=0.01;
  fg$alpha1.opt=0.6;
  fg$tau2.opt=fg$tau1.opt/10;
  fg$alpha2.opt=1-fg$alpha1.opt;
  fg$tau3.opt=fg$tau1.opt*2;

    
ub<-c(fo.opt=Inf,fm.opt=Inf,tau1.opt=Inf,alpha1.opt=1,tau2.opt=Inf,alpha2.opt=1,tau3.opt=Inf);
lb<-c(fo.opt=-Inf,fm.opt=-Inf,tau1.opt=-Inf,alpha1.opt=0,tau2.opt=-Inf,alpha2.opt=0,tau3.opt=-Inf);

#fitting str function 
fit_str(x,y,fo.opt = fg$fo.opt,fm.opt=fg$fm.opt,wm,d,fg=fg,ub=ub,lb=lb)

#Output object needed for the plotting function
#it is stored in the Environment e2
#it would be nice to change the name of this object so that it's different
#from the sti function
res<-local(res, env = e2)
    
    
res$model_str=wm;
names(str.fit)<-names(res)
#add res into str.fit data frame in appropriate row
str.fit[d-1,]<-res


#make str plots that will write to already opened pdf file.
#define parameters again here to make certain they are correct
x<-data.str$time/1e6; # change time to seconds
y<-data.str[,d]
wplot<-'str'
res<-str.fit[d-1,]
titl<-names(data.str)[d]

#plot str function
psi_plots_str(x,y, fit1,fit1u,fit1l,titl,res)
      
      
#remove objects to prevent carry over
#rm(x)
#rm(y)
#rm(res)

} #end of str loop through data.str columns

#turn off pdf catcher
dev.off() # puts plots in .pdf created, adobe cannot be open


################################################################################
##Screen version
#temporary fix to plot the STR on the screen, this is far from being ideal 
#we are repeating the same calculations and fits that were previously done for the pdf
#maybe it would be better to first create an object with all the data needed for the
#plottings and then do one pdf and one screen version. The non-fitted would be
#marked as NA

#create the layout using the number of light steps
dev.new()
nf<-layout(matrix(c(1:(ceiling(no_light_steps/4)*4)),nrow=ceiling(no_light_steps/4),ncol=4,byrow=TRUE))
par(oma=c(4,4,3,4),mar=c(0,0,0,0), xpd=NA,tcl=-0.3,bg="white",cex=0.8,cex.axis=0.9,cex.lab=0.9,bty="o",las=1,mgp=c(3,0.5,0),adj=0.5)


#loop through columns (starting at second to avoid the time column) of data.str i.e. the relaxation curves
for(d in 2:ncol(data.str)){
  #define our new x, y etc for str model inputs
  x<-data.str$time
  y<-data.str[,d]
  wm<-settings$str_model
  
  # change time to seconds
  x<-x/1e6;
  
  # First guesses and limits for relaxation fit
  fg=list();
  fg$fo.opt=min(y);
  fg$fm.opt=max(y);
  fg$tau1.opt=0.01;
  fg$alpha1.opt=0.6;
  fg$tau2.opt=fg$tau1.opt/10;
  fg$alpha2.opt=1-fg$alpha1.opt;
  fg$tau3.opt=fg$tau1.opt*2;
  
  
  ub<-c(fo.opt=Inf,fm.opt=Inf,tau1.opt=Inf,alpha1.opt=1,tau2.opt=Inf,alpha2.opt=1,tau3.opt=Inf);
  lb<-c(fo.opt=-Inf,fm.opt=-Inf,tau1.opt=-Inf,alpha1.opt=0,tau2.opt=-Inf,alpha2.opt=0,tau3.opt=-Inf);
  
  #fitting str function 
  fit_str(x,y,fo.opt = fg$fo.opt,fm.opt=fg$fm.opt,wm,d,fg=fg,ub=ub,lb=lb)
  
  #Output object needed for the plotting function
  #it is stored in the Environment e2
  #it would be nice to change the name of this object so that it's different
  #from the sti function
  res<-local(res, env = e2)
  
  
  res$model_str=wm;
  names(str.fit)<-names(res)
  #add res into str.fit data frame in appropriate row
  str.fit[d-1,]<-res
  
  
  #make str plots that will write to already opened pdf file.
  #define parameters again here to make certain they are correct
  x<-data.str$time/1e6; # change time to seconds
  y<-data.str[,d]
  wplot<-'str'
  res<-str.fit[d-1,]
  titl<-names(data.str)[d]
  
  #plot str function
  psi_plots_str_screen(x,y, fit1,fit1u,fit1l,titl,res)
  
  

} #end of str loop through data.str columns
dev.set()


#end of if STR==TRUE
}

###############Output the data generated with names of rows included.


if(rep_names[1]!=F){
    rep_names_long<-rep(rep_names, each=no_light_steps)
}

row.names <- names(data.sti[2:ncol(data.sti)])
lc_in_file<-as.character(t(as.data.frame(strsplit(row.names, split='_')))[,1])
step_in_lc<-as.character(t(as.data.frame(strsplit(row.names, split='_')))[,2])


#a selection step in case the str fit was set to false
#and a decision step in the case the replicate name was defined

if(settings$run_str==TRUE){

if(rep_names[1]!=F){
    results<-cbind(rep_names_long, lc_in_file, step_in_lc, formatted.time, sti.fit, str.fit)
    names(results)[1:4]<-c('replicate','lc_in_file','step_in_lc','measurement_time')
} else {
    results<-cbind(lc_in_file, step_in_lc, formatted.time, sti.fit, str.fit)
    names(results)[1:3]<-c('lc_in_file','step_in_lc','measurement_time')
}
}else{
  if(rep_names[1]!=F){
    results<-cbind(rep_names_long, lc_in_file, step_in_lc, formatted.time, sti.fit)
    names(results)[1:4]<-c('replicate','lc_in_file','step_in_lc','measurement_time')
  } else {
    results<-cbind(lc_in_file, step_in_lc, formatted.time, sti.fit)
    names(results)[1:3]<-c('lc_in_file','step_in_lc','measurement_time')
  }  
  
}




#write out the sti and str fitting results
write.csv(results, file=paste(nam, '_fitted_sti_and_str.csv', sep=''),row.names=F)

#BJ: creating the object on the workspace for further processing, eg. RLC model fitting
fit_fluo<<-results




#write out the meta.data for the run
if(rep_names[1]!=F){
    meta.data<-cbind(rep_names,meta.data)}
write.csv(meta.data, file=paste(nam, '_meta_data.csv',sep=''), row.names=F)



###########################################################################
#calculating light intensity at each light step and returning that object
par_intensity(protocol,no_light_steps,cali,nam)


#creating a ETR and rETR dataframe and exporting it to the workspace
etr_calculations(fit_fluo,no_light_steps)

#creating a NPQ dataframe and exporting it to the workspace
npq_calculations(fit_fluo,no_light_steps)

}

