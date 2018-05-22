#' @title Initialize constants
#' @description Creates default values for the constants
#' @return list with the constants to use in the calculations
#' @keywords internal
#' @export
initialize_constants <- function() {
  constants<-list()
  constants$photons.in.umol <- 6.02214129e17
  constants$meter.sq.in.angstrom.sq <- 1.0e-20
  return(as.list(constants))
}


#' @title Calculate power per flashlet
#' @description Calculate power per flashlet
#' @return power levels
#' @keywords internal
#' @export
pw_level<- function(protocol){
  constants<-initialize_constants()
  cali<-local(cali, env = e2)

  #it looks for the type of actinic light used to extract the right info from
  #the calibration file

  blue<-meta.data[1,"Blue_light_curves_enable"]
  red<-meta.data[1,"Red_light_curves_enable"]

  if (protocol=='RLC'){
    power.level <- meta.data[1,"F_Voltage"]
    #print(power.level)
    power <- cali[which(cali[,"Power.Level"]==power.level),"Red.Flash"]*meta.data[1,"ST1"]*constants$photons.in.umol*constants$meter.sq.in.angstrom.sq
  }else{
    if (red==1){
      power.level <- meta.data[1,"F_Voltage"]
      #print(power.level)
      power <- cali[which(cali[,"Power.Level"]==power.level),"Red.Flash"]*meta.data[1,"ST1_light"]*constants$photons.in.umol*constants$meter.sq.in.angstrom.sq
    }
    if (blue==1){
      power.level <- meta.data[1,"F2_Voltage"]
      #print(power.level)
      power <- cali[which(cali[,"Power.Level"]==power.level),"Blue.Flash"]*meta.data[1,"ST1_light"]*constants$photons.in.umol*constants$meter.sq.in.angstrom.sq
    }
  }
  #print (power)
  #power<-assign("power",power, envir = e2)
  return(power)
}

#' @title Calculates the light intensity at each light level
#' @description Calculates the light intensity at each light level using the voltage values
#' @return PAR levels
#' @keywords internal
#' @export
par_intensity<-function(protocol,no_light_steps,cali,nam){

  #no_light_steps<-local(no_light_steps,env = e2)

  blue<-meta.data[1,"Blue_light_curves_enable"]
  red<-meta.data[1,"Red_light_curves_enable"]

  if (protocol=='RLC'){
    lc_voltage<-names(meta.data)[which(grepl('Amber',names(meta.data)))]
    lc_voltage<-lc_voltage[which(lc_voltage!='Amber_Flashlet_Voltage')]
    lc_voltage<-lc_voltage[which(lc_voltage!='Amber_RLC_Levels')]
    lc_voltage<-lc_voltage[which(lc_voltage!='Amber_Actinic')]
    lc_voltage<-lc_voltage[which(lc_voltage!='Amber_Flashlet')]
    voltages<-data.frame(step=rep(NA,no_light_steps), voltage=NA, par=NA)
    for(i in 1:no_light_steps){
      voltages$step[i]<-lc_voltage[which(lc_voltage==paste('Amber',i-1,sep='_'))]
      voltages$voltage[i]<-meta.data[1, voltages$step[i]]
      #look up PAR using voltage in calibration file - must already be loaded under 'cali'
      voltages$par[i]<-cali[cali$Power.Level==voltages$voltage[i], paste('Red','.Actinic',sep='')]
      #print(voltages)
    }
  }else{
    if (blue==1){
      lc_voltage<-names(meta.data)[which(grepl('Blue',names(meta.data)))]
      lc_voltage<-lc_voltage[which(lc_voltage!='Blue_light_curves_enable')]
      lc_voltage<-lc_voltage[which(lc_voltage!='Blue_light_curves')]

      voltages<-data.frame(step=rep(NA,no_light_steps/2), voltage=NA, par=NA)
      for(i in 1:(no_light_steps/2)){
        voltages$step[i]<-lc_voltage[which(lc_voltage==paste('Blue_light_curves',i-1,sep='_'))]
        voltages$voltage[i]<-meta.data[1, voltages$step[i]]
        #look up PAR using voltage in calibration file - must already be loaded under 'cali'
        voltages$par[i]<-cali[cali$Power.Level==voltages$voltage[i], paste('Blue','.Actinic',sep='')]
        #print(voltages)
      }
    }

    if (red==1){
      lc_voltage<-names(meta.data)[which(grepl('Red',names(meta.data)))]
      lc_voltage<-lc_voltage[which(lc_voltage!='Red_light_curves_enable')]
      lc_voltage<-lc_voltage[which(lc_voltage!='Red_light_curves')]

      voltages<-data.frame(step=rep(NA,no_light_steps/2), voltage=NA, par=NA)
      for(i in 1:(no_light_steps/2)){
        voltages$step[i]<-lc_voltage[which(lc_voltage==paste('Red_light_curves',i-1,sep='_'))]
        voltages$voltage[i]<-meta.data[1, voltages$step[i]]
        #look up PAR using voltage in calibration file - must already be loaded under 'cali'
        voltages$par[i]<-cali[cali$Power.Level==voltages$voltage[i], paste('Red','.Actinic',sep='')]
        #print(voltages)
      }
    }

  }



  #write out voltages data frame for use in other analyses later on
  write.csv(voltages, file=paste(nam,'_voltages_par.csv',sep=''), row.names=F)

  #BJ: creating the object on the workspace for further processing, eg. RLC model fitting
  volt<<-voltages



}

#' @title Calculates ETR and rETR for each PAR
#' @description Calculates ETR and rETR for each PAR
#' @return object with ETR and rETR
#' @keywords internal
#' @export
etr_calculations<-function(fit_fluo,no_light_steps){

  # initialize data frame
  etr_names<-c('PAR','sigma','eff','rETR','ETR')
  etr_fit <- as.data.frame(matrix(nrow=(nrow(fit_fluo)),ncol=length(etr_names)))
  names(etr_fit)=etr_names
  rm(etr_names)

  #store the Fo values
  f<-numeric()
  fo<-numeric()
  a<-1
  if (nrow(fit_fluo)>no_light_steps){
    for (i in 1:(nrow(fit_fluo)/no_light_steps)){
      f<-fit_fluo$fo_sti[a]
      fo<-c(fo,f)
      a<-a+no_light_steps
    }
  }else{
    fo<-fit_fluo$fo_sti[1]
  }

  #store the Fm values
  f<-numeric()
  fm<-numeric()
  a<-1
  if (nrow(fit_fluo)>no_light_steps){
    for (i in 1:(nrow(fit_fluo)/no_light_steps)){
      f<-fit_fluo$fm_sti[a]
      fm<-c(fm,f)
      a<-a+no_light_steps
    }
  }else{
    fm<-fit_fluo$fm_sti[1]
  }

  #storing the light values in the etr.fit object
  etr_fit$PAR<-volt$par

  #storing the sigma values in the etr.fit object
  etr_fit$sigma<-fit_fluo$sigma_sti

  #PSII efficiency
  etr_fit$eff<-(fit_fluo$fm_sti - fit_fluo$fo_sti)/fit_fluo$fm_sti

  #rETR
  etr_fit$rETR<-etr_fit$eff*0.5*etr_fit$PAR

  #ETR
  #ETR = PAR * sigma * ((F'm-F)/F'm)*(Fm/(Fm-F0))*6.022/1000
  #check this equation with Doug
  #for the moment this loop only applies if there's no replicates
  #i.e. 1 set of light and 1 set of dark measurements
  #need to test it with a file with several replicates to see how to develop this
  #to incorporate several replicates.

  for (i in 1:nrow(fit_fluo)){
    if (i <= no_light_steps){
      a<-1
      etr_fit$ETR[i]<- etr_fit$PAR[i] * etr_fit$sigma[i] * ((fit_fluo$fm_sti[i]-fit_fluo$fo_sti[i])/fit_fluo$fm_sti[i])*(fm[a]/(fm[a]-fo[a]))*6.022/1000
    }else{
      a<-2
      etr_fit$ETR[i]<- etr_fit$PAR[i] * etr_fit$sigma[i] * ((fit_fluo$fm_sti[i]-fit_fluo$fo_sti[i])/fit_fluo$fm_sti[i])*(fm[a]/(fm[a]-fo[a]))*6.022/1000
    }
  }
  etr_fit<<-etr_fit
}




#npq_calculations


#' @title Calculate NPQ
#' @description Calculate several parameters related to NPQ
#' @return object with NPQ
#' @keywords internal
#' @export
#'
npq_calculations<-function(fit_fluo,no_light_steps){

  # initialize data frame
  npq_names<-c('PAR','NPQ','YNPQ','YNO')
  npq_fit <- as.data.frame(matrix(nrow=(nrow(fit_fluo)),ncol=length(npq_names)))
  names(npq_fit)=npq_names
  rm(npq_names)

  #store the Fm values, in Doug's protocol there are always two Fm values
  #one for light and another for 2 seconds after dark
  f<-numeric()
  fm<-numeric()
  a<-1
  if (nrow(fit_fluo)>no_light_steps){
    for (i in 1:(nrow(fit_fluo)/no_light_steps)){
      f<-fit_fluo$fm_sti[a]
      fm<-c(fm,f)
      a<-a+no_light_steps
    }
  }else{
    fm<-fit_fluo$fm_sti[1]
  }


  #store the F values
  f<-numeric()
  fo<-numeric()
  a<-1
  if (nrow(fit_fluo)>no_light_steps){
    for (i in 1:(nrow(fit_fluo)/no_light_steps)){
      f<-fit_fluo$fo_sti[a]
      fo<-c(fo,f)
      a<-a+no_light_steps
    }
  }else{
    fo<-fit_fluo$fo_sti[1]
  }


  #storing the light values in the npq.fit object
  npq_fit$PAR<-volt$par

  #storing the NPQ values in the npq.fit object
  for (i in 1:nrow(fit_fluo)){
    if (i <= no_light_steps){
      a<-1
      npq_fit$NPQ[i]<-(fm[a]-fit_fluo$fm_sti[i])/fit_fluo$fm_sti[i]
    }else{
      a<-2
      npq_fit$NPQ[i]<-(fm[a]-fit_fluo$fm_sti[i])/fit_fluo$fm_sti[i]
    }
  }

  #storing the Y(NPQ) values in the npq.fit object
  #Y(NPQ)=(F/Fm')-(F/Fm)
  for (i in 1:nrow(fit_fluo)){
    if (i <= no_light_steps){
      a<-1
      npq_fit$YNPQ[i]<-(fit_fluo$fo_sti[i]/fit_fluo$fm_sti[i])-(fit_fluo$fo_sti[i]/fm[a])
    }else{
      a<-2
      npq_fit$YNPQ[i]<-(fit_fluo$fo_sti[i]/fit_fluo$fm_sti[i])-(fit_fluo$fo_sti[i]/fm[a])
    }
  }

  #storing the Y(NO) values in the npq.fit object
  #Y(NO)=F/Fm
  for (i in 1:nrow(fit_fluo)){
    if (i <= no_light_steps){
      a<-1
      npq_fit$YNO[i]<-(fit_fluo$fo_sti[i]/fm[a])
    }else{
      a<-2
      npq_fit$YNO[i]<-(fit_fluo$fo_sti[i]/fm[a])
    }
  }


  #exporting the final table to the workspace
  npq_fit<<-npq_fit
}

