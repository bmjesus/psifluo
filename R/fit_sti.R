#' @title Fitting function for a FRRf induction curve
#' @description This will fit a model to the induction phase using the equations of Zbigniew S. Kolber, Ondrej Prasil, Paul G. Falkowski (1998) Measurements of variable chlorophyll Fuorescence using fast repetition rate techniques: defining methodology and experimental protocols. Biochimica et Biophysica Acta, 1367: 88-106.
#' @param x Flashlet energy (photons . A^-2 .  flashlet{^-1}) used to generate a x axis of cumulative energy by multiplying by flashlet count. A typical value from a flashlet calibration file of Doug Campbell's lab is 0.0002167971 (photons . A^-2 . flashlet^-1). This value is wrong if you are not using his machine and you should supply the correct value for your machine if you want to obtain correct sigma values (A^2 . PSII^{-1}).
#' @param y Fluorescence data vector.
#' @param fo Initial value for fo, default (NULL) is the 1st point of the fluorescence data vector y.
#' @param fm Initial value for fm, default (NULL) is the average of the maximum 5 points of the fluorescence data vector y.
#' @param sigma.opt Initial value for sigma, default = 300 A^2 PSII^-1.
#' @param  rho.opt Initial value for rho excitonic connectivity among PSII units, default = 0.1.
#' @param  fit_model How many parameters to fit to the curve. Default is "all". Options are: "all", "sigma_rho" and "sigma". "All" - fo, fm, sigma and rho are fitted; "sigma_rho" -  sigma and rho are fitted, fo and fm are fixed at user-defined levels, defaults are fo =  1st point of the fluorescence data vector y, fm = average of the maximum 5 points of the fluorescence data vector y;  "sigma" - fo, fm and rho are fixed, only sigma is fitted.
#' @param  plots Turns on and off plotting of data and fitted curves (TRUE, FALSE). Default is TRUE.
#' @return The function returns a list with: 1- fitted parameters and respective SE, p-values (list); 2 - flashlet energy (single element numeric vector); 3- fitted values (numeric vector); 4 - cumulative energy (numeric vector).
#' @keywords external
#' @examples
#'
#' fit_sti(x = 0.0002167971, y = sti_data$light_1,fit_model = 'all')
#'
#' fit_sti(x = 0.0002167971, y = sti_data$light_1,fit_model = 'sigma_rho')
#' @export
#'
fit_sti<-function(x
                  ,y
                  ,fo = NULL
                  ,fm = NULL
                  ,sigma.opt = 300
                  ,rho.opt = 0.1
                  ,fit_model = 'all'
                  ,plots = TRUE){

#setting up default Fo and Fm values if the user defines none
  if(is.null(fo)){
    fo.opt=y[1]
  }else{fo.opt=fo}

#Fm using the last 5 points in the curve
#if(is.null(fm)){
#fm.opt=mean(y[(length(y)-5):length(y)])
#}else{fm.opt=fm}

#Fm using the highest 5 points in the curve
if(is.null(fm)){
fm.opt=mean(utils::tail(sort(y),5))
}else{fm.opt=fm}

#print(fm.opt)



  #storing the value of flashlet energy for the output list
  flashlet_energy<-x

  #creating a variable that is the flashlet energy at each fluorescence measurement
  x<-rep(x,times=length(y))

  #storing the value of the cumulative energy for the output list
  cumulative_energy<-numeric()
  for (i in 1:length(x)){
    if (i==1)
    {
      cumulative_energy[i]<-flashlet_energy
    }else{
      cumulative_energy[i]<-cumulative_energy[i-1]+flashlet_energy
    }
  }
  #cumulative_energy<<-cumulative_energy


  #starting values for the model fitting all parameters (fo, fm, sigma and rho)
  fg<-list(fo.opt=fo.opt,fm.opt=fm.opt,sigma.opt=sigma.opt,rho.opt=0);
  #boundaries
  lb<-c(fo.opt=fo.opt*0.5,fm.opt=0.6*fm.opt,sigma.opt=50,rho.opt=0);
  ub<-c(fo.opt=fo.opt*1.5,fm.opt=fo.opt/0.2,sigma.opt=1500,rho.opt=1);

  #starting values for the model fitting sigma and rho
  fg2<-list(sigma.opt=sigma.opt,rho.opt=rho.opt);
  #boundaries
  lb2<-c(sigma.opt=50,rho.opt=0);
  ub2<-c(sigma.opt=1500,rho.opt=1);

  #starting values for the model fitting sigma
  fg3<-list(sigma.opt=sigma.opt);
  lb3<-c(sigma.opt=50);
  ub3<-c(sigma.opt=1500);

  ##########################################################################
  #####section to fit the Kolber model using different options


  switch(fit_model,all={
  #first option
  #Fo,Fm,Sigma and Rho are fitted

  lsfit.light <- tryCatch({
    minpack.lm::nlsLM(y ~ kolber_sti(incident=x, sigma=sigma.opt, fo=fo.opt, fm=fm.opt, rho=rho.opt),  start=fg, algorithm="port", upper=ub, lower=lb, trace=F,
                      control=stats::nls.control(maxiter=1024) )
  },error=function(e){
    print(paste("Error: could not fit one induction dataset with model ",fit_model,".",sep=""));
    return(NULL)
  });

  }
  ,sigma_rho={

  #Fo and Fm fixed, Sigma and rho fitted
  lsfit.light <- tryCatch({
    minpack.lm::nlsLM(y ~ kolber_sti(incident=x, sigma=sigma.opt, fo=fo.opt, fm=fm.opt, rho=rho.opt),  start=fg2, algorithm="port", upper=ub2, lower=lb2, trace=F,control=stats::nls.control(maxiter=1024) )
  },error=function(e){
    print(paste("Error: could not fit one induction dataset with model ",fit_model,".",sep=""));
    return(NULL)
  });
  }
  ,sigma={
    #Fo, Fm and rho fixed, Sigma fitted
    lsfit.light <- tryCatch({
      minpack.lm::nlsLM(y ~ kolber_sti(incident=x, sigma=sigma.opt, fo=fo.opt, fm=fm.opt, rho=rho.opt),  start=fg3, algorithm="port", upper=ub3, lower=lb3, trace=F,control=stats::nls.control(maxiter=1024) )
    },error=function(e){
      print(paste("Error: could not fit one induction dataset with model ",fit_model,".",sep=""));
      return(NULL)
    });
  }

  #end of switch section
  )


  ##########################################################################
  #this section stores the fitted parameters and respective errors in a list
  #it needs 1 section for each fitting model (all, sigma_rho, sigma)

  switch(fit_model,all={
  if(!is.null(lsfit.light)){
    res=as.list(summary(lsfit.light)$coefficients);
  } else{
    res<-as.list(numeric(4*4)*NaN);
  }
  names(res)<-c('fo_sti','fm_sti','sigma_sti','rho_sti','fo_se_sti','fm_se_sti','sigma_se_sti','rho_se_sti','fo_tvalue_sti','fm_tvalue_sti','sigma_tvalue_sti','rho_tvalue_sti','fo_p_sti','fm_p_sti','sigma_p_sti','rho_p_sti');
  }
  ,sigma_rho={
    if(!is.null(lsfit.light)){
      res=as.list(c(fo.opt,fm.opt,
        summary(lsfit.light)$coefficients[1],
        summary(lsfit.light)$coefficients[2],
        NaN,NaN,
        summary(lsfit.light)$coefficients[3],
        summary(lsfit.light)$coefficients[4],
        NaN,NaN,
        summary(lsfit.light)$coefficients[5],
        summary(lsfit.light)$coefficients[6],
        NaN,NaN,
        summary(lsfit.light)$coefficients[7],
        summary(lsfit.light)$coefficients[8]
        ))
      #print(res)
    } else{
      res<-as.list(numeric(4*4)*NaN);
    }
    names(res)<-c('fo_sti','fm_sti','sigma_sti','rho_sti','fo_se_sti','fm_se_sti','sigma_se_sti','rho_se_sti','fo_tvalue_sti','fm_tvalue_sti','sigma_tvalue_sti','rho_tvalue_sti','fo_p_sti','fm_p_sti','sigma_p_sti','rho_p_sti');
  }
  ,sigma={
    if(!is.null(lsfit.light)){
      res=as.list(c(fo.opt,fm.opt,
                    summary(lsfit.light)$coefficients[1],
                    rho.opt,
                    NaN,NaN,
                    summary(lsfit.light)$coefficients[2],
                    NaN,
                    NaN,NaN,
                    summary(lsfit.light)$coefficients[3],
                    NaN,
                    NaN,NaN,
                    summary(lsfit.light)$coefficients[4],
                    NaN
      ))
      #print(res)
    } else{
      res<-as.list(numeric(4*4)*NaN);
    }
    names(res)<-c('fo_sti','fm_sti','sigma_sti','rho_sti','fo_se_sti','fm_se_sti','sigma_se_sti','rho_se_sti','fo_tvalue_sti','fm_tvalue_sti','sigma_tvalue_sti','rho_tvalue_sti','fo_p_sti','fm_p_sti','sigma_p_sti','rho_p_sti');
  }

  )
  ##########################################################################
#print(res)

#this calculates the fitted values using the parameters that were estimated above
fitted_values<-tryCatch({fit1<-kolber_sti(incident=x, sigma=res$sigma_sti, fo=res$fo_sti, fm=res$fm_sti, rho=res$rho_sti)},error=function(Fit.Fail){return("Fit Fail")})

##########################################################################
#Optional plot
if (plots==TRUE){
graphics::plot(cumulative_energy,y,ylab="Fluorescence",xlab=expression(paste("Cumulative energy (photons ",ring(A)^{-2},")")),las=1)
  graphics::points(cumulative_energy,fitted_values,col=2,type='l')



#this is for the legend at the bottom right
# Light: Print optimized parameter values on the plot ####
leg.L.model<-paste("Model: ",fit_model,sep="");
leg.L.main<-"Parameter \u00B1 SE" # \u00B1 is the unicode character for the plus/minus symbol
leg.L.Fo<-tryCatch({paste("Fo =", round(res$fo_sti, digits=3), "\u00B1", round(res$fo_se_sti, digits=3), sep=" ")},error=function(Fit.Fail){return("Fit Fail")})
leg.L.Fm<-tryCatch({paste("Fm =", round(res$fm_sti, digits=3), "\u00B1", round(res$fm_se_sti, digits=3), sep=" ")},error=function(Fit.Fail){return("Fit Fail")})
leg.L.sigma<-tryCatch({paste("sigma =", round(res$sigma_sti, digits=0), "\u00B1", round(res$sigma_se_sti, digits=0), sep=" ")},error=function(Fit.Fail){return("Fit Fail")})
leg.L.rho<-tryCatch({paste("rho =", round(res$rho_sti, digits=3), "\u00B1", round(res$rho_se_sti, digits=3), sep=" ")},error=function(Fit.Fail){return("Fit Fail")})
graphics::legend("bottomright",legend=c(leg.L.model,leg.L.main,leg.L.Fo,leg.L.Fm,leg.L.sigma,leg.L.rho),bty="n",text.font=c(2,2,1,1,1,1),adj=c(0, 0.5))
}

##########################################################################
##########################################################################

#calculating PSII quantum efficiency and adding it to the fitted_parameters dataframe
#print(fitted_values$fo_sti)

psII_eff<-(res$fm_sti-res$fo_sti)/res$fm_sti

res$psII_eff_sti<-psII_eff


##########################################################################
##########################################################################

#Output list
output<-list(flashlet_energy,res,fitted_values,cumulative_energy)
names(output)<-c('flashlet_energy','fitted_parameters','fitted_values','cumulative_energy')
return(output)



}
