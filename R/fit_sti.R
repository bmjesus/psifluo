#' @title Fitting function for the induction stage
#' @description this will fit the Kolber et al. 1998 model to the induction step
#' @param x time
#' @param y fluorescence
#' @param fo.opt 
#' @param fm.opt
#' @param sigma.opt
#' @param  rho.opt
#' @param  sigma_max
#' @param  wm
#' @param  sti.fit
#' @param c is the object that has the sti data stored, it can have several columns
#' if there are several measurement
#' @return a value that will be used in the fitting function
#' @keywords external
#' @export
#this is the fitting function for the sti model
fit_sti<-function(x,y,fo.opt,fm.opt,sigma.opt,rho.opt,sigma_max,wm,sti.fit,c){
  
  
#Set starting, lower bound, and upper bound for sigma, fo, fm, rho.  To be used when solving for optimal values
fg<-list(fo.opt=y[1],fm.opt=mean(y[(length(x)-5):length(x)]),sigma.opt=300,rho.opt=0);
lb<-c(fo.opt=-Inf,fm.opt=-Inf,sigma.opt=0,rho.opt=0);
ub<-c(fo.opt=Inf,fm.opt=Inf,sigma.opt=sigma_max,rho.opt=1);

# Try to fit fofmsigp model data and catch errors
lsfit.light <- tryCatch({
  minpack.lm::nlsLM(y ~ psiworx.model.sti(incident=x, sigma=sigma.opt, fo=fo.opt, fm=fm.opt, rho=rho.opt),  start=fg, algorithm="port", upper=ub, lower=lb, trace=F,
        control=nls.control(maxiter=1024) )
},error=function(e){
  ##BJ:it would be nice to have the light step and light curve
  #at which the model was not fitted
  print(paste("Error: could not fit with model ",wm,".",sep=""));
  return(NULL)
});

if(!is.null(lsfit.light)){
  res=as.list(summary(lsfit.light)$coefficients);
} else{
  res<-as.list(numeric(4*4)*NaN);
}
names(res)<-c('fo_sti','fm_sti','sigma_sti','rho_sti','fo_se_sti','fm_se_sti','sigma_se_sti','rho_se_sti','fo_tvalue_sti','fm_tvalue_sti','sigma_tvalue_sti','rho_tvalue_sti','fo_p_sti','fm_p_sti','sigma_p_sti','rho_p_sti');


res$model_sti<-wm;
names(sti.fit)<-names(res)
sti.fit[c-1,]<-res


# Light: Calculate fitted values using optimized parameters
#BJ seems like this also fits the upper and lower intervals using the SE for the different parameters
tryCatch({fit1<-psiworx.model.sti(incident=x, sigma=res$sigma_sti, fo=res$fo_sti, fm=res$fm_sti, rho=res$rho_sti)},error=function(Fit.Fail){return("Fit Fail")})
tryCatch({fit1u<-psiworx.model.sti(incident=x, sigma=res$sigma_sti+res$sigma_se_sti, fo=res$fo_sti+res$fo_se_sti, fm=res$fm_sti+res$fm_se_sti, rho=res$rho_sti+res$rho_se_sti)},error=function(Fit.Fail){return("Fit Fail")})
tryCatch({fit1l<-psiworx.model.sti(incident=x, sigma=res$sigma_sti-res$sigma_se_sti, fo=res$fo_sti-res$fo_se_sti, fm=res$fm_sti-res$fm_se_sti, rho=res$rho_sti-res$rho_se_sti)},error=function(Fit.Fail){return("Fit Fail")})


res<-assign("res",res, envir = e2)
fit1<-assign("fit1",fit1, envir = e2)
fit1u<-assign("fit1u",fit1u, envir = e2)
fit1l<-assign("fit1l",fit1l, envir = e2)


}