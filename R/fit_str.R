#this is the fitting function for the str model
fit_str<-function(x,y,fo.opt,fm.opt,wm,d,fg,ub,lb){
  
# Try to fit str model to data and catch errors
switch(wm,
       tau1={
         rlfit.light <- tryCatch({
          minpack.lm::nlsLM(y ~psiworx.model.str(x=x, fo=fo.opt, fm=fm.opt, tau1=tau1.opt, wm=wm),  start=fg[1:3], algorithm="port", upper=ub[1:3], lower=lb[1:3], trace=F,
                 control=nls.control(maxiter=1024) )
         },error=function(e){
           print(paste("Error: could not fit with model ",wm," at dataset = ",names(data.str)[d],sep=""));
           return(NULL)
         })
         # parse fit results
         if(!is.null(rlfit.light)){
           res=as.list(summary(rlfit.light)$coefficients);
         } else{
           res<-as.list(numeric(3*4)*NaN);
         }
         names(res)<-c('fo_str','fm_str','tau1_str','fo_se_str','fm_se_str','tau1_se_str','fo_tvalue_str','fm_tvalue_str','tau1_tvalue_str','fo_p_str','fm_p_str','tau1_p_str');
         res<-c(res,alpha1_str=NaN,tau2_str=NaN,alpha2_str=NaN,tau3_str=NaN,alpha1_se_str=NaN,tau2_se_str=NaN,alpha2_se_str=NaN,tau3_se_str=NaN,alpha1_tvalue_str=NaN,tau2_tvalue_str=NaN,alpha2_tvalue_str=NaN,tau3_tvalue_str=NaN,alpha1_p_str=NaN,tau2_p_str=NaN,alpha2_p_str=NaN,tau3_p_str=NaN);
       },
       
       tau2={
         rlfit.light <- tryCatch({
           minpack.lm::nlsLM(y ~psiworx.model.str(x=x, fo=fo.opt, fm=fm.opt, tau1=tau1.opt, alpha1=alpha1.opt, tau2=tau2.opt, wm=wm),  start=fg[1:5], algorithm="port", upper=ub[1:5], lower=lb[1:5], trace=F,
                 control=nls.control(maxiter=1024) )
         },error=function(e){
           print(paste("Error: could not fit with model ",wm," at dataset = ",names(data.str)[d],sep=""));
           return(NULL)
         })
         # parse fit results
         if(!is.null(rlfit.light)){
           res=as.list(summary(rlfit.light)$coefficients);
         } else{
           res<-as.list(numeric(5*4)*NaN);
         }
         names(res)<-c('fo_str','fm_str','tau1_str','alpha1_str','tau2_str','fo_se_str','fm_se_str','tau1_se_str','alpha1_se_str','tau2_se_str','fo_tvalue_str','fm_tvalue_str','tau1_tvalue_str','alpha1_tvalue_str','tau2_tvalue_str','fo_p_str','fm_p_str','tau1_p_str','alpha1_p_str','tau2_p_str');
         res<-c(res,alpha2_str=NaN,tau3_str=NaN,alpha2_se_str=NaN,tau3_se_str=NaN,alpha2_tvalue_str=NaN,tau3_tvalue_str=NaN,alpha2_p_str=NaN,tau3_p_str=NaN);
       },
       
       tau3={
         rlfit.light <- tryCatch({
           minpack.lm::nlsLM(y ~psiworx.model.str(x=x, fo=fo.opt, fm=fm.opt, tau1=tau1.opt, alpha1=alpha1.opt, tau2=tau2.opt, alpha2=alpha2.opt, tau3=tau3.opt, wm=wm),  start=fg, algorithm="port", upper=ub, lower=lb, trace=F,
                 control=nls.control(maxiter=1024) )
         },error=function(e){
           print(paste("Error: could not fit with model ",wm," at dataset = ",names(data.str)[d],sep=""));
           return(NULL)
         })
         # parse fit results
         if(!is.null(rlfit.light)){
           res=as.list(summary(rlfit.light)$coefficients);
         } else{
           res<-as.list(numeric(7*4)*NaN);
         }
         names(res)<-c('fo_str','fm_str','tau1_str','alpha1_str','tau2_str','alpha2_str','tau3_str','fo_se_str','fm_se_str','tau1_se_str','alpha1_se_str','tau2_se_str','alpha2_se_str','tau3_se_str','fo_tvalue_str','fm_tvalue_str','tau1_tvalue_str','alpha1_tvalue_str','tau2_tvalue_str','alpha2_tvalue_str','tau3_tvalue_str','fo_p_str','fm_p_str','tau1_p_str','alpha1_p_str','tau2_p_str','alpha2_p_str','tau3_p_str');
       }
); # end switch
  

res<-assign("res",res, envir = e2)
  
}