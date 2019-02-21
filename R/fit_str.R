#' @title Fitting function for a FRRf re-opening curve
#' @description  This will fit a model to the re-opening phase using the equations of Zbigniew S. Kolber, Ondrej Prasil, Paul G. Falkowski (1998) Measurements of variable chlorophyll Fuorescence using fast repetition rate techniques: defining methodology and experimental protocols. Biochimica et Biophysica Acta, 1367: 88-106.
#' @param x Vector with cumulative time in microseconds
#' @param y Fluorescence vector
#' @param fo Initial value for fo, if NULL then it defaults to the last point of the relaxation curve
#' @param fm Initial value for fm, if NULL then it defaults to the first point of the relaxation curve
#' @param  tau_model Options: "tau1", "tau2", "tau3" and "tau1_subset" (the user selects a subset of the data starting from the first point of the re-opening curve). Default is tau1, fitting re-opening with a single phase exponential decay with lifetime of tau1.
#' @param  subset_time Amount of time to fit the tau1_subset model (default = 500 us)
#' @param  plots Turns on and off the plotting (TRUE, FALSE). Default is TRUE
#' @return The function returns a list with: 1 - fitted parameters (list); 2 - fitted values (numeric vector); 3 - time in us (numeric vector)
#' @keywords external
#' @examples
#'
#' fit_str(x = str_data$time, y = str_data$light_1,tau_model = 'tau1')
#'
#' fit_str(x = str_data$time, y = str_data$light_1,tau_model = 'tau3')
#' @export

#this is the fitting function for the str model
fit_str<-function(x,
                  y,
                  fo = NULL,
                  fm = NULL,
                  tau_model = "tau1",
                  subset_time = 500,
                  plots = TRUE){


  #setting up default values if the user defines none
  if(is.null(fo)){
    fo.opt=y[length(y)]
  }else{fo.opt=fo}

  if(is.null(fm)){
    fm.opt=y[1]
  }else{fm.opt=fm}

#creating subset for fitting 'tau1_subset' model
#based on time
dt_str<-as.data.frame(cbind(x,y))
#print(dt_str)
x_subset<-dt_str$x[dt_str[,1]<=subset_time]
y_subset<-y[1:length(x_subset)]

# change time to seconds
x<-x/1e6;
x_subset<-x_subset/1e6;


# starting values for tau1
  fg=list();
  fg$fo.opt=min(y);
  fg$fm.opt=max(y);
  fg$tau1.opt=0.001;
#boundaries for tau1
lb<-c(fo.opt=fg$fo.opt*0.6,fm.opt=fg$fm.opt*0.7,
      tau1.opt=0)
ub<-c(fo.opt=fg$fo.opt*1.2,fm.opt=fg$fo.opt/0.2,
      tau1.opt=Inf)


  # starting values for tau2
  fg2=list();
  fg2$fo.opt=min(y);
  fg2$fm.opt=max(y);
  fg2$tau1.opt=0.001;
  fg2$alpha1.opt=0.6;
  fg2$tau2.opt=fg$tau1.opt*10;
  #boundaries for tau2
  lb2<-c(fo.opt=fg2$fo.opt*0.6,fm.opt=fg2$fm.opt*0.7,
        tau1.opt=0,alpha1.opt=0,tau2.opt=0)
  ub2<-c(fo.opt=fg2$fo.opt*1.2,fm.opt=fg2$fo.opt/0.2,
        tau1.opt=Inf,alpha1.opt=Inf, tau2.opt=Inf)


  # starting values for tau3
  fg3=list();
  fg3$fo.opt=min(y);
  fg3$fm.opt=max(y);
  fg3$tau1.opt=0.001;
  fg3$alpha1.opt=0.6;
  fg3$tau2.opt=fg3$tau1.opt*10;
  fg3$alpha2.opt=1-fg3$alpha1.opt;
  fg3$tau3.opt=fg3$tau1.opt*100;
  #boundaries for tau3
  lb3<-c(fo.opt=fg3$fo.opt*0.6,fm.opt=fg3$fm.opt*0.7,
         tau1.opt=0,alpha1.opt=-Inf,tau2.opt=0,
         alpha2.opt=0,tau3.opt=0)
  ub3<-c(fo.opt=fg3$fo.opt*1.2,fm.opt=fg3$fo.opt/0.2,
         tau1.opt=Inf,alpha1.opt=Inf, tau2.opt=Inf,
         alpha2.opt=Inf,tau3.opt=Inf)


  # starting values for tau1_subset
  fg_subset=list();
  fg_subset$tau1.opt=0.0005;
  #boundaries for tau1_subset
  lb_subset<-c(tau1.opt=0)
  ub_subset<-c(tau1.opt=Inf)




# Try to fit str model to data and catch errors
#outputs a list named "res" with the fitted parameters
#the "res" order is the same for all models (values,SE,t-values,p-values)
switch(tau_model,
       tau1={
         rlfit.light <- tryCatch({
          minpack.lm::nlsLM(y ~kolber_str(x=x, fo=fo.opt, fm=fm.opt,tau1 = tau1.opt, tau_model="tau1"),  start=fg, algorithm="port",upper=ub,lower=lb,trace=F,
                 control = stats::nls.control(maxiter=1024) )
         },error=function(e){
           print(paste("Error: could not fit one re-opening dataset with model ",tau_model,sep=""));
           return(NULL)
         })
         # parse fit results
         if(!is.null(rlfit.light)){
           res=as.list(summary(rlfit.light)$coefficients);
         } else{
           res<-as.list(numeric(3*4)*NaN);
         }
         #name the list elements to facilitate selection
         names(res)<-c('fo_str','fm_str','tau1_str','fo_se_str','fm_se_str','tau1_se_str','fo_tvalue_str','fm_tvalue_str','tau1_tvalue_str','fo_p_str','fm_p_str','tau1_p_str');
         #reorder to be the same in all models
         res<-as.list(c(res$fo_str,res$fm_str,res$tau1_str,alpha1_str=NaN,tau2_str=NaN,alpha2_str=NaN,tau3_str=NaN,
           res$fo_se_str,res$fm_se_str,res$tau1_se_str,alpha1_se_str=NaN,tau2_se_str=NaN,alpha2_se_str=NaN,tau3_se_str=NaN,res$fo_tvalue_str,res$fm_tvalue_str,res$tau1_tvalue_str,alpha1_tvalue_str=NaN,tau2_tvalue_str=NaN,alpha2_tvalue_str=NaN,tau3_tvalue_str=NaN,res$fo_p_str,res$fm_p_str,res$tau1_p_str,alpha1_p_str=NaN,tau2_p_str=NaN,alpha2_p_str=NaN,tau3_p_str=NaN))
         #rename the new order
         names(res)[c(1,2,3,8,9,10,15,16,17,22,23,24)]<-c('fo_str','fm_str','tau1_str','fo_se_str','fm_se_str','tau1_se_str','fo_tvalue_str','fm_tvalue_str','tau1_tvalue_str','fo_p_str','fm_p_str','tau1_p_str')
         #adding alpha3  to the dataframe
         res$alpha3_str=NaN
         res$alpha3_se_str=NaN
         },

       tau2={
         rlfit.light <- tryCatch({
           minpack.lm::nlsLM(y ~kolber_str(x=x, fo=fo.opt, fm=fm.opt, tau1=tau1.opt, alpha1=alpha1.opt, tau2=tau2.opt,tau_model="tau2"),  start=fg2, algorithm="port",upper=ub2,lower=lb2 ,trace=F,
                 control = stats::nls.control(maxiter=1024) )
         },error=function(e){
           print(paste("Error: could not fit one re-opening dataset with model ",tau_model,sep=""));
           return(NULL)
         })
         # parse fit results
         if(!is.null(rlfit.light)){
           res=as.list(summary(rlfit.light)$coefficients);
         } else{
           res<-as.list(numeric(5*4)*NaN);
         }
         #name the list elements to facilitate selection
         names(res)<-c('fo_str','fm_str','tau1_str','alpha1_str','tau2_str','fo_se_str','fm_se_str','tau1_se_str','alpha1_se_str','tau2_se_str','fo_tvalue_str','fm_tvalue_str','tau1_tvalue_str','alpha1_tvalue_str','tau2_tvalue_str','fo_p_str','fm_p_str','tau1_p_str','alpha1_p_str','tau2_p_str');
         #reorder to be the same in all models
         res<-as.list(c(res$fo_str,res$fm_str,res$tau1_str,alpha1_str=res$alpha1_str,tau2_str=res$tau2_str,alpha2_str=NaN,tau3_str=NaN,
res$fo_se_str,res$fm_se_str,res$tau1_se_str,res$alpha1_se_str,res$tau2_se_str,alpha2_se_str=NaN,tau3_se_str=NaN,res$fo_tvalue_str,res$fm_tvalue_str,res$tau1_tvalue_str,res$alpha1_tvalue_str,res$tau2_tvalue_str,alpha2_tvalue_str=NaN,tau3_tvalue_str=NaN,res$fo_p_str,res$fm_p_str,res$tau1_p_str,res$alpha1_p_str,res$tau2_p_str,alpha2_p_str=NaN,tau3_p_str=NaN))
         #rename the new order
         names(res)[c(1,2,3,4,5,8,9,10,11,12,15,16,17,18,19,22,23,24,25,26)]<-c('fo_str','fm_str','tau1_str','alpha1_str','tau2_str','fo_se_str','fm_se_str','tau1_se_str','alpha1_se_str','tau2_se_str','fo_tvalue_str','fm_tvalue_str','tau1_tvalue_str','alpha1_tvalue_str','tau2_tvalue_str','fo_p_str','fm_p_str','tau1_p_str','alpha1_p_str','tau2_p_str')
         #adding alpha2 and SE to the dataframe
         res$alpha2_str=1-res$alpha1_str
         res$alpha2_se_str=res$alpha1_se_str
         #adding alpha3  to the dataframe
         res$alpha3_str=NaN
         res$alpha3_se_str=NaN
       },

       tau3={
         rlfit.light <- tryCatch({
           minpack.lm::nlsLM(y ~kolber_str(x=x, fo=fo.opt, fm=fm.opt, tau1=tau1.opt, alpha1=alpha1.opt, tau2=tau2.opt, alpha2=alpha2.opt, tau3=tau3.opt,tau_model="tau3"),  start=fg3, upper=ub3,lower = lb3,
                             algorithm="port", trace=F,
                 control = stats::nls.control(maxiter=1024) )
         },error=function(e){
           print(paste("Error: could not fit one re-opening dataset with model ",tau_model,names(y),sep=""));
           return(NULL)
         })
         # parse fit results
         if(!is.null(rlfit.light)){
           res=as.list(summary(rlfit.light)$coefficients);
         } else{
           res<-as.list(numeric(7*4)*NaN);
         }
         #name the list elements to facilitate selection
         names(res)<-c('fo_str','fm_str','tau1_str','alpha1_str','tau2_str','alpha2_str','tau3_str','fo_se_str','fm_se_str','tau1_se_str','alpha1_se_str','tau2_se_str','alpha2_se_str','tau3_se_str','fo_tvalue_str','fm_tvalue_str','tau1_tvalue_str','alpha1_tvalue_str','tau2_tvalue_str','alpha2_tvalue_str','tau3_tvalue_str','fo_p_str','fm_p_str','tau1_p_str','alpha1_p_str','tau2_p_str','alpha2_p_str','tau3_p_str');
         #reorder to be the same in all models
         res<-as.list(c(res$fo_str,res$fm_str,res$tau1_str,alpha1_str=res$alpha1_str,tau2_str=res$tau2_str,res$alpha2_str,res$tau3_str,res$fo_se_str,res$fm_se_str,res$tau1_se_str,res$alpha1_se_str,res$tau2_se_str,res$alpha2_se_str,res$tau3_se_str,res$fo_tvalue_str,res$fm_tvalue_str,res$tau1_tvalue_str,res$alpha1_tvalue_str,res$tau2_tvalue_str,res$alpha2_tvalue_str,res$tau3_tvalue_str,res$fo_p_str,res$fm_p_str,res$tau1_p_str,res$alpha1_p_str,res$tau2_p_str,res$alpha2_p_str,res$tau3_p_str))
         #rename the new order
         names(res)<-c('fo_str','fm_str','tau1_str','alpha1_str','tau2_str','alpha2_str','tau3_str','fo_se_str','fm_se_str','tau1_se_str','alpha1_se_str','tau2_se_str','alpha2_se_str','tau3_se_str','fo_tvalue_str','fm_tvalue_str','tau1_tvalue_str','alpha1_tvalue_str','tau2_tvalue_str','alpha2_tvalue_str','tau3_tvalue_str','fo_p_str','fm_p_str','tau1_p_str','alpha1_p_str','tau2_p_str','alpha2_p_str','tau3_p_str')
         #adding alpha3  to the dataframe
         res$alpha3_str=1-(res$alpha1_str+res$alpha2_str)
         res$alpha3_se_str=NaN
       },
       tau1_subset={
         rlfit.light <- tryCatch({
           minpack.lm::nlsLM(y_subset ~subset_str(x=x_subset, fo=fg$fo.opt, fm=fg$fm.opt,tau1 = tau1.opt),  start=fg_subset,upper = ub_subset,lower=lb_subset, algorithm="port", trace=F,
                             control = stats::nls.control(maxiter=1024) )
         },error=function(e){
           print(paste("Error: could not fit one re-opening dataset with model ",tau_model,sep=""));
           return(NULL)
         })
         # parse fit results
         if(!is.null(rlfit.light)){
           res=as.list(summary(rlfit.light)$coefficients);
         } else{
           res<-as.list(numeric(3*4)*NaN);
         }
         #name the list elements to facilitate selection
         names(res)<-c('tau1_str','tau1_se_str','tau1_tvalue_str','tau1_p_str');
         #reorder to be the same in all models
         res<-as.list(c(fg$fo.opt,fg$fm.opt,res$tau1_str,alpha1_str=NaN,tau2_str=NaN,alpha2_str=NaN,tau3_str=NaN,fo_se_str=NaN,fm_se_str=NaN,res$tau1_se_str,alpha1_se_str=NaN,tau2_se_str=NaN,alpha2_se_str=NaN,tau3_se_str=NaN,fo_tvalue_str=NaN,fm_tvalue_str=NaN,res$tau1_tvalue_str,alpha1_tvalue_str=NaN,tau2_tvalue_str=NaN,alpha2_tvalue_str=NaN,tau3_tvalue_str=NaN,fo_p_str=NaN,fm_p_str=NaN,res$tau1_p_str,alpha1_p_str=NaN,tau2_p_str=NaN,alpha2_p_str=NaN,tau3_p_str=NaN));
         #rename the new order
         names(res)[c(1,2,3,10,17,24)]<-c('fo_str','fm_str','tau1_str','tau1_se_str','tau1_tvalue_str','tau1_p_str')
         #adding alpha3  to the dataframe
         res$alpha3_str=NaN
         res$alpha3_se_str=NaN
                        }
); # end switch

#res<<-res
#print(res)

##########################################################################
#Calculate fitted values using optimized parameters

if (tau_model=='tau1_subset'){
  fitted_values<-tryCatch(fit1<-{subset_str(x=x_subset, fo=res$fo_str, fm=res$fm_str, tau1=res$tau1_str)},error=function(Fit.Fail){return("Fit Fail")})
}else{
  fitted_values<-tryCatch(fit1<-{kolber_str(x=x, fo=res$fo_str, fm=res$fm_str, tau1=res$tau1_str, alpha1=res$alpha1_str, tau2=res$tau2_str, alpha2=res$alpha2_str, tau3=res$tau3_str, tau_model=tau_model)},error=function(Fit.Fail){return("Fit Fail")})
}

#########################################################################


##########################################################################
#Optional plot
#print(fitted_values)

    if (plots==TRUE){
    graphics::plot(x=x,y=y, col="blue", xlab="Time (log(s))", ylab="Fluorescence yield", bty="l",log="x",las=1);

  if (tau_model=='tau1_subset'){
    graphics::points(x_subset,fitted_values,col=2,type='l')
  }else{
    graphics::points(x,fitted_values,col=2,type='l')
  }



  # Light: Print optimized parameter values on the plot ####
  leg.L.model<-paste("Model: ",tau_model,sep="");
  leg.L.main<-"Parameter \u00B1 SE" # \u00B1 is the unicode character for the plus/minus symbol
  leg.L.Fo<-tryCatch({paste("Fo =", round(res$fo_str, digits=3), "\u00B1", round(res$fo_se_str, digits=3), sep=" ")},error=function(Fit.Fail){return("Fit Fail")})
  leg.L.Fm<-tryCatch({paste("Fm =", round(res$fm_str, digits=3), "\u00B1", round(res$fm_se_str, digits=3), sep=" ")},error=function(Fit.Fail){return("Fit Fail")})
  leg.L.tau1<-tryCatch({paste("tau1 =", round(res$tau1_str, digits=5), "\u00B1", round(res$tau1_se_str, digits=5), sep=" ")},error=function(Fit.Fail){return("Fit Fail")})
  leg.L.alpha1<-tryCatch({paste("alpha1 =", round(res$alpha1_str, digits=2), "\u00B1", round(res$alpha1_se_str, digits=2), sep=" ")},error=function(Fit.Fail){return("Fit Fail")})
  leg.L.tau2<-tryCatch({paste("tau2 =", round(res$tau2_str, digits=5), "\u00B1", round(res$tau2_se_str, digits=5), sep=" ")},error=function(Fit.Fail){return("Fit Fail")})
  leg.L.alpha2<-tryCatch({paste("alpha2 =", round(res$alpha2_str, digits=2), "\u00B1", round(res$alpha2_se_str, digits=2), sep=" ")},error=function(Fit.Fail){return("Fit Fail")})
  leg.L.tau3<-tryCatch({paste("tau3 =", round(res$tau3_str, digits=5), "\u00B1", round(res$tau3_se_str, digits=5), sep=" ")},error=function(Fit.Fail){return("Fit Fail")})
  leg.L.alpha3<-tryCatch({paste("alpha3 =", round(res$alpha3_str, digits=2), "\u00B1", round(res$alpha3_se_str, digits=2), sep=" ")},error=function(Fit.Fail){return("Fit Fail")})


  if(is.na(res$alpha1_str) && is.na(res$alpha2_str)){
    graphics::legend("topright",legend=c(leg.L.model,leg.L.main,leg.L.Fo,leg.L.Fm,leg.L.tau1),bty="n",text.font=c(2,2,1,1,1),adj=c(0, 0.5))
  }
  if(!is.na(res$alpha1_str)&is.na(res$tau3_str)){

    graphics::legend("topright",legend=c(leg.L.model,leg.L.main,leg.L.Fo,leg.L.Fm,leg.L.tau1,leg.L.alpha1,leg.L.tau2,leg.L.alpha2),bty="n",text.font=c(2,2,1,1,1,1,1,1),adj=c(0, 0.5))
  }
  if(!is.na(res$alpha3_str)){

    graphics::legend("topright",legend=c(leg.L.model,leg.L.main,leg.L.Fo,leg.L.Fm,leg.L.tau1,leg.L.alpha1,leg.L.tau2,leg.L.alpha2,leg.L.tau3,leg.L.alpha3),bty="n",text.font=c(2,2,1,1,1,1,1,1,1,1),adj=c(0, 0.5))
    }


}
##########################################################################


##########################################################################
#Output list
#res<-assign("res",res, envir = e2)
output<-list(res,fitted_values,x)
names(output)<-c('fitted_parameters','fitted_values','time')
return(output)


}
