#this function plots all the datasets in the the PSI file and should be directed
#to a pdf file in the fitting function
#the input for the plotting function are the raw data at the fitting results
#each function bellow NEEDS documentation

#to fit the sti data
psi_plots<- function(y, fit1,fit1u,fit1l,titl,res){

fit1<-local(fit1, env = e2)
fit1u<-local(fit1u, env = e2)
fit1l<-local(fit1l, env = e2)

  
plot(y, col="blue", xlab="Flash number", ylab="Fluorescence yield", bty="l",main=titl);

#Draw black line for the fitted model and red for the CI

tryCatch({points(fit1, type="l", col="black")},error=function(Fit.Fail){return("Fit Fail")})
tryCatch({points(fit1u, type="l", lty=2, col="red")},error=function(Fit.Fail){return("Fit Fail")})
tryCatch({points(fit1l, type="l", lty=2, col="red")},error=function(Fit.Fail){return("Fit Fail")})


#this is for the legend at the bottom right
# Light: Print optimized parameter values on the plot ####
  leg.L.model<-paste("Model: ",res$model_sti,sep="");
  leg.L.main<-"Parameter \u00B1 SE" # \u00B1 is the unicode character for the plus/minus symbol
  leg.L.Fo<-tryCatch({paste("Fo =", round(res$fo_sti, digits=3), "\u00B1", round(res$fo_se_sti, digits=3), sep=" ")},error=function(Fit.Fail){return("Fit Fail")})
  leg.L.Fm<-tryCatch({paste("Fm =", round(res$fm_sti, digits=3), "\u00B1", round(res$fm_se_sti, digits=3), sep=" ")},error=function(Fit.Fail){return("Fit Fail")})
  leg.L.sigma<-tryCatch({paste("sigma =", round(res$sigma_sti, digits=0), "\u00B1", round(res$sigma_se_sti, digits=0), sep=" ")},error=function(Fit.Fail){return("Fit Fail")})
  leg.L.rho<-tryCatch({paste("rho =", round(res$rho_sti, digits=3), "\u00B1", round(res$rho_se_sti, digits=3), sep=" ")},error=function(Fit.Fail){return("Fit Fail")})
  legend("bottomright",legend=c(leg.L.model,leg.L.main,leg.L.Fo,leg.L.Fm,leg.L.sigma,leg.L.rho),bty="n",text.font=c(2,2,1,1,1,1),adj=c(0, 0.5))
  

} 



#to fit the str data
psi_plots_str<- function(x,y, fit1,fit1u,fit1l,titl,res){

  tryCatch({fit1l<-psiworx.model.str(x=x, fo=res$fo_str-res$fo_se_str, fm=res$fm_str-res$fm_se_str, tau1=res$tau1_str-res$tau1_se_str, alpha1=res$alpha1_str-res$alpha1_se_str, tau2=res$tau2_str-res$tau2_se_str, alpha2=res$alpha2_str-res$alpha2_se_str, tau3=res$tau3_str-res$tau3_se_str, wm=res$model_str)},error=function(Fit.Fail){return("Fit Fail")})
  # Plot fit of relaxation period
  plot(x=x,y=y, col="blue", xlab="Time (s)", ylab="Fluorescence yield", bty="l",main=titl,log="x");
  
  # Light: Calculate fitted values using optimized parameters
  tryCatch({fit1<-psiworx.model.str(x=x, fo=res$fo_str, fm=res$fm_str, tau1=res$tau1_str, alpha1=res$alpha1_str, tau2=res$tau2_str, alpha2=res$alpha2_str, tau3=res$tau3_str, wm=res$model_str)},error=function(Fit.Fail){return("Fit Fail")})
  tryCatch({fit1u<-psiworx.model.str(x=x, fo=res$fo_str+res$fo_se_str, fm=res$fm_str+res$fm_se_str, tau1=res$tau1_str+res$tau1_se_str, alpha1=res$alpha1_str+res$alpha1_se_str, tau2=res$tau2_str+res$tau2_se_str, alpha2=res$alpha2_str+res$alpha2_se_str, tau3=res$tau3_str+res$tau3_se_str, wm=res$model_str)},error=function(Fit.Fail){return("Fit Fail")})
  
  #Draw red line for the fitted model
  clr<-"red";
  tryCatch({points(x=x,y=fit1, type="l", col="black")},error=function(Fit.Fail){return("Fit Fail")})
  tryCatch({points(x=x,y=fit1u, type="l", lty=2, col="red")},error=function(Fit.Fail){return("Fit Fail")})
  tryCatch({points(x=x,y=fit1l, type="l", lty=2, col="red")},error=function(Fit.Fail){return("Fit Fail")})
  
  # Light: Print optimized parameter values on the plot ####
  leg.L.model<-paste("Model: ",res$model_str,sep="");
  leg.L.main<-"Parameter \u00B1 SE" # \u00B1 is the unicode character for the plus/minus symbol
  leg.L.Fo<-tryCatch({paste("Fo =", round(res$fo_str, digits=3), "\u00B1", round(res$fo_se_str, digits=3), sep=" ")},error=function(Fit.Fail){return("Fit Fail")})
  leg.L.Fm<-tryCatch({paste("Fm =", round(res$fm_str, digits=3), "\u00B1", round(res$fm_se_str, digits=3), sep=" ")},error=function(Fit.Fail){return("Fit Fail")})
  leg.L.tau1<-tryCatch({paste("tau1 =", round(res$tau1_str, digits=5), "\u00B1", round(res$tau1_se_str, digits=5), sep=" ")},error=function(Fit.Fail){return("Fit Fail")})
  if(is.na(res$alpha1_str) && is.na(res$alpha2_str)){
    legend("topright",legend=c(leg.L.model,leg.L.main,leg.L.Fo,leg.L.Fm,leg.L.tau1),bty="n",text.font=c(2,2,1,1,1),adj=c(0, 0.5))
  }
  if(!is.na(res$alpha1_str)){
    leg.L.alpha1<-tryCatch({paste("alpha1 =", round(res$alpha1_str, digits=2), "\u00B1", round(res$alpha1_se_str, digits=2), sep=" ")},error=function(Fit.Fail){return("Fit Fail")})
    leg.L.tau2<-tryCatch({paste("tau2 =", round(res$tau2_str, digits=5), "\u00B1", round(res$tau2_se_str, digits=5), sep=" ")},error=function(Fit.Fail){return("Fit Fail")})
    legend("topright",legend=c(leg.L.model,leg.L.main,leg.L.Fo,leg.L.Fm,leg.L.tau1,leg.L.alpha1,leg.L.tau2),bty="n",text.font=c(2,2,1,1,1,1,1),adj=c(0, 0.5))
  }
  if(!is.na(res$alpha2_str)){
    leg.L.alpha2<-tryCatch({paste("alpha2 =", round(res$alpha2_str, digits=2), "\u00B1", round(res$alpha2_se_str, digits=2), sep=" ")},error=function(Fit.Fail){return("Fit Fail")})
    leg.L.tau3<-tryCatch({paste("tau3 =", round(res$tau3_str, digits=5), "\u00B1", round(res$tau3_se_str, digits=5), sep=" ")},error=function(Fit.Fail){return("Fit Fail")})
    legend("topright",legend=c(leg.L.model,leg.L.main,leg.L.Fo,leg.L.Fm,leg.L.tau1,leg.L.alpha1,leg.L.tau2,leg.L.alpha2,leg.L.tau3),bty="n",text.font=c(2,2,1,1,1,1,1,1,1),adj=c(0, 0.5))
  }
  
  
}



###############################################################################
##Functions for the plots on the screen, the above were for the pdfs
###############################################################################

#to fit the STI data
psi_plots_screen<- function(y, fit1,fit1u,fit1l,titl,res){
  

  fit1<-local(fit1, env = e2)
  fit1u<-local(fit1u, env = e2)
  fit1l<-local(fit1l, env = e2)
  
  
  plot(y, col="blue", xlab="", ylab="", bty="l",xaxt='n',yaxt='n',cex=0.6,pch=21,bg='blue')
  legend("topleft",legend = titl,cex=0.6,bty='n')
  #Draw black line for the fitted model and red for the CI
  
  tryCatch({points(fit1, type="l", col="black")},error=function(Fit.Fail){return("Fit Fail")})
  #tryCatch({points(fit1u, type="l", lty=2, col="red")},error=function(Fit.Fail){return("Fit Fail")})
  #tryCatch({points(fit1l, type="l", lty=2, col="red")},error=function(Fit.Fail){return("Fit Fail")})
  

  #this is for the legend at the bottom right
  leg.L.model<-paste("Model: ",res$model_sti,sep="");
  leg.L.main<-"Parameter \u00B1 SE" # \u00B1 is the unicode character for the plus/minus symbol
  leg.L.Fo<-tryCatch({paste("Fo =", round(res$fo_sti, digits=3), "\u00B1", round(res$fo_se_sti, digits=3), sep=" ")},error=function(Fit.Fail){return("Fit Fail")})
  leg.L.Fm<-tryCatch({paste("Fm =", round(res$fm_sti, digits=3), "\u00B1", round(res$fm_se_sti, digits=3), sep=" ")},error=function(Fit.Fail){return("Fit Fail")})
  leg.L.sigma<-tryCatch({paste("sigma =", round(res$sigma_sti, digits=0), "\u00B1", round(res$sigma_se_sti, digits=0), sep=" ")},error=function(Fit.Fail){return("Fit Fail")})
  leg.L.rho<-tryCatch({paste("rho =", round(res$rho_sti, digits=3), "\u00B1", round(res$rho_se_sti, digits=3), sep=" ")},error=function(Fit.Fail){return("Fit Fail")})
  legend("bottomright",legend=c(leg.L.model,leg.L.main,leg.L.Fo,leg.L.Fm,leg.L.sigma,leg.L.rho),bty="n",text.font=c(2,2,1,1,1,1),adj=c(0, 0.5),cex=0.6)
  
  
} 

###############################################################################
###plots STR on the screen
###############################################################################


psi_plots_str_screen<- function(x,y, fit1,fit1u,fit1l,titl,res){
  
  tryCatch({fit1l<-psiworx.model.str(x=x, fo=res$fo_str-res$fo_se_str, fm=res$fm_str-res$fm_se_str, tau1=res$tau1_str-res$tau1_se_str, alpha1=res$alpha1_str-res$alpha1_se_str, tau2=res$tau2_str-res$tau2_se_str, alpha2=res$alpha2_str-res$alpha2_se_str, tau3=res$tau3_str-res$tau3_se_str, wm=res$model_str)},error=function(Fit.Fail){return("Fit Fail")})
  # Plot fit of relaxation period
  plot(x=x,y=y, col="blue", xlab="", ylab="", bty="l",log="x",xaxt='n',yaxt='n',cex=0.6,pch=21,bg='blue')
legend("bottomleft",legend = titl,cex=0.6,bty='n')
  
  # Light: Calculate fitted values using optimized parameters
  tryCatch({fit1<-psiworx.model.str(x=x, fo=res$fo_str, fm=res$fm_str, tau1=res$tau1_str, alpha1=res$alpha1_str, tau2=res$tau2_str, alpha2=res$alpha2_str, tau3=res$tau3_str, wm=res$model_str)},error=function(Fit.Fail){return("Fit Fail")})
  tryCatch({fit1u<-psiworx.model.str(x=x, fo=res$fo_str+res$fo_se_str, fm=res$fm_str+res$fm_se_str, tau1=res$tau1_str+res$tau1_se_str, alpha1=res$alpha1_str+res$alpha1_se_str, tau2=res$tau2_str+res$tau2_se_str, alpha2=res$alpha2_str+res$alpha2_se_str, tau3=res$tau3_str+res$tau3_se_str, wm=res$model_str)},error=function(Fit.Fail){return("Fit Fail")})
  
  #Draw red line for the fitted model
  clr<-"red";
  tryCatch({points(x=x,y=fit1, type="l", col="black")},error=function(Fit.Fail){return("Fit Fail")})
  #tryCatch({points(x=x,y=fit1u, type="l", lty=2, col="red")},error=function(Fit.Fail){return("Fit Fail")})
  #tryCatch({points(x=x,y=fit1l, type="l", lty=2, col="red")},error=function(Fit.Fail){return("Fit Fail")})
  
  
  #BJ: had to remove the SE from the plots due to lack of space
  leg.L.model<-paste("Model: ",res$model_str,sep="");
  #leg.L.main<-"Parameter" # \u00B1 is the unicode character for the plus/minus symbol
  leg.L.Fo<-tryCatch({paste("Fo =", round(res$fo_str, digits=3), sep=" ")},error=function(Fit.Fail){return("Fit Fail")})
  leg.L.Fm<-tryCatch({paste("Fm =", round(res$fm_str, digits=3), sep=" ")},error=function(Fit.Fail){return("Fit Fail")})
  leg.L.tau1<-tryCatch({paste("tau1 =", round(res$tau1_str, digits=5), sep=" ")},error=function(Fit.Fail){return("Fit Fail")})
  if(is.na(res$alpha1_str) && is.na(res$alpha2_str)){
    legend("topright",legend=c(leg.L.model,leg.L.Fo,leg.L.Fm,leg.L.tau1),bty="n",text.font=c(2,2,1,1,1),adj=c(0, 0.5),cex=0.6)
  }
  if(!is.na(res$alpha1_str)){
    leg.L.alpha1<-tryCatch({paste("alpha1 =", round(res$alpha1_str, digits=2), sep=" ")},error=function(Fit.Fail){return("Fit Fail")})
    leg.L.tau2<-tryCatch({paste("tau2 =", round(res$tau2_str, digits=5), sep=" ")},error=function(Fit.Fail){return("Fit Fail")})
    legend("topright",legend=c(leg.L.model,leg.L.Fo,leg.L.Fm,leg.L.tau1,leg.L.alpha1,leg.L.tau2),bty="n",text.font=c(2,2,1,1,1,1,1),adj=c(0, 0.5),cex=0.6)
  }
  if(!is.na(res$alpha2_str)){
    leg.L.alpha2<-tryCatch({paste("alpha2 =", round(res$alpha2_str, digits=2), sep=" ")},error=function(Fit.Fail){return("Fit Fail")})
    leg.L.tau3<-tryCatch({paste("tau3 =", round(res$tau3_str, digits=5), sep=" ")},error=function(Fit.Fail){return("Fit Fail")})
    legend("topright",legend=c(leg.L.model,leg.L.Fo,leg.L.Fm,leg.L.tau1,leg.L.alpha1,leg.L.tau2,leg.L.alpha2,leg.L.tau3),bty="n",text.font=c(2,2,1,1,1,1,1,1,1),adj=c(0, 0.5),cex=0.6)
  }
  
  
}




