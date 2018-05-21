#' @title Fitting function for the EP model
#' @description  Fitting function for the EP model
#' @param light light values
#' @param etr etr values
#' @param a1
#' @param b1 
#' @param c1 
#' @return the fitted values for the EP model
#' @export 
#' @keywords external
ep_model<-function(light,etr,a1=0.000001,b1=0.001,c1=3)
{
  
  #set up the dataframe  to be processed
  pe3<<-as.data.frame(cbind(light,etr))
  
 # names(pe3)<-c("light","etr")
  
  min.func<-function(params,data=pe3)
  {  a<-params[1]
  b<-params[2]
  c<-params[3]
  return( sum((data$etr-data$light/(a*data$light^2+b*data$light+c))^2))
  }
  
  etr_sim<-optim(par=c(a1, b1 , c1),fn=min.func,control=list("maxit"=1000000,pgtol=0.01) )
  
  print(etr_sim)
  
  a2<-etr_sim$par[1]
  b2<-etr_sim$par[2]
  c2<-etr_sim$par[3]
  
  alpha2<-(1/c2)
  etrmax2<-1/(b2+2*(a2*c2)^0.5)
  Eopt2<-(c2/a2)^0.5

  #if the model does not fit the ETRmax then an alternative is proposed 
  #in the form of the maximum observed value, this will allow to estimate Ek
  if (etrmax2=="NaN")
  {
    etrmax2<-max(data$etr)
    Eopt2<-"NA"
  }
  Ek2<-etrmax2/alpha2
  
  x2<-seq(0,max(pe3$light),by=5)
  
  predicted<-x2/(a2*x2^2+b2*x2+c2)

  print(paste("Alpha=",round(alpha2,2)))
  print(paste("ETR=",round(etrmax2,2)))
  print(paste("Ek=",round(Ek2,2)))
  print(paste("Eopt=",round(Eopt2,2)))

  
  plot(pe3$light,pe3$etr,pch=21,xlab="Light",ylab="ETR")
  lines(x2,predicted,col="blue")
  
}

#######Simplified version without the parameters a, b and c

#' @title Fitting function for the EP model (simplified)
#' @description  Fitting function for the EP model (simplified)
#' @param light light values
#' @param etr etr values
#' @param Eopt
#' @param alpha 
#' @param Ps
#' @return the fitted values for the EP model
#' @export 
#' @keywords external
ep_model_simple<-function(light,etr,alpha=0.1,Eopt=300,Ps=200)
{
  
  #set up the dataframe  to be processed
  pe3<<-as.data.frame(cbind(light,etr))
  
  # names(pe3)<-c("light","etr")
  
  min.func<-function(params,data=pe3)
  {  alpha<-params[1]
  Eopt<-params[2]
  Ps<-params[3]
  return( sum((  data$etr-data$light/(data$light^2*(1/(alpha*Eopt^2))+(data$light/Ps)-((2*data$light)/(alpha*Eopt))+(1/alpha)))^2))
  }
  

  etr_sim<-optim(par=c(alpha, Eopt , Ps),fn=min.func,control=list("maxit"=1000000,pgtol=0.01) )
  
  print(etr_sim)
  
  alpha2<-etr_sim$par[1]
  Eopt2<-etr_sim$par[2]
  Ps2<-etr_sim$par[3]
  
  #alpha2<-(1/c2)
  #etrmax2<-1/(b2+2*(a2*c2)^0.5)
  #Eopt2<-(c2/a2)^0.5
  
  #if the model does not fit the ETRmax then an alternative is proposed 
  #in the form of the maximum observed value, this will allow to estimate Ek
  if (Ps2=="NaN")
  {
    Ps2<-max(data$etr)
    Eopt2<-"NA"
  }
  Ek2<-Ps2/alpha2
  
  x2<-seq(0,max(pe3$light),by=5)
  
  predicted<-x2/(x2^2*(1/(alpha2*Eopt2^2))+(x2/Ps2)-((2*x2)/(alpha2*Eopt2))+(1/alpha2))
  
  
  print(paste("Alpha=",round(alpha2,2)))
  print(paste("ETR=",round(Ps2,2)))
  print(paste("Ek=",round(Ek2,2)))
  print(paste("Eopt=",round(Eopt2,2)))
  
  
  plot(pe3$light,pe3$etr,pch=21,xlab="Light",ylab="ETR")
  lines(x2,predicted,col="blue")
  
}


