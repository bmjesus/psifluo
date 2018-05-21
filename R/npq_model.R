#' @title Fitting function for the EP model
#' @description  Fitting function for the EP model
#' @param light light values
#' @param npq npq values
#' @return the parameters for the fitted model
#' @export 
#' @keywords external
npq_model<-function(light,npq)

{


#this constructs the function that will be optmized in the next step. the criteria is to minimize the sum of the square of the deviations
#the fitting follows the npq model


pe3<-as.data.frame(cbind(light,npq))
#pe3<<-teste
names(pe3)<-c("light","npq")
min.func<-function(params,data=pe3)
#Secao para estimar os parametros da curva luminosa.
pe3<<-teste
min.func<-function(params,data=pe3)
{ NPQmax<-params[1]
  E50<-params[2]
  hill<-params[3]
  return(sum( (data$npq-NPQmax*(data$light^hill/((E50^hill + data$light^hill))))^2))
}

npq_sim<<-optim(par=c(2,20,1),fn=min.func,hessian=TRUE)

#print(boot(teste,min.func$hill,900))
NPQmax2<<-npq_sim$par[1]
E502<<-npq_sim$par[2]
hill2<<-npq_sim$par[3]
 
x2<-seq(0,max(light),by=5)

plot(pe3$light,pe3$npq,xlab="Light",ylab="NPQ",las=1,pch=21,bg=1,cex=1.5,cex.lab=1.5,cex.axis=1.5)
lines(x2,NPQmax2*(x2^hill2)/(E502^hill2 + x2^hill2),col="blue",lwd=2)

predicted<-NPQmax2*(pe3$light^hill2)/(E502^hill2 + pe3$light^hill2)
WSSR<<-sum((pe3$npq-predicted)^2)
S2<<-WSSR/(length(pe3$npq)-3)
hessian<-npq_sim$hessian
var_parameters<<-sqrt(abs(diag(solve(-1*npq_sim$hessian)))*S2)
#summary(npq_sim)
print(NPQmax2)
print(E502)
print(hill2)
#print(var_parameters)
#print(WSSR)
#print(S2)
}

