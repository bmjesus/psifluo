#' @title Model (Kolber et al. 1998) - induction
#' @description definition of the mode that will be fitted to the induction step
#' @param incident flashlets, numeric
#' @param sigma value of sigma PSII
#' @param fo Fo value
#' @param fm Fm value
#' @param rho connectivity value (rho)
#' @return a value that will be used in the fitting function
#' @keywords internal
#' @export
kolber_sti<-function(incident, sigma, fo, fm, rho){
  c<-rep(0, times=length(incident)) # initialize cumulative excitation
  c[2]<-incident[1]*sigma           # get the first excitation value, using c[2] to properly fit curve to "first point"
  # calculate fraction of closed PSII centres, with connectivity
  for (ii in 3:length(c)){
    c[ii]<-c[ii-1]+((incident[ii]*sigma*(1-c[ii-1]))/(1-(rho*c[ii-1])));
  }
  return(fo + (((fm - fo)*c*(1-rho))/(1-(c*rho))))
}


#' @title Model (Kolber et al. 1998) - relaxation
#' @description definition of the model that will be fitted to the relaxation step
#' @param tau_model the model type defined by the user (i.e. tau1, tau2 or tau3)
#' @param fo Fo value
#' @param fm Fm value
#' @param tau1 parameter that will store the value of tau1
#' @param tau2 parameter that will store the value of tau2
#' @param tau3 parameter that will store the value of tau3
#' @param alpha1 no idea what is it
#' @param alpha2 no idea what is it
#' @return a value that will be used in the fitting function
#' @keywords internal
#' @export
kolber_str<-function(x,fo,fm,tau1,alpha1,tau2,alpha2,tau3,tau_model){
switch(tau_model,
         tau1={
           y <- fo + (fm-fo) * exp(-x/tau1);
         },
         tau2={
           y <- fo + (fm-fo) * (alpha1 * exp(-x/tau1) + (1-alpha1)*exp(-x/tau2));
         },
         tau3={
           y <- fo + (fm-fo) * (alpha1 * exp(-x/tau1) + alpha2*exp(-x/tau2) + (1-alpha1-alpha2)*exp(-x/tau3));
         }
  )
  return(y);
}

#' @title Simple negative exponential to a subset of the data in the relaxation curve
#' @description Definition of the model that will be fitted to the relaxation step
#' @param x  Time
#' @param fo Fo value
#' @param fm Fm value
#' @param tau1 Parameter that will store the value of tau1
#' @return A numerical value that will be used in the fitting function
#' @keywords internal
#' @export
subset_str<-function(x,fo,fm,tau1){
           y <- fo + (fm-fo) * exp(-x/tau1);
           return(y);
}

