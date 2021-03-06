% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/models.R
\name{kolber_str}
\alias{kolber_str}
\title{Model (Kolber et al. 1998) - relaxation}
\usage{
kolber_str(x, fo, fm, tau1, alpha1, tau2, alpha2, tau3, tau_model)
}
\arguments{
\item{fo}{Fo value}

\item{fm}{Fm value}

\item{tau1}{parameter that will store the value of tau1}

\item{alpha1}{no idea what is it}

\item{tau2}{parameter that will store the value of tau2}

\item{alpha2}{no idea what is it}

\item{tau3}{parameter that will store the value of tau3}

\item{tau_model}{the model type defined by the user (i.e. tau1, tau2 or tau3)}
}
\value{
a value that will be used in the fitting function
}
\description{
definition of the model that will be fitted to the relaxation step
}
\keyword{internal}
