\name{MS_Var_Likelihood}
\alias{MS_Var_Likelihood}
\title{ Likelihood Function for Markov Switching VAR }
\description{ 
  \code{MS_Var} This is the likelihood function for the Markov Switching 
VAR which is minimized by nlminb.
}

\usage{
MS_Var_Likelihood(params)
}

\arguments{
  \item{ params}{ A Vector of Initial Parameter estimates }
}

\value{
  a scalar: 
  \item{LL}{sum of log likelihood for the model}
}

\author{ James Eustace \email{james.eustace@omam.co.uk}}}
