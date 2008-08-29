\name{MS_Var}
\alias{MS_Var}
\title{ Estimate Markov Switching VAR }
\description{ 
  \code{MS_Var} Estimates the Markov Switching VAR. Set-up optimization 
parameters and reports output.
}

\usage{
MS_Var(dep)
}

\arguments{
  \item{ dep}{ Matrix of endogenous variables }
}

\value{
  a list containing: 

  \item{transitionmatrix}{Markov transition probabilities for each state}
  \item{coefficients}{constant coefficients for each state}
  \item{Loglikelihood}{sum of log likelihood for the model}
  \item{densities}{vector state probability densities}
  \item{residuals}{vector of residuals for each state}
  \item{FilteredProbs}{vector of filtered probabilities for each state}
  \item{smoothed Probs }{vector of filtered probabilities for each state  }
}

\author{ James Eustace \email{james.eustace@omam.co.uk}}}

\examples{
#  Use the data-set 'Hamilton_Data.csv' to reproduce Hamilton's results
#  from "Long Swings in the Dollar: Are They in the Data and 
#  Do Markets Know It?", American Economic Review, Sept. 1990

data("EX_R")
data("INT_D")
dep<-cbind(EX_R,INT_D)
res<-MS_Var(dep)

}
