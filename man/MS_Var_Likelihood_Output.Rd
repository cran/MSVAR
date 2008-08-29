\name{MS_Var_Likelihood_Output}
\alias{MS_Var_Likelihood_Output}
\title{ Likelihood Function Output for Markov Switching VAR }
\description{ 
  \code{MS_Var} This function gathers the output from the likelihood function using 
the optimal parameters found by nlminb.
}

\usage{
MS_Var_Likelihood_Output(params)
}

\arguments{
  \item{ params}{ A Vector of Initial Parameter estimates }
}

\value{
  a list containing
  \item{transition matrix}{Markov transition probabilities for each state}
  \item{coefficients}{constant coefficients for each state}
  \item{Log-likelihood}{sum of log likelihood for the model}
  \item{densities}{vector state probability densities}
  \item{residuals}{vector of residuals for each state}
  \item{Filtered Probs}{vector of filtered probabilities for each state}
}

\author{ James Eustace \email{james.eustace@omam.co.uk}}}
