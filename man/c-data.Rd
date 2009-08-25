\name{data_gau}

\alias{gen_bayesgau}

\alias{order_features}

\title{Data generation and processing}

\description{\code{gendata_bayesgau} generates data from a Bayesian Gaussian classification model. Model specification please refer to the paper.

\code{order_features} orders features by the F-statistic.
}

\usage{
gen_bayesgau (n,p,G,tau_nu=1,tau_mu=1,p_tau_x=c(2,1) )
order_features (features, response)
}

\arguments{

\item{n}{the number of cases.}

\item{p}{the number of features.}

\item{G}{the number of groups.}

\item{tau_nu}{the inverse of variance of Gaussian distribution for generating \eqn{\nu}.}

\item{tau_mu}{the inverse of variance of Gaussian distribution for generating \eqn{\mu}.}

\item{p_tau_x}{a vector of 2 numbers, specifying the Gamma distribution as prior for the inverse of the variance for the distribution of features, the first number is shape, the second is rate.}

\item{features}{the features, with the rows for the cases.}

\item{response}{the response values.}

}

\value{

The function \code{gendata_bayesgau} returns the following items:
\item{X}{the features, with the row standing for the cases.}

\item{y}{the response values.}

\item{mu}{the values of mu, with columns for groups.}

\item{sigmas}{the variance generating features.}

\item{nu}{the value of nu for features.}

The function \code{order_features} returns the following iterms:

\item{vars}{the indice of features, in decreasing order of F-statistic.}

\item{fstats}{the values of F-statistic for features indexed by \code{vars}.}

}

\examples{

# generate data for Gaussian classification models

data <- gen_bayesgau (n = 100,p = 100,G = 2,
        tau_nu=100,tau_mu = 100, p_tau_x = c(4,1))
i_sel <- order_features (data$X, data$y)
}


\keyword{datagen}
