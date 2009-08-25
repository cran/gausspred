\name{train_pred_gau}

\alias{train_pred_gau}

\alias{training_gau}

\alias{predict_gau}

\alias{crossvalid_gau}

\title{Training with Markov chain sampling, predicting for test cases, and evaluating performance with cross-validation}

\description{
\code{training_gau} trains the Gaussian classification models with Markov chain Monte Carlo.

\code{predict_gau} uses the posterior samples returned by \code{training_gau} to predict the response values of test cases. 

\code{crossvalid_gau} uses cross-validation to evaluate the prediction performance, which can be used to evaluate the goodness of the prediction methods, also used to tune parameters used in predictions, such as of prior distributions and Markov chain sampling.
}

\usage{
training_gau (
    ## arguments specifying data sets
    G,features,response,
    ## arguments specifying priors
    prior_y=rep(1,G),
    p_tau_nu =  c(alpha=2,w=0.5),
    p_tau_mu =  c(alpha=2,w=0.5),
    p_tau_x  =  c(alpha=2,w=0.5),
    ## arguments specifying Gibbs sampling
    nos_super = 100, nos_trans = 5, ini_taus=rep(1,3),
    ## arguments specifying correcting for bias 
    cor=0, p=ncol(features),cutoff=1, min_qf=exp(-10),
    nos_lambda=1000, stepsize_log_tau=-2, no_steps=10 
    )
predict_gau (features, out_tr, pred_range, thin)
crossvalid_gau ( 
    ## arguments specifying data sets and crossvalidation
    no_fold, k, G, features,response, lmd_trf,
    ## arguments specifying priors
    prior_y=rep(1,G),
    p_tau_nu =  c(alpha=2,w=0.5),
    p_tau_mu =  c(alpha=2,w=0.5),
    p_tau_x  =  c(alpha=2,w=0.5),
    ## arguments specifying Gibbs sampling
    nos_super = 100, nos_trans = 5, ini_taus=rep(1,3),
    ## arguments specifying correcting for bias 
    cor=0, min_qf=exp(-10),nos_lambda=1000, 
    stepsize_log_tau=-2, no_steps=10 ,
    ## arguments specifying prediction
    pred_range, thin =1
    )
}

\arguments{

\item{}{\bold{Arguments of \code{training_gau} and \code{crossvalid_gau}:}}


\item{G}{the number of groups, ie, number of possibilities of response}

\item{features}{the features, with the rows for the cases.}

\item{response}{the response values.}

\item{prior_y}{a vector of length 'G', specifying the Dirichlet prior distribution for probabilities of the response.}

\item{p_tau_nu, p_tau_mu, p_tau_x}{vectors of 2 numbers, specifying the Gamma distribution as prior for the inverse of the variance of the distribution of \eqn{\nu}, \eqn{\mu}, and features respectively; the first number is shape, the second is rate.}

\item{nos_super, nos_trans}{\code{nos_super} of super Markov chain transitions are run, with \code{nos_trans} Markov chain iterations for each. Only the last state of each super transition is saved. This is used to avoid saving Markov chain state for each iteration. }

\item{ini_taus}{a vector of length 3, specifying initial values for \eqn{\tau^\nu}, \eqn{\tau^\mu}, \eqn{\tau^x}.}

\item{cor}{taking value 0 or 1, indicating whether bias-correction is to be applied.}

\item{p}{the number of total features before selection. This number needs to be supplied by users other than inferred from other arguments.}

\item{cutoff}{the cutoff of F-statistic used to select features. This number needs to be supplied by users other than inferred from other arguments.}

\item{min_qf}{the minimum value of "f" used to cut the infinite summation in calculating correction factor. Details see the paper.}

\item{nos_lambda}{the number of random numbers for \eqn{\Lambda} in approximating the correction factor. Details see the paper.}

\item{stepsize_log_tau}{the stepsize of Gaussian proposal used in sampling log \code{tau_mu} when bias-correction is applied.}

\item{no_steps}{iterations of Metropolis sampling for log \code{tau_mu}.}

\item{}{\bold{Arguments only of \code{predict_gau}:}}

\item{out_tr}{output of Markov chain sampling returned by \code{training_gau}}

\item{pred_range}{the range of super Markov chain transitions used to predict the response values of test cases.}

\item{thin}{only 1 sample for every \code{thin} samples are used in predicting, chosen evenly.}

\item{}{\bold{Arguments only of \code{crossvalid_gau}:}}

\item{no_fold}{the number of subsets of the data in making cross-validation assessment.}

\item{k}{the number of features selected.}

\item{lmd_trf}{the value of lambda used to estimate covariance matrix, which is used to transform data with Choleski decomposition. The larger this number is the estimated covariance matrix is closer to diagonal.}

}

\value{

The function \code{training_gau} returns the following values:

\item{mu}{an array of three dimensions, storing Markov chain samples of \eqn{\mu}, with the first dimension for different features, the 2nd dimension for different groups, the third dimension for different Markov chain super transitions.}

\item{nu}{a matrix, storing Markov chain samples of \eqn{\nu}, with rows for features, the columns for Markov chain iterations.}

\item{tau_x}{a vector, storing Markov chain samples of \eqn{\tau^x}.}

\item{tau_mu}{a vector, storing Markov chain samples of \eqn{\tau^\mu}.}

\item{tau_nu}{a vector, storing Markov chain samples of \eqn{\tau^\nu}.}

\item{freq_y}{the posterior mean of the probabilities for response.}

Both \code{predict_gau} and \code{crossvalid_gau} return a matrix of the predictive probabilities, with rows for cases, columns for different groups (different values of response).

}

\examples{

##### this is a full demonstration of using this package ######
###############################################################

## parameter setting
n <- 200+400
p <- 400
G <- 6
p_tau_x  <- c(4,1)
p_tau_mu <- c(1.5,0.01)
p_tau_nu <- c(1.5,0.01)
tau_mu <- 100 

## generate a data set
data <- gen_bayesgau (n,p,G,tau_nu=100,tau_mu,p_tau_x )

## specifying cases as training set
ix_tr <- 1:200 

## ordering features by F-statistic
i_sel <- order_features(data$X[ix_tr,],data$y[ix_tr])
vars <- i_sel$vars[1:10]
cutoff <- i_sel$fstat[10]

## training model with bias-corrected method
out_tr_cor <- training_gau(
    G = G, data$X[ix_tr,vars,drop=FALSE], data$y[ix_tr], 
    prior_y = rep(1,G),
    p_tau_nu, p_tau_mu, p_tau_x ,
    nos_super = 400, nos_trans = 1, ini_taus=rep(1,3),
    ## information on correcting for bias 
    cor=1, p=p,cutoff=cutoff, min_qf=exp(-10),nos_lambda=100, 
    stepsize_log_tau=0.5, no_steps=5 
    )

## make prediction
out_pred_cor <- predict_gau( 
    data$X[-(ix_tr),vars,drop=FALSE], out_tr_cor, 
    pred_range=c(50,400), thin = 1) 

## define 0-1 loss function
Mlosser <- matrix (1,G,G)
diag(Mlosser) <- 0

## randomly generate a loss function
Mloss <- matrix(1,G,G)
Mloss <- matrix(exp(rnorm(G^2,0,2)),G,G)
diag(Mloss) <- 0

## evaluate prediction with test cases 

## calculating average minus log probabilities
amlp_cor <- comp_amlp(out_pred_cor,data$y[-ix_tr])

## calculating error rate
er_cor <- comp_loss(out_pred_cor,data$y[-ix_tr],Mlosser)

## calculating average loss from the randomly generated loss function
l_cor <- comp_loss(out_pred_cor,data$y[-ix_tr],Mloss)

}

\keyword{classif}
