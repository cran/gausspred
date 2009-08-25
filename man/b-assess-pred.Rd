\name{assess_prediction}

\alias{assess_prediction}

\alias{comp_amlp}

\alias{comp_er}

\alias{comp_loss}

\alias{cal_tab}

\title{Functions for evaluating predictive probabilities given true values of response of test cases}

\description{These functions evaluate predictive probabilities with average minus log probabilities, error rate, and average loss for a defined loss function, or calculate calibration table. 
}

\usage{
comp_amlp (probs_pred, responses)
comp_er (probs_pred, responses)
comp_loss (probs_pred, y_true, Mloss)
cal_tab (probs_pred, true_y, ix_y, no_cat=10)
}

\arguments{

\item{probs_pred}{a matrix of the predictive probabilities, with rows for cases, columns for groups (different values of response).}

\item{Mloss}{a matrix defining a loss function, with rows for true values, and columns for predicted values.}

\item{responses, y_true, true_y}{a vector of true values of response in test cases.}

\item{ix_y}{the index of column used to produce calibration table.}

\item{no_cat}{number of categories in producing calibration table.}



}

\value{
\code{comp_amlp} returns average minus log probabilities, \code{comp_er} returns error rate, \code{comp_loss} returns average loss, and expected loss, \code{cal_tab} returns a calibration data frame.
}

\examples{
## See train_pred_gau
}

\seealso{

\code{\link{train_pred_gau}}

}

\keyword{classif}
