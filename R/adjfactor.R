gen_QF_lmd <- function(cutoff, nos_g, p_tau_x=c(1.5,0.01),
                       min_qf=exp(-10),nos_lambda=1000)
{   
    G <- length(nos_g)
    n <- sum(nos_g)

    comp_lambda <- function(z)
    { z.bar <- sum(z * nos_g) / n  
      sum ( (z - z.bar)^2 * nos_g )
    }

    QF <- c()
    k = 1
    while(TRUE){
      QF[k] <- 
          pf(cutoff*(G-1)/(G-1 + 2*(k-1)), G-1 + 2*(k-1), n-G)
      if( QF[k] < min_qf ) break
      k = k + 1
    }
    #print(exp(cent_F_quan))
    lambdas <- 
       0.5 * rgamma(nos_lambda,p_tau_x[1]/2, p_tau_x[1] * p_tau_x[2] /2 ) * 
       apply(matrix(rnorm(nos_lambda*G),G,nos_lambda),2,comp_lambda)
       
    list(lambdas = lambdas, QF=QF)
}

comp_adjfactor <- function(tau_mu, QF_lambdas)
{
   lambdas <- QF_lambdas$lambdas / tau_mu
   QF <- QF_lambdas$QF
   adjfactor <- 0
   .C("comp_adjfactor",length(QF),length(lambdas),QF,lambdas,
       adjfactor=adjfactor,PACKAGE = "gausspred")$adjfactor

}

