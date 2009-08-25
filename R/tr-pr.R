training_gau <- function
      (## arguments specifying data sets
       G,features,response,
       ## arguments specifying priors
       prior_y=rep(1,G),
       p_tau_nu =  c(alpha=2,w=0.5),
       p_tau_mu =  c(alpha=2,w=0.5),
       p_tau_x  =  c(alpha=2,w=0.5),
       ## arguments specifying Gibbs sampling
       nos_super = 100, nos_trans = 5, ini_taus=rep(1,3),
       ## arguments specifying correcting for bias 
       cor=0, p=ncol(features),cutoff=1, min_qf=exp(-10),nos_lambda=1000, 
       stepsize_log_tau=-2, no_steps=10 
      )
{   n <- nrow(features)
    k <- ncol(features)
    if( p < k) stop("In 'training_gau', argument 'p' is wrong")    
    
    ## Compute sufficient statistic
    nos_g <- tapply(rep(1,n),INDEX = response, sum)    
    if(any(nos_g < 2)) stop("Less than 2 cases in some group") 
    if(length(nos_g)!=G) stop("group number is wrong")
    NOS_G <- matrix( nos_g, k, G, byrow = TRUE )
    x_bar <- matrix(0, k, G)    
    for(g in 1:G){
       x_bar[,g] <- apply( features[response == g,],2, mean )
    }
    X_ssq <- apply(features^2,2,sum )
    
    ## preliminary comp for correcting bias
    if(cor == 1 & k < p){
        QF_lmd <- gen_QF_lmd(cutoff, nos_g, p_tau_x,min_qf,nos_lambda)
    }
    ## Markov chain storage
    mu <- array(0, dim = c(k,G,nos_super) )
    nu <- array(0, dim = c(k,nos_super) )
    tau_x <- array(0, dim = c(k,nos_super)) 
    tau_mu <- rep(0, nos_super)
    tau_nu <- rep(0, nos_super)
    
    ## Markov chain state
    pred_means_smu <- pred_taus_smu <- sq_smu <- smu <- matrix(0,k,G) 
    snu <- rep(0,k)
    stau_x <- rep(ini_taus[3],k)
    stau_mu <- ini_taus[2]
    log_tau_mu <- log(stau_mu)
    stau_nu <- ini_taus[1]
    
    ## Start Markov chain super-transition
    i_super <- 1
    while(i_super <= nos_super)
    {   
       ## start Gibbs sampling
       i_mc <- 1
       while(i_mc <= nos_trans)
       {   
           ## update mu 
           tau_ng <- stau_x %*% t(nos_g)
           pred_taus_mu <- stau_mu + tau_ng
           pred_means_mu <- (x_bar * tau_ng + stau_mu * snu) / pred_taus_mu
           smu <- matrix(rnorm(G*k),k,G) / sqrt(pred_taus_mu) + pred_means_mu
           sq_smu <- smu^2
           
           ## update tau_x
           rates_tau_x <- 0.5 * p_tau_x [1]*p_tau_x [2] + 0.5 * X_ssq + 
                          0.5 * apply(sq_smu * NOS_G,1,sum) - 
                          apply(smu * x_bar * NOS_G, 1, sum)
           stau_x <- rgamma(k, (p_tau_x[1] + n)/2) / rates_tau_x           
           
           ## update nu
           snu <- rnorm(k) / sqrt(stau_nu + G*stau_mu) + 
                  apply(smu,1,sum) / (stau_nu / stau_mu + G)
           ssq_snu <- sum(snu ^ 2)
           
           ## update tau_mu
           shape_tau_mu <- (p_tau_mu[1] +  G * k) / 2
           rate_tau_mu <- 0.5 * p_tau_mu[1] * p_tau_mu[2] + 0.5 * sum(sq_smu) + 
                          0.5 * G * ssq_snu - sum( apply(smu,1,sum)*snu)
           if(cor == 0 || p == k){ 
               ## do not correct for bias
               stau_mu <- rgamma(1, shape_tau_mu , rate_tau_mu )
           }
           else { 
               ## do correct for bias
               ## define log posterior density
               log_post_tau_mu <- function(log_tau_mu){
                 ## term from selected features
                 dgamma(exp(log_tau_mu),shape_tau_mu,rate_tau_mu,
                 log = TRUE ) + 
                 ## term from discarded features
                 (p-k) * log (comp_adjfactor(exp(log_tau_mu),QF_lmd)) + 
                 ## Jocobian from log transformation
                 log_tau_mu
               }
               
               ## start Metrhopolis sampling
               i_mh <- 1
               log_post <- log_post_tau_mu(log_tau_mu)
               while(i_mh <= no_steps){
                  new_log_tau_mu <- rnorm(1,log_tau_mu, stepsize_log_tau)
                  new_log_post <- log_post_tau_mu(new_log_tau_mu)
                  if(log(runif(1)) < new_log_post - log_post)
                  {   log_tau_mu <- new_log_tau_mu
                      log_post <- new_log_post
                  }
                  i_mh <- i_mh + 1
               }
               stau_mu <- exp(log_tau_mu)
           }
           
           ## update tau_nu
           stau_nu <- rgamma ( 1, (p_tau_nu[1] + k)/2, 
                               0.5 * p_tau_nu[1] * p_tau_nu[2] + 0.5 * ssq_snu ) 
           i_mc <- i_mc + 1
       }
       
       ## write states into Markov chain arrays
       mu[,,i_super] <- smu
       nu[,i_super] <- snu
       tau_x[,i_super] <- stau_x
       tau_mu[i_super] <- stau_mu
       tau_nu[i_super] <- stau_nu
       i_super <- i_super + 1 
    }
    
    ## posterior mean of y
    nos_g + prior_y -> post_y
    freq_y <- post_y / sum( post_y) 
    
    list(mu=mu,nu=nu,tau_x=tau_x,tau_mu=tau_mu,tau_nu = tau_nu,
         freq_y = freq_y)  
}

predict_gau <- function (features, out_tr, pred_range, thin)
{
    mu_dim <- dim(out_tr$mu) 
    splsize <- mu_dim[3]
    G <- mu_dim[2]
    k <- mu_dim[1]
    n <- nrow(features)
    
    if( k != ncol(features) ) 
       stop("In 'predict_gau': Test and training data NOT match") 
    
    if(pred_range[1] < 1) pred_range[1] <- 1
    if(pred_range[2] > splsize) pred_range[2] <- splsize    
    
    ## prepare indice of samples used to predict 
    ix_pred <- pred_range[1] + thin * 
       seq(0, floor( ( pred_range[2] - pred_range[1] ) / thin ) ) 
    
    probs_pred <- matrix(0,n,G)    
    
    .C ( "pred_gau", n,k,G, length(ix_pred), features,
         out_tr$mu[,,ix_pred], 
         1/sqrt(out_tr$tau_x[,ix_pred]), log(out_tr$freq_y),
         probs_pred=probs_pred, PACKAGE = "gausspred" )$ probs_pred
}



crossvalid_gau <- function 
      (## arguments specifying data sets and crossvalidation
       no_fold,k, G,features,response, lmd_trf,
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
{
  if(! is.matrix(features) )
    stop("'features' must be a matrix with rows for cases")
        
  n <- nrow(features)
  p <- ncol(features)
  #folding the data sets
  test_sets <- partition_even( 1:n, no_fold )
  probs_pred <- matrix(0,n,G)
  for(fd in seq(1,no_fold)) {
      id_test <- na.omit( test_sets[fd,] )
      id_trn <- (1:n)[-id_test]
      
      if(lmd_trf > 0) {
        if(lmd_trf < 1000) {
          ## transforming data
          L <- gen_transformer(features[id_trn,], response[id_trn], lmd_trf)
          features_L <- features %*% chol2inv(L)
        }
        else 
        {
          sigmas_x <- sqrt( apply(features[id_trn,], 2, pooled_var, 
                                  y = response[id_trn], G = G) )
          features_L <- t( t(features) / sigmas_x )
        }
      }
      else {
        features_L <- features
      }
      
      ## ordering features by F-statistic for selecting features
      if(k < p){
          info_sel <- order_features(features_L[id_trn,],response[id_trn])
          vars_sel <- info_sel$vars[1:k]
          cutoff <- info_sel$fstat[k]
      }
      else {
          vars_sel <- 1:p
          cutoff <- 0 ## here cutoff is only a dummy
      }
      
      ## sampling for lsigmas with slice sampling
      out_sampling <- training_gau ( 
          G, features_L[id_trn,vars_sel], response[id_trn],
          rep(1,G), p_tau_nu, p_tau_mu, p_tau_x,
          nos_super, nos_trans, ini_taus,
          cor, p, cutoff, min_qf, nos_lambda, 
          stepsize_log_tau, no_steps 
      )

      ## making prediction for test cases
      probs_pred[id_test,] <-
      predict_gau(features_L[id_test, vars_sel],out_sampling, pred_range, thin)
  }  
  out_sampling$probs_pred <- probs_pred
  out_sampling
}

gen_transformer <- function(features_tr, y_tr, lmd_trf = 1)
{  
  # estimate pooled covariance matrix
  n <- length(y_tr)
  nos_g <- tapply(rep(1,n),INDEX=y_tr,sum)
  G <- length(nos_g)
  
  sum_cov <- 0
  for(g in seq(1,G) ) {
      sum_cov <- sum_cov + 
      cov(features_tr[which(y_tr==g),]) * (nos_g[g] - 1)
  }
  pooled_cov <- sum_cov / (n - G)
  
  #shrinking pooled covariance
  
  shrinked_cov <- pooled_cov
  
  diag(shrinked_cov) <- diag(shrinked_cov) *(1 + lmd_trf)    
  
  shrinked_cov <- shrinked_cov / ( 1 + lmd_trf )
  
  #computing the choleski decomposition of shrunken pooled covariance
  chol(shrinked_cov) # this linear mulplier is returned.
}

log_sum_exp <- function(lx)
{  mlx <- max(lx)
   log(sum(exp(lx - mlx))) + mlx
}

partition_even <- function(items, g)
{ n <- length(items)
  m <- ceiling(n/g)
  partition_index <- matrix(c(1:n, rep(NA, m * g - n)), nrow=g,ncol=m)
  partition <- matrix(NA,g,m)
  for(i in 1:g)
     partition[i,] <- items[ partition_index[i,] ]

  partition
}

pooled_var <- function(x,y, G)
{   n <- length(y)
    ssq <- function(x) sum( (x-mean(x))^2 )
    ssq_g <- tapply(x,INDEX=y,ssq)
    sum(ssq_g) / (n-G)
}

