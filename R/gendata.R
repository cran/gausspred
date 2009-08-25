gen_bayesgau <- function(n,p,G,tau_nu=1,tau_mu=1,p_tau_x=c(2,1) )
{
    y <- sort(sample(1:G, n, replace=TRUE))   
    nos_g <- tapply(rep(1,n), INDEX = y, sum)
    
    nu <- rnorm(p) / sqrt(tau_nu)
    mu <- matrix(rnorm(G*p)/sqrt(tau_mu)+nu,G,p,byrow=TRUE)
    sigmas <- 1 / sqrt( rgamma( p, p_tau_x[1]/2, p_tau_x[1]*p_tau_x[2]/2 ) ) 
    
    X <- c()
    for(g in 1:G){
        X <- rbind(X, 
                   matrix( rnorm(nos_g[g]*p) * sigmas + mu[g,], 
                           nos_g[g],p, byrow = TRUE ) 
             )
    }
    
    ix_rn <- sample(1:n)
    
    list( X=X[ix_rn,],y=y[ix_rn], mu = t(mu), sigmas = sigmas, nu = nu )
}

order_features <- function(features, response)
{   fstats <- apply(features,2,comp_fstat,y=response)
    vars <- order(fstats,decreasing=TRUE)
    list(vars=vars, fstats=fstats[vars])
}

comp_fstat <- function(x,y)
{   ssq <- function(x) sum( (x-mean(x))^2 )
    means_g <- tapply(x,INDEX=y,mean)
    ssq_g <- tapply(x,INDEX=y,ssq)
    x_bar <- mean(x)
    n <- length(y)
    nos_g <- tapply(rep(1,n),INDEX=y,sum)
    G <- length(nos_g)
    
    sum(nos_g*(means_g - x_bar)^2) / sum(ssq_g) * (n-G) / (G-1)
}



