
comp_amlp <- function(probs_pred, responses)
{   mlp <- rep(0, nrow(probs_pred))
    for(i in 1:nrow(probs_pred)) 
        mlp[i] <- log(probs_pred[i,responses[i]])
    - mean(mlp)
}

comp_er <- function(probs_pred,responses)
{  
    ypred <- apply(probs_pred,1,which.max)
    mean(ypred!=responses)
}


comp_loss <- function(probs_pred, y_true, Mloss)
{
     loss_pred <- probs_pred %*% Mloss
     y_pred <- apply(loss_pred,1,which.min)
     loss_exp <- mean( apply(loss_pred,1,min) )
     loss <- 0
     for(i in 1:nrow(probs_pred)) {
         loss <- loss + Mloss[y_true[i],y_pred[i]]
     } 
     c(loss=loss / nrow(probs_pred), eloss=loss_exp)
}

cal_tab <- function(probs_pred, true_y, ix_y, no_cat=10)
{
    category <- character(no_cat)
    no <- meanprob <- prop1 <- rep(0,no_cat)
    n <- length(true_y)
    for( i in seq(1,no_cat) )
    {
        category[i]=paste((i-1)/no_cat,"-",i/no_cat,sep="")
        
        label <- ( (probs_pred[,ix_y] <= i/no_cat) & 
                   (probs_pred[,ix_y] > (i-1)/no_cat) )
        
        meanprob[i] <- mean(probs_pred[label,ix_y])
        no[i] <- sum(label)
        prop1[i]= mean(true_y[label] == ix_y)
    }
    data.frame(category,no,meanprob,prop1)    
}
