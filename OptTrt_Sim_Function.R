library(gtools)
library(MASS)
library(np)


expit <- function(x){
    return( exp(x)/(1+exp(x)) )
}


IND.OPT.RULE <- function(index,S,OR.EST){
    expectedY <- OR.EST[index,1:(M[index]+1)]
    weight <- function(xx){ return( as.numeric(choose(M[index],0:M[index])*rep(xx,M[index]+1)^(0:M[index])*rep(1-xx,M[index]+1)^(M[index] - 0:M[index])) ) }
    value <- function(p){ (sum( expectedY*weight(p) ) - S)^2 }
    thr <- rep(0,11)
    thr[1] <- optimize(value, tol = 10^(-10), lower=0, upper=1)$minimum
    for(jjj in 2:11){
        thr[jjj] <- optim(runif(1),value,method="Brent",lower=0,upper=1,control=list(maxit=10000))$par
    }
    values <- rep(0,11)
    for(jjj in 1:jjj){
        values[jjj] <- value(thr[jjj])
    }
    
    return(thr[which.min(values)])
}


IND.OPT.RULE.start <- function(ii,S,k.loss,type="IPW",Data.Cluster,PS.EST,OR.EST){
    test <- matrix(0,5,2)
    for(kk in 1:dim(test)[1]){
        temp <- optimize(function(xxx){ return(Loss.Entire.scalar(xxx,
                                                                  Data.Cluster$Y[ii],
                                                                  Data.Cluster$A[ii],
                                                                  NULL,
                                                                  M[ii],
                                                                  PS.EST[ii],
                                                                  OR.EST[ii,],
                                                                  type=type,S,k.loss)$loss[1]) }, c(kk/dim(test)[1]-1/dim(test)[1],kk/dim(test)[1]+1/dim(test)[1]))
        test[kk,] <- c(temp$minimum,temp$objective)
    }
    
    return( test[which.min(test[,2]),1] ) 
    
}




IND.OPT.RULE.start.Ext <- function(ii,S,k.loss,type="IPW",Data.Cluster,PS.EST,OR.EST){
    test <- matrix(0,51,2)
    for(kk in 1:dim(test)[1]){
        rrr <- c((kk-26)/5,(kk-25)/5)
        temp <- optimize(function(xxx){ return(Loss.Entire.scalar.Ext(xxx,
                                                                      Data.Cluster$Y[ii],
                                                                      Data.Cluster$A[ii],
                                                                      NULL,
                                                                      M[ii],
                                                                      PS.EST[ii],
                                                                      OR.EST[ii,],
                                                                      type=type,S,k.loss)) }, interval=rrr)
        test[kk,] <- c(temp$minimum,temp$objective)
    }
    
    return( test[which.min(test[,2]),1] ) 
    
}


Loss.Individual <- function(aa,psi,S,k.loss=0){
    
    eta <- 0.001
    
    M <- length(psi)-1
    
    if( abs(aa-0.5)<=0.5 ){
        
        temp.out.neg <- rep(0,M+1)
        temp.out.pos <- rep(0,M+1)
        deri.out.neg <- rep(0,M+1)
        deri.out.pos <- rep(0,M+1)
        
        if(aa==1){
            temp.out.neg.grad <- rep(0,M+1)
            temp.out.pos.grad <- rep(0,M+1)
        }
        
        
        for(tt in 0:M){
            temp.neg <- rep(0,M-tt+1)
            temp.pos <- rep(0,M-tt+1)
            deri.neg <- rep(0,M-tt+1)
            deri.pos <- rep(0,M-tt+1)
            if(aa==1){
                temp.neg.grad <- rep(0,M-tt+1)
                temp.pos.grad <- rep(0,M-tt+1)
            }
            
            for(ll in 0:(M-tt)){
                temp.neg[ll+1] <- choose(M-tt,ll) * max(-psi[tt+1] * (-1)^(ll)*((aa^(ll+tt+1)-k.loss^(ll+tt+1))/(ll+tt+1)),0)
                temp.pos[ll+1] <- choose(M-tt,ll) * max(psi[tt+1] * (-1)^(ll)*((aa^(ll+tt+1)-k.loss^(ll+tt+1))/(ll+tt+1)),0)
                deri.neg[ll+1] <- choose(M-tt,ll) * max(-psi[tt+1] * (-1)^(ll),0) * aa^(ll+tt)
                deri.pos[ll+1] <- choose(M-tt,ll) * max(psi[tt+1] * (-1)^(ll),0) * aa^(ll+tt)
                if(aa==1){
                    temp.neg.grad[ll+1] <- choose(M-tt,ll) * max(-psi[tt+1] * (-1)^(ll),0)
                    temp.pos.grad[ll+1] <- choose(M-tt,ll) * max(psi[tt+1] * (-1)^(ll),0)
                }
            }
            temp.out.neg[tt+1] <- sum(temp.neg)
            temp.out.pos[tt+1] <- sum(temp.pos)
            deri.out.neg[tt+1] <- sum(deri.neg)
            deri.out.pos[tt+1] <- sum(deri.pos)
            if(aa==1){
                temp.out.neg.grad[tt+1] <- sum(temp.neg.grad)
                temp.out.pos.grad[tt+1] <- sum(temp.pos.grad)
            }
        }
        
        output <- list()
        output$neg <- sum(temp.out.neg) + max(S,0)*(aa-k.loss)
        output$pos <- sum(temp.out.pos) + max(-S,0)*(aa-k.loss)
        output$loss <- sum(temp.out.pos-temp.out.neg) - S*(aa-k.loss)
        
        if(abs(aa-0.5)<0.5){
            
            output$deri.neg <- sum(deri.out.neg) + max(S,0)
            output$deri.pos <- sum(deri.out.pos) + max(-S,0)
            output$deri.loss <- output$deri.pos - output$deri.neg
            
        } else if (aa==0) {
            
            output$deri.neg <- 0.5*( sum(deri.out.neg) + max(S,0) - eta )
            output$deri.pos <- 0.5*( sum(deri.out.pos) + max(-S,0) -2*eta )
            output$deri.loss <- output$deri.pos - output$deri.neg
            
        } else if (aa==1) {
            
            Neg1.grad <- sum(temp.out.neg.grad) + max(S,0)
            Pos1.grad <- sum(temp.out.pos.grad) + max(-S,0)
            
            output$deri.neg <- 0.5*( sum(deri.out.neg) + max(S,0) + sum(deri.out.neg) + max(S,0) + eta )
            output$deri.pos <- 0.5*( sum(deri.out.pos) + max(-S,0) + sum(deri.out.pos) + max(-S,0) + 2*eta )
            output$deri.loss <- output$deri.pos - output$deri.neg
            
        }
        
        
        
    } else if (aa<0) {
        
        temp.out.neg <- rep(0,M+1)
        temp.out.pos <- rep(0,M+1)
        temp.out.neg.grad <- rep(0,M+1)
        temp.out.pos.grad <- rep(0,M+1)
        
        for(tt in 0:M){
            temp.neg <- rep(0,M-tt+1)
            temp.pos <- rep(0,M-tt+1)
            temp.neg.grad <- rep(0,M-tt+1)
            temp.pos.grad <- rep(0,M-tt+1)
            for(ll in 0:(M-tt)){
                temp.neg[ll+1] <- choose(M-tt,ll) * max(-psi[tt+1] * (-1)^(ll)*(-k.loss^(ll+tt+1)/(ll+tt+1)),0)
                temp.pos[ll+1] <- choose(M-tt,ll) * max(psi[tt+1] * (-1)^(ll)*(-k.loss^(ll+tt+1)/(ll+tt+1)),0)
                temp.neg.grad[ll+1] <- 0
                temp.pos.grad[ll+1] <- 0
            }
            temp.out.neg[tt+1] <- sum(temp.neg)
            temp.out.pos[tt+1] <- sum(temp.pos)
            temp.out.neg.grad[tt+1] <- sum(temp.neg.grad)
            temp.out.pos.grad[tt+1] <- sum(temp.pos.grad)
        }
        
        Neg1 <- sum(temp.out.neg) + max(S,0)*(-k.loss)
        Pos1 <- sum(temp.out.pos) + max(-S,0)*(-k.loss)
        
        Neg1.grad <- sum(temp.out.neg.grad) + max(S,0)
        Pos1.grad <- sum(temp.out.pos.grad) + max(-S,0)
        
        output <- list()
        output$pos <- Pos1 + max(Pos1.grad,Neg1.grad)*(-aa) + 2*eta*(-aa) 
        output$neg <- Neg1 + max(Pos1.grad,Neg1.grad)*(-aa) + 2*eta*(-aa) - eta + eta*exp(aa)
        output$loss <- output$pos - output$neg
        
        output$deri.neg <- -max(Pos1.grad,Neg1.grad) - 2*eta + eta*exp(aa)
        output$deri.pos <- -max(Pos1.grad,Neg1.grad) - 2*eta 
        output$deri.loss <- output$deri.pos - output$deri.neg
        
    } else if (aa>1) {
        
        temp.out.neg <- rep(0,M+1)
        temp.out.pos <- rep(0,M+1)
        temp.out.neg.grad <- rep(0,M+1)
        temp.out.pos.grad <- rep(0,M+1)
        
        for(tt in 0:M){
            temp.neg <- rep(0,M-tt+1)
            temp.pos <- rep(0,M-tt+1)
            temp.neg.grad <- rep(0,M-tt+1)
            temp.pos.grad <- rep(0,M-tt+1)
            for(ll in 0:(M-tt)){
                temp.neg[ll+1] <- choose(M-tt,ll) * max(-psi[tt+1] * (-1)^(ll)*((1-k.loss^(ll+tt+1))/(ll+tt+1)),0)
                temp.pos[ll+1] <- choose(M-tt,ll) * max(psi[tt+1] * (-1)^(ll)*((1-k.loss^(ll+tt+1))/(ll+tt+1)),0)
                temp.neg.grad[ll+1] <- choose(M-tt,ll) * max(-psi[tt+1] * (-1)^(ll),0)
                temp.pos.grad[ll+1] <- choose(M-tt,ll) * max(psi[tt+1] * (-1)^(ll),0)
            }
            temp.out.neg[tt+1] <- sum(temp.neg)
            temp.out.pos[tt+1] <- sum(temp.pos)
            temp.out.neg.grad[tt+1] <- sum(temp.neg.grad)
            temp.out.pos.grad[tt+1] <- sum(temp.pos.grad)
        }
        
        Neg1 <- sum(temp.out.neg) + max(S,0)*(1-k.loss)
        Pos1 <- sum(temp.out.pos) + max(-S,0)*(1-k.loss)
        
        Neg1.grad <- sum(temp.out.neg.grad) + max(S,0)
        Pos1.grad <- sum(temp.out.pos.grad) + max(-S,0)
        
        output <- list()
        output$pos <- Pos1 + max(Pos1.grad,Neg1.grad)*(aa-1) + 2*eta*(aa-1) 
        output$neg <- Neg1 + max(Pos1.grad,Neg1.grad)*(aa-1) + 2*eta*(aa-1) - eta + eta*exp(1-aa)
        output$loss <- output$pos - output$neg
        
        output$deri.neg <- max(Pos1.grad,Neg1.grad) + 2*eta - eta*exp(1-aa)
        output$deri.pos <- max(Pos1.grad,Neg1.grad) + 2*eta 
        output$deri.loss <- output$deri.pos - output$deri.neg
    }
    
    return( output )
    
}


Loss.Individual.Ext <- function(aa,psi,S,k.loss=0){
    
    M <- length(psi)-1
    
    temp.out.neg <- rep(0,M+1)
    temp.out.pos <- rep(0,M+1)
    deri.out.neg <- rep(0,M+1)
    deri.out.pos <- rep(0,M+1)
    
    
    for(tt in 0:M){
        temp.neg <- rep(0,M-tt+1)
        temp.pos <- rep(0,M-tt+1)
        
        for(ll in 0:(M-tt)){
            temp.neg[ll+1] <- choose(M-tt,ll) * max(-psi[tt+1] * (-1)^(ll)*((aa^(ll+tt+1)-k.loss^(ll+tt+1))/(ll+tt+1)),0)
            temp.pos[ll+1] <- choose(M-tt,ll) * max(psi[tt+1] * (-1)^(ll)*((aa^(ll+tt+1)-k.loss^(ll+tt+1))/(ll+tt+1)),0)
        }
        temp.out.neg[tt+1] <- sum(temp.neg)
        temp.out.pos[tt+1] <- sum(temp.pos)
    }
    
    output <- list()
    output$neg <- sum(temp.out.neg) + max(S,0)*(aa-k.loss)
    output$pos <- sum(temp.out.pos) + max(-S,0)*(aa-k.loss)
    output$loss <- sum(temp.out.pos-temp.out.neg) - S*(aa-k.loss)
    return( output )
    
}




Psi <- function(y,a,x,ps.EST,or.EST,M.temp){               ## A is proportion!
    a <- round(a*M.temp)
    psi.ipw <- rep(0,M.temp+1)
    psi.or <- rep(0,M.temp+1)
    psi.dr <- rep(0,M.temp+1)
    for(aaa in 0:M.temp){
        psi.ipw[aaa+1] <- as.numeric( choose(M.temp,aaa)*y*as.numeric(a==aaa)/ps.EST )
        psi.or[aaa+1] <- as.numeric( choose(M.temp,aaa)*or.EST[aaa+1] )
        psi.dr[aaa+1] <- as.numeric( choose(M.temp,aaa)*((y-or.EST[aaa+1])*as.numeric(a==aaa)/ps.EST + or.EST[aaa+1]) )
    }
    
    psi <- list()
    psi$ipw <- psi.ipw
    psi$or <- psi.or
    psi$dr <- psi.dr
    
    return(psi)
}

GaussianKernel <- function(X,gamma=10^(-0)){
    
    K <- matrix(0,dim(X)[1],dim(X)[1])
    
    X <- matrix(as.numeric(data.matrix(X)),dim(X)[1],dim(X)[2])
    
    for(ii in 1:dim(X)[1]){
        # for(jj in 1:dim(X)[1]){
        #     K[ii,jj] <- sum( (X[ii,]-X[jj,])^2 )
        # }
        
        K[ii,] <- apply( (X - matrix( rep( as.numeric(X[ii,]), dim(X)[1] ), dim(X)[1], dim(X)[2], byrow=T))^2, 1, sum )
    }
    
    K <- exp( -gamma*K )
    
    return(K)
}

GaussianKernel.Cross <- function(X,x2,gamma=10^(-0)){
    
    K <- matrix(0,dim(x2)[1],dim(X)[1])
    
    X <- matrix(as.numeric(data.matrix(X)),dim(X)[1],dim(X)[2])
    
    for(ii in 1:dim(x2)[1]){
        
        # for(jj in 1:dim(X)[1]){
        #     K[ii,jj] <- sum( (x2[ii,]-X[jj,])^2 )
        # }
        
        K[ii,] <- apply( (X - matrix( rep( as.numeric(x2[ii,]), dim(X)[1] ), dim(X)[1], dim(X)[2], byrow=T))^2, 1, sum )
    }
    
    K <- exp( -gamma*K )
    
    return(K)
}



Loss.Entire <- function(alpha,Y,A,X,M,PS.EST,OR.EST,type="IPW",K.type="Gaussian",K.matrix,S,lambda,k.loss){
    if(K.type=="Gaussian") {
        alpha <- matrix(alpha,length(Y),1)
    } else if (K.type=="Linear") {
        alpha <- matrix(alpha,dim(X)[2],1)
    }
    
    result <- matrix(0,length(Y),6)
    colnames(result) <- c("loss","pos","neg","deri.loss","deri.pos","deri.neg")
    for(ii in 1:length(Y)){
        psi <- Psi(Y[ii],A[ii],X[ii,],PS.EST[ii],OR.EST[ii,],M[ii])                ## A is proportion!
        if(K.type=="Gaussian") {
            aa <- K.matrix[ii,]%*%alpha
        } else if (K.type=="Linear") {
            aa <- X[ii,]%*%alpha
        }
        if(type=="IPW"){
            loss.result <- Loss.Individual(aa,psi$ipw,S,k.loss)
        } else if (type=="OR") {
            loss.result <- Loss.Individual(aa,psi$or,S,k.loss)
        } else if (type=="DR") {
            loss.result <- Loss.Individual(aa,psi$dr,S,k.loss)
        }
        
        result[ii,] <- c(loss.result$loss,
                         loss.result$pos,
                         loss.result$neg,
                         loss.result$deri.loss,
                         loss.result$deri.pos,
                         loss.result$deri.neg)
    }
    
    output <- list()
    if(K.type=="Gaussian"){
        penalty <- lambda/2*(t(alpha)%*%K.matrix%*%(alpha))
    } else if (K.type=="Linear"){
        penalty <- lambda/2*(t(alpha)%*%(t(X)%*%X)%*%(alpha))
    }
    output$loss <- apply(result,2,mean)[1:3] + c( penalty , penalty , 0)
    if(K.type=="Gaussian"){
        output$grad.neg <- K.matrix%*%result[,6]/dim(K.matrix)[1]
        output$grad.pos <- K.matrix%*%result[,5]/dim(K.matrix)[1] + lambda*K.matrix%*%(alpha)
        output$grad <- output$grad.pos - output$grad.neg
    } else if (K.type=="Linear"){
        output$grad.neg <- t(X)%*%result[,6]/length(Y)
        output$grad.pos <- t(X)%*%result[,5]/length(Y) + lambda*((t(X)%*%X)%*%alpha)
        output$grad <- output$grad.pos - output$grad.neg
    }
    
    
    return(output)
}




Loss.Entire.intercept <- function(alpha.previous,b.intercept,Y,A,X,M,PS.EST,OR.EST,type="IPW",K.type="Gaussian",K.matrix,S,lambda,k.loss){
    if(K.type=="Gaussian") {
        alpha <- matrix( c(b.intercept,alpha.previous),length(alpha.previous)+1,1)
    } else if (K.type=="Linear") {
        alpha <- matrix( c(b.intercept,alpha.previous),length(alpha.previous)+1,1)
    }
    
    result <- matrix(0,length(Y),6)
    colnames(result) <- c("loss","pos","neg","deri.loss","deri.pos","deri.neg")
    for(ii in 1:length(Y)){
        psi <- Psi(Y[ii],A[ii],X[ii,],PS.EST[ii],OR.EST[ii,],M[ii])                ## A is proportion!
        if(K.type=="Gaussian") {
            aa <- alpha[1] + K.matrix[ii,]%*%alpha[-1]
        } else if (K.type=="Linear") {
            aa <- alpha[1] + X[ii,]%*%alpha[-1]
        }
        if(type=="IPW"){
            loss.result <- Loss.Individual(aa,psi$ipw,S,k.loss)
        } else if (type=="OR") {
            loss.result <- Loss.Individual(aa,psi$or,S,k.loss)
        } else if (type=="DR") {
            loss.result <- Loss.Individual(aa,psi$dr,S,k.loss)
        }
        
        result[ii,] <- c(loss.result$loss,
                         loss.result$pos,
                         loss.result$neg,
                         loss.result$deri.loss,
                         loss.result$deri.pos,
                         loss.result$deri.neg)
    }
    
    output <- list()
    if(K.type=="Gaussian"){
        penalty <- lambda/2*(t(alpha[-1])%*%K.matrix%*%(alpha[-1]))
    } else if (K.type=="Linear"){
        penalty <- lambda/2*(t(alpha[-1])%*%(t(X)%*%X)%*%(alpha[-1]))
    }
    output$loss <- apply(result,2,mean)[1:3] + c( penalty , penalty , 0)
    if(K.type=="Gaussian"){
        output$grad.neg <- rbind(1,K.matrix)%*%result[,6]/dim(K.matrix)[1]
        output$grad.pos <- rbind(1,K.matrix)%*%result[,5]/dim(K.matrix)[1] + lambda*c(0, K.matrix%*%(alpha[-1]))
        output$grad <- output$grad.pos - output$grad.neg
    } else if (K.type=="Linear"){
        output$grad.neg <- rbind(1,t(X))%*%result[,6]/length(Y)
        output$grad.pos <- rbind(1,t(X))%*%result[,5]/length(Y) + lambda*c(0, ((t(X)%*%X)%*%alpha[-1]))
        output$grad <- output$grad.pos - output$grad.neg
    }
    
    
    return(output)
}






Loss.Entire.scalar <- function(aa.vec,Y,A,X,M,PS.EST,OR.EST,type="IPW",S,k.loss){      
    result <- matrix(0,length(Y),6)
    colnames(result) <- c("loss","pos","neg","deri.loss","deri.pos","deri.neg")
    for(ii in 1:length(Y)){
        if(!is.null(dim(OR.EST)) ){
            psi <- Psi(Y[ii],A[ii],X[ii,],PS.EST[ii],OR.EST[ii,],M[ii])     ## A is proportion!
        } else {
            psi <- Psi(Y[ii],A[ii],X[ii,],PS.EST[ii],OR.EST,M[ii])
        }
        aa <- aa.vec[ii]
        if(type=="IPW"){
            loss.result <- Loss.Individual(aa,psi$ipw,S,k.loss)
        } else if (type=="OR") {
            loss.result <- Loss.Individual(aa,psi$or,S,k.loss)
        } else if (type=="DR") {
            loss.result <- Loss.Individual(aa,psi$dr,S,k.loss)
        }
        
        result[ii,] <- c(loss.result$loss,
                         loss.result$pos,
                         loss.result$neg,
                         loss.result$deri.loss,
                         loss.result$deri.pos,
                         loss.result$deri.neg)
    }
    
    output <- list()
    output$loss <- apply(result,2,mean)[1:3]
    output$grad.neg <- apply(result,2,mean)[6]
    output$grad.pos <- apply(result,2,mean)[5]
    output$grad <- apply(result,2,mean)[4]
    
    return(output)
}


Loss.Entire.scalar.Ext <- function(aa.vec,Y,A,X,M,PS.EST,OR.EST,type="IPW",S,k.loss){      
    result <- rep(0,length(Y))
    for(ii in 1:length(Y)){
        if(!is.null(dim(OR.EST)) ){
            psi <- Psi(Y[ii],A[ii],X[ii,],PS.EST[ii],OR.EST[ii,],M[ii])     ## A is proportion!
        } else {
            psi <- Psi(Y[ii],A[ii],X[ii,],PS.EST[ii],OR.EST,M[ii])
        }
        aa <- aa.vec[ii]
        if(type=="IPW"){
            loss.result <- Loss.Individual.Ext(aa,psi$ipw,S,k.loss)
        } else if (type=="OR") {
            loss.result <- Loss.Individual.Ext(aa,psi$or,S,k.loss)
        } else if (type=="DR") {
            loss.result <- Loss.Individual.Ext(aa,psi$dr,S,k.loss)
        }
        
        result[ii] <- c(loss.result$loss)
    }
    
    output <- result
    
    return(output)
}


Winsorization <- function(x){
    return( x*as.numeric(abs(x-0.5)<=0.5) + 0*as.numeric(x<0) + 1*as.numeric(x>1) )
}


Evaluation <- function(RuleA, RuleY){
    TP <- sum(RuleA*RuleY)
    TN <- sum((1-RuleA)*(1-RuleY))
    FP <- sum(RuleA*(1-RuleY))
    FN <- sum((1-RuleA)*RuleY)
    Accuracy <- (TP+TN)/(TP+FP+FN+TN)
    Precision <- TP/(TP+FP)
    Recall <- TP/(TP+FN)
    F1 <- 2*Recall*Precision/(Recall+Precision)
    if(is.nan(F1)){ F1 <- 0 }
    MCC <- (TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
    if(is.nan(MCC)){ MCC <- 0 }
    
    TP <- sum((1-RuleA)*(1-RuleY))
    TN <- sum((RuleA)*(RuleY))
    FP <- sum((1-RuleA)*(RuleY))
    FN <- sum((RuleA)*(1-RuleY))
    Accuracy2 <- (TP+TN)/(TP+FP+FN+TN)
    Precision2 <- TP/(TP+FP)
    Recall2 <- TP/(TP+FN)
    F2 <- 2*Recall2*Precision2/(Recall2+Precision2)
    if(is.nan(F2)){ F2 <- 0 }
    
    R <- matrix(c(TP,TN,FP,FN,Accuracy,Precision,Recall,F1,F2,MCC),1,10)
    R <- data.frame(R)
    colnames(R) <- c("TP (A>OMP,Y>S)","TN (A<OMP,Y<S)","FP (A>OMP,Y<S)","FN  (A<OMP,Y>S)",
                     "Accuracy","Precision","Recall","F1","F2","MCC")
    return(R)
}
