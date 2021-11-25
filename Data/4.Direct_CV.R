############################################
# By running this R-file, the empirical risks under each parameter are evaluated.
# The results are saved as "Summary_CV_S###.csv" files in "CV" folder.
# The description is written in Section A.4 of the supplementary material.
# The code will take a lot of time, so it is recommended to use parallel computing.
# To implement parallel computing, split ss.iter=1,...,10 and BATCH=1,...,3300 into multiple jobs.
############################################

source("../OptTrt_Sim_Function.R")

library(gtools)
library(MASS)

Cl.Data <- read.csv("Sen_Training_1417.csv")
Test.Data <- read.csv("Sen_Test_18.csv")
s.grid <- (64:73)/100

N <- dim(Cl.Data)[1]
M <- Cl.Data$C_size

tot.cv <- 5
NF.split <- 5
pos.X.svm <- 3:dim(Cl.Data)[2]

Covariate <- matrix( as.numeric(data.matrix( Cl.Data[,pos.X.svm] )), N, length(pos.X.svm))
Covariate <- scale(Covariate)

for(ss.iter in 1:10){
    
    S <- s.grid[ss.iter]
    
    g.grid <- seq(-2, 0.5, by=0.25)
    l.grid <- c(-7,-3,+1)
    P.mat <-  expand.grid(1:100,10^g.grid, 10^l.grid)
    RULE.NAME <- paste(c("Adj2.IPW.","Adj2.OR.","Adj2.DR."),S*1000,sep="")
    TYPE <- "DR"
    
    Train.Error <- matrix(0,length(P.mat)[1],2)
    Test.Error <- matrix(0,length(P.mat)[1],2)
    
    CV.RISK.VALUE <- data.frame(cbind(P.mat[,1],3,P.mat[,2],P.mat[,3],0,0,0,0))
    colnames(CV.RISK.VALUE) <- c("SS.No","Type","gamma","lambda",
                                 "Train1.Adj","Train2.Adj","Test1.Adj","Test2.Adj")
    
    for(BATCH in 1:3300){
        
        nf.iter <- P.mat[BATCH,1]
        gamma <- P.mat[BATCH,2]
        lambda <- P.mat[BATCH,3]
        
        
        for(SS.iter in 1:2){
            
            
            col.NF.mat <- colnames(read.csv(sprintf("NF/NF_Train_Adj_B%0.4d.csv",nf.iter)))
            NF.mat <- as.matrix(read.csv(sprintf("NF/NF_Train_Adj_B%0.4d.csv",nf.iter)))
            SS.index1 <- which( NF.mat[,col.NF.mat=="SS.ind"]==SS.iter )
            SS.index2 <- which( NF.mat[,col.NF.mat=="SS.ind"]!=SS.iter )
            
            PS.EST <- NF.mat[,which(col.NF.mat=="PS.OL.Est")]
            OR.EST <- NF.mat[,which(col.NF.mat=="OR0.Est"):which(col.NF.mat=="OR22.Est")]
            OR.EST.VEC <- NF.mat[,which(col.NF.mat=="OR.Est")]
            
            Ind.Rule.IPW <- NF.mat[,which(col.NF.mat==RULE.NAME[1])]
            Ind.Rule.OR <-  NF.mat[,which(col.NF.mat==RULE.NAME[2])]
            Ind.Rule.DR <-  NF.mat[,which(col.NF.mat==RULE.NAME[3])]
            
            
            CV.train <- rep(0,tot.cv)
            CV.result <- rep(0,tot.cv)
            
            fold.list <- sample(SS.index1,length(SS.index1))
            fold.grid <- round(seq(1,length(SS.index1),length=tot.cv+1))
            
            valid <- list()
            valid[[1]] <- sort( fold.list[fold.grid[1]:fold.grid[2]] )
            for(ii in 2:tot.cv){
                valid[[ii]] <- sort( fold.list[(fold.grid[ii]+1):fold.grid[ii+1]] )
            }
            
            train <- list()
            for(ii in 1:tot.cv){
                train[[ii]] <- setdiff(SS.index1, valid[[ii]])
            }
            
            for(train.index in 1:tot.cv){
                
                K.matrix <- GaussianKernel.Cross(Covariate[ train[[train.index]], ], 
                                                 Covariate[ train[[train.index]], ],gamma=gamma)
                K.matrix.approx <- K.matrix
                rc.value <- rcond(K.matrix)
                INDICATOR <- 0
                while(rc.value<.Machine$double.eps){
                    K.matrix.approx <- K.matrix.approx + diag( rep(.Machine$double.eps*10^3,dim(K.matrix)[1]) )
                    rc.value <- rcond(K.matrix.approx)
                    INDICATOR <- INDICATOR+1
                    if(INDICATOR>0){ print(paste(c("Approximated K.matrix ",INDICATOR," times"),collapse="")) }
                }
                
                if(TYPE=="IPW"){
                    b.previous <- mean( Ind.Rule.IPW[ train[[train.index]] ] )
                    alpha.previous <- solve(K.matrix.approx,Ind.Rule.IPW[ train[[train.index]] ]-b.previous)
                } else if (TYPE=="OR") {
                    b.previous <- mean( Ind.Rule.OR[ train[[train.index]] ] )
                    alpha.previous <- solve(K.matrix.approx,Ind.Rule.OR[ train[[train.index]] ]-b.previous)
                } else if (TYPE=="DR") {
                    b.previous <- mean( Ind.Rule.DR[ train[[train.index]] ] )
                    alpha.previous <- solve(K.matrix.approx,Ind.Rule.DR[ train[[train.index]] ]-b.previous)
                }
                
                Y.train <- Cl.Data$Y[ train[[train.index]] ]
                A.train <- Cl.Data$A[ train[[train.index]] ]
                X.train <- Covariate[ train[[train.index]], ]
                M.train <- M[ train[[train.index]] ]
                PS.EST.train <- PS.EST[ train[[train.index]] ]
                OR.EST.train <- OR.EST[ train[[train.index]], ]
                
                Y.valid <- Cl.Data$Y[ valid[[train.index]] ]
                A.valid <- Cl.Data$A[ valid[[train.index]] ]
                X.valid <- Covariate[ valid[[train.index]], ]
                M.valid <- M[ valid[[train.index]] ]
                PS.EST.valid <- PS.EST[ valid[[train.index]] ]
                OR.EST.valid <- OR.EST[ valid[[train.index]], ]
                
                
                error.alpha <- 1
                error.value <- 1
                iter <- 0
                convergence <- FALSE
                
                K.TYPE <- "Gaussian"
                
                print("Ready for DC Algorithm")
                
                while(!convergence){
                    
                    loss.entire <- Loss.Entire.intercept(alpha.previous,
                                                         b.previous,
                                                         Y.train,
                                                         A.train,
                                                         X.train,
                                                         M.train,
                                                         PS.EST.train,
                                                         OR.EST.train,
                                                         type=TYPE,K.type=K.TYPE,K.matrix,S,lambda=lambda,k.loss=0)
                    beta.next <- loss.entire$grad.neg
                    
                    Step2 <- optim(c(b.previous,alpha.previous), 
                                   fn=function(aaa){ return( Loss.Entire.intercept(aaa[-1],
                                                                                   aaa[1],
                                                                                   Y.train,
                                                                                   A.train,
                                                                                   X.train,
                                                                                   M.train,
                                                                                   PS.EST.train,
                                                                                   OR.EST.train,
                                                                                   type=TYPE,K.type=K.TYPE,K.matrix,S,lambda=lambda,k.loss=0 )$loss[2] - t(beta.next)%*%matrix(aaa,length(aaa),1) ) },
                                   gr=function(aaa){ return( Loss.Entire.intercept(aaa[-1],
                                                                                   aaa[1],
                                                                                   Y.train,
                                                                                   A.train,
                                                                                   X.train,
                                                                                   M.train,
                                                                                   PS.EST.train,
                                                                                   OR.EST.train,
                                                                                   type=TYPE,K.type=K.TYPE,K.matrix,S,lambda=lambda,k.loss=0 )$grad.pos - beta.next ) },
                                   method="L-BFGS-B")
                    
                    b.next <- as.numeric(Step2$par)[1]
                    alpha.next <- as.numeric(Step2$par)[-1]
                    
                    temp.loss <- Loss.Entire.intercept(alpha.next,
                                                       b.next,
                                                       Y.train,
                                                       A.train,
                                                       X.train,
                                                       M.train,
                                                       PS.EST.train,
                                                       OR.EST.train,
                                                       type=TYPE,K.matrix,K.type=K.TYPE,S,lambda=lambda,k.loss=0)
                    
                    error.alpha <- sqrt( (b.previous - b.next)^2 + sum((alpha.previous-alpha.next)^2) )
                    error.value <- loss.entire$loss[1] - temp.loss$loss[1]
                    
                    b.previous <- b.next
                    alpha.previous <- alpha.next
                    
                    iter <- iter+1
                    
                    # print(paste("Iter: ",iter,", log10(gamma): ",log(gamma,10),", log10(lambda): ",log(lambda,10), 
                    #             ", Error.alpha (10^6): ",round(error.alpha*10^6,4),", Error.value (10^6): ",round(error.value*10^6,4),
                    #             ", Loss: ",temp.loss$loss[1], sep=""))
                    
                    if( (error.value<10^(-6))|iter>100 ){
                        convergence <- TRUE
                    }
                    
                    
                    
                }
                
                
                SVM.Rule.train <- Winsorization( b.next + K.matrix%*%alpha.next )
                
                K.matrix.valid <- GaussianKernel.Cross(Covariate[ train[[train.index]], ], Covariate[ valid[[train.index]], ],gamma=gamma)
                
                CV.train[train.index] <-  Loss.Entire.scalar(SVM.Rule.train,
                                                             Y.train,
                                                             A.train,
                                                             X.train,
                                                             M.train,
                                                             PS.EST.train,
                                                             OR.EST.train,
                                                             type=TYPE,S,k.loss=0)$loss[1]
                
                SVM.Rule.valid <- Winsorization(b.next + K.matrix.valid%*%alpha.next)
                
                CV.result[train.index] <-  Loss.Entire.scalar(SVM.Rule.valid,
                                                              Y.valid,
                                                              A.valid,
                                                              X.valid,
                                                              M.valid,
                                                              PS.EST.valid,
                                                              OR.EST.valid,
                                                              type=TYPE,S,k.loss=0)$loss[1]
                
                
            }
            
            CV.RISK.VALUE[BATCH,4+SS.iter] <- mean(CV.train)
            CV.RISK.VALUE[BATCH,6+SS.iter] <- mean(CV.result)
        }
        
    }
    
    
    write.csv(CV.RISK.VALUE,sprintf("CV/Summary_CV_S%0.3d.csv",S*1000),row.names=FALSE)
    
}



