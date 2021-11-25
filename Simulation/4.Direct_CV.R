############################################
# By running this R-file, the empirical risks under each parameter are evaluated.
# The results are saved as "Summary_CV_B##.csv" files in "CV" folder.
# The description is written in Section A.4 of the supplementary material.
# The code will take a lot of time, so it is recommended to use parallel computing.
# To implement parallel computing, split DATA.GRID=1,...,5 and BATCH=1,...,825 into multiple jobs.
############################################


source("../OptTrt_Sim_Function.R")


for(DATA.GRID in 1:5){
    
    ss.grid <- 1:25
    loss.grid <- c("DR")
    g.grid <- 10^seq(-1.25,1.25,by=0.25)
    l.grid <- 10^seq(-7,-1,length=3)
    S.grid <- (45:55)/100
    
    R.seed.mat <- expand.grid(l.grid,g.grid,3,ss.grid, DATA.GRID)
    
    
    pos.AX <- 2:12
    pos.X.svm <- 3:12
    
    tot.cv <- 5
    
    Cl.Data <- read.csv(sprintf("Data/Cl_Data_B%0.4d.csv",DATA.GRID))
    N <- dim(Cl.Data)[1]
    M <- Cl.Data$X5M
    
    Train.Error <- matrix(0,dim(R.seed.mat)[1]*length(S.grid),2)
    Test.Error <- matrix(0, dim(R.seed.mat)[1]*length(S.grid),2)
    
    CV.RISK.VALUE <- cbind( rep(R.seed.mat[,5],each=11),
                            rep(R.seed.mat[,4],each=11),
                            rep(R.seed.mat[,3],each=11),
                            rep(R.seed.mat[,1],each=11),
                            rep(R.seed.mat[,2],each=11),
                            rep(S.grid,dim(R.seed.mat)[1]),
                            Train.Error, Test.Error )
                            
    colnames(CV.RISK.VALUE) <- c("Data.No","SS.No","Type","lambda","gamma","S","Train1.Adj","Train2.Adj","Test1.Adj","Test2.Adj")
    
    
    for(BATCH in 1:825){ ## This takes a lot of time. Recommend to use parallel computing.
        
        SS.GRID <- R.seed.mat[BATCH,4]
        TYPE <- if(R.seed.mat[BATCH,3]==1){"IPW"} else if (R.seed.mat[BATCH,3]==2){"OR"} else {"DR"}
        gamma <- R.seed.mat[BATCH,2]
        lambda <- R.seed.mat[BATCH,1]
        
        NF.Data <- read.csv(sprintf("NF/NF_Train_Adj_B%0.4d_SS%0.4d.csv",DATA.GRID,SS.GRID))
        
        PS.EST <- NF.Data$PS.Est
        OR.EST <- cbind( NF.Data$OR.Est.0, NF.Data$OR.Est.1, NF.Data$OR.Est.2,
                         NF.Data$OR.Est.3, NF.Data$OR.Est.4, NF.Data$OR.Est.5,
                         NF.Data$OR.Est.6, NF.Data$OR.Est.7, NF.Data$OR.Est.8, 
                         NF.Data$OR.Est.9, NF.Data$OR.Est.10)
        
        SS <- list()
        SS$MS <- which(NF.Data$SSID==1)
        SS$AS <- which(NF.Data$SSID==2)
        
        
        
        for(P.iter in 1:length(S.grid)){
            
            S <- S.grid[P.iter]
            Ind.Rule.DR <- NF.Data[,paste("Adj2.DR.",S*1000,sep="")]
            
            for(SS.iter in 1:2){
                
                CV.train <- rep(0,tot.cv)
                CV.result <- rep(0,tot.cv)
                if(SS.iter==1){ SS.index1 <- SS$MS ; SS.index2 <- SS$AS } else { SS.index1 <- SS$AS ; SS.index2 <- SS$MS }
                
                Covariate <- matrix( as.numeric(data.matrix( Cl.Data[,pos.X.svm] )), N, length(pos.X.svm))
                Covariate <- scale(Covariate)
                
                fold.list <- sample(SS.index1,length(SS.index1))
                fold.grid <- c(0,50,100,150,200,250)
                
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
                        # loss.entire$loss
                        
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
                        
                        if( (error.value<10^(-6))|iter>5 ){
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
                
                CV.RISK.VALUE[(BATCH-1)*dim(R.seed.mat)[1]+P.iter,6+SS.iter] <- mean(CV.train)
                CV.RISK.VALUE[(BATCH-1)*dim(R.seed.mat)[1]+P.iter,8+SS.iter] <- mean(CV.result)
                
            }
            
            
            
            
        }
        
        
    }
    
    write.csv(CV.RISK.VALUE,sprintf("CV/Summary_CV_B%0.2d.csv",DATA.GRID),row.names=FALSE)
    
}
