############################################
# By running this R-file, the direct OMARs for the training sets are obtained where parameters are chosen from CV.
# The results are saved as "RULE_S###_Test.csv" files in "RULE" folder.
# The code will take a lot of time, so it is recommended to use parallel computing.
# To implement parallel computing, split ss.iter=1,...,10, nf.iter=1,...,100 into multiple jobs.
############################################

source("../OptTrt_Sim_Function.R")

Cl.Data <- read.csv("Sen_Training_1417.csv")

N <- dim(Cl.Data)[1]
M <- Cl.Data$C_size

s.grid <- (64:73)/100

for(ss.iter in 1:10){
    
    S <- s.grid[ss.iter]
    
    TYPE <- "DR"
    
    pos.X.svm <- 3:dim(Cl.Data)[2]
    
    CV.RISK <- read.csv(sprintf("CV/Summary_CV_S%0.3d.csv",S*1000))
    
    CV.RISK$Train.Adj.avg <- 0.5*(CV.RISK$Train1.Adj+CV.RISK$Train2.Adj)
    CV.RISK$Test.Adj.avg <- 0.5*(CV.RISK$Test1.Adj+CV.RISK$Test2.Adj)
    
    
    RESULT <- data.frame(matrix(0,100,3*N+3))
    
    colnames(RESULT) <- c("NF.iter",
                          paste("Avg_",1:N,sep=""),
                          paste("Sep_",1:N,sep=""),
                          paste("MS_alpha_",0:(514),sep=""),
                          paste("AS_alpha_",0:(513),sep=""))
    
    
    for(nf.iter in 1:100){
        
        Agg.CV <- CV.RISK[CV.RISK$SS.No==nf.iter,]
        
        
        Agg.CV.DR <- Agg.CV[Agg.CV$Type==3,]
        CV.PARA <- matrix(c(Agg.CV.DR[which.min(Agg.CV.DR$Test1.Adj),]$lambda,
                            Agg.CV.DR[which.min(Agg.CV.DR$Test1.Adj),]$gamma,
                            Agg.CV.DR[which.min(Agg.CV.DR$Test2.Adj),]$lambda,
                            Agg.CV.DR[which.min(Agg.CV.DR$Test2.Adj),]$gamma),1,4,byrow=T)
        
        gamma.v <- (CV.PARA[1,c(2,4)])
        lambda.v <- (CV.PARA[1,c(1,3)])
        min.CV <- c(min(Agg.CV.DR$Test1.Adj),
                    min(Agg.CV.DR$Test2.Adj))
        
        RULE.NAME <- paste(c("Adj2.IPW.","Adj2.OR.","Adj2.DR."),S*1000,sep="")
        
        
        
        Covariate <- matrix( as.numeric(data.matrix( Cl.Data[,pos.X.svm] )), N, length(pos.X.svm))
        Covariate <- scale(Covariate)
        
        SVM.Rule.W <- matrix(0,N,2)
        alpha.vec <- list()
        
        col.NF.mat <- colnames(read.csv(sprintf("NF/NF_Train_Adj_B%0.4d.csv",nf.iter)))
        NF.mat <- as.matrix(read.csv(sprintf("NF/NF_Train_Adj_B%0.4d.csv",nf.iter)))
        
        for(SS.iter in 1:2){
            
            gamma <- gamma.v[SS.iter]
            lambda <- lambda.v[SS.iter]
            ZeroRule <- as.numeric(min.CV[SS.iter]>0)
            
            SS.index1 <- which( NF.mat[,col.NF.mat=="SS.ind"]==SS.iter )
            SS.index2 <- which( NF.mat[,col.NF.mat=="SS.ind"]!=SS.iter )
            
            PS.EST <- NF.mat[,which(col.NF.mat=="PS.OL.Est")]
            OR.EST <- NF.mat[,which(col.NF.mat=="OR0.Est"):which(col.NF.mat=="OR22.Est")]
            OR.EST.VEC <- NF.mat[,which(col.NF.mat=="OR.Est")]
            
            Ind.Rule.IPW <- NF.mat[,which(col.NF.mat==RULE.NAME[1])]
            Ind.Rule.OR <-  NF.mat[,which(col.NF.mat==RULE.NAME[2])]
            Ind.Rule.DR <-  NF.mat[,which(col.NF.mat==RULE.NAME[3])]
            
            K.matrix <- GaussianKernel.Cross(Covariate[SS.index1,], Covariate[SS.index1,], gamma=gamma)
            K.matrix.approx <- K.matrix
            rc.value <- rcond(K.matrix)
            INDICATOR <- 0
            while(rc.value<.Machine$double.eps){
                K.matrix.approx <- K.matrix.approx + diag( rep(.Machine$double.eps*10^3,dim(K.matrix)[1]) )
                rc.value <- rcond(K.matrix.approx)
                INDICATOR <- INDICATOR+1
                if(INDICATOR>0){ print(paste(c("Approximated K.matrix ",INDICATOR," times"),collapse="")) }
            }
            
            
            b.previous <- mean( Ind.Rule.DR[SS.index1] )
            alpha.previous <- solve(K.matrix.approx,Ind.Rule.DR[SS.index1]-b.previous)
            
            Y <- Cl.Data$Y
            A <- Cl.Data$A
            
            
            error.alpha <- 1
            error.value <- 1
            iter <- 0
            convergence <- FALSE
            
            K.TYPE <- "Gaussian"
            
            print("Ready for DC Algorithm")
            
            while(!convergence){
                
                loss.entire <- Loss.Entire.intercept(alpha.previous,
                                                     b.previous,
                                                     Y[SS.index1],
                                                     A[SS.index1],
                                                     Covariate[SS.index1,],
                                                     M[SS.index1],
                                                     PS.EST[SS.index1],
                                                     OR.EST[SS.index1,],
                                                     type=TYPE,K.type=K.TYPE,K.matrix,S,lambda=lambda,k.loss=0)
                # loss.entire$loss
                
                beta.next <- loss.entire$grad.neg
                
                Step2 <- optim(c(b.previous,alpha.previous), 
                               fn=function(aaa){ return( Loss.Entire.intercept(aaa[-1],
                                                                               aaa[1],
                                                                               Y[SS.index1],
                                                                               A[SS.index1],
                                                                               Covariate[SS.index1,],
                                                                               M[SS.index1],
                                                                               PS.EST[SS.index1],
                                                                               OR.EST[SS.index1,],
                                                                               type=TYPE,K.type=K.TYPE,K.matrix,S,lambda=lambda,k.loss=0 )$loss[2] - t(beta.next)%*%matrix(aaa,length(aaa),1) ) },
                               gr=function(aaa){ return( Loss.Entire.intercept(aaa[-1],
                                                                               aaa[1],
                                                                               Y[SS.index1],
                                                                               A[SS.index1],
                                                                               Covariate[SS.index1,],
                                                                               M[SS.index1],
                                                                               PS.EST[SS.index1],
                                                                               OR.EST[SS.index1,],
                                                                               type=TYPE,K.type=K.TYPE,K.matrix,S,lambda=lambda,k.loss=0 )$grad.pos - beta.next ) },
                               method="L-BFGS-B")
                
                b.next <- as.numeric(Step2$par)[1]
                alpha.next <- as.numeric(Step2$par)[-1]
                
                temp.loss <- Loss.Entire.intercept(alpha.next,
                                                   b.next,
                                                   Y[SS.index1],
                                                   A[SS.index1],
                                                   Covariate[SS.index1,],
                                                   M[SS.index1],
                                                   PS.EST[SS.index1],
                                                   OR.EST[SS.index1,],
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
            
            K.Full <- GaussianKernel.Cross(Covariate[SS.index1,], Covariate, gamma=gamma)
            SVM.Rule.W[,SS.iter] <- b.next + K.Full%*%alpha.next
            alpha.vec[[SS.iter]] <- c(b.next,alpha.next)
            
        }
        
        SS.MS <- which(NF.mat[,col.NF.mat=="SS.ind"]==1)
        SS.AS <- which(NF.mat[,col.NF.mat=="SS.ind"]==2)
        
        SVM.Rule.Avg <- apply(SVM.Rule.W,1,mean)
        pos.SS <- rep(0,N)
        pos.SS[SS.MS] <- 1
        SVM.Rule.Sep <- apply( SVM.Rule.W*cbind(pos.SS,1-pos.SS), 1, sum)
        RESULT[nf.iter,] <- c(nf.iter, SVM.Rule.Avg, SVM.Rule.Sep, alpha.vec[[1]], alpha.vec[[2]])
        
    }
    
    write.csv(RESULT,sprintf("RULE/RULE_S%0.3d.csv",S*1000),row.names=FALSE)
    
    
}

