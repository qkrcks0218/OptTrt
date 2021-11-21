############################################
# By running this R-file, the direct OMARs for the training sets are obtained where parameters are chosen from CV.
# The results are saved as "RULE_B##_S###_DR.csv" files in "RULE" folder.
# The code will take a lot of time, so it is recommended to use parallel computing.
# To implement parallel computing, split DATA.GRID=1,...,5, thr.grid=1,...,11, and BATCH=1,...,25 into multiple jobs.
############################################

source("../OptTrt_Sim_Function.R")

for(DATA.GRID in 1:5){
    
    for(thr.grid in 1:11){
        
        ss.grid <- 1:25
        loss.grid <- c(3)
        S.grid <- ((45:55)/100)[thr.grid]
        R.seed.mat <- expand.grid(S.grid, loss.grid, ss.grid, DATA.GRID)
        
        RESULT <- data.frame(matrix(0,25,1503))
        
        
        Cl.Data <- read.csv(sprintf("Data/Cl_Data_B%0.4d.csv",DATA.GRID))
        N <- dim(Cl.Data)[1]
        M <- Cl.Data$X5M
        
        pos.AX <- 2:12
        pos.X.svm <- 3:12
        
        CV.RISK <- read.csv(sprintf("CV/Summary_CV_B%0.2d.csv",DATA.GRID))
        
        colnames(RESULT) <- c("NF.iter",
                              paste("Avg_",1:N,sep=""),
                              paste("Sep_",1:N,sep=""),
                              paste("MS_alpha_",0:(250),sep=""),
                              paste("AS_alpha_",0:(250),sep=""))
        
        for(BATCH in 1:25){
            
            DATA.GRID <- R.seed.mat[BATCH,4]
            SS.GRID <- R.seed.mat[BATCH,3]
            TYPE <- if(R.seed.mat[BATCH,2]==1){"IPW"} else if(R.seed.mat[BATCH,2]==2){"OR"} else {"DR"}
            S <- R.seed.mat[BATCH,1]
            
            Agg.CV <- CV.RISK[CV.RISK$SS.No==SS.GRID,]
            Agg.CV.DR <- Agg.CV[Agg.CV$Type==3 & Agg.CV$S==S,]
            CV.PARA <- matrix(c(Agg.CV.DR[which.min(Agg.CV.DR$Test1.Adj),]$lambda,
                                Agg.CV.DR[which.min(Agg.CV.DR$Test1.Adj),]$gamma,
                                Agg.CV.DR[which.min(Agg.CV.DR$Test2.Adj),]$lambda,
                                Agg.CV.DR[which.min(Agg.CV.DR$Test2.Adj),]$gamma),1,4,byrow=T)
            
            gamma.v <- (CV.PARA[1,c(2,4)])
            lambda.v <- (CV.PARA[1,c(1,3)])
            min.CV <- c(min(Agg.CV.DR$Test1.Adj),
                        min(Agg.CV.DR$Test2.Adj))
            
            Covariate <- matrix( as.numeric(data.matrix( Cl.Data[,pos.X.svm] )), N, length(pos.X.svm))
            Covariate <- scale(Covariate)
            
            SVM.Rule.W <- matrix(0,N,2)
            alpha.vec <- list()
            
            col.NF.mat <- colnames(read.csv(sprintf("NF/NF_Train_Adj_B%0.4d_SS%0.4d.csv",DATA.GRID,SS.GRID)))
            NF.mat <- as.matrix(read.csv(sprintf("NF/NF_Train_Adj_B%0.4d_SS%0.4d.csv",DATA.GRID,SS.GRID)))
            
            for(SS.iter in 1:2){
                
                gamma <- gamma.v[SS.iter]
                lambda <- lambda.v[SS.iter]
                ZeroRule <- as.numeric(min.CV[SS.iter]>0)
                
                SS.index1 <- which( NF.mat[,col.NF.mat=="SSID"]==SS.iter )
                SS.index2 <- which( NF.mat[,col.NF.mat=="SSID"]!=SS.iter )
                
                PS.EST <- NF.mat[,which(col.NF.mat=="PS.Est")]
                OR.EST <- NF.mat[,which(col.NF.mat=="OR.Est.0"):which(col.NF.mat=="OR.Est.10")]
                OR.EST.VEC <- NF.mat[,which(col.NF.mat=="OR.Est")]
                
                Ind.Rule.DR <-  NF.mat[,which(col.NF.mat==paste(c("Adj2.DR."),S*1000,sep=""))]
                
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
                
                if(TYPE=="IPW"){
                    b.previous <- mean( Ind.Rule.IPW[SS.index1] )
                    alpha.previous <- solve(K.matrix.approx,Ind.Rule.IPW[SS.index1]-b.previous)
                } else if (TYPE=="OR") {
                    b.previous <- mean( Ind.Rule.OR[SS.index1] )
                    alpha.previous <- solve(K.matrix.approx,Ind.Rule.OR[SS.index1]-b.previous)
                } else if (TYPE=="DR") {
                    b.previous <- mean( Ind.Rule.DR[SS.index1] )
                    alpha.previous <- solve(K.matrix.approx,Ind.Rule.DR[SS.index1]-b.previous)
                }
                
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
            
            SS.MS <- which(NF.mat[,col.NF.mat=="SSID"]==1)
            SS.AS <- which(NF.mat[,col.NF.mat=="SSID"]==2)
            
            SVM.Rule.Avg <- apply(SVM.Rule.W,1,mean)
            pos.SS <- rep(0,N)
            pos.SS[SS.MS] <- 1
            SVM.Rule.Sep <- apply( SVM.Rule.W*cbind(pos.SS,1-pos.SS), 1, sum)
            RESULT[BATCH,] <- c(SS.GRID,
                                SVM.Rule.Avg,
                                SVM.Rule.Sep, 
                                alpha.vec[[1]], 
                                alpha.vec[[2]])
            
        }
        
        write.csv(RESULT,sprintf("RULE/RULE_B%0.2d_S%0.3d_%s.csv",DATA.GRID,S*1000,TYPE),row.names=FALSE)
        
        
    }
}


