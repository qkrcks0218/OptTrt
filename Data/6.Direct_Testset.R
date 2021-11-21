############################################
# By running this R-file, the direct OMARs for the test sets are obtained where parameters are chosen from CV.
# The results are saved as "RULE_S###_Test.csv" files in "RULE" folder.
############################################

source("../OptTrt_Sim_Function.R")

Cl.Data <- read.csv("Sen_Training_1417.csv")

N <- dim(Cl.Data)[1]
M <- Cl.Data$C_size

Test.Data <- read.csv("Sen_Test_18.csv")
N.test <- dim(Test.Data)[1]
M.test <- Test.Data$C_size


s.grid <- (64:73)/100

for(ss.iter in 1:10){
    
    S <- s.grid[ss.iter]
    
    pos.X.svm <- 3:dim(Cl.Data)[2]
    
    
    Covariate <- matrix( as.numeric(data.matrix( Cl.Data[,pos.X.svm] )), N, length(pos.X.svm))
    meanC <- apply(Covariate,2,mean)
    sdC <- apply(Covariate,2,sd)
    Covariate <- scale(Covariate)
    
    Test.Covariate <- matrix( as.numeric(data.matrix( Test.Data[,pos.X.svm] )), N.test, length(pos.X.svm))
    Test.Covariate <- (Test.Covariate - matrix(meanC,N.test,length(pos.X.svm),byrow=T))/matrix(sdC,N.test,length(pos.X.svm),byrow=T)
    
    CV.RISK <- read.csv(sprintf("CV/Summary_CV_S%0.3d.csv",S*1000))

    CV.RISK$Train.Adj.avg <- 0.5*(CV.RISK$Train1.Adj+CV.RISK$Train2.Adj)
    CV.RISK$Test.Adj.avg <- 0.5*(CV.RISK$Test1.Adj+CV.RISK$Test2.Adj)


    RULE.DR <-  read.csv(sprintf("RULE/RULE_S%0.3d.csv", S*1000))[,c(1,(2*N+2):(3*N+3))]

    RULE.IPW.Test <- RULE.OR.Test <- RULE.DR.Test <- matrix(0,N.test,200)


    for(iter in 1:100){

        Agg.CV <- CV.RISK[CV.RISK$SS.No==iter,]
        Agg.CV.DR <- Agg.CV[Agg.CV$Type==3,]
        CV.PARA <- matrix(c(Agg.CV.DR[which.min(Agg.CV.DR$Test1.Adj),]$lambda,
                            Agg.CV.DR[which.min(Agg.CV.DR$Test1.Adj),]$gamma,
                            Agg.CV.DR[which.min(Agg.CV.DR$Test2.Adj),]$lambda,
                            Agg.CV.DR[which.min(Agg.CV.DR$Test2.Adj),]$gamma),1,4,byrow=T)
        
        gamma.DR <-   (CV.PARA[1,c(2,4)])
        lambda.DR <-  (CV.PARA[1,c(1,3)])

        col.NF.mat <- colnames(read.csv(sprintf("NF/NF_Train_Adj_B%0.4d.csv",iter)))
        NF.mat <- as.matrix(read.csv(sprintf("NF/NF_Train_Adj_B%0.4d.csv",iter)))
        SS.index1 <- which( NF.mat[,col.NF.mat=="SS.ind"] ==1 )
        SS.index2 <- which( NF.mat[,col.NF.mat=="SS.ind"] ==2 )


        if(sum(RULE.DR[,1]==iter)>0){

            pos <- which(RULE.DR[,1]==iter)

            K.matrix1 <- GaussianKernel.Cross(Covariate[SS.index1,], Test.Covariate, gamma=gamma.DR[1])
            K.matrix2 <- GaussianKernel.Cross(Covariate[SS.index2,], Test.Covariate, gamma=gamma.DR[2])

            alpha0.1 <- RULE.DR[pos,2]
            alphaK.1 <- RULE.DR[pos,2+(1:length(SS.index1))]
            alpha0.2 <- RULE.DR[pos,length(SS.index1)+3]
            alphaK.2 <- RULE.DR[pos,length(SS.index1)+3+(1:length(SS.index2))]

            RULE.DR.Test[,iter] <- alpha0.1 + K.matrix1%*%matrix(as.numeric(alphaK.1),length(SS.index1),1)
            RULE.DR.Test[,iter+100] <- alpha0.2 + K.matrix2%*%matrix(as.numeric(alphaK.2),length(SS.index2),1)

        } else {
            RULE.DR.Test[,iter] <- RULE.DR.Test[,iter+100] <- NA
        }

        print(iter)

    }

    RULE.IPW.Test <- 0
    RULE.OR.Test  <- 0
    RULE.DR.Test  <- RULE.DR.Test[,which(!is.na(apply(RULE.DR.Test,2,sum)))]

    write.csv(RULE.DR.Test,sprintf( "RULE/RULE_S%0.3d_Test.csv", S*1000),row.names=FALSE)
    
}

