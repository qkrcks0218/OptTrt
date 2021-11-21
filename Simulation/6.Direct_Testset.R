############################################
# By running this R-file, the direct OMARs for the test sets are obtained where parameters are chosen from CV.
# The results are saved as "RULE_Test_B##_S###_DR.csv" files in "RULE" folder.
############################################

source("../OptTrt_Sim_Function.R")

DG <- 1:5
TI <- 1:25

for(DATA.GRID in DG){
    
    Data.Test <- read.csv(sprintf("Data/Cl_Test_Data_B%0.4d.csv",DATA.GRID))
    N.test <- dim(Data.Test)[1]
    Data.Training <- read.csv(sprintf("Data/Cl_Data_B%0.4d.csv",DATA.GRID))
    N <- dim(Data.Training)[1]
    pos.X.svm <- 3:12
    
    S.grid <- (45:55)/100
    
    for(S.iter in 1:length(S.grid)){
        
        S <- S.grid[S.iter]
        
        RULE.DR <-  read.csv(sprintf("RULE/RULE_B%0.2d_S%0.3d_DR.csv", DATA.GRID,S.grid[S.iter]*1000))[,-(2:(2*N+1))]
        
        RULE.DR.Test <- matrix(0,N.test,50)
        
        for(iter in 1:25){
            
            CV.RISK <- read.csv(sprintf("CV/Summary_CV_B%0.2d.csv",DATA.GRID))
            Agg.CV <- CV.RISK[CV.RISK$SS.No==iter,]
            
            Agg.CV.DR <- Agg.CV[Agg.CV$Type==3 & Agg.CV$S==S,]
            CV.PARA <- matrix(c(Agg.CV.DR[which.min(Agg.CV.DR$Test1.Adj),]$lambda,
                                Agg.CV.DR[which.min(Agg.CV.DR$Test1.Adj),]$gamma,
                                Agg.CV.DR[which.min(Agg.CV.DR$Test2.Adj),]$lambda,
                                Agg.CV.DR[which.min(Agg.CV.DR$Test2.Adj),]$gamma),1,4,byrow=T)
            
            gamma.DR <- (CV.PARA[1,c(2,4)])
            lambda.DR <- (CV.PARA[1,c(1,3)])
            
            
            col.NF.mat <- colnames(read.csv(sprintf("NF/NF_Train_Adj_B%0.4d_SS%0.4d.csv",DATA.GRID,iter)))[-1]
            NF.mat <- as.matrix(read.csv(sprintf("NF/NF_Train_Adj_B%0.4d_SS%0.4d.csv",DATA.GRID,iter)))[,-1]
            SS.index1 <- which( NF.mat[,col.NF.mat=="SSID"] ==1 )
            SS.index2 <- which( NF.mat[,col.NF.mat=="SSID"] ==2 )
            
            
            Covariate <- NF.mat[,pos.X.svm]
            
            MEAN <- apply(Covariate,2,mean)
            SD <- apply(Covariate,2,sd)
            
            
            Covariate <- scale(Covariate)
            
            Test.Covariate <- Data.Test[,pos.X.svm] 
            Test.Covariate <- (Test.Covariate - matrix(MEAN,dim(Test.Covariate)[1],dim(Test.Covariate)[2],byrow=T))/
                matrix(SD,dim(Test.Covariate)[1],dim(Test.Covariate)[2],byrow=T)
            
              if(sum(RULE.DR[,1]==iter)>0){
                
                pos <- which(RULE.DR[,1]==iter)
                
                K.matrix1 <- GaussianKernel.Cross(Covariate[SS.index1,], Test.Covariate, gamma=gamma.DR[1])
                K.matrix2 <- GaussianKernel.Cross(Covariate[SS.index2,], Test.Covariate, gamma=gamma.DR[2])
                
                alpha0.1 <- RULE.DR[pos,2]
                alphaK.1 <- RULE.DR[pos,2+(1:length(SS.index1))]
                alpha0.2 <- RULE.DR[pos,length(SS.index1)+3]
                alphaK.2 <- RULE.DR[pos,length(SS.index1)+3+(1:length(SS.index2))]
                
                RULE.DR.Test[,iter] <- alpha0.1 + K.matrix1%*%matrix(as.numeric(alphaK.1),length(SS.index1),1)
                RULE.DR.Test[,iter+25] <- alpha0.2 + K.matrix2%*%matrix(as.numeric(alphaK.2),length(SS.index2),1)
                
            } else {
                RULE.DR.Test[,iter] <- RULE.DR.Test[,iter+25] <- NA
            }
            
            print(iter)
            
        }
        
        write.csv(RULE.DR.Test,sprintf( "RULE/RULE_Test_B%0.2d_S%0.3d_DR.csv", DATA.GRID,S*1000),row.names=FALSE)
        
    }
    
}





