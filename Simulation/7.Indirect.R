############################################
# By running this R-file, the indirect OMARs for the test sets using the Lasso and random forest are obtained.
# The results are saved as "NF_Test_LR_Data_B##.csv" files in "NF_LR" folder.
# The code will take a lot of time, so it is recommended to use parallel computing.
# To implement parallel computing, split data.grid=1,...,5 into multiple jobs.
############################################

source("../OptTrt_Sim_Function.R")
source("../MySL.R")

library(glmnet)
library(caret)

pos.AX <- 2:12
pos.X <- 3:12

for(data.grid in 1:5){
    
    S.grid <- (45:55)/100
    Sample.Index <- data.grid
    
    Cl.Data <- read.csv(sprintf("Data/Cl_Data_B%0.4d.csv",Sample.Index))[,1:max(pos.AX)]
    N <- dim(Cl.Data)[1]
    M <- Cl.Data$X5M
    
    # Nuisance function estimation
    
    Lasso.Fit <- RF.Fit <- list()
    FORM <- as.formula("~(A+X1+X2+X3+X4+X5M+U1+U2+U3+U4+U5)^2+
                       I(A^2)+I(X1^2)+I(X2^2)+I(X3^2)+I(X4^2)+I(X5M^2)+
                       I(U1^2)+I(U2^2)+I(U3^3)+I(U4^2)+I(U5^2)")
    FORM.Y <- as.formula("Y~A+X1+X2+X3+X4+X5M+U1+U2+U3+U4+U5")
    
    Lasso.Fit <- cv.glmnet(model.matrix(FORM,data=Cl.Data),Cl.Data$Y)
    RF.Fit <- caret::train( FORM.Y,
                            data=Cl.Data,
                            method="rf",
                            trControl=caret::trainControl(method = "cv",number = 5),
                            verbose=FALSE)
    
    # Test Set
    
    Cl.Data.Test <- read.csv(sprintf("Data/Cl_Test_Data_B%0.4d.csv",Sample.Index))[,1:max(pos.AX)]
    N.Test <- dim(Cl.Data.Test)[1]
    M.Test <- Cl.Data.Test$X5M
    
    Lasso.EST.Test <- RF.EST.Test <- matrix(0,N.Test,max(M.Test)+1)
    Lasso.EST.VEC.Test <- RF.EST.VEC.Test <- rep(0,N.Test)
    
    for(ii in 1:N.Test){
        TempD <- data.frame( matrix( as.numeric(Cl.Data.Test[ii,pos.AX]), M.Test[ii]+1, length(pos.AX), byrow=T ) )
        colnames(TempD) <- colnames(Cl.Data.Test)[pos.AX]
        TempD$A <- seq(0,M.Test[ii],by=1)/M.Test[ii]
        Lasso.EST.Test[ii,1:(M.Test[ii]+1)] <- predict(Lasso.Fit,newx=model.matrix(FORM,data=TempD),s="lambda.min")
        RF.EST.Test[ii,1:(M.Test[ii]+1)] <- predict(RF.Fit,newdata=TempD)
    }
    
    for(ii in 1:N.Test){
        Lasso.EST.VEC.Test[ii] <- Lasso.EST.Test[ii,round(Cl.Data.Test$A[ii]*M.Test[ii]+1)]
        RF.EST.VEC.Test[ii] <- RF.EST.Test[ii,round(Cl.Data.Test$A[ii]*M.Test[ii]+1)]
    }
    
    # Grid Search
    
    Ind.Rule.Lasso <- Ind.Rule.RF <- matrix(0,N.Test,length(S.grid))
    M <- M.Test
    
    for(S.iter in 1:length(S.grid)){
        
        S <- S.grid[S.iter]
        
        for(ii in 1:N.Test){
            Ind.Rule.Lasso[ii,S.iter] <- IND.OPT.RULE.start(ii,S,k.loss=0,type="OR",Data.Cluster=Cl.Data.Test,PS.EST=rep(0,N.Test),OR.EST=Lasso.EST.Test)
            Ind.Rule.RF[ii,S.iter] <- IND.OPT.RULE.start(ii,S,k.loss=0,type="OR",Data.Cluster=Cl.Data.Test,PS.EST=rep(0,N.Test),OR.EST=RF.EST.Test)
        }
        
    }
    
    
    ############################################
    # Save
    ############################################
    
    SSvec <- rep(1,N.Test)
    NF.Data <- data.frame(cbind(1:N.Test,
                                Ind.Rule.Lasso,
                                Ind.Rule.RF))
    colnames(NF.Data) <- c("GP",
                           paste("Lasso.Rule.",S.grid),
                           paste("RF.Rule.",S.grid))
    
    write.csv(NF.Data,sprintf("NF_LR/NF_Test_LR_Data_B%0.2d.csv",Sample.Index),row.names=FALSE)
    
    
    
}


