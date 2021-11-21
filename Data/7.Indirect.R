############################################
# By running this R-file, the indirect OMARs for the test sets using the Lasso and random forest are obtained.
# The results are saved as "RULE_LR_Test.csv" files in "NF_LR" folder.
############################################

library(glmnet)
library(caret)

source("../OptTrt_Sim_Function.R")
source("../MySL.R")

Data.Cluster <- read.csv("Sen_Training_1417.csv")
Data.test <- read.csv("Sen_Test_18.csv")

S <- (64:73)/100

FORM.Lasso <- as.formula("~(A+C_size+C_rt_C+X_HHsize+X_HHchildsize+X_Nojob+X_head_edu+X_resp_edu+X_resp_age+X_age)^2+
                         I(A^2)+I(C_size^2)+I(C_rt_C^2)+I(X_HHsize^2)+I(X_HHchildsize^2)+
                             I(X_Nojob^2)+I(X_head_edu^2)+I(X_resp_edu^2)+I(X_resp_age^2)+I(X_age^2)")
FORM.RF <- as.formula("Y~A+C_size+C_rt_C+X_HHsize+X_HHchildsize+X_Nojob+X_head_edu+X_resp_edu+X_resp_age+X_age")


N <- dim(Data.Cluster)[1]
M <- Data.Cluster$C_size
N.test <- dim(Data.test)[1]
M.test <- Data.test$C_size

OR.EST <- matrix(0,N,max(M)+1)
OR.EST.VEC <- rep(0,N)

RR <- read.csv(sprintf("NF/NF_Train_Adj_B%0.4d.csv",1))


Lasso.Fit <- cv.glmnet(model.matrix(FORM.Lasso,data=Data.Cluster),Data.Cluster$Y)

RF.Fit <- caret::train( FORM.RF ,
                        data=Data.Cluster,
                        method="rf",
                        trControl=caret::trainControl(method = "cv",number = 5),
                        verbose=FALSE)

Lasso.EST <- RF.EST <- matrix(0,N.test,max(M.test)+1)
Lasso.EST.VEC.test <- RF.EST.VEC.test <- rep(0,N.test)

for(ii in 1:N.test){
    TempD <- data.frame( matrix( as.numeric(Data.test[ii,-1]), M.test[ii]+1, dim(Data.test[ii,-1])[2], byrow=T ) )
    colnames(TempD) <- colnames(Data.test[ii,-1])
    TempD$A <- seq(0,M.test[ii],by=1)/M.test[ii]
    Lasso.EST[ii,1:(M.test[ii]+1)] <- predict(Lasso.Fit,newx=model.matrix(FORM.Lasso,data=TempD),s="lambda.min")
    RF.EST[ii,1:(M.test[ii]+1)] <- predict(RF.Fit,newdata=TempD)
}

for(ii in 1:N.test){
    Lasso.EST.VEC.test[ii] <- Lasso.EST[ii,round(Data.test$A[ii]*M.test[ii]+1)]
    RF.EST.VEC.test[ii] <- RF.EST[ii,round(Data.test$A[ii]*M.test[ii]+1)]
}

M <- M.test

Ind.Rule.Lasso.test <- Ind.Rule.RF.test <- matrix(0,N.test,10)
for(ii in 1:N.test){
    
    for(ss in 1:10){
        Ind.Rule.Lasso.test[ii,ss] <- IND.OPT.RULE.start(ii,S[ss],type="OR",Data.Cluster=Data.test,PS.EST=rep(0,N.test),OR.EST=Lasso.EST,k.loss=0)
        Ind.Rule.RF.test[ii,ss] <- IND.OPT.RULE.start(ii,S[ss],type="OR",Data.Cluster=Data.test,PS.EST=rep(0,N.test),OR.EST=RF.EST,k.loss=0)
    }
}

NF.Data.test <- data.frame(cbind(1:N.test, 
                                 Ind.Rule.Lasso.test,
                                 Ind.Rule.RF.test,
                                 Lasso.EST,
                                 RF.EST))
colnames(NF.Data.test) <- c("GP",
                            paste("Lasso.Rule.",S*1000,sep=""),
                            paste("RF.Rule.",S*1000,sep=""),
                            paste("Lasso.Est.",0:max(M.test),sep=""),
                            paste("RF.Est.",0:max(M.test),sep=""))
write.csv(NF.Data.test,sprintf("NF_LR/RULE_LR_Test.csv"),row.names=FALSE)



