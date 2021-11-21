############################################
# By running this R-file, the outcome regressions and propensity scores are estimated. 
# The results are saved as "NF_Train_Q_B####.csv" files in "NF" folder.
# The description is written in Section 3.3 of the main paper.
# The code will take a lot of time, so it is recommended to use parallel computing.
# To implement parallel computing, split BATCH=1,...,100 into multiple jobs.
############################################

source("../MySL.R")
source("../OptTrt_Sim_Function.R")

Data.Cluster <- read.csv("Sen_Training_1417.csv")
Test.Data <- read.csv("Sen_Test_18.csv")

SL.hpara <- list()
SL.hpara$SLL <- c(1,2,3,4,5,6,7,9,10)
# Superlearner basic learning algorithms: 
# 1: GLM
# 2: lasso/ridge
# 3: earth
# 4: GAM
# 5: xgboost
# 6: polynomial spline
# 7: random forest
# 9: gbm
# 10: 1-layer MLP
SL.hpara$MLPL <- c(2,4)
SL.hpara$MTRY <- c(2,4)
SL.hpara$NMN <- 50
SL.hpara$MLPdecay <- 10^c(-1,-3)


for(BATCH in 1:100){
    
    
    S <- c(64:73)/100
    
    Aq <- c(0,quantile(Data.Cluster$A,c(1:2)/3),1)
    
    N <- dim(Data.Cluster)[1]
    M <- Data.Cluster$C_size
    N.test <- dim(Test.Data)[1]
    M.test <- Test.Data$C_size
    
    DA <- rep(0,N)
    M1 <- max(M)+1
    for(ii in 1:N){
        DA[ii] <- floor(Data.Cluster$A[ii]*M1)
        if(DA[ii]==M1){
            DA[ii] <- max(M)
        }
    }
    
    pos.AX <- 2:dim(Data.Cluster)[2]
    pos.X <- 3:dim(Data.Cluster)[2]
    pos.conti <- 4:dim(Data.Cluster)[2]
    
    Train.X <- Data.Cluster[,pos.conti]
    meanX <- apply(Train.X,2,mean)
    sdX <- apply(Train.X,2,sd)
    
    Data.Cluster[,pos.conti] <- matrix(as.numeric( scale(Train.X) ),N,length(pos.conti))
    Test.Data[,pos.conti] <- (Test.Data[,pos.conti] - matrix(meanX,dim(Test.Data)[1],length(pos.conti),byrow=T))/
        matrix(sdX,dim(Test.Data)[1],length(pos.conti),byrow=T)
    
    
    
    
    PS.EST <- rep(0,N)
    OR.EST <- matrix(0,N,max(M)+1)
    OR.EST.VEC <- rep(0,N)
    
    
    TEST <- 0
    while(TEST==0){
        MS.Cl <- sort( sample(1:N,round(N/2)) )
        AS.Cl <- (1:N)[-MS.Cl]
        
        TEST <- as.numeric(length( unique(DA[MS.Cl]) ) == length( unique(DA[AS.Cl]) ))
        
    }
    
    #############################
    # OR Estimation
    #############################
    
    Fit.MS.Y <- MySL(Data.Cluster[MS.Cl,], 1, pos.AX, Ydist=gaussian(), 
                     SL.list=SL.hpara$SLL, MTRY=SL.hpara$MTRY, MLPL=SL.hpara$MLPL, NMN=SL.hpara$NMN,MLPdecay=SL.hpara$MLPdecay)
    Fit.AS.Y <- MySL(Data.Cluster[AS.Cl,], 1, pos.AX, Ydist=gaussian(), 
                     SL.list=SL.hpara$SLL, MTRY=SL.hpara$MTRY, MLPL=SL.hpara$MLPL, NMN=SL.hpara$NMN,MLPdecay=SL.hpara$MLPdecay)
    
    for(ii in 1:N){
        if( sum(MS.Cl==ii)>0 ){
            TempD <- data.frame( matrix( as.numeric(Data.Cluster[ii,pos.AX]), M[ii]+1, dim(Data.Cluster[ii,pos.AX])[2], byrow=T ) )
            colnames(TempD) <- colnames(Data.Cluster[ii,pos.AX])
            TempD$A <- seq(0,M[ii],by=1)/M[ii]
            OR.EST[ii,1:(M[ii]+1)] <- predict(Fit.AS.Y,newdata=TempD,onlySL=TRUE)$pred
        } else {
            TempD <- data.frame( matrix( as.numeric(Data.Cluster[ii,pos.AX]), M[ii]+1, dim(Data.Cluster[ii,pos.AX])[2], byrow=T ) )
            colnames(TempD) <- colnames(Data.Cluster[ii,pos.AX])
            TempD$A <- seq(0,M[ii],by=1)/M[ii]
            OR.EST[ii,1:(M[ii]+1)] <- predict(Fit.MS.Y,newdata=TempD,onlySL=TRUE)$pred
        }
    }
    
    
    for(ii in 1:N){
        OR.EST.VEC[ii] <- OR.EST[ii,Data.Cluster$A[ii]*M[ii]+1]
    }
    
    #############################
    # PS Estimation
    #############################
    
    PS.Fit <- list()
    
    Data.Cluster.A <- Data.Cluster
    
    DA <- rep(0,N)
    M1 <- max(M)+1
    for(ii in 1:N){
        DA[ii] <- floor(Data.Cluster.A$A[ii]*M1)
        if(DA[ii]==M1){
            DA[ii] <- max(M)
        }
    }
    Data.Cluster.A$DA <- DA
    
    PS.formula <- as.formula(paste("as.factor(DA)~(",
                                   paste(colnames(Data.Cluster)[pos.X],collapse="+"),")",
                                   sep=""))
    
    PS.Fit$MS.OL <- polr(PS.formula,data=Data.Cluster.A[MS.Cl,])
    PS.Fit$AS.OL <- polr(PS.formula,data=Data.Cluster.A[AS.Cl,])
    
    PMAT.OL <- matrix(0,N,M1)
    PMAT.OL[MS.Cl,] <- predict(PS.Fit$AS.OL,type="probs",newdata=Data.Cluster.A[MS.Cl,])
    PMAT.OL[AS.Cl,] <- predict(PS.Fit$MS.OL,type="probs",newdata=Data.Cluster.A[AS.Cl,])
    
    PS.EST.OL <- rep(0,N)
    for(ii in 1:N){
        Ab <- round((Data.Cluster.A$A[ii]*M[ii]-1))/M[ii]
        if(Ab>=0){
            pu <- floor(Ab*M1)
            bt <- (pu+2):(DA[ii]+1)
            PS.EST.OL[ii] <- sum(PMAT.OL[ii,(pu+2):(DA[ii]+1)])
        } else {
            PS.EST.OL[ii] <- sum(PMAT.OL[ii,1:(DA[ii]+1)])
        }
    }
    
    PS.EST.grid.Mat.F.OL <- matrix(0,N,3)
    for(ii in 1:N){
        aa <- floor(Aq[2]*23)
        PS.EST.grid.Mat.F.OL[ii,1] <- sum(PMAT.OL[ii,1:(aa+1)])
        aa <- floor(Aq[3]*23)
        PS.EST.grid.Mat.F.OL[ii,2] <- sum(PMAT.OL[ii,1:(aa+1)])-PS.EST.grid.Mat.F.OL[ii,1]
        PS.EST.grid.Mat.F.OL[ii,3] <- 1-sum(PS.EST.grid.Mat.F.OL[ii,1:2])
    }
    
    
    
    CL.ind <- rep(1,N)
    CL.ind[AS.Cl] <- 2
    
    #############################
    # Rule
    #############################
    
    Ind.Rule.IPW <- Ind.Rule.OR <- Ind.Rule.DR <- matrix(0,N,10)
    
    for(ss in 1:10){
        
        S.temp <- S[ss]
        
        for(ii in 1:N){
            Ind.Rule.DR[ii,ss]  <- round( IND.OPT.RULE.start(ii,S.temp,k.loss=0,type="DR", Data.Cluster=Data.Cluster,PS.EST=PS.EST.OL,OR.EST=OR.EST), 3)
            print(c(ss,ii))
        }
        
    }
    
    IPW.S.name <- paste("IPW.Start.",S*1000 ,sep="")
    OR.S.name  <- paste("OR.Start." ,S*1000 ,sep="")
    DR.S.name  <- paste("DR.Start." ,S*1000 ,sep="")
    
    
    RRR <- data.frame( cbind(Data.Cluster, CL.ind,
                             PS.EST.OL,
                             PS.EST.grid.Mat.F.OL,
                             OR.EST, OR.EST.VEC,
                             Ind.Rule.IPW, Ind.Rule.OR, Ind.Rule.DR) )
    colnames(RRR) <- c(colnames(Data.Cluster),"SS.ind",
                       "PS.OL.Est",paste("PS.OL.Bin",1:dim(PS.EST.grid.Mat.F.OL)[2],".Est",sep=""),
                       paste("OR",0:max(M),".Est",sep=""),"OR.Est",
                       IPW.S.name, OR.S.name, DR.S.name)
    
    write.csv(RRR,sprintf("NF/NF_Train_dB%0.4d.csv",BATCH),row.names=FALSE)
    
}
