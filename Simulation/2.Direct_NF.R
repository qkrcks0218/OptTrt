############################################
# By running this R-file, the outcome regressions and propensity scores are estimated. 
# The results are saved as "NF_Data_B####_SS####.csv" files in "NF" folder.
# The description is written in Section 3.3 of the main paper.
# The code will take a lot of time, so it is recommended to use parallel computing.
# To implement parallel computing, split BATCH=1,...,125 into multiple jobs.
############################################

source("../OptTrt_Sim_Function.R")
source("../MySL.R")

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
SL.hpara$MTRY <- c(2,4)                # random forest parameters
SL.hpara$MLPL <- c(2,4)                # number of nodes in MLP
SL.hpara$NMN <- 25                     # gbm parameter
SL.hpara$MLPdecay <- 10^c(-3,-4,-5)    # MLP decay parameter

pos.AX <- 2:12
pos.X <- 3:12

S.grid <- c(45:55)/100

data.grid <- 1:5
ss.grid <- 1:25
R.seed.mat <- expand.grid(ss.grid, data.grid)

for(BATCH in 1:125){                            ## This will take a lot of time. Recommend to use parallel computing.
    
    data.index <- R.seed.mat[BATCH,2]
    ss.index <- R.seed.mat[BATCH,1]
    
    ############################################
    # Read Data
    ############################################
    
    Data.Cluster <- read.csv(sprintf("Data/Cl_Data_B%0.4d.csv",data.index))
    N <- dim(Data.Cluster)[1]
    M <- Data.Cluster$X5M
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
    
    ############################################
    # SS
    ############################################
    
    check <- 0
    
    while(check==0){
        
        MS.Cl <- sort( sample(1:N,round(N/2)) )
        AS.Cl <- (1:N)[-MS.Cl]
        
        if(length(table(DA[MS.Cl])) == length(table(DA[AS.Cl]))){
            check <- 1
        }
    }
    
    ############################################
    # Nuisance Function Estimation : Takes long time
    ############################################
    
    OR.Fit <- PS.Fit <- list()
    
    OR.Fit$MS <- MySL(Data.Cluster[MS.Cl,], 1, pos.AX, Ydist=gaussian(), 
                      SL.list=SL.hpara$SLL, MTRY=SL.hpara$MTRY, MLPL=SL.hpara$MLPL, NMN=SL.hpara$NMN,MLPdecay=SL.hpara$MLPdecay)
    OR.Fit$AS <- MySL(Data.Cluster[AS.Cl,], 1, pos.AX, Ydist=gaussian(), 
                      SL.list=SL.hpara$SLL, MTRY=SL.hpara$MTRY, MLPL=SL.hpara$MLPL, NMN=SL.hpara$NMN,MLPdecay=SL.hpara$MLPdecay)
    
    PS.formula <- as.formula(paste("as.factor(DA)~(",
                                   paste(colnames(Data.Cluster)[pos.X],collapse="+"),")",
                                   sep=""))
    
    PS.Fit$MS.OL <- polr(PS.formula,data=Data.Cluster.A[MS.Cl,])
    PS.Fit$AS.OL <- polr(PS.formula,data=Data.Cluster.A[AS.Cl,])
    
    
    ############################################
    # Nuisance Function Evaluation
    ############################################
    
    OR.EST <- matrix(0,N,max(M)+1)
    OR.EST.VEC <- rep(0,N)
    
    for(ii in 1:N){
        if( sum(MS.Cl==ii)>0 ){
            TempD <- data.frame( matrix( as.numeric(Data.Cluster[ii,pos.AX]), M[ii]+1, length(pos.AX), byrow=T ) )
            colnames(TempD) <- colnames(Data.Cluster)[pos.AX]
            TempD$A <- seq(0,M[ii],by=1)/M[ii]
            OR.EST[ii,1:(M[ii]+1)] <- predict(OR.Fit$AS,newdata=TempD,onlySL=TRUE)$pred
        } else {
            TempD <- data.frame( matrix( as.numeric(Data.Cluster[ii,pos.AX]), M[ii]+1, length(pos.AX), byrow=T ) )
            colnames(TempD) <- colnames(Data.Cluster)[pos.AX]
            TempD$A <- seq(0,M[ii],by=1)/M[ii]
            OR.EST[ii,1:(M[ii]+1)] <- predict(OR.Fit$MS,newdata=TempD,onlySL=TRUE)$pred
        }
    }
    
    for(ii in 1:N){
        OR.EST.VEC[ii] <- OR.EST[ii,round(Data.Cluster$A[ii]*M[ii]+1)]
    }
    
    
    PMAT <- matrix(0,N,M1)
    PMAT[MS.Cl,] <- predict(PS.Fit$AS.OL,type="probs",newdata=Data.Cluster.A[MS.Cl,])
    PMAT[AS.Cl,]  <- predict(PS.Fit$MS.OL,type="probs",newdata=Data.Cluster.A[AS.Cl,])
    
    PS.EST <- rep(0,N)
    for(ii in 1:N){
        Ab <- round((Data.Cluster.A$A[ii]*M[ii]-1))/M[ii]
        if(Ab>=0){
            pu <- floor(Ab*M1)
            bt <- (pu+2):(DA[ii]+1)
            PS.EST[ii] <- sum(PMAT[ii,(pu+2):(DA[ii]+1)])
        } else {
            PS.EST[ii] <- sum(PMAT[ii,1:(DA[ii]+1)])
        }
    }
    
    ############################################
    # Obtain Starting Point
    ############################################
    
    Ind.Rule.IPW <- Ind.Rule.OR <- Ind.Rule.DR <- Ind.Rule.OR.original <- matrix(0,N,length(S.grid))
    for(ii in 1:N){
        for(ss in 1:length(S.grid)){
            Ind.Rule.IPW[ii,ss] <- IND.OPT.RULE.start(ii,S.grid[ss],k.loss=0,type="IPW",Data.Cluster=Data.Cluster,PS.EST=PS.EST,OR.EST=OR.EST)
            Ind.Rule.OR[ii,ss] <- IND.OPT.RULE.start(ii,S.grid[ss],k.loss=0,type="OR",Data.Cluster=Data.Cluster,PS.EST=PS.EST,OR.EST=OR.EST)
            Ind.Rule.DR[ii,ss] <- IND.OPT.RULE.start(ii,S.grid[ss],k.loss=0,type="DR",Data.Cluster=Data.Cluster,PS.EST=PS.EST,OR.EST=OR.EST)
            Ind.Rule.OR.original[ii,ss] <- IND.OPT.RULE(ii,S.grid[ss],OR.EST=OR.EST)
        }
        print(ii)
    }
    
    ############################################
    # Save
    ############################################
    
    SSvec <- rep(1,N)
    SSvec[AS.Cl] <- 2
    NF.Data <- data.frame(cbind(Data.Cluster,
                                1:N,SSvec,PS.EST,OR.EST.VEC,OR.EST,
                                Ind.Rule.OR.original,Ind.Rule.IPW,Ind.Rule.OR,Ind.Rule.DR))
    colnames(NF.Data) <- c(colnames(Data.Cluster),
                           "GP","SSID","PS.Est","OR.Est",paste("OR.Est.",0:max(M),sep=""),
                           paste("Indirect.OR.Rule",S.grid,sep=""),
                           paste("IPW.Starting",S.grid,sep=""),
                           paste("OR.Starting",S.grid,sep=""),
                           paste("DR.Starting",S.grid,sep=""))
    write.csv(NF.Data,sprintf("NF/NF_Data_B%0.4d_SS%0.4d.csv",R.seed.mat[BATCH,2],R.seed.mat[BATCH,1]),row.names=FALSE)
    
    
}
