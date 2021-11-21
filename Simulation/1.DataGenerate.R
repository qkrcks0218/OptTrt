############################################
# By running this R-file, the data sets used for the simulation are generated. 
# The training and test data sets are saved as 
# "Cl_Data_B####.csv" and "Cl_Test_Data_B####.csv" in "Data" folder, respectively.
############################################

source("../OptTrt_Sim_Function.R")
source("../MySL.R")


for(BATCH in 1:5){
    
    TRUE.OR.Function <- function(a,x){
        0.1*( 4 + 2*(1+0.25*x[,4]^2)*a + 0.25*as.numeric(x[,1]*x[,2]>=1) + 0.25*x[,2]*x[,3] )
    }
    
    N <- 1000 ; N.train <- 500
    M <- sample(c(2:10),N,prob=rep(1,9),replace=TRUE)
    GP <- rep(1:N,M)
    X <- cbind(rnorm(sum(M)),rnorm(sum(M)),rbinom(sum(M),1,0.5),rep(rnorm(N),M),rep(M,M),
               rnorm(sum(M)),rnorm(sum(M)),rbinom(sum(M),1,0.5),rep(rnorm(N),M),rep(rnorm(N),M))
    mX <- aggregate(X~GP,FUN="mean")[,-1]
    
    Basic.GP.PS <- expit( -2 + 1.5*as.numeric(mX[,1]>=0.5)+2.5*mX[,2]^2+2*mX[,3]-2*mX[,4] )
    tA <- rbinom(N,M,Basic.GP.PS)
    A <- tA/M
    
    True.Cl.PS <- choose(M,tA)*Basic.GP.PS^tA*(1-Basic.GP.PS)^(M-tA)
    True.Cl.OR <- TRUE.OR.Function(A,mX)
    
    oY <- True.Cl.OR + rnorm(N)*0.25/sqrt( mX[,5] )
    Data.Cluster <- cbind(oY, A, mX, 1:N)
    colnames(Data.Cluster) <- c("Y", "A", "X1", "X2" , "X3", "X4", "X5M", 
                                "U1","U2","U3","U4","U5",
                                "GP")
    
    ############################################
    # True Rule
    ############################################
    
    pos.AX <- 2:12
    S.grid <- c(45:55)/100
    
    TRUE.OR.EST <- matrix(NA,N,max(M)+1)
    for(ii in 1:N){
        for(jj in 0:M[ii]){
            TRUE.OR.EST[ii,jj+1] <- as.numeric(TRUE.OR.Function(jj/M[ii],Data.Cluster[ii,pos.AX[-1] ]))
        }
    }
    
    TRUE.OPT.RULE <- matrix(0,N,length(S.grid))
    for(ii in 1:N){
        for(ss in 1:length(S.grid)){
            TRUE.OPT.RULE[ii,ss] <- IND.OPT.RULE(index=ii,S=S.grid[ss],OR.EST=TRUE.OR.EST)
        }
        print(c(BATCH,ii))
    }
    
    
    ############################################
    # Save
    ############################################
    
    Data.Cluster.Comb <- cbind(Data.Cluster, True.Cl.PS, True.Cl.OR, TRUE.OR.EST, TRUE.OPT.RULE)
    colnames(Data.Cluster.Comb) <- c( colnames(Data.Cluster), "True.PS", "True.OR" ,
                                      paste("True.OR.",0:max(M),sep=""),
                                      paste("True.Rule.",S.grid,sep=""))
    write.csv(Data.Cluster.Comb[1:N.train,],sprintf("Data/Cl_Data_B%0.4d.csv",BATCH),row.names=FALSE)
    write.csv(Data.Cluster.Comb[-(1:N.train),],sprintf("Data/Cl_Test_Data_B%0.4d.csv",BATCH),row.names=FALSE)
    
}


