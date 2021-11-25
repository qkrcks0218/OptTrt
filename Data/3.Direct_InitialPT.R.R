############################################
# By running this R-file, the initial values for the DC algorithm are evaluated,
# and are saved as "NF_Train_Adj_B####.csv" files in "NF" folder.
# The description is written in Section A.3 of the supplementary material.
############################################

source("../OptTrt_Sim_Function.R")

Data.Cluster <- read.csv("Sen_Training_1417.csv")
Test.Data <- read.csv("Sen_Test_18.csv")

s.grid <- (64:73)/100

N <- dim(Data.Cluster)[1]
M <- Data.Cluster$C_size

NF.result <- list()
NF.new.result <- list()


for(ii in 1:100){
    
    NF.result[[ii]] <- read.csv(sprintf("NF/NF_Train_B%0.4d.csv",ii))
    
    CL <- colnames(NF.result[[ii]])
    PS.EST <- NF.result[[ii]]$PS.Est
    OR.EST <- NF.result[[ii]][,which(CL=="OR0.Est"):which(CL=="OR22.Est")]
    IPW.score <- NF.result[[ii]]$Y/NF.result[[ii]]$PS.OL.Est
    OR.score <- NF.result[[ii]]$OR.Est
    DR.score <- (NF.result[[ii]]$Y - NF.result[[ii]]$OR.Est)/NF.result[[ii]]$PS.OL.Est + NF.result[[ii]]$OR.Est
    
    NewStart1 <- NewStart2 <- NewStart3 <- NewStart4 <- matrix(0,N,30)
    NewStart6 <- matrix(0,N,30)
    
    TYPE <- "DR"
    SCORE <- DR.score
    
    for(ss in 1:10){
        
        S <- s.grid[ss]*1000
        
        RULE.name <- paste(TYPE,".Start.",S,sep="")
        Test.Rule <- NF.result[[ii]][,RULE.name]
        
        b01 <- which(abs(Test.Rule-0.5)<0.5)
        Test.Rule1 <- Test.Rule2 <- Test.Rule3 <- Test.Rule4 <- Test.Rule
        Test.Rule6 <- Test.Rule
        
        if(length(b01)>1){
            
            LM1 <- lm( Test.Rule[b01]~SCORE[b01] )
            LM2 <- lm( SCORE[b01]~Test.Rule[b01] )
            
            
            U.expand <- which(Test.Rule==1)
            Test.Rule1[U.expand] <- apply(cbind(1,LM1$coefficients[1]+LM1$coefficients[2]*SCORE[U.expand]),1,max)
            
            L.expand <- which(Test.Rule==0)
            Test.Rule1[L.expand] <- apply(cbind(0,LM1$coefficients[1]+LM1$coefficients[2]*SCORE[L.expand]),1,min)
            
            U.expand <- which(Test.Rule==1)
            Test.Rule2[U.expand] <- apply(cbind(1,(SCORE[U.expand]-LM2$coefficients[1])/LM2$coefficients[2]),1,max)
            
            L.expand <- which(Test.Rule==0)
            Test.Rule2[L.expand] <- apply(cbind(0,(SCORE[L.expand]-LM2$coefficients[1])/LM2$coefficients[2]),1,min)

        } 
        
        NewStart1[,ss+20] <- Test.Rule1
        NewStart2[,ss+20] <- Test.Rule2

    }
    
    New.NF <- cbind(NF.result[[ii]], NewStart1, NewStart2)
    New.NF <- data.frame(New.NF)
    colnames(New.NF) <- c(colnames(NF.result[[ii]]),
                          paste("Adj1.IPW.",s.grid*1000,sep=""),
                          paste("Adj1.OR.",s.grid*1000,sep=""),
                          paste("Adj1.DR.",s.grid*1000,sep=""),
                          paste("Adj2.IPW.",s.grid*1000,sep=""),
                          paste("Adj2.OR.",s.grid*1000,sep=""),
                          paste("Adj2.DR.",s.grid*1000,sep=""))
    write.csv(New.NF, sprintf("NF/NF_Train_aAdj_B%0.4d.csv",ii))
    
    print(ii)
    
}   
