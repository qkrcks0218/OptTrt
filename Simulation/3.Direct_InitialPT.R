############################################
# By running this R-file, the initial values for the DC algorithm are evaluated,
# and are saved as "NF_Train_Adj_B####_SS####.csv" files in "NF" folder.
# The description is written in Section A.3 of the supplementary material.
############################################

source("../OptTrt_Sim_Function.R")

for(Data.Num in 1:5){
    
    Data.Cluster <- read.csv(sprintf("Data/Cl_Data_B%0.4d.csv",Data.Num))
    
    S <- (45:55)/100
    
    N <- dim(Data.Cluster)[1]
    M <- Data.Cluster$X5M
    
    NF.result <- list()
    NF.new.result <- list()
    
    for(ii in 1:25){
        
        NF.result[[ii]] <- read.csv(sprintf("NF/NF_Data_B%0.4d_SS%0.4d.csv",Data.Num,ii))
        
        CL <- colnames(NF.result[[ii]])
        PS.EST <- NF.result[[ii]]$PS.Est
        OR.EST <- NF.result[[ii]][,which(CL=="OR.Est.0"):which(CL=="OR.Est.5")]
        IPW.score <- Data.Cluster$Y/NF.result[[ii]]$PS.Est
        OR.score <- NF.result[[ii]]$OR.Est
        DR.score <- (Data.Cluster$Y - NF.result[[ii]]$OR.Est)/NF.result[[ii]]$PS.Est + NF.result[[ii]]$OR.Est
        
        NewStart1 <- NewStart2 <- NewStart3 <- NewStart4 <- matrix(0,N,3*length(S))
        
        for(tt in 1:3){
            
            if(tt==1){
                TYPE <- "IPW"
                SCORE <- IPW.score
            } else if (tt==2){
                TYPE <- "OR"
                SCORE <- OR.score
            } else {
                TYPE <- "DR"
                SCORE <- DR.score
            }
            
            for(ss in 1:length(S)){
                
                SS <- S[ss]
                
                RULE.name <- paste(TYPE,".Starting",SS,sep="")
                Test.Rule <- round(NF.result[[ii]][,RULE.name],3)
                
                b01 <- which(abs(Test.Rule-0.5)<0.5)
                Test.Rule1 <- Test.Rule2 <- Test.Rule
                
                
                if(length(b01)>1){
                    
                    LM1 <- lm( Test.Rule[b01]~SCORE[b01] )
                    LM2 <- lm( SCORE[b01]~Test.Rule[b01] )
                    
                    
                    U.expand <- which(Test.Rule>=1)
                    Test.Rule1[U.expand] <- apply(cbind(1,LM1$coefficients[1]+LM1$coefficients[2]*SCORE[U.expand]),1,max)
                    
                    L.expand <- which(Test.Rule<=0)
                    Test.Rule1[L.expand] <- apply(cbind(0,LM1$coefficients[1]+LM1$coefficients[2]*SCORE[L.expand]),1,min)
                    
                    U.expand <- which(Test.Rule>=1)
                    Test.Rule2[U.expand] <- apply(cbind(1,(SCORE[U.expand]-LM2$coefficients[1])/LM2$coefficients[2]),1,max)
                    
                    L.expand <- which(Test.Rule<=0)
                    Test.Rule2[L.expand] <- apply(cbind(0,(SCORE[L.expand]-LM2$coefficients[1])/LM2$coefficients[2]),1,min)
                    
                    
                } 
                
                NewStart1[,ss+(tt-1)*length(S)] <- Test.Rule1
                NewStart2[,ss+(tt-1)*length(S)] <- Test.Rule2
                
            }
        }
        
        New.NF <- cbind(NF.result[[ii]], NewStart1, NewStart2)
        New.NF <- data.frame(New.NF)
        colnames(New.NF) <- c(colnames(NF.result[[ii]]),
                              paste("Adj1.IPW.",S*1000,sep=""),
                              paste("Adj1.OR.",S*1000,sep=""),
                              paste("Adj1.DR.",S*1000,sep=""),
                              paste("Adj2.IPW.",S*1000,sep=""),
                              paste("Adj2.OR.",S*1000,sep=""),
                              paste("Adj2.DR.",S*1000,sep=""))
        write.csv(New.NF, sprintf("NF/NF_Train_Adj_B%0.4d_SS%0.4d.csv",Data.Num,ii))
        
        
        print(ii)
        
    }   
    
}

