############################################
# By running this R-file, the simulation results are summarized and the plots are drawn.
############################################

source("../OptTrt_Sim_Function.R")

Rule.Our.Ind <- Rule.Our.DR <- True.Rule <- Rule.Our.Lasso <- Rule.Our.RF <- list()

DG <- 1:5
TI <- 1:25

for(DATA.GRID in DG){
  
  print(DATA.GRID)
  Cl.Data <- read.csv(sprintf("Data/Cl_Test_Data_B%0.4d.csv",DATA.GRID))
  N <- dim(Cl.Data)[1]
  M <- Cl.Data$X5M
  
  S.grid <- (45:55)/100
  
  True.Rule[[DATA.GRID]] <- matrix(0,N,length(S.grid))
  
  for(S.iter in 1:length(S.grid)){
    True.Rule[[DATA.GRID]][,S.iter] <- Cl.Data[,paste("True.Rule.",S.grid[S.iter],sep="")]
  }
  
  
  Rule.Our.DR[[DATA.GRID]] <- matrix(0,N,length(S.grid))
  
  
  for(S.iter in 1:length(S.grid)){
    
    T3 <- read.csv(sprintf("RULE/RULE_Test_B%0.2d_S%0.3d_DR.csv", DATA.GRID,S.grid[S.iter]*1000))
    
    Rule.Our.DR[[DATA.GRID]][,S.iter] <-  apply(T3[,c(TI,TI+25)],1,median)
    
  }
  
  LR.Result <- read.csv(sprintf("NF_LR/NF_Test_LR_Data_B%0.2d.csv",DATA.GRID))[,-1]
  Rule.Our.Lasso[[DATA.GRID]] <- Winsorization(round(LR.Result[,1:length(S.grid)],3))
  Rule.Our.RF[[DATA.GRID]] <- Winsorization(round(LR.Result[,length(S.grid)+1:length(S.grid)],3))
  
}

Our.DR <- Our.Lasso <- Our.RF <- matrix(0,length(DG)*N,length(S.grid))
for(DATA.GRID in DG){
  for(ss in 1:length(S.grid)){
    Our.DR[(DATA.GRID-1)*N+1:N,ss]      <- abs( Winsorization(Rule.Our.DR[[DATA.GRID]][,ss]) - Winsorization(True.Rule[[DATA.GRID]][,ss]) )
    Our.Lasso[(DATA.GRID-1)*N+1:N,ss]   <- abs( Winsorization(Rule.Our.Lasso[[DATA.GRID]][,ss]) - Winsorization(True.Rule[[DATA.GRID]][,ss]) )
    Our.RF[(DATA.GRID-1)*N+1:N,ss]      <- abs( Winsorization(Rule.Our.RF[[DATA.GRID]][,ss]) - Winsorization(True.Rule[[DATA.GRID]][,ss]) )
    
  }
}



## Plot 


png("Plot/Figure1_3.png",width=11,height=2.75,unit="in",res=500)

layout(matrix(c(1,2),1,2),widths = c(4,1))
par(mar=c(3,3.5,1,0))

MGR <- matrix(0,dim(Our.DR)[1],dim(Our.DR)[2]*3)
XL <- rep(0,dim(Our.DR)[2]*3)
for(jj in 1:dim(Our.DR)[2]){
  MGR[,3*jj-2] <- abs(Our.DR)[,jj]
  MGR[,3*jj-1] <- abs(Our.Lasso)[,jj]
  MGR[,3*jj-0] <- abs(Our.RF)[,jj]
  XL[3*jj+(-2):0] <- c(-1,0,1) + 6*jj
}

COL <- c("red",rgb(0,206/255,206/255,alpha=0.8),	rgb(255/255,165/255,0,alpha=0.8))

boxplot(MGR, at=XL, width=rep(0.02,length(XL)), boxlty=0,
        col=COL,
        whiskcol=COL,
        staplecol=COL,
        outcol=COL,
        medlwd=1,
        cex=0.1,pch=19,
        names=rep("",33),ylab="",line=2,ylim=c(0,1), xaxt='n')
axis(1,at=(1:11)*6,label=(45:55)/100)
title(xlab="T",line=2)
title(ylab=expression("|True OMAR - Estimated OMAR|"),line=2.5,cex.lab=0.9)
abline(v=(2:11)*6-3,col=rgb(0,0,1,0.25),lty=3)

par(mar=c(3,1,1,0))
plot.new()
points(rep(0.05,3),seq(9,5,by=-2)/10,pch=15,cex=2,col=COL)
text(rep(0.15,3),seq(9,5,by=-2)/10,pos=4,
     c("Direct","Indirect-Lasso","Indirect-RF"))

dev.off()








Value <- list()

S.grid <- (45:55)/100
DG <- 1:5

for(ss in 1:11){
  
  RR <- matrix(0,4,10)
  
  for(DATA.GRID in DG){
    
    S <- S.grid[ss]
    
    T.DR <-    rep(0,N)
    T.Lasso <- rep(0,N)
    T.RF <-    rep(0,N)
    T.True <-  rep(0,N)
    T.Y <-     rep(0,N)
    
    
    
    Cl.Data <- read.csv(sprintf("Data/Cl_Test_Data_B%0.4d.csv",DATA.GRID))
    M <- Cl.Data$X5M
    
    T.DR  <-    as.numeric(Winsorization(Rule.Our.DR[[DATA.GRID]][,ss]) < Cl.Data$A)
    T.Lasso <- as.numeric(Winsorization(Rule.Our.Lasso[[DATA.GRID]][,ss])   < Cl.Data$A)
    T.RF <-    as.numeric(Winsorization(Rule.Our.RF[[DATA.GRID]][,ss])      < Cl.Data$A)
    T.True <-  as.numeric(Winsorization(True.Rule[[DATA.GRID]][,ss])    < Cl.Data$A)
    T.Y <- as.numeric(S < Cl.Data$Y)
    
    
    
    RR <- RR + matrix(as.numeric(rbind(c( Evaluation(T.DR,    T.Y)),
                                       c( Evaluation(T.Lasso, T.Y)),
                                       c( Evaluation(T.RF,    T.Y)),
                                       c( Evaluation(T.True,  T.Y)) )),4,10)
    
  }
  RR <- RR/length(DG)
  Value[[ss]] <- matrix(as.numeric(RR),4,10)
}


RESULT <- data.frame(cbind(rep(1:4,11),rbind(Value[[1]],Value[[2]],Value[[3]],
                                             Value[[4]],Value[[5]],Value[[6]],
                                             Value[[7]],Value[[8]],Value[[9]],
                                             Value[[10]],Value[[11]])))
colnames(RESULT) <- c("Method","","","","","Accuracy","","","F1","F2","MCC")

COL2 <- c("red",rgb(0,206/255,206/255,alpha=0.8),	rgb(255/255,165/255,0,alpha=0.8),rgb(0,0,0,0.5))
LTY <- c(1,5,5,1)
LWD <- c(1.75,1.75,1.75,1.75)

png("Plot/Figure2_3.png",width=8,height=2.0,unit="in",res=500)

layout(matrix(c(1,2,3,4,5),1,5,byrow=T),widths = c(2,2,2,2,1.5))
par(mar=c(4,2.5,2,0.5))

plot(S.grid,RESULT$Accuracy[RESULT$Method==1],type='l',ylim=c(0.725,0.8),
     col=COL2[1],lty=LTY[1],lwd=LWD[1],
     xlab="",ylab="",main="")
title(xlab="",ylab="",line=2)
title(main="Accuracy",font.main=1)
par(new=T)
points(S.grid,RESULT$Accuracy[RESULT$Method==2],type='l',col=COL2[2],lty=LTY[2],lwd=LWD[2]); par(new=T)
points(S.grid,RESULT$Accuracy[RESULT$Method==3],type='l',col=COL2[3],lty=LTY[3],lwd=LWD[3]); par(new=T)
points(S.grid,RESULT$Accuracy[RESULT$Method==4],type='l',col=COL2[4],lty=LTY[4],lwd=LWD[4])


plot(S.grid,RESULT$F1[RESULT$Method==1],type='l',ylim=c(0.675,0.8),
     col=COL2[1],lty=LTY[1],lwd=LWD[1],
     xlab="",ylab="",main="")
title(xlab="",ylab="",line=2)
title(main="F1",font.main=1)
par(new=T)
points(S.grid,RESULT$F1[RESULT$Method==2],type='l',col=COL2[2],lty=LTY[2],lwd=LWD[2]); par(new=T)
points(S.grid,RESULT$F1[RESULT$Method==3],type='l',col=COL2[3],lty=LTY[3],lwd=LWD[3]); par(new=T)
points(S.grid,RESULT$F1[RESULT$Method==4],type='l',col=COL2[4],lty=LTY[4],lwd=LWD[4])


plot(S.grid,RESULT$F2[RESULT$Method==1],type='l',ylim=c(0.6,0.85),
     col=COL2[1],lty=LTY[1],lwd=LWD[1],
     xlab="",ylab="",main="")
title(xlab="",ylab="",line=2)
title(main="Inverse F1",font.main=1)
par(new=T)
points(S.grid,RESULT$F2[RESULT$Method==2],type='l',col=COL2[2],lty=LTY[2],lwd=LWD[2]); par(new=T)
points(S.grid,RESULT$F2[RESULT$Method==3],type='l',col=COL2[3],lty=LTY[3],lwd=LWD[3]); par(new=T)
points(S.grid,RESULT$F2[RESULT$Method==4],type='l',col=COL2[4],lty=LTY[4],lwd=LWD[4])

plot(S.grid,RESULT$MCC[RESULT$Method==1],type='l',ylim=c(0.425,0.55),
     col=COL2[1],lty=LTY[1],lwd=LWD[1],
     xlab="",ylab="",main="")
title(xlab="",ylab="",line=2)
title(main="MCC",font.main=1)
par(new=T)
points(S.grid,RESULT$MCC[RESULT$Method==2],type='l',col=COL2[2],lty=LTY[2],lwd=LWD[2]); par(new=T)
points(S.grid,RESULT$MCC[RESULT$Method==3],type='l',col=COL2[3],lty=LTY[3],lwd=LWD[3]); par(new=T)
points(S.grid,RESULT$MCC[RESULT$Method==4],type='l',col=COL2[4],lty=LTY[4],lwd=LWD[4])

par(mar=c(4,0,2,0))
plot.new()
segments(0,0.90,0.3,0.90,col=COL2[1],lty=LTY[1],lwd=LWD[1]+0.25)
segments(0,0.75,0.3,0.75,col=COL2[2],lty=LTY[2],lwd=LWD[2]+0.25)
segments(0,0.60,0.3,0.60,col=COL2[3],lty=LTY[3],lwd=LWD[3]+0.25)
segments(0,0.45,0.3,0.45,col=COL2[4],lty=LTY[4],lwd=LWD[4]+0.25)

text(rep(0.35,4),c(90,75,60,45)/100,
     c("Direct","Indirect-Lasso","Indirect-RF","True"),
     pos=4)

mtext("T",side=1,line=-1,outer=TRUE)

dev.off()
