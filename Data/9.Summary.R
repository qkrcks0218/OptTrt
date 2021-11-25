############################################
# By running this R-file, the simulation results are summarized and the plots are drawn.
############################################

library(ggplot2)
library(dplyr)
library(maps)
library(GADMTools)
library(ptinpoly)
library(ggpubr)
library(RColorBrewer)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)

source("../OptTrt_Sim_Function.R")

# Classification

Test.Data <- read.csv("Sen_Test_18.csv")
N.test <- dim(Test.Data)[1]
M.test <- Test.Data$C_size

DATA <- read.csv("RULE_Merged.csv")


Value <- RRR <- list()

for(ss in 1:10){
  
  S.quart <- S <- c((64:73)/100)[ss]
  
  RULE.DR.Test.Median <- DATA[,which(colnames(DATA)==sprintf("RULE.DR.Test.Median.%0.2d",S*100))]
  RULE.Lasso.Test.Median <- DATA[,which(colnames(DATA)==sprintf("RULE.Lasso.Test.Median.%0.2d",S*100))]
  RULE.RF.Test.Median <- DATA[,which(colnames(DATA)==sprintf("RULE.RF.Test.Median.%0.2d",S*100))]
  
  T.DR <- as.numeric(RULE.DR.Test.Median       < Test.Data$A)
  T.Lasso <- as.numeric(RULE.Lasso.Test.Median < Test.Data$A)
  T.RF <- as.numeric(RULE.RF.Test.Median       < Test.Data$A)
  
  OUTCOME <- Test.Data$Y
  
  RR <- as.matrix(rbind(c( Evaluation(T.DR, as.numeric(S    < OUTCOME)),sum(RULE.DR.Test.Median*M.test)/sum(M.test) ),
                        c( Evaluation(T.Lasso, as.numeric(S < OUTCOME)),sum(RULE.Lasso.Test.Median*M.test)/sum(M.test) ),
                        c( Evaluation(T.RF, as.numeric(S    < OUTCOME)),sum(RULE.RF.Test.Median*M.test)/sum(M.test) )))
  Value[[ss]] <- matrix(as.numeric(RR),3,11)
  
}



RESULT <- data.frame(cbind(rep(1:3,10),rbind(Value[[1]],Value[[2]],Value[[3]],
                                             Value[[4]],Value[[5]],Value[[6]],
                                             Value[[7]],Value[[8]],Value[[9]],
                                             Value[[10]])))
colnames(RESULT) <- c("Method","","","","","Accuracy","","","F1","F2","MCC","Prop")

COL2 <- c("red",rgb(0,206/255,206/255,alpha=0.8),	rgb(255/255,165/255,0,alpha=0.8),rgb(0,0,0,0.5))
LTY <- c(1,5,5,1)
LWD <- c(1.75,1.75,1.75,1.75)

S.grid <- (64:73)/100

png("Plot/CPMeasure.png",width=8,height=2,unit="in",res=500)

layout(matrix(c(1,2,3,4,5),1,5,byrow=T),widths = c(2,2,2,2,1.5))
par(mar=c(4,2.5,2,0.5))

plot(S.grid,RESULT$Accuracy[RESULT$Method==1],type='l',ylim=c(0.5,0.65),
     col=COL2[1],lty=LTY[1],lwd=LWD[1],
     xlab="",ylab="",main="",xaxt='n')
title(xlab="",ylab="",line=2)
title(main="Accuracy",font.main=1)
axis(1,at=S.grid)
par(new=T)
points(S.grid,RESULT$Accuracy[RESULT$Method==2],type='l',col=COL2[2],lty=LTY[2],lwd=LWD[2]); par(new=T)
points(S.grid,RESULT$Accuracy[RESULT$Method==3],type='l',col=COL2[3],lty=LTY[3],lwd=LWD[3])


plot(S.grid,RESULT$F1[RESULT$Method==1],type='l',ylim=c(0.35,0.75),
     col=COL2[1],lty=LTY[1],lwd=LWD[1],
     xlab="",ylab="",main="",xaxt='n')
title(xlab="",ylab="",line=2)
axis(1,at=S.grid)
title(main="F1",font.main=1)
par(new=T)
points(S.grid,RESULT$F1[RESULT$Method==2],type='l',col=COL2[2],lty=LTY[2],lwd=LWD[2]); par(new=T)
points(S.grid,RESULT$F1[RESULT$Method==3],type='l',col=COL2[3],lty=LTY[3],lwd=LWD[3])

plot(S.grid,RESULT$F2[RESULT$Method==1],type='l',ylim=c(0.35,0.7),
     col=COL2[1],lty=LTY[1],lwd=LWD[1],
     xlab="",ylab="",main="",xaxt='n')
title(xlab="",ylab="",line=2)
axis(1,at=S.grid)
title(main="Inverse F1",font.main=1)
par(new=T)
points(S.grid,RESULT$F2[RESULT$Method==2],type='l',col=COL2[2],lty=LTY[2],lwd=LWD[2]); par(new=T)
points(S.grid,RESULT$F2[RESULT$Method==3],type='l',col=COL2[3],lty=LTY[3],lwd=LWD[3])

plot(S.grid,RESULT$MCC[RESULT$Method==1],type='l',ylim=c(0.0,0.25),
     col=COL2[1],lty=LTY[1],lwd=LWD[1],
     xlab="",ylab="",main="",xaxt='n')
axis(1,at=S.grid)
title(xlab="",ylab="",line=2)
title(main="MCC",font.main=1)
par(new=T)
points(S.grid,RESULT$MCC[RESULT$Method==2],type='l',col=COL2[2],lty=LTY[2],lwd=LWD[2]); par(new=T)
points(S.grid,RESULT$MCC[RESULT$Method==3],type='l',col=COL2[3],lty=LTY[3],lwd=LWD[3])


par(mar=c(4,0,2,0))
plot.new()
segments(0,0.90,0.3,0.90,col=COL2[1],lty=LTY[1],lwd=LWD[1])
segments(0,0.75,0.3,0.75,col=COL2[2],lty=LTY[2],lwd=LWD[2])
segments(0,0.60,0.3,0.60,col=COL2[3],lty=LTY[3],lwd=LWD[3])

text(rep(0.35,3),c(90,75,60)/100,
     c("Direct","Indirect-Lasso","Indirect-RF"),
     pos=4)

mtext("T",side=1,line=-1,outer=TRUE)

dev.off()


# Amount of Trt, Heatplots

R1 <- DATA$R1

SEN2 <- gadm_sf_loadCountries(c("SEN"), level=2, basefile = "./gadm36_SEN_shp")

POS <- list()
POS[[ 1]] <- which( SEN2$sf$NAME_2==sort(unique((SEN2$sf$NAME_2)))[14] )
POS[[ 2]] <- which( SEN2$sf$NAME_2==sort(unique((SEN2$sf$NAME_2)))[15] )
POS[[ 3]] <- which( SEN2$sf$NAME_2==sort(unique((SEN2$sf$NAME_2)))[19] )
POS[[ 4]] <- which( SEN2$sf$NAME_2==sort(unique((SEN2$sf$NAME_2)))[20] )
POS[[ 5]] <- which( SEN2$sf$NAME_2==sort(unique((SEN2$sf$NAME_2)))[24] )
POS[[ 6]] <- which( SEN2$sf$NAME_2==sort(unique((SEN2$sf$NAME_2)))[26] )
POS[[ 7]] <- which( SEN2$sf$NAME_2==sort(unique((SEN2$sf$NAME_2)))[28] )
POS[[ 8]] <- which( SEN2$sf$NAME_2==sort(unique((SEN2$sf$NAME_2)))[30] )
POS[[ 9]] <- which( SEN2$sf$NAME_2==sort(unique((SEN2$sf$NAME_2)))[35] )
POS[[10]] <- which( SEN2$sf$NAME_2==sort(unique((SEN2$sf$NAME_2)))[38] )
POS[[11]] <- which( SEN2$sf$NAME_2==sort(unique((SEN2$sf$NAME_2)))[40] )
POS[[12]] <- which( SEN2$sf$NAME_2==sort(unique((SEN2$sf$NAME_2)))[42] )
POS[[13]] <- which( SEN2$sf$NAME_2==sort(unique((SEN2$sf$NAME_2)))[44] )

SEN2$sf$NAME_2[ POS[[ 1]] ] <- "Guediawaye"
SEN2$sf$NAME_2[ POS[[ 2]] ] <- "Guinguineo"
SEN2$sf$NAME_2[ POS[[ 3]] ] <- "Kebemer"
SEN2$sf$NAME_2[ POS[[ 4]] ] <- "Kedougou"
SEN2$sf$NAME_2[ POS[[ 5]] ] <- "Linguere"
SEN2$sf$NAME_2[ POS[[ 6]] ] <- "Maleme Hodar"
SEN2$sf$NAME_2[ POS[[ 7]] ] <- "Mbacke"
SEN2$sf$NAME_2[ POS[[ 8]] ] <- "Medina Yoro Foula"
SEN2$sf$NAME_2[ POS[[ 9]] ] <- "Ranerou Ferlo"
SEN2$sf$NAME_2[ POS[[10]] ] <- "Salemata"
SEN2$sf$NAME_2[ POS[[11]] ] <- "Sedhiou"
SEN2$sf$NAME_2[ POS[[12]] ] <- "Thies"
SEN2$sf$NAME_2[ POS[[13]] ] <- "Velingara"

NNM <- SEN2$sf$NAME_2
SEN2$sf$District15 <- NNM

DATA.M <- aggregate(cbind(DATA$RULE.Lasso.Test.Median.64,
                          DATA$RULE.RF.Test.Median.64,
                          DATA$RULE.DR.Test.Median.64,
                          
                          DATA$RULE.Lasso.Test.Median.67,
                          DATA$RULE.RF.Test.Median.67,
                          DATA$RULE.DR.Test.Median.67,
                          
                          DATA$RULE.Lasso.Test.Median.70,
                          DATA$RULE.RF.Test.Median.70,
                          DATA$RULE.DR.Test.Median.70,
                          
                          DATA$RULE.Lasso.Test.Median.73,
                          DATA$RULE.RF.Test.Median.73,
                          DATA$RULE.DR.Test.Median.73,
                          
                          1,
                          (1-DATA$C_rt_C)/DATA$C_size,
                          DATA$X_HHsize,
                          DATA$X_HHchildsize,
                          -DATA$X_age,
                          -DATA$X_head_edu,
                          DATA$X_Nojob,
                          -DATA$X_resp_age,
                          -DATA$X_resp_edu)*DATA$C_size~DATA$R1,FUN="sum")

COUNT <- aggregate(rep(1,dim(DATA)[1])~DATA$R1,FUN="sum")
DATA.M <- cbind(DATA.M,COUNT[,2])

colN <- c("NAME_2",
          "Lasso.64","RF.64","DR.64",
          "Lasso.67","RF.67","DR.67",
          "Lasso.70","RF.70","DR.70",
          "Lasso.73","RF.73","DR.73",
          "C_size","C_rt_C","X_HHsize",
          "X_HHchildsize","X_age","X_head_edu","X_Nojob",
          "X_resp_age","X_resp_edu","COUNT")
colnames(DATA.M) <- colN
DATA.M <- cbind(DATA.M[,1],DATA.M[,2:13]/DATA.M$C_size,DATA.M[,14:15]/DATA.M[,23],DATA.M[,16:22]/DATA.M$C_size,DATA.M[,23])
colnames(DATA.M) <- colN

mycol <- brewer.pal(3,"Blues")[c(1,3)]
SEN2$sf$Lasso.64 <- SEN2$sf$RF.64 <- SEN2$sf$DR.64 <- 0
SEN2$sf$Lasso.67 <- SEN2$sf$RF.67 <- SEN2$sf$DR.67 <- 0
SEN2$sf$Lasso.70 <- SEN2$sf$RF.70 <- SEN2$sf$DR.70 <- 0
SEN2$sf$Lasso.73 <- SEN2$sf$RF.73 <- SEN2$sf$DR.73 <- 0


SEN2$sf$C_size <- 0
SEN2$sf$C_rt_C <- 0
SEN2$sf$X_HHsize <- 0
SEN2$sf$X_HHchildsize <- 0
SEN2$sf$X_age <- 0
SEN2$sf$X_head_edu <- 0
SEN2$sf$X_Nojob <- 0
SEN2$sf$X_resp_age <- 0
SEN2$sf$X_resp_edu <- 0
SEN2$sf$COUNT <- 0
SEN2$sf$X <- 0
SEN2$sf$Y <- 0


COORDINATE <- cbind(SEN2$sf$NAME_2,st_coordinates(st_centroid(SEN2$sf)))
for(ii in 1:length(SEN2$sf$NAME_2)){
  
  SEN2$sf$Lasso.64         [ii] <- ( (DATA.M$Lasso.64          ))[DATA.M[,1]==SEN2$sf$NAME_2[ii]]
  SEN2$sf$RF.64            [ii] <- ( (DATA.M$RF.64             ))[DATA.M[,1]==SEN2$sf$NAME_2[ii]]
  SEN2$sf$DR.64            [ii] <- ( (DATA.M$DR.64             ))[DATA.M[,1]==SEN2$sf$NAME_2[ii]]
  
  SEN2$sf$Lasso.67         [ii] <- ( (DATA.M$Lasso.67          ))[DATA.M[,1]==SEN2$sf$NAME_2[ii]]
  SEN2$sf$RF.67            [ii] <- ( (DATA.M$RF.67             ))[DATA.M[,1]==SEN2$sf$NAME_2[ii]]
  SEN2$sf$DR.67            [ii] <- ( (DATA.M$DR.67             ))[DATA.M[,1]==SEN2$sf$NAME_2[ii]]
  
  SEN2$sf$Lasso.70         [ii] <- ( (DATA.M$Lasso.70          ))[DATA.M[,1]==SEN2$sf$NAME_2[ii]]
  SEN2$sf$RF.70            [ii] <- ( (DATA.M$RF.70             ))[DATA.M[,1]==SEN2$sf$NAME_2[ii]]
  SEN2$sf$DR.70            [ii] <- ( (DATA.M$DR.70             ))[DATA.M[,1]==SEN2$sf$NAME_2[ii]]
  
  SEN2$sf$Lasso.73         [ii] <- ( (DATA.M$Lasso.73          ))[DATA.M[,1]==SEN2$sf$NAME_2[ii]]
  SEN2$sf$RF.73            [ii] <- ( (DATA.M$RF.73             ))[DATA.M[,1]==SEN2$sf$NAME_2[ii]]
  SEN2$sf$DR.73            [ii] <- ( (DATA.M$DR.73             ))[DATA.M[,1]==SEN2$sf$NAME_2[ii]]
  
  SEN2$sf$C_size          [ii] <- ( (DATA.M$C_size          ))[DATA.M[,1]==SEN2$sf$NAME_2[ii]]
  SEN2$sf$C_rt_C          [ii] <- ( (DATA.M$C_rt_C          ))[DATA.M[,1]==SEN2$sf$NAME_2[ii]]
  SEN2$sf$X_HHsize        [ii] <- ( (DATA.M$X_HHsize        ))[DATA.M[,1]==SEN2$sf$NAME_2[ii]]
  SEN2$sf$X_HHchildsize   [ii] <- ( (DATA.M$X_HHchildsize   ))[DATA.M[,1]==SEN2$sf$NAME_2[ii]]
  SEN2$sf$X_age           [ii] <- ( (DATA.M$X_age           ))[DATA.M[,1]==SEN2$sf$NAME_2[ii]]
  SEN2$sf$X_head_edu      [ii] <- ( (DATA.M$X_head_edu      ))[DATA.M[,1]==SEN2$sf$NAME_2[ii]]
  SEN2$sf$X_Nojob         [ii] <- ( (DATA.M$X_Nojob         ))[DATA.M[,1]==SEN2$sf$NAME_2[ii]]
  SEN2$sf$X_resp_age      [ii] <- ( (DATA.M$X_resp_age      ))[DATA.M[,1]==SEN2$sf$NAME_2[ii]]
  SEN2$sf$X_resp_edu      [ii] <- ( (DATA.M$X_resp_edu      ))[DATA.M[,1]==SEN2$sf$NAME_2[ii]]
  SEN2$sf$COUNT           [ii] <- ( (DATA.M$COUNT           ))[DATA.M[,1]==SEN2$sf$NAME_2[ii]]
  
  SEN2$sf$X[ii] <- as.numeric(COORDINATE[,2][COORDINATE[,1]==SEN2$sf$NAME_2[ii]])
  SEN2$sf$Y[ii] <- as.numeric(COORDINATE[,3][COORDINATE[,1]==SEN2$sf$NAME_2[ii]])
}


WINS <- function(x,a,b){
  if(is.null(a)){
    XX <- as.numeric(x<=b)*x + as.numeric(x>b)*b
  } else if(is.null(b)){
    XX <- as.numeric(x>=a)*x + as.numeric(x<a)*a
  } else {
    XX <- as.numeric(x<=a)*a + as.numeric(a<x & x<b)*x + as.numeric(x>=b)*b
  }
  return(XX)
}


SEN2$sf$DR.64 <- WINS(SEN2$sf$DR.64,0.2,0.8)
SEN2$sf$Lasso.64 <- WINS(SEN2$sf$Lasso.64,0.2,0.8)
SEN2$sf$RF.64 <- WINS(SEN2$sf$RF.64,0.2,0.8)

SEN2$sf$DR.67 <- WINS(SEN2$sf$DR.67,0.2,0.8)
SEN2$sf$Lasso.67 <- WINS(SEN2$sf$Lasso.67,0.2,0.8)
SEN2$sf$RF.67 <- WINS(SEN2$sf$RF.67,0.2,0.8)

SEN2$sf$DR.70 <- WINS(SEN2$sf$DR.70,0.2,0.8)
SEN2$sf$Lasso.70 <- WINS(SEN2$sf$Lasso.70,0.2,0.8)
SEN2$sf$RF.70 <- WINS(SEN2$sf$RF.70,0.2,0.8)

SEN2$sf$DR.73 <- WINS(SEN2$sf$DR.73,0.2,0.8)
SEN2$sf$Lasso.73 <- WINS(SEN2$sf$Lasso.73,0.2,0.8)
SEN2$sf$RF.73 <- WINS(SEN2$sf$RF.73,0.2,0.8)



G1 <- G2 <- G3 <- list()
LB <- c(expression(phantom(x)<="0.2"),
        "0.3","0.4","0.5","0.6","0.7",
        expression(phantom(x)>="0.8"))

G1[[1]] <- ggplot() + geom_sf(data=SEN2$sf,aes(fill=(DR.64)),lwd=0.1)  + 
  scale_fill_distiller(palette = "Blues", direction=1, limits=c(0.2,0.8),
                       labels=LB) + theme_void() +
  labs(title = "", fill = "") + theme(legend.position = "none")
G2[[1]] <- ggplot() + geom_sf(data=SEN2$sf,aes(fill=(Lasso.64)),lwd=0.1)  + 
  scale_fill_distiller(palette = "Blues", direction=1, limits=c(0.2,0.8),
                       labels=LB) + theme_void() +
  labs(title = "", fill = "") + theme(legend.position = "none")
G3[[1]] <- ggplot() + geom_sf(data=SEN2$sf,aes(fill=(RF.64)),lwd=0.1)  + 
  scale_fill_distiller(palette = "Blues", direction=1, limits=c(0.2,0.8),
                       labels=LB) + theme_void() +
  labs(title = "", fill = "") + theme(legend.position = "none")

G1[[2]] <- ggplot() + geom_sf(data=SEN2$sf,aes(fill=(DR.67)),lwd=0.1)  + 
  scale_fill_distiller(palette = "Blues", direction=1, limits=c(0.2,0.8),
                       labels=LB) + theme_void() +
  labs(title = "", fill = "")  + theme(legend.position = "none") 
G2[[2]] <- ggplot() + geom_sf(data=SEN2$sf,aes(fill=(Lasso.67)),lwd=0.1)  + 
  scale_fill_distiller(palette = "Blues", direction=1, limits=c(0.2,0.8),
                       labels=LB) + theme_void() +
  labs(title = "", fill = "") + theme(legend.position = "none") 
G3[[2]] <- ggplot() + geom_sf(data=SEN2$sf,aes(fill=(RF.67)),lwd=0.1)  + 
  scale_fill_distiller(palette = "Blues", direction=1, limits=c(0.2,0.8),
                       labels=LB) + theme_void() +
  labs(title = "", fill = "") + theme(legend.position = "none") 

G1[[3]] <- ggplot() + geom_sf(data=SEN2$sf,aes(fill=(DR.70)),lwd=0.1)  + 
  scale_fill_distiller(palette = "Blues", direction=1, limits=c(0.2,0.8),
                       labels=LB) + theme_void() +
  labs(title = "", fill = "")  + theme(legend.position = "none") 
G2[[3]] <- ggplot() + geom_sf(data=SEN2$sf,aes(fill=(Lasso.70)),lwd=0.1)  + 
  scale_fill_distiller(palette = "Blues", direction=1, limits=c(0.2,0.8),
                       labels=LB) + theme_void() +
  labs(title = "", fill = "") + theme(legend.position = "none") 
G3[[3]] <- ggplot() + geom_sf(data=SEN2$sf,aes(fill=(RF.70)),lwd=0.1)  + 
  scale_fill_distiller(palette = "Blues", direction=1, limits=c(0.2,0.8),
                       labels=LB) + theme_void() +
  labs(title = "", fill = "") + theme(legend.position = "none") 

G1[[4]] <- ggplot() + geom_sf(data=SEN2$sf,aes(fill=(DR.73)),lwd=0.1)  + 
  scale_fill_distiller(palette = "Blues", direction=1, limits=c(0.2,0.8),
                       labels=LB) + theme_void() +
  labs(title = "", fill = "")  + theme(legend.position = "none") 
G2[[4]] <- ggplot() + geom_sf(data=SEN2$sf,aes(fill=(Lasso.73)),lwd=0.1)  + 
  scale_fill_distiller(palette = "Blues", direction=1, limits=c(0.2,0.8),
                       labels=LB) + theme_void() +
  labs(title = "", fill = "") + theme(legend.position = "none") 
G3[[4]] <- ggplot() + geom_sf(data=SEN2$sf,aes(fill=(RF.73)),lwd=0.1)  + 
  scale_fill_distiller(palette = "Blues", direction=1, limits=c(0.2,0.8),
                       labels=LB) + theme_void() +
  labs(title = "", fill = "") + theme(legend.position = "none") 

G4 <- ggplot() + geom_sf(data=SEN2$sf,aes(fill=(DR.70)),lwd=0.1)  + 
  scale_fill_distiller(palette = "Blues", direction=1, limits=c(0.2,0.8),
                       labels=LB) + theme_void() +
  labs(title = "Direct (T=0.70)", fill = "") 

LEGEND <- get_legend(G4)

gg51 <- LEGEND



png("Plot/TotalAmount.png",
    height=2.5,width=10,res=500,unit="in")

COL2 <- c("white","white","red",rgb(0,206/255,206/255,alpha=0.8),	rgb(255/255,165/255,0,alpha=0.8),rgb(0,0,0,0.5))
LTY <- c(0,0,1,5,5,1)
LWD <- c(0,0,1.75,1.75,1.75,1.75)*1.5
Yaxis <-    ggplotGrob(ggplot() + annotate("text", x = 0.5, y = 0.5, size=6,
                                           label ="Proportion of Households \n That Require \n WASH Facilities",srt=90) + theme_void())
S.grid <- 64:73
Trt.Prop2 <- matrix(0,10,4)

for(gg in 1:10){
  
  CN <- c( which( colnames(DATA)==paste("RULE.DR.Test.Median.",   S.grid[gg],sep="") ),
           which( colnames(DATA)==paste("RULE.Lasso.Test.Median.",S.grid[gg],sep="") ),
           which( colnames(DATA)==paste("RULE.RF.Test.Median.",   S.grid[gg],sep="") ) )
  Trt.Prop2[gg,] <- c(S.grid[gg],apply(DATA[,CN]*DATA$C_size,2,sum)/sum(DATA$C_size))
  
}

layout(matrix(1:3,1,3),widths=c(1,6,2))

par(mar=c(3,0,1,0))

plot.new()
text(0.5,0.5,"Average of \nthe Estimated OMARs",
     srt=90,cex=1.75)

par(mar=c(3,1,1,2))

plot(S.grid,Trt.Prop2[,2],type='l',col=COL2[3],lty=LTY[3],lwd=LWD[3],ylim=c(0.2,0.8),
     xaxt="n",xlab="",ylab="",yaxt="n")
axis(1,at=64:73,label=c("T=0.64","","","T=0.67","","","T=0.70","","","T=0.73"),
     cex.axis=1.5)
axis(2,at=2:8/10,label=2:8/10,cex.axis=1.5)
lines(S.grid,Trt.Prop2[,3],col=COL2[4],lty=LTY[4],lwd=LWD[4])
lines(S.grid,Trt.Prop2[,4],col=COL2[5],lty=LTY[5],lwd=LWD[5])
polygon(c(69.5,0,0,69.5,69.5),c(2,2,-2,-2,2),col=rgb(0,0,0,0.1),border=NA)
text(64,0.75,"Direct rule suggests", pos=4,cex=1.5)
text(64,0.68,"more WASH facilities",pos=4,cex=1.5)
text(70,0.32,"Direct rule suggests", pos=4,cex=1.5)
text(70,0.25,"less WASH facilities",pos=4,cex=1.5)

par(mar=c(3,0,1,0))
plot.new()
segments(0,0.90,0.3,0.90,col=COL2[3],lty=LTY[3],lwd=LWD[3])
segments(0,0.75,0.3,0.75,col=COL2[4],lty=LTY[4],lwd=LWD[4])
segments(0,0.60,0.3,0.60,col=COL2[5],lty=LTY[5],lwd=LWD[5])

text(rep(0.35,3),c(90,75,60)/100,
     c("Direct","Indirect-Lasso","Indirect-RF"),pos=4,cex=1.5)

dev.off()



gDR <-    ggplotGrob(ggplot() + annotate("text", x = 0.5, y = 0.5, size=6, label ="Direct",srt=90) + theme_void())
gLasso <- ggplotGrob(ggplot() + annotate("text", x = 0.5, y = 0.5, size=6, label ="Indirect-Lasso",srt=90) + theme_void() )
gRF <-    ggplotGrob(ggplot() + annotate("text", x = 0.5, y = 0.5, size=6, label ="Indirect-RF",srt=90) + theme_void() )

gT64 <-    ggplotGrob(ggplot() + annotate("text", x = 0.5, y = 0.5, size=6, label ="T=0.64") + theme_void())
gT67 <-    ggplotGrob(ggplot() + annotate("text", x = 0.5, y = 0.5, size=6, label ="T=0.67") + theme_void() )
gT70 <-    ggplotGrob(ggplot() + annotate("text", x = 0.5, y = 0.5, size=6, label ="T=0.70") + theme_void() )
gT73 <-    ggplotGrob(ggplot() + annotate("text", x = 0.5, y = 0.5, size=6, label ="T=0.73") + theme_void() )

gg11 <- ggplotGrob(G1[[1]]); gg12 <- ggplotGrob(G2[[1]]); gg13 <- ggplotGrob(G3[[1]])
gg21 <- ggplotGrob(G1[[2]]); gg22 <- ggplotGrob(G2[[2]]); gg23 <- ggplotGrob(G3[[2]])
gg31 <- ggplotGrob(G1[[3]]); gg32 <- ggplotGrob(G2[[3]]); gg33 <- ggplotGrob(G3[[3]])
gg41 <- ggplotGrob(G1[[4]]); gg42 <- ggplotGrob(G2[[4]]); gg43 <- ggplotGrob(G3[[4]])




png("Plot/Grad1.png",
    height=6,width=10,res=500,unit="in")

lay <- rbind(c(24,1:4,23), cbind(matrix(5:19,3,5),c(22,20,21)))

grid.arrange(
  gT64,gT67,gT70,gT73,
  gDR,gLasso,gRF,
  gg11,gg12,gg13,
  gg21,gg22,gg23,
  gg31,gg32,gg33,
  gg41,gg42,gg43,gg51,layout_matrix=lay,widths=c(1,3,3,3,3,1),heights=c(1,5,5,5))

dev.off()



DDD <- matrix(0,10,3)

for(gg in 64:73){
  
  CN <- c( which( colnames(DATA)==paste("RULE.DR.Test.Median.",gg,sep="") ),
           which( colnames(DATA)==paste("RULE.Lasso.Test.Median.",gg,sep="") ),
           which( colnames(DATA)==paste("RULE.RF.Test.Median.",gg,sep="") ) )
  
  DDD[gg-63,] <- apply(cbind( DATA[,CN]*DATA$C_size ),2,mean)
  
  
}

add.coor.X <- add.coor.Y <- rep(0,dim(SEN2$sf)[1])
add.coor.X[c(which(SEN2$sf$District15=="Dakar"))] <- -0.1
add.coor.Y[c(which(SEN2$sf$District15=="Dakar"),
             which(SEN2$sf$District15=="Pikine"))] <- - 0.1
add.coor.Y[c(which(SEN2$sf$District15=="Guediawaye"))] <- 0.1
add.coor.Y[19] <- add.coor.Y[19] -0.05

POSA <- c(29,19)

G4 <- ggplot() + geom_sf(data=SEN2$sf,aes(fill=(DR.70)),lwd=0.1)  + 
  scale_fill_distiller(palette = "Blues", direction=1, limits=c(0.2,0.8),
                       labels=LB) + theme_void() +
  labs(title = "Average OMPR Under T=0.70", fill = "") + 
  geom_text(aes( SEN2$sf$X[POSA]+add.coor.X[POSA], 
                 SEN2$sf$Y[POSA]+add.coor.Y[POSA], 
                 label=rep("*",2)),col="black")
DR70.L <- get_legend(G4)

DR70.noL <- ggplot() + geom_sf(data=SEN2$sf,aes(fill=(DR.70)),lwd=0.1)  +
  scale_fill_distiller(palette = "Blues", direction=1, limits=c(0.2,0.8),
                       labels=LB) + theme_void() +
  labs(title = "Average OMPR Under T=0.70", fill = "") + 
  geom_text(aes( SEN2$sf$X[POSA]+add.coor.X[POSA], 
                 SEN2$sf$Y[POSA]+add.coor.Y[POSA], 
                 label=rep("*",2)),col="black") + theme(legend.position = "none")


RURAL <- ggplot() + geom_sf(data=SEN2$sf,aes(fill=WINS(C_rt_C,0.2,1.0)),lwd=0.1)  +
  scale_fill_distiller(palette = "Blues", direction=1,
                       labels=c(expression(phantom(x)<="0.2"),
                                "0.4","0.6","0.8","1.0")) + theme_void() +
  labs(title = "Proportion of Rural Areas", fill = "") + 
  geom_text(aes( SEN2$sf$X[POSA]+add.coor.X[POSA], 
                 SEN2$sf$Y[POSA]+add.coor.Y[POSA], 
                 label=rep("*",2)),col="white")
RURAL.L <- get_legend(RURAL)

RURAL.noL <-  ggplot() + geom_sf(data=SEN2$sf,aes(fill=WINS(C_rt_C,0.2,1.0)),lwd=0.1)  +
  scale_fill_distiller(palette = "Blues", direction=1,
                       labels=c(expression(phantom(x)<="0.2"),
                                "0.4","0.6","0.8","1.0")) + theme_void() +
  labs(title = "Proportion of Rural Areas", fill = "") + 
  geom_text(aes( SEN2$sf$X[POSA]+add.coor.X[POSA], 
                 SEN2$sf$Y[POSA]+add.coor.Y[POSA], 
                 label=rep("*",2)),col="white") + theme(legend.position = "none")

HS <- ggplot() + geom_sf(data=SEN2$sf,aes(fill=WINS(X_HHsize,8,14.4)),lwd=0.1)  +
  scale_fill_distiller(palette = "Blues", direction=1) + theme_void() +
  labs(title = "Average Household Size", fill = "") + 
  geom_text(aes( SEN2$sf$X[POSA]+add.coor.X[POSA], 
                 SEN2$sf$Y[POSA]+add.coor.Y[POSA], 
                 label=rep("*",2)),col="black")

HS.L <- get_legend(HS)

HS.noL <- ggplot() + geom_sf(data=SEN2$sf,aes(fill=WINS(X_HHsize,8,14.4)),lwd=0.1)  +
  scale_fill_distiller(palette = "Blues", direction=1) + theme_void() +
  labs(title = "Average Household Size", fill = "") + 
  geom_text(aes( SEN2$sf$X[POSA]+add.coor.X[POSA], 
                 SEN2$sf$Y[POSA]+add.coor.Y[POSA], 
                 label=rep("*",2)),col="black") + theme(legend.position = "none")




png("Plot/Grad2.png",
    height=2.5,width=10,res=500,unit="in")


grid.arrange(DR70.noL,DR70.L,
             HS.noL,HS.L,
             RURAL.noL,RURAL.L,
             layout_matrix=matrix(1:6,1,6),
             widths=c(3,1,3,1,3,1))

dev.off()




