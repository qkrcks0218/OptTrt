############################################
# By running this R-file, the assumptions are assessed as in Section 5.2.
############################################

library(geosphere)
library(KernSmooth)

Data.Cluster <- read.csv("Sen_Training_1417.csv")
Test.Data <- read.csv("Sen_Test_18.csv")

N <- dim(Data.Cluster)[1]
M <- Data.Cluster$C_size

NF.result <- list()
for(ii in 1:100){
    NF.result[[ii]] <- read.csv(sprintf("NF/NF_Train_B%0.4d.csv",ii))
    print(ii)
}

CL.name <- colnames(NF.result[[1]])

NF <- array(unlist(NF.result), c(dim(NF.result[[1]])[1], dim(NF.result[[1]])[2], 100))
NF.median <- apply(NF,c(1:2),median) # 20-42
NF.mean <- apply(NF,c(1:2),mean)


############# Assessment of Monotonicity Assumption

Index <- matrix(0,dim(NF.result[[1]])[1],max(M))
for(ii in 1:dim(NF.result[[1]])[1]){
    Index[ii,1:M[ii]] <- 1
}

Diff <- matrix(0,dim(NF.result[[1]])[1],max(M))
for(jj in 1:max(M)){
    Diff[,jj] <- NF.median[,which(CL.name=="OR0.Est")+jj]-NF.median[,which(CL.name=="OR0.Est")-1+jj]
}


sum(apply(Diff*Index<0,1,sum)) ; sum(apply(Diff*Index<0,1,sum))/sum(M)
sum( apply(abs(Diff)*Index*(Diff*Index<0),1,sum) ) / sum( apply(abs(Diff)*Index,1,sum) )


Diff <- matrix(0,dim(NF.result[[1]])[1],max(M))
for(jj in 1:max(M)){
    Diff[,jj] <- NF.mean[,which(CL.name=="OR0.Est")+jj]-NF.mean[,which(CL.name=="OR0.Est")-1+jj]
}
sum(apply(Diff*Index<0,1,sum)) ; sum(apply(Diff*Index<0,1,sum))/sum(M)
sum( apply(abs(Diff)*Index*(Diff*Index<0),1,sum) ) / sum( apply(abs(Diff)*Index,1,sum) ) #3%full, #2% 5% # 0.3% 10%

NF.median <- data.frame(NF.median)
colnames(NF.median) <- CL.name

NF.median$OR.Max.Est <- apply(NF.median[,17:39],1,max)


############# Assessment of Covariate Balance

PP <- dim(Data.Cluster)[2]-2

A.group <- list()
A.q <- c(0,quantile( NF.median[,2], (1:2)/3 ),1)
A.group[[1]] <- which( NF.median[,2] <= A.q[2])
A.group[[2]] <- which( A.q[2] < NF.median[,2] & NF.median[,2] <= A.q[3])
A.group[[3]] <- which( A.q[3] < NF.median[,2])

PS.loc <- c(which(CL.name=="PS.OL.Bin1.Est"), which(CL.name=="PS.OL.Bin2.Est"),
            which(CL.name=="PS.OL.Bin3.Est"))

TTT <- SI <- list()
LL <- 10

for(tt in 1:100){
    
    NF.Temp <- NF.result[[tt]]
    
    Tstat.NA <- Tstat.A <- matrix(0,PP,3)
    Size <- matrix(0,PP,6)
    
    for(kk in 1:3){
        PS <- NF.Temp[,PS.loc[kk]]
        PS.q <- c(0,quantile(PS[A.group[[kk]]],c(1:(LL-1))/LL),1)
        
        w <- SE <- Est <- rep(0,LL)
        
        for(pp in 1:PP){
            
            for(ag in 1:LL){
                GPSgp <- which( PS.q[ag] < PS & PS <= PS.q[ag+1])
                
                p1 <- intersect(GPSgp,A.group[[kk]])
                p2 <- intersect(GPSgp,(1:N)[-A.group[[kk]]])
                Size[pp,kk] <- length(p1) 
                Size[pp,kk+3] <- length(p2)
                
                v1 <- if(length(p1)>1){ var(NF.median[p1,pp+2]) } else { 0 }
                v2 <- if(length(p2)>1){ var(NF.median[p2,pp+2]) } else { 0 }
                
                Est[ag] <- mean(NF.median[p1,pp+2]) - mean(NF.median[p2,pp+2])
                SE[ag] <- sqrt( v1/length(p1) + v2/length(p2) )
                w[ag] <- length(GPSgp)/N
                
            }
            Tstat.A[pp,kk] <- sum( Est*w ) / sqrt( sum( SE^2*w^2 ) )
            
            Tstat.NA[pp,kk] <- t.test( NF.result[[tt]][A.group[[kk]],pp+2], 
                                       NF.result[[tt]][-A.group[[kk]],pp+2] )$statistic
            
        }
    }
    
    
    TTT[[tt]] <- data.frame( round(cbind( Tstat.NA, Tstat.A),3), row.names=CL.name[3:(PP+2)] )
    SI[[tt]] <- Size
    
    
}

png("Plot/Balance.png",height=4,width=12,unit="in",res=500)

RN <- c("Cluster size",
        "Cluster location",
        "HH size",
        "HH's num of children",
        "HH adults' employment",
        "HH head's education",
        "Respondent's education",
        "Respondent's age",
        "HH children's age")

NF.TTT <- array(unlist(TTT), c(dim(TTT[[1]])[1], dim(TTT[[1]])[2], 100))
NF.TTT.median <- apply(NF.TTT,c(1:2),median)
print( data.frame(NF.TTT.median) , 
       row.names=)

NF.SI <- array(unlist(SI), c(dim(SI[[1]])[1], dim(SI[[1]])[2], 100))
NF.SI.median <- apply(NF.SI,c(1,2),median)

layout(matrix(c(1,2),1,2),widths=c(5,1))
par(mar=c(4,3,2,0.5))

plot(1:9-0.3,NF.TTT.median[,1],
     ylim=c(-14,12),xlim=c(0.8,9.2),
     ylab="",xlab="",yaxt="n",xaxt="n",pch=16)
points(1:9-0.2,NF.TTT.median[,2],pch=17)
points(1:9-0.1,NF.TTT.median[,3],pch=18)

points(1:9+0.1,NF.TTT.median[,4],pch=16,col=2)
points(1:9+0.2,NF.TTT.median[,5],pch=17,col=2)
points(1:9+0.3,NF.TTT.median[,6],pch=18,col=2)

axis(1,at=c(1,3,5,7,9),tck=-0.05,labels=RN[c(1,3,5,7,9)])
axis(1,at=c(2,4,6,8),tck=-0.1,labels=rep("",4))
axis(1,at=c(2,4,6,8),labels=RN[c(2,4,6,8)],line=1.5,tick=F)

axis(2,at=c(-10,qnorm(0.025),0,qnorm(0.975),10),
     label=c(-10,-1.96,0,1.96,10),las=2)
abline(h=qnorm(0.975),lty=3)
abline(h=0,lty=1)
abline(h=qnorm(0.025),lty=3)

gg <- 1
polygon(c(gg-10,gg+0.5,gg+0.5,gg-10,gg-10),
        c(20,20,-20,-20,20),col=rgb(0,0,0,0.1),border=NA)

gg <- 9
polygon(c(gg-0.5,gg+10,gg+10,gg-0.5,gg-0.5),
        c(20,20,-20,-20,20),col=rgb(0,0,0,0.1),border=NA)

for(gg in c(3,5,7)){
    polygon(c(gg-0.5,gg+0.5,gg+0.5,gg-0.5,gg-0.5),
            c(20,20,-20,-20,20),col=rgb(0,0,0,0.1),border=NA)
}

par(mar=c(4,0,2,0))

plot.new()
points(0.05,0.95,pch=16,col=1)
points(0.05,0.85,pch=17,col=1)
points(0.05,0.75,pch=18,col=1)
points(0.05,0.65,pch=16,col=2)
points(0.05,0.55,pch=17,col=2)
points(0.05,0.45,pch=18,col=2)

text(0.1,0.95,expression("Unadjusted, "~A["1"]),pos=4)
text(0.1,0.85,expression("Unadjusted, "~A["2"]),pos=4)
text(0.1,0.75,expression("Unadjusted, "~A["3"]),pos=4)
text(0.1,0.65,expression("Adjusted, "~A["1"]),pos=4)
text(0.1,0.55,expression("Adjusted, "~A["2"]),pos=4)
text(0.1,0.45,expression("Adjusted, "~A["3"]),pos=4)

mtext("t-statistics for Covariate Balance Assessment",line=-1.5,outer=T,cex=1.2)

dev.off()



############# Assessment of Overlap Assumption

png("Plot/Overlap.png",height=2.5,width=8,unit="in",res=500)

PP <- dim(Data.Cluster)[2]-2

A.group <- list()
A.q <- c(0,quantile( NF.median[,2], (1:2)/3 ),1)
A.group[[1]] <- which( NF.median[,2] <= A.q[2])
A.group[[2]] <- which( A.q[2] < NF.median[,2] & NF.median[,2] <= A.q[3])
A.group[[3]] <- which( A.q[3] < NF.median[,2])

PS.loc <- c(which(CL.name=="PS.OL.Bin1.Est"), which(CL.name=="PS.OL.Bin2.Est"),
            which(CL.name=="PS.OL.Bin3.Est"))

TTT <- SI <- list()
LL <- 10

NF.Temp <- NF.median


Tstat.NA <- Tstat.A <- matrix(0,PP,3)
Size <- matrix(0,PP,6)

par(mar=c(3.5,3.5,2,1),mfrow=c(1,3))

kk <- 1
PS <- NF.Temp[,PS.loc[kk]]
hist(PS[A.group[[kk]]],breaks=seq(0,1,length=26),col=rgb(1,0,0,0.2),border=NA,
     prob=T,ylim=c(0,4),main="",xlab="",ylab="",xaxt="n",yaxt="n'")
par(new=T)
hist(PS[-A.group[[kk]]],breaks=seq(0,1,length=26),col=rgb(0,0,1,0.2),border=NA,
     prob=T,ylim=c(0,4),main="",xlab="",ylab="")
legend("topright",legend=c(sprintf("Group 1: [%0.3f,%0.3f]",min(PS[A.group[[kk]]]),max(PS[A.group[[kk]]])),
                           sprintf("Group 2,3: [%0.3f,%0.3f]",min(PS[-A.group[[kk]]]),max(PS[-A.group[[kk]]]))),
       col=c(rgb(1,0,0,0.2),rgb(0,0,1,0.2)),bty="n",
       pch=15,pt.cex=2)
title(main="Group 1 vs Group 2,3")

kk <- 2
PS <- NF.Temp[,PS.loc[kk]]
hist(PS[A.group[[kk]]],breaks=seq(0.1,0.45,length=26),col=rgb(1,0,0,0.2),border=NA,
     prob=T,ylim=c(0,15),main="",xlab="",ylab="",xaxt="n",yaxt="n'")
par(new=T)
hist(PS[-A.group[[kk]]],breaks=seq(0.1,0.45,length=26),col=rgb(0,0,1,0.2),border=NA,
     prob=T,ylim=c(0,15),main="",xlab="",ylab="")
legend("topleft",legend=c(sprintf("Group 2: [%0.3f,%0.3f]",min(PS[A.group[[kk]]]),max(PS[A.group[[kk]]])),
                          sprintf("Group 1,3: [%0.3f,%0.3f]",min(PS[-A.group[[kk]]]),max(PS[-A.group[[kk]]]))),
       col=c(rgb(1,0,0,0.2),rgb(0,0,1,0.2)),bty="n",
       pch=15,pt.cex=2)
title(main="Group 2 vs Group 1,3")

kk <- 3
PS <- NF.Temp[,PS.loc[kk]]
hist(PS[A.group[[kk]]],breaks=seq(0,1,length=26),col=rgb(1,0,0,0.2),border=NA,
     prob=T,ylim=c(0,6),main="",xlab="",ylab="",xaxt="n",yaxt="n'")
par(new=T)
hist(PS[-A.group[[kk]]],breaks=seq(0,1,length=26),col=rgb(0,0,1,0.2),border=NA,
     prob=T,ylim=c(0,6),main="",xlab="",ylab="")
legend("topright",legend=c(sprintf("Group 3: [%0.3f,%0.3f]",min(PS[A.group[[kk]]]),max(PS[A.group[[kk]]])),
                           sprintf("Group 1,2: [%0.3f,%0.3f]",min(PS[-A.group[[kk]]]),max(PS[-A.group[[kk]]]))),
       col=c(rgb(1,0,0,0.2),rgb(0,0,1,0.2)),bty="n",
       pch=15,pt.cex=2)
title(main="Group 3 vs Group 1,2")

mtext("PS",side=1,outer=TRUE,line=-1)
mtext("Density",side=2,outer=TRUE,line=-1.25)

dev.off()



############# Assessment of Stratified Inference

Res <- NF.median$Y-NF.median$OR.Est

PP <- function(jj,RN){
    jj <- jj+1
    if(length(unique(NF.median[,jj]))>2){
        plot( NF.median[, jj], Res,pch=16, cex=0.5, col=rgb(0,0,0,0.15), bg=rgb(1,0,0))
        title(main=RN)
        lines(locpoly(NF.median[, jj], Res,bandwidth = dpill(NF.median[, jj],Res),degree=0) ,col=2, lwd=2)
    } else {
        plot( as.numeric(NF.median[,jj]==max(NF.median[,jj])), Res,pch=16, cex=0.5, col=rgb(0,0,0,0.15), bg=rgb(1,0,0))
        title(main=RN)
        x0 <- 0 ; x1 <- 1
        y0 <- mean(Res[NF.median[,jj]==min(NF.median[,jj])])
        y1 <- mean(Res[NF.median[,jj]==max(NF.median[,jj])])
        slope <- (y1-y0)/(x1-x0)
        intercept <- y1-slope*x1
        abline(a=intercept,b=slope,col=2, lwd=2)
    }
    abline(h=0,col=4,lty=2)
    
}

RN <- c("Prop. of treatment",
        "Cluster size","Urban/Rural","Avg. household size",
        "Avg. num. of children",
        "Prop. of parents having no job",
        "Prop. of educated HH heads",
        "Prop. of educated HH respondents",
        "Avg. age of respondents",
        "Avg. age of children")

png("Plot/ResidualPlot.png",height=5,width=12,unit="in",res=500)

par(mar=c(3,3,2,0.5), oma=c(2,2,0.5,0.5))
layout(cbind(c(1,12),matrix(2:11,2,5,byrow=T)))

plot( NF.median$OR.Est, Res, pch=16, cex=0.5, col=rgb(0,0,0,0.15), bg=rgb(1,0,0))
title(main="OR Estimate")
lines(locpoly(NF.median$OR.Est, Res,bandwidth = dpill(NF.median$OR.Est,Res),degree=0) ,col=2, lwd=2)
abline(h=0,col=4,lty=2)

for(jj in 1:10){
    PP(jj,RN[jj])
}
mtext("OR Estimate / Regressor",side=1,outer=T,line=0.5)
mtext("Residual",side=2,outer=T,line=0.5)

dev.off()


############# Assessment of Partial Interference Using Geographical Distance

Coordinate <- read.csv("RULE_Merged.csv")
DIST <- distm(cbind(Coordinate$X,Coordinate$Y))
DIST <- c(DIST[upper.tri(DIST)])
mean(DIST<10000)

