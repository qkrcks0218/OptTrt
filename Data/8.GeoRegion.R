############################################
# By running this R-file, the simulation results are merged with GPS file "Data/GPS_2018/SNGE81FL.shp".
# The summary is saved as "RULE_Mergged.csv".
############################################

library(ggplot2)
library(dplyr)
library(maps)

source("../OptTrt_Sim_Function.R")

RULE.IPW.Test.Median <- RULE.OR.Test.Median <- RULE.DR.Test.Median <- RULE.Lasso.Test.Median <- RULE.RF.Test.Median <- matrix(0,213,10)
S.grid <- (64:73)/100

Data.Cluster <- read.csv("Sen_Training_1417.csv")
Test.Data <- read.csv("Sen_Test_18.csv")

for(jj in 1:10){
  
  S.quart <- S <- S.grid[jj]
  RULE.DR.Test <-  read.csv(sprintf("RULE/RULE_S%0.3d_Test.csv", S*1000))
  RULE.DR.Test.Median[,jj] <- Winsorization(apply(RULE.DR.Test,1,median))
  
  LRdata <- read.csv(sprintf("NF_LR/RULE_LR_Test.csv"))
  
  L1.name <- sprintf("Lasso.Rule.%0.3d",round(1000*S.quart))
  R1.name <- sprintf("RF.Rule.%0.3d",round(1000*S.quart))
  
  RULE.Lasso.Test.Median[,jj] <- round(LRdata[,L1.name],3)
  RULE.RF.Test.Median[,jj] <-    round(LRdata[,R1.name],3)
}

GPS <- rgdal::readOGR("GPS_2018/SNGE81FL.shp") 
R1 <- c("Dakar","Pikine","Mbour","Mbour","Mbour","Mbour","Thies","Thies","Thies","Thies","Thies","Thies","Pikine","Thies","Tivaouane","Tivaouane","Tivaouane","Tivaouane","Kebemer","Kebemer","Kebemer","Linguere","Linguere","Pikine","Linguere","Linguere","Louga","Louga","Louga","Louga","Louga","Louga","Louga","Fatick","Rufisque","Fatick","Fatick","Fatick","Fatick","Fatick","Fatick","Foundiougne","Foundiougne","Foundiougne","Foundiougne","Rufisque","Gossas","Gossas","Gossas","Kolda","Kolda","Kolda","Kolda","Kolda","Kolda","Velingara","Rufisque","Velingara","Velingara","Velingara","Velingara","Velingara","Velingara","Velingara","Medina Yoro Foula","Medina Yoro Foula","Matam","Rufisque","Matam","Matam","Matam","Matam","Matam","Matam","Kanel","Kanel","Kanel","Kanel","Rufisque","Kanel","Ranerou Ferlo","Ranerou Ferlo","Kaffrine","Kaffrine","Kaffrine","Kaffrine","Kaffrine","Birkilane","Birkilane","Rufisque","Koungheul","Koungheul","Koungheul","Koungheul","Maleme Hodar","Maleme Hodar","Maleme Hodar","Kedougou","Kedougou","Kedougou","Guediawaye","Kedougou","Kedougou","Kedougou","Kedougou","Salemata","Salemata","Saraya","Saraya","Saraya","Saraya","Dakar","Guediawaye","Saraya","Sedhiou","Sedhiou","Sedhiou","Sedhiou","Sedhiou","Bounkiling","Bounkiling","Bounkiling","Bounkiling","Bignona","Goudomp","Goudomp","Goudomp","Goudomp","Goudomp","Bignona","Bignona","Bignona","Bignona","Bignona","Oussouye","Ziguinchor","Dakar","Ziguinchor","Ziguinchor","Ziguinchor","Ziguinchor","Ziguinchor","Bambey","Bambey","Bambey","Diourbel","Diourbel","Dakar","Diourbel","Diourbel","Mbacke","Mbacke","Mbacke","Mbacke","Mbacke","Mbacke","Mbacke","Mbacke","Dakar","Mbacke","Dagana","Dagana","Dagana","Dagana","Podor","Podor","Podor","Podor","Podor","Dakar","Podor","Saint-Louis","Saint-Louis","Saint-Louis","Saint-Louis","Saint-Louis","Saint-Louis","Bakel","Tambacounda","Tambacounda","Dakar","Tambacounda","Tambacounda","Tambacounda","Tambacounda","Tambacounda","Goudiry","Goudiry","Koupentoum","Koupentoum","Koupentoum","Pikine","Koupentoum","Kaolack","Kaolack","Kaolack","Kaolack","Kaolack","Kaolack","Kaolack","Kaolack","Kaolack","Pikine","Nioro du Rip","Nioro du Rip","Nioro du Rip","Nioro du Rip","Nioro du Rip","Nioro du Rip","Guinguineo","Mbour","Mbour","Mbour")
R2 <- c("Dakar","Dakar","Thies","Thies","Thies","Thies","Thies","Thies","Thies","Thies","Thies","Thies","Dakar","Thies","Thies","Thies","Thies","Thies","Louga","Louga","Louga","Louga","Louga","Dakar","Louga","Louga","Louga","Louga","Louga","Louga","Louga","Louga","Louga","Fatick","Dakar","Fatick","Fatick","Fatick","Fatick","Fatick","Fatick","Fatick","Fatick","Fatick","Fatick","Dakar","Fatick","Fatick","Fatick","Kolda","Kolda","Kolda","Kolda","Kolda","Kolda","Kolda","Dakar","Kolda","Kolda","Kolda","Kolda","Kolda","Kolda","Kolda","Kolda","Kolda","Matam","Dakar","Matam","Matam","Matam","Matam","Matam","Matam","Matam","Matam","Matam","Matam","Dakar","Matam","Matam","Matam","Kaffrine","Kaffrine","Kaffrine","Kaffrine","Kaffrine","Kaffrine","Kaffrine","Dakar","Kaffrine","Kaffrine","Kaffrine","Kaffrine","Kaffrine","Kaffrine","Kaffrine","Kedougou","Kedougou","Kedougou","Dakar","Kedougou","Kedougou","Kedougou","Kedougou","Kedougou","Kedougou","Kedougou","Kedougou","Kedougou","Kedougou","Dakar","Dakar","Kedougou","Sedhiou","Sedhiou","Sedhiou","Sedhiou","Sedhiou","Sedhiou","Sedhiou","Sedhiou","Sedhiou","Ziguinchor","Sedhiou","Sedhiou","Sedhiou","Sedhiou","Sedhiou","Ziguinchor","Ziguinchor","Ziguinchor","Ziguinchor","Ziguinchor","Ziguinchor","Ziguinchor","Dakar","Ziguinchor","Ziguinchor","Ziguinchor","Ziguinchor","Ziguinchor","Diourbel","Diourbel","Diourbel","Diourbel","Diourbel","Dakar","Diourbel","Diourbel","Diourbel","Diourbel","Diourbel","Diourbel","Diourbel","Diourbel","Diourbel","Diourbel","Dakar","Diourbel","Saint-Louis","Saint-Louis","Saint-Louis","Saint-Louis","Saint-Louis","Saint-Louis","Saint-Louis","Saint-Louis","Saint-Louis","Dakar","Saint-Louis","Saint-Louis","Saint-Louis","Saint-Louis","Saint-Louis","Saint-Louis","Saint-Louis","Tambacounda","Tambacounda","Tambacounda","Dakar","Tambacounda","Tambacounda","Tambacounda","Tambacounda","Tambacounda","Tambacounda","Tambacounda","Tambacounda","Tambacounda","Tambacounda","Dakar","Tambacounda","Kaolack","Kaolack","Kaolack","Kaolack","Kaolack","Kaolack","Kaolack","Kaolack","Kaolack","Dakar","Kaolack","Kaolack","Kaolack","Kaolack","Kaolack","Kaolack","Kaolack","Thies","Thies","Thies")

coordinate <- cbind(GPS$DHSCLUST,GPS$LONGNUM,GPS$LATNUM)
TESTID <- read.csv("TestSetID.csv")

CL.ID <- TESTID[,1]

Sort.coordinate <- matrix(0,length(CL.ID),14)
for(ii in 1:dim(Sort.coordinate)[1]){
  Sort.coordinate[ii,] <- c(coordinate[coordinate[,1]==CL.ID[ii],],
                            as.numeric(TESTID[ii,-1]))
}

AA <- cbind(Sort.coordinate,R1,R2,
            RULE.Lasso.Test.Median, RULE.RF.Test.Median,RULE.DR.Test.Median)
colnames(AA)[1:16] <- c("ID","X","Y",
                        colnames(TESTID)[-1],"R1","R2")
colnames(AA)[17:26]    <- paste("RULE.Lasso.Test.Median.",64:73,sep="")
colnames(AA)[17:26+10] <- paste("RULE.RF.Test.Median.",64:73,sep="")
colnames(AA)[17:26+20] <- paste("RULE.DR.Test.Median.",64:73,sep="")

write.csv(AA,"RULE_Merged.csv",row.names=F)




