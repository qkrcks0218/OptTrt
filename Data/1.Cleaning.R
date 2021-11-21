############################################
# By running this R-file, the DHS data sets are cleaned.
# The training and test data sets are saved as "Sen_Training_1417.csv" and "Sen_Test_18.csv", respectively.
# The test set including the identifier is saved as "TestSetID.csv".
############################################

##########################################
# Download the following data sets from https://dhsprogram.com/data/available-datasets.cfm
# SNKR70FL.DTA , SNKR7HFL.DTA , SNKR7IFL.DTA , SNKR7ZFL.DTA , SNKR81FL.DTA
# Place the data sets in the following paths.
# Data/DHS/2014/SNKR70FL.DTA
# Data/DHS/2015/SNKR7HFL.DTA
# Data/DHS/2016/SNKR7IFL.DTA
# Data/DHS/2017/SNKR7ZFL.DTA
# Data/DHS/2018/SNKR81FL.DTA
# 
# Next, download "SNGE81FL.ZIP" from https://dhsprogram.com/data/dataset/Senegal_Continuous-DHS_2018.cfm?flag=1
# Unzip the zip file and place the 10 files including "SNGE81FL.shp" are in "Data/GPS_2018" folder.
##########################################


library(readstata13)

LData <- list()
Data.List <- c("2014/SNKR70FL.DTA",
               "2015/SNKR7HFL.DTA",
               "2016/SNKR7IFL.DTA",
               "2017/SNKR7ZFL.DTA",
               "2018/SNKR81FL.DTA")

for(jj in 1:5){
    
    ChildRaw <- read.dta13(sprintf("DHS/%s",Data.List[jj]))
    ShortChildRaw <- ChildRaw[,c("v001","v002","v003","v008",
                                 "v012","v024","v025","v113","v116",
                                 "v133","v136","v137","v138","v152",
                                 "v160","v702","v705","v717",
                                 "b3","b5","b9","h11")]
    colnames(ShortChildRaw) <- c("clusterid","HHid","indid","InterviewDate",
                                 "X_resp_age","C_region","C_residencetype","X_water","X_toilet",
                                 "X_resp_edu","X_HHsize","X_HHchildsize","X_HHWsize","X_head_age",
                                 "X_share_t","X_head_edu","X_head_job","X_resp_job",
                                 "X_birth","X_alive","X_res","Y_D")
    ShortChildRaw$X_age <- ShortChildRaw$InterviewDate-ShortChildRaw$X_birth
    ShortChildRaw <- ShortChildRaw[,-which(colnames(ShortChildRaw)=="InterviewDate")]
    ShortChildRaw <- ShortChildRaw[,-which(colnames(ShortChildRaw)=="X_birth")] 
    
    LData[[jj]] <- ShortChildRaw
    
    LData[[jj]]$clusterid <- paste(substring(Data.List[jj],1,4),"_",LData[[jj]]$clusterid,sep="")

}

OData <- rbind(LData[[1]],LData[[2]],LData[[3]],LData[[4]],LData[[5]])
OData$X_share_t[is.na(OData$X_share_t)] <- "yes"
OData$X_head_edu[is.na(OData$X_head_edu) | OData$X_head_edu==98 ] <- 0
OData$X_resp_edu[is.na(OData$X_resp_edu) | OData$X_resp_edu==98 ] <- 0
OData$X_head_job <- as.character(OData$X_head_job)
OData$X_resp_job <- as.character(OData$X_resp_job)
OData$X_head_job[ is.na(OData$X_head_job) ] <- "did not work"
OData$X_resp_job[ is.na(OData$X_resp_job) ] <- "did not work"

Data <- OData[!is.na(OData$Y_D) & OData$Y_D!="don't know" & OData$X_alive=="yes" &           # Remove invalid data
                  (OData$X_water!="not a dejure resident" |
                       OData$X_water!="not de jure") & !is.na(OData$X_water) &
                  OData$X_toilet!="not a dejure resident" & !is.na(OData$X_toilet) &
                  OData$X_share_t!="not a dejure resident" & 
                  !is.na(OData$X_res) & OData$X_res=="respondent" ,]



###### Outcome ####### 

Data$Y <- 1-as.numeric(Data$Y_D=="yes, last two weeks")                                      # Indicator whether child had diarrhea in the past two weeks

###### Treatment ####### 

Data$X_water_ind <- as.numeric(Data$X_water=="piped into dwelling"|
                                   Data$X_water=="piped to yard/plot") # 2.7%

Data$X_share_t <- as.numeric(Data$X_share_t=="no") # even, weak
Data$X_toilet <- as.numeric(Data$X_toilet=="flush toilet"|
                                Data$X_toilet=="flush to piped sewer system"|
                                Data$X_toilet=="flush to septic tank"|
                                Data$X_toilet=="flush to pit latrine"|
                                Data$X_toilet=="flush to somewhere else"|
                                Data$X_toilet=="flush, don't know where ")
Data$A <- as.numeric( Data$X_toilet*Data$X_share_t+Data$X_water_ind > 0 )

###### Covariate #######

Data$C_rt_C <- as.numeric(Data$C_residencetype=="urban")                                     # Indicator of urban area
Data$C_size <- 0                                                                             # Make empty vector to contain cluster size later
Data$X_resp_edu <- as.numeric(Data$X_resp_edu>0)
Data$X_head_edu <- as.numeric(Data$X_head_edu>0)
Data$X_edu <- Data$X_resp_edu*Data$X_head_edu
J1 <- as.numeric(Data$X_resp_job=="professional/technical/managerial")
J2 <- as.numeric(Data$X_head_job=="professional/technical/managerial")
NJ1 <- as.numeric(Data$X_resp_job=="did not work"|Data$X_resp_job=="not working")
NJ2 <- as.numeric(Data$X_head_job=="did not work")

Data$X_HEjob <- as.numeric(J1+J2>1)
Data$X_Nojob <- NJ1*NJ2 

###### Child-level Data #######

Child.Data <- Data[,c("Y","A","clusterid","C_size",
                      "C_rt_C",
                      "X_HHsize",
                      "X_HHchildsize",
                      "X_Nojob",
                      "X_head_edu","X_resp_edu",
                      "X_resp_age","X_age",
                      "HHid")]
Child.Data <- data.frame(Child.Data)
colnames(Child.Data) <- c("Y","A","clusterid","C_size",
                          "C_rt_C",
                          "X_HHsize",
                          "X_HHchildsize",
                          "X_Nojob",
                          "X_head_edu","X_resp_edu",
                          "X_resp_age","X_age",
                          "HHid")

###### Household-level Data #######

HH.Data <- aggregate(.~HHid+clusterid,data=Child.Data,FUN=mean)


# Outcome

HH.Data$Y <- floor(HH.Data$Y) # Household having at least one sick child is coded as Y=0

# Covariate

HH.Data$C_size <- rep(0,dim(HH.Data)[1])
CID <- sort(unique(Child.Data$clusterid))
dim(HH.Data)
s <- 0
for(ii in 1:length(CID)){
    M.temp <- sum( HH.Data$clusterid==CID[ii] )
    HH.Data$C_size[(s+1):(s+M.temp)] <- M.temp
    s <- s+M.temp
}


###### Cluster-level Data #######


Cl.Data <- aggregate(.~clusterid,data=HH.Data,FUN="mean")[,-2]                          # Remove clusterid and HHid

Cl.Data.Training <- Cl.Data[as.numeric(substring(Cl.Data$clusterid,1,4))<2018,-1]
Cl.Data.Test <- Cl.Data[substring(Cl.Data$clusterid,1,4)=="2018",-1]

write.csv(Cl.Data.Training,"Sen_Training_1417.csv",row.names=FALSE)
write.csv(Cl.Data.Test,"Sen_Test_18.csv",row.names=FALSE)

TEST.ID <- Cl.Data[substring(Cl.Data$clusterid,1,4)=="2018",]
TEST.ID[,1] <- apply(matrix(TEST.ID[,1],213,1),1,
                     function(x){substring(x,6,nchar(x))})
write.csv(TEST.ID,"TestSetID.csv",row.names=FALSE)






