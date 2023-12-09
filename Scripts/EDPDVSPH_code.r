#Upload csv files 
Merced_Species_List<-read.csv("~/Data/csv_files/Vernal_Pool_CommunityDiversity/Merced_Species_List.csv")
FULL_COMMUNITY_csv<-read.csv("~/Data/csv_files/Vernal_Pool_CommunityDiversity/FULL_COMMUNITY_csv.csv")
COMPETITIVETRAITTEST <- read_csv("~/Data/csv_files/Vernal_Pool_CommunityDiversity/COMPETITIVETRAITTEST_REDONE.csv")
ENV <- read_csv("~/Data/csv_files/Vernal_Pool_CommunityDiversity/ENVIRONMENTAL_VAR.csv")
#Install Packages
install.packages(c("data.table","readr","caper","picante","V.PhyloMaker2","sjPlot","ggplot","Biostrings","lmtest","vegan","ecodist","spaa"))
lapply(c("data.table","readr","caper","picante","V.PhyloMaker2","sjPlot","ggplot","lmtest","Biostrings","vegan","ecodist","spaa"), require, character.only = TRUE)
list.of.packages <- c("data.table","readr","caper","picante","V.PhyloMaker2","sjPlot","ggplot","lmtest","Biostrings","vegan","ecodist","spaa")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
BiocManager::install("ggtree")
library(ggtree)
#Make Phylogeny using Vphylo.maker2
SynthPhylo<-V.PhyloMaker2::phylo.maker(Merced_Species_List,tree=GBOTB.extended.TPL,nodes=nodes.info.1.TPL,output.sp.list = TRUE,output.tree = TRUE,scenarios = "S3")
SynthTree<-SynthPhylo$scenario.3

#Save tree in newick format
write.tree(phy,"~/VP_synthesis_Phylogeny.tre")

#Pull out only Merced Clusters
cluster<-FULL_COMMUNITY_csv[,9:51]
rownames(cluster)<-FULL_COMMUNITY_csv$Clusters

###############################
##Community Turnover Analysis##
###############################

#Calculate Jaccard and Sorenson's Dissimilarity Index
pzy<-FULL_COMMUNITY_csv[493:519,]
pyz_f<-pzy[,c(-1,-2,-3,-4,-5,-6,-7,-8)]
rownames(pyz_f)<-pzy$Clusters
sorensen<-distance(pzy_f,method="sorensen")
jaccard<-vegdist(pzy_f,method="jaccard")
jaccard_list<-dist2list(jaccard)
sorensen_list<-dist2list(sorensen)
beta<-merge(sorensen_list,jaccard_list,by=c("col","row"))
for(i in 1:nrow(ju)){
    ju[i,"Pool"]<-ifelse((substring(ju[i,1],first=1,last=1))==(substring(ju[i,2],first=1,last=1)),"SAME","NO")
    ju[i,"Zone"]<-ifelse((substring(ju[i,1],first=4,last=4))==(substring(ju[i,2],first=4,last=4)),"SAME","NO")
    ju[i,"Year"]<-ifelse((substring(ju[i,1],first=2,last=2))==(substring(ju[i,2],first=2,last=2)),"SAME","NO")
}
for(i in 1:nrow(ju)){
    if (ju[i,"Pool"]=="SAME" & ju[i,"Zone"]=="SAME" & ju[i,"Year"]=="NO"){
        ju[i,"type"]<-"HT"
    } else if (ju[i,"Pool"]=="SAME" & ju[i,"Zone"]=="NO" & ju[i,"Year"]=="NO") {
        ju[i,"type"]<-"VT"
    } else if (ju[i,"Pool"]=="SAME" & ju[i,"Zone"]=="NO" & ju[i,"Year"]=="SAME") {
        ju[i,"type"]<-"V"
    } else if (ju[i,"Pool"]=="NO" & ju[i,"Zone"]=="NO" & ju[i,"Year"]=="SAME") {
        ju[i,"type"]<-"HV"
    } else if (ju[i,"Pool"]=="NO" & ju[i,"Zone"]=="SAME" & ju[i,"Year"]=="SAME"){
        ju[i,"type"]<-"H"
    } else if (ju[i,"Pool"]=="NO" & ju[i,"Zone"]=="SAME" & ju[i,"Year"]=="NO"){
        ju[i,"type"]<-"HHT"
    } else if (ju[i,"Pool"]=="NO" & ju[i,"Zone"]=="NO" & ju[i,"Year"]=="NO"){
        ju[i,"type"]<-"VVT"
    } else {
        ju[i,"type"]<-"NONE"
    }
}

################################
###Phylogenetic Beta Analysis###
################################
clusterZone_Pool<-as.matrix(FULL_COMMUNITY_csv[17:25,9:51])
row.names(clusterZone_Pool)<-<-FULL_COMMUNITY_csv[17:25,"Clusters"]
beta_zone.year <- comdist(clusterZone_Pool,phydist)
beta_zone.year.clusters <- hclust(beta_zone.year)

clusterSeason_Pool<-as.matrix(FULL_COMMUNITY_csv[458:466,9:51])
row.names(clusterSeason_Pool)<-FULL_COMMUNITY_csv[458:466,"Clusters"]
beta_zone.year <- comdist(clusterSeason_Pool,phydist)
beta_season.year.clusters <- hclust(beta_zone.year)

################################
#####MPD and MNTD Analysis######
################################

#Calculate MPD and MNTD for Each Zone
clusterZone<-as.matrix(FULL_COMMUNITY_csv[4:9,9:51])
row.names(clusterZone)<-FULL_COMMUNITY_csv[4:9,"Clusters"]
mpd.result_A <- ses.mpd(clusterZone,phydist,null.model="taxa.labels",abundance.weighted=FALSE,runs=1000)
mntd.result_A <- ses.mntd(clusterZone,phydist,null.model="taxa.labels",abundance.weighted=FALSE,runs=1000)
#Merge MPD and MNTD
mpd.result<-cbind(cluster=clusterZone$Clusters,mpd.result)
mntd.result<-cbind(cluster=clusterZone$Clusters,mntd.result)
MPD_MNTD_result <- merge(mpd.result,mntd.result,by="cluster",sort=FALSE)

#Create Column with interpretation of Results
for (i in 1:nrow(MPD_MNTD_result)){
    MPD_MNTD_result$Tree_Wide[i] <- if (is.na(MPD_MNTD_result$mpd.obs.p[i])){
        print("Missing")
    } else if (MPD_MNTD_result$mpd.obs.p[i] < 0.05 && MPD_MNTD_result$mpd.obs.z[i] < 0){
    print("Clustered")
    } else if (MPD_MNTD_result$mpd.obs.p[i] > 0.95 && MPD_MNTD_result$mpd.obs.z[i] > 0){
    print("Over-Dispersed")
    } else {
    print("Random")
    }
}
for (i in 1:nrow(MPD_MNTD_result)){
    MPD_MNTD_result$Tips[i] <- if (is.na(MPD_MNTD_result$mntd.obs.p[i])){
        print("Missing")
    } else if (MPD_MNTD_result$mntd.obs.p[i] < 0.05 && MPD_MNTD_result$mntd.obs.z[i] < 0){
    print("Clustered")
    } else if (MPD_MNTD_result$mntd.obs.p[i] > 0.95 && MPD_MNTD_result$mntd.obs.z[i] > 0){
    print("Over-Dispersed")
    } else {
    print("Random")
    }
}

#Make table of only SES annd p values of MPD and MNTD results
Pool.Zone_cluster <- MPD_MNTD_result[1:15,c(1,2,7,8,15,16,17,18,19)]
tab_df(Pool.Zone_cluster,title="Table: MPD and MNTD Results *Modify*",file="MPD_MNTD_result.docx")

###########################################
#MNTD,MPD,PDfaith for Competition Analysis#
###########################################
#Calculate MNTD,MPD, and PD standardized effect sizes for clusters of zone,pool,season, and year (n=81)
clusterComp<-as.matrix(FULL_COMMUNITY_csv[377:457,9:51])
row.names(clusterComp)<-<-FULL_COMMUNITY_csv[377:457,"Clusters"]
mpd.result_A <- ses.mpd(clusterComp,phydist,null.model="taxa.labels",abundance.weighted=FALSE,runs=1000)
mntd.result_A <- ses.mntd(clusterComp,phydist,null.model="taxa.labels",abundance.weighted=FALSE,runs=1000)
pd.result_A <- ses.pd(clus,SynthTree,null.model="richness",runs=1000)
Phylo_Frame<-data.frame(row.names(clusterComp),mpd.result_A$mpd.obs.z,mntd.result_A$mntd.obs.z,pd.result_A$pd.obs.z)
CompFrame<-merge(Phylo_Frame,COMPETITIVETRAITTEST_REDONE,by="Clusters")

#Perform Linear Models for Competition by Zone
mntd_B<-lm(mntd.obs.z~PC_mean,data=CompFrame[CompFrame$Zone=="B",])
mntd_E<-lm(mntd.obs.z~PC_mean,data=CompFrame[CompFrame$Zone=="E",])
mntd_U<-lm(mntd.obs.z~PC_mean,data=CompFrame[CompFrame$Zone=="U",])
mpd_B<-lm(mpd.obs.z~PC_mean,data=CompFrame[CompFrame$Zone=="B",])
mpd_E<-lm(mpd.obs.z~PC_mean,data=CompFrame[CompFrame$Zone=="E",])
mpd_U<-lm(mpd.obs.z~PC_mean,data=CompFrame[CompFrame$Zone=="U",])
mntd_B<-lm(mntd.obs.z~PC_variance,data=CompFrame[CompFrame$Zone=="B",])
mntd_E<-lm(mntd.obs.z~PC_variance,data=CompFrame[CompFrame$Zone=="E",])
mntd_U<-lm(mntd.obs.z~PC_variance,data=CompFrame[CompFrame$Zone=="U",])
mpd_B<-lm(mpd.obs.z~PC_variance,data=CompFrame[CompFrame$Zone=="B",])
mpd_E<-lm(mpd.obs.z~PC_variance,data=CompFrame[CompFrame$Zone=="E",])
mpd_U<-lm(mpd.obs.z~PC_variance,data=CompFrame[CompFrame$Zone=="U",])
#Perform Linear Model for Innundation by Zone (only Bottom and Edge given Upland doesn't flood)
mntd_B_I<-lm(mntd.obs.z~innundation_length,data=CompFrame[CompFrame$Zone=="B",])
mntd_E_I<-lm(mntd.obs.z~innundation_length,data=CompFrame[CompFrame$Zone=="E",])
mpd_B_I<-lm(mpd.obs.z~innundation_length,data=CompFrame[CompFrame$Zone=="B",])
mpd_E_I<-lm(mpd.obs.z~innundation_length,data=CompFrame[CompFrame$Zone=="E",])


###############################################################
#Mantel Tests of Community Turnover and Phylogenetic Diversity#
###############################################################
ENV <- read_csv("Documents/ENVIRONMENTAL_VAR.csv")
#Create Environmental Distance Matrix for Seasonal and Annual Temperature and Precipitation
season_mean<-ecodist::distance(ENV$SEASON_MEANTEMP, method = "euclidean")
season_max<-ecodist::distance(ENV$SEASON_MAXTEMP, method = "euclidean")
season_precip<-ecodist::distance(ENV$SEASON_PRECIP, method = "euclidean")
annual_mean<-ecodist::distance(ENV$YEAR_MEANTEMP, method = "euclidean")
annual_max<-ecodist::distance(ENV$YEAR_MAXTEMP, method = "euclidean")
annual_precip<-ecodist::distance(ENV$YEAR_PRECIP, method = "euclidean")

#Calculate Jaccard Dissimilarity Index for Whole Pool Community
clusSeason_Year<-FULL_COMMUNITY_csv[467:493,9:51]
rownames(clusSeason_Year)<-FULL_COMMUNITY_csv[467:493,1]
jaccardClus<-as.matrix(vegdist(clusSeason_Year, method="jaccard"))

#Calculate Faith's PD for Whole Pool Community
Faithpd<-pd(clusSeason,SynthTree,include.root=TRUE)
faithClus<-as.matrix(ecodist::distance(FaithClus$PD, method = "euclidean"))

#Perform Mantel Test for Whole Pool Community
s_mean_mantel<-mantel_test(season_mean,jaccardClus,nboot=1000)
s_max_mantel<-mantel_test(season_max,jaccardClus,nboot=1000)
s_precip_mantel<-mantel_test(season_precip,jaccardClus,nboot=1000)
a_mean_mantel<-mantel_test(annual_mean,jaccardClus,nboot=1000)
a_max_mantel<-mantel_test(annual_max,jaccardClus,nboot=1000)
a_precip_mantel<-mantel_test(annual_precip,jaccardClus,nboot=1000)

s_mean_mantel<-mantel_test(season_mean,faithClus,nboot=1000)
s_max_mantel<-mantel_test(season_max,faithClus,nboot=1000)
s_precip_mantel<-mantel_test(season_precip,faithClus,nboot=1000)
a_mean_mantel<-mantel_test(annual_mean,faithClus,nboot=1000)
a_max_mantel<-mantel_test(annual_max,faithClus,nboot=1000)
a_precip_mantel<-mantel_test(annual_precip,faithClus,nboot=1000)

#######################################################################
#Mantel Tests of Community Turnover and Phylogenetic Diversity by Zone#
#######################################################################

#Calculate Jaccard Dissimilarity Index for each Zone
clusSeason_Year<-FULL_COMMUNITY_csv[467:493,9:51]
rownames(clusSeason_Year)<-FULL_COMMUNITY_csv[467:493,1]
jaccardClus_B<-as.matrix(vegdist(clusSeason_Year[clusSeason_Year$Zone="B",], method="jaccard"))
jaccardClus_E<-as.matrix(vegdist(clusSeason_Year[clusSeason_Year$Zone="E",], method="jaccard"))
jaccardClus_U<-as.matrix(vegdist(clusSeason_Year[clusSeason_Year$Zone="U",], method="jaccard"))

#Calculate Faith's PD for Whole Pool Community
Faithpd<-pd(clusSeason,SynthTree,include.root=TRUE)
faithClus_B<-as.matrix(ecodist::distance(FaithClus[Faithpd$Zone=="B","PD"], method = "euclidean"))
faithClus_E<-as.matrix(ecodist::distance(FaithClus[Faithpd$Zone=="E","PD"], method = "euclidean"))
faithClus_U<-as.matrix(ecodist::distance(FaithClus[Faithpd$Zone=="U","PD"], method = "euclidean"))

#Perform Mantel Test for Each Zonal Community
s_mean_mantel_B<-mantel_test(season_mean,jaccardClus_B,nboot=1000)
s_max_mantel_B<-mantel_test(season_max,jaccardClus_B,nboot=1000)
s_precip_mantel_B<-mantel_test(season_precip,jaccardClus_B,nboot=1000)
a_mean_mantel_B<-mantel_test(annual_mean,jaccardClus_B,nboot=1000)
a_max_mantel_B<-mantel_test(annual_max,jaccardClus_B,nboot=1000)
a_precip_mantel_B<-mantel_test(annual_precip,jaccardClus_B,nboot=1000)

s_mean_mantel_B<-mantel_test(season_mean,faithClus_B,nboot=1000)
s_max_mantel_B<-mantel_test(season_max,faithClus_B,nboot=1000)
s_precip_mantel_B<-mantel_test(season_precip,faithClus_B,nboot=1000)
a_mean_mantel_B<-mantel_test(annual_mean,faithClus_B,nboot=1000)
a_max_mantel_B<-mantel_test(annual_max,faithClus_B,nboot=1000)
a_precip_mantel_B<-mantel_test(annual_precip,faithClus_B,nboot=1000)

s_mean_mantel_E<-mantel_test(season_mean,jaccardClus_E,nboot=1000)
s_max_mantel_E<-mantel_test(season_max,jaccardClus_E,nboot=1000)
s_precip_mantel_E<-mantel_test(season_precip,jaccardClus_E,nboot=1000)
a_mean_mantel_E<-mantel_test(annual_mean,jaccardClus_E,nboot=1000)
a_max_mantel_E<-mantel_test(annual_max,jaccardClus_E,nboot=1000)
a_precip_mantel_E<-mantel_test(annual_precip,jaccardClus_E,nboot=1000)

s_mean_mantel_E<-mantel_test(season_mean,faithClus_E,nboot=1000)
s_max_mantel_E<-mantel_test(season_max,faithClus_E,nboot=1000)
s_precip_mantel_E<-mantel_test(season_precip,faithClus_E,nboot=1000)
a_mean_mantel_E<-mantel_test(annual_mean,faithClus_E,nboot=1000)
a_max_mantel_E<-mantel_test(annual_max,faithClus_E,nboot=1000)
a_precip_mantel_E<-mantel_test(annual_precip,faithClus_E,nboot=1000)

s_mean_mantel_U<-mantel_test(season_mean,jaccardClus_U,nboot=1000)
s_max_mantel_U<-mantel_test(season_max,jaccardClus_U,nboot=1000)
s_precip_mantel_U<-mantel_test(season_precip,jaccardClus_U,nboot=1000)
a_mean_mantel_U<-mantel_test(annual_mean,jaccardClus_U,nboot=1000)
a_max_mantel_U<-mantel_test(annual_max,jaccardClus_U,nboot=1000)
a_precip_mantel_U<-mantel_test(annual_precip,jaccardClus_U,nboot=1000)

s_mean_mantel_U<-mantel_test(season_mean,faithClus_U,nboot=1000)
s_max_mantel_U<-mantel_test(season_max,faithClus_U,nboot=1000)
s_precip_mantel_U<-mantel_test(season_precip,faithClus_U,nboot=1000)
a_mean_mantel_U<-mantel_test(annual_mean,faithClus_U,nboot=1000)
a_max_mantel_U<-mantel_test(annual_max,faithClus_U,nboot=1000)
a_precip_mantel_U<-mantel_test(annual_precip,faithClus_U,nboot=1000)

##########################################################
#Faith's PD of observation week per pool and year by zone#
##########################################################

#Find best fitting polynomial function

polyreg_pd_week <- lm(pd.obs.z ~ poly(week,k,raw=TRUE),data=ses.pd.info)
degree = 3
K = 10
k_list=list()
regress_list=list()
Breusch_Pagan_list=list()
mse=matrix(data=NA,nrow=K,ncol=degree)
for (x in PD_shuff_LIST){
    for(i in 1:K){
    folds <- cut(seq(1,nrow(x)),breaks=K,labels=FALSE)
    index <- which(folds==i,arr.ind=TRUE)
    test <- x[index,]
    train <- x[-index,]

    for(j in 1:degree){
        fittrain = lm(pd.obs.z ~ poly(week,j),data = train)
        fittest = predict(fittrain, new = test)
        mse[i,j]=mean((fittest-test$pd.obs.z)^2)
     }
    }
    mse_5 <- na.omit(mse)
    min_mse <- min(colMeans(mse_5))
    v <- colMeans(mse_5)
    k = which(v==min_mse)
    k_list[[x[1,3]]]<-k
    #Run Polynomial regression using resulting best fit polynomial for
    #each zone
    for (y in PD_ZONE){
        q = y[1,3]
        d = x[1,3]
        if(q==d){
            regress <- lm(pd.obs.z~poly(week,k,raw=TRUE),data = y)
            #Perform Breusch-Pagan test to determine if model residuals exhibit heteroscedasticity
            bpmodel <- bptest(regress)
            #Plot the Regression
            ggplot(y,aes(x=week,y=pd.obs.z))+
            geom_point()+
            stat_smooth(method="lm", formula = y ~ poly(x,k,raw=TRUE),linewidth=1)+
            xlab('Week')+
            ylab('Phylogenetic Diversity')
            #Save Regressions as pdf
            name=paste(d,"Polynomial_regression.pdf",sep="_")
            ggsave(
                filename=name,
                plot=last_plot(),
                device=pdf,
                dpi=300,
                limitsize=TRUE
            )
            }
        regress_list[[d]]<-regress
        Breusch_Pagan_list[[d]]<- bpmodel
    }
}

#Get Summaries of Each Regression
for (i in regression_list){
    print(summary(i))
}

