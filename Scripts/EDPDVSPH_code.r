#Upload csv files 
Merced_Species_List<-read_csv("~/Documents/Github_Projects/Vernal_Pool_CommunityDiversity/Data/csv_files/Merced_Species_List.csv")
FULL_COMMUNITY_csv<-read_csv("~/Documents/Github_Projects/Vernal_Pool_CommunityDiversity/Data/csv_files/FULL_COMMUNITY_csv.csv")
COMPETITIVETRAITTEST <- read_csv("~/Documents/Github_Projects/Vernal_Pool_CommunityDiversity/Data/csv_files/COMPETITIVETRAITTEST.csv")
ENV <- read_csv("~/Documents/Github_Projects/Vernal_Pool_CommunityDiversity/Data/csv_files/ENVIRONMENTAL_VAR.csv")
#Install Packages
install.packages(c("matrixcalc","data.table","readr","caper","picante","V.PhyloMaker2","sjPlot","ggplot","Biostrings","lmtest","vegan","ecodist","spaa"))
lapply(c("matrixcalc","data.table","readr","caper","picante","V.PhyloMaker2","sjPlot","ggplot","lmtest","Biostrings","vegan","ecodist","spaa"), require, character.only = TRUE)
list.of.packages <- c("matrixcalc","data.table","readr","caper","picante","V.PhyloMaker2","sjPlot","ggplot","lmtest","Biostrings","vegan","ecodist","spaa")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
BiocManager::install("ggtree")
library(ggtree)
#Make Phylogeny using Vphylo.maker2
SynthPhylo<-V.PhyloMaker2::phylo.maker(Merced_Species_List,tree=GBOTB.extended.TPL,nodes=nodes.info.1.TPL,output.sp.list = TRUE,output.tree = TRUE,scenarios = "S3")
SynthTree<-SynthPhylo$scenario.3

#Save tree in newick format
write.tree(SynthTree,"~/VP_synthesis_Phylogeny.tre")

###############################
##Community Turnover Analysis##
###############################

#Calculate Jaccard and Sorenson's Dissimilarity Index
pzy<-FULL_COMMUNITY_csv[493:519,]
pyz_f<-pzy[,9:50]
rownames(pyz_f)<-pzy$Clusters
sorensen<-ecodist::distance(pyz_f,method="sorensen")
jaccard<-vegan::vegdist(pyz_f,method="jaccard")
dist2list <- function (dist, tri=TRUE) {
  if (!class(dist) == "dist") { stop("Error: The input data must be a dist object.\n") }

  dat <- as.data.frame(as.matrix(dist))
  if (is.null(names(dat))) {
    rownames(dat) <- paste(1:nrow(dat))
  }
  value <- stack(dat)$values
  rnames <- rownames(dat)
  namecol <- expand.grid(rnames, rnames)
  colnames(namecol) <- c("col", "row")
  res <- data.frame(namecol, value)

  if(tri == TRUE){    # return only lower triangular part of dist
    res <- res[-which(upper.tri(as.matrix(dist), diag = T)), ]
  }

  return(res)
}
jaccard_list<-dist2list(jaccard,tri=T)
sorensen_list<-dist2list(sorensen,tri=T)
beta<-merge(sorensen_list,jaccard_list,by=c("col","row"))
for(i in 1:nrow(beta)){
    beta[i,"Pool"]<-ifelse((substring(beta[i,1],first=1,last=1))==(substring(beta[i,2],first=1,last=1)),"SAME","NO")
    beta[i,"Zone"]<-ifelse((substring(beta[i,1],first=4,last=4))==(substring(beta[i,2],first=4,last=4)),"SAME","NO")
    beta[i,"Year"]<-ifelse((substring(beta[i,1],first=2,last=2))==(substring(beta[i,2],first=2,last=2)),"SAME","NO")
}
for(i in 1:nrow(beta)){
    if (beta[i,"Pool"]=="SAME" & beta[i,"Zone"]=="SAME" & beta[i,"Year"]=="NO"){
        beta[i,"type"]<-"HT"
    } else if (beta[i,"Pool"]=="SAME" & beta[i,"Zone"]=="NO" & beta[i,"Year"]=="NO") {
        beta[i,"type"]<-"VT"
    } else if (beta[i,"Pool"]=="SAME" & beta[i,"Zone"]=="NO" & beta[i,"Year"]=="SAME") {
        beta[i,"type"]<-"V"
    } else if (beta[i,"Pool"]=="NO" & beta[i,"Zone"]=="NO" & beta[i,"Year"]=="SAME") {
        beta[i,"type"]<-"HV"
    } else if (beta[i,"Pool"]=="NO" & beta[i,"Zone"]=="SAME" & beta[i,"Year"]=="SAME"){
        beta[i,"type"]<-"H"
    } else if (beta[i,"Pool"]=="NO" & beta[i,"Zone"]=="SAME" & beta[i,"Year"]=="NO"){
        beta[i,"type"]<-"HHT"
    } else if (beta[i,"Pool"]=="NO" & beta[i,"Zone"]=="NO" & beta[i,"Year"]=="NO"){
        beta[i,"type"]<-"VVT"
    } else {
        beta[i,"type"]<-"NONE"
    }
}
beta<-beta[!(beta$type=="NONE"),]
colnames(beta)[3]<-"sorensen"
colnames(beta)[4]<-"jaccard"
beta_sorensen<-lm(sorensen~type,data=beta)
TukeyHSD(aov(beta_sorensen))
beta_jaccard<-lm(jaccard~type,data=beta)
TukeyHSD(aov(beta_jaccard))
################################
###Phylogenetic Beta Analysis###
################################
phydist<-cophenetic(SynthTree)
clusterZone_Pool<-as.matrix(FULL_COMMUNITY_csv[16:24,9:50])
names<-FULL_COMMUNITY_csv[16:24,1]
rownames(clusterZone_Pool)<-names$Clusters
beta_zone.year <- comdist(clusterZone_Pool,phydist)
beta_zone.year.clusters <- hclust(beta_zone.year)

clusterSeason_Pool<-as.matrix(FULL_COMMUNITY_csv[457:465,9:50])
names<-FULL_COMMUNITY_csv[457:465,1]
rownames(clusterSeason_Pool)<-names$Clusters
beta_zone.year <- comdist(clusterSeason_Pool,phydist)
beta_season.year.clusters <- hclust(beta_zone.year)

################################
#####MPD and MNTD Analysis######
################################

#Calculate MPD and MNTD for Each Zone
phydist<-cophenetic(SynthTree)
clusterZone<-as.matrix(FULL_COMMUNITY_csv[4:9,9:50])
names<-FULL_COMMUNITY_csv[4:9,1]
rownames(clusterZone)<-names$Clusters
mpd.result_A <- ses.mpd(clusterZone,phydist,null.model="taxa.labels",abundance.weighted=FALSE,runs=1000)
mntd.result_A <- ses.mntd(clusterZone,phydist,null.model="taxa.labels",abundance.weighted=FALSE,runs=1000)
#Merge MPD and MNTD
MPD_MNTD_result<-data.frame("Clusters"=rownames(mntd.result_A),"Richness"=mpd.result_A$ntaxa,"mntd.obs.z"=mntd.result_A$mntd.obs.z,"mntd.obs.p"=mntd.result_A$mntd.obs.p,"mpd.obs.z"=mpd.result_A$mpd.obs.z,"mpd.obs.p"=mpd.result_A$mpd.obs.p)
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

###########################################
#MNTD,MPD,PDfaith for Competition Analysis#
###########################################
#Calculate MNTD,MPD, and PD standardized effect sizes for clusters of zone,pool,season, and year (n=81)
clusterComp<-as.matrix(FULL_COMMUNITY_csv[376:456,9:50])
names<-FULL_COMMUNITY_csv[376:456,1]
rownames(clusterComp)<-names$Clusters
mpd.result_A <- ses.mpd(clusterComp,phydist,null.model="taxa.labels",abundance.weighted=FALSE,runs=1000)
mntd.result_A <- ses.mntd(clusterComp,phydist,null.model="taxa.labels",abundance.weighted=FALSE,runs=1000)
pd.result_A <- ses.pd(clusterComp,SynthTree,null.model="richness",runs=1000)
Phylo_Frame<-data.frame("Clusters"=rownames(clusterComp),"MPD"=mpd.result_A$mpd.obs.z,"MNTD"=mntd.result_A$mntd.obs.z,"PD"=pd.result_A$pd.obs.z)
CompFrame<-merge(Phylo_Frame,COMPETITIVETRAITTEST,by="Clusters")

#Multiple Regression Model of Phylogenetic, Taxonmic, and Functional Trait responnses to Zone, Season, and Zone by Season Interactions. 
CWM_H<-anova(lm(CWM.H~ZONE+SEASON+ZONE*SEASON,data=CompFrame))
CWM_LA<-anova(lm(CWM.LA~ZONE+SEASON+ZONE*SEASON,data=CompFrame))
PD_P<-anova(lm(PD~ZONE+SEASON+ZONE*SEASON,data=CompFrame))
MPD_P<-anova(lm(MPD~ZONE+SEASON+ZONE*SEASON,data=CompFrame))
MNTD_P<-anova(lm(MNTD~ZONE+SEASON+ZONE*SEASON,data=CompFrame))
Richness_T<-anova(lm(Richness~ZONE+SEASON+ZONE*SEASON,data=CompFrame))

#Perform Linear Models for Competition by Zone
mntd_B_M<-lm(MNTD~PC1_Mean,data=CompFrame[CompFrame$ZONE=="B",])
mntd_E_M<-lm(MNTD~PC1_Mean,data=CompFrame[CompFrame$ZONE=="E",])
mntd_U_M<-lm(MNTD~PC1_Mean,data=CompFrame[CompFrame$ZONE=="U",])
mpd_B_M<-lm(MPD~PC1_Mean,data=CompFrame[CompFrame$ZONE=="B",])
mpd_E_M<-lm(MPD~PC1_Mean,data=CompFrame[CompFrame$ZONE=="E",])
mpd_U_M<-lm(MPD~PC1_Mean,data=CompFrame[CompFrame$ZONE=="U",])
mntd_B_V<-lm(MNTD~PC1_Variance,data=CompFrame[CompFrame$ZONE=="B",])
mntd_E_V<-lm(MNTD~PC1_Variance,data=CompFrame[CompFrame$ZONE=="E",])
mntd_U_V<-lm(MNTD~PC1_Variance,data=CompFrame[CompFrame$ZONE=="U",])
mpd_B_V<-lm(MPD~PC1_Variance,data=CompFrame[CompFrame$ZONE=="B",])
mpd_E_V<-lm(MPD~PC1_Variance,data=CompFrame[CompFrame$ZONE=="E",])
mpd_U_V<-lm(MPD~PC1_Variance,data=CompFrame[CompFrame$ZONE=="U",])
#Perform Linear Model for Innundation by Zone (only Bottom and Edge given Upland doesn't flood)
mntd_B_I<-lm(MNTD~innundation_length,data=CompFrame[CompFrame$ZONE=="B",])
mntd_E_I<-lm(MNTD~innundation_length,data=CompFrame[CompFrame$ZONE=="E",])
mpd_B_I<-lm(MPD~innundation_length,data=CompFrame[CompFrame$ZONE=="B",])
mpd_E_I<-lm(MPD~innundation_length,data=CompFrame[CompFrame$ZONE=="E",])


###############################################################
#Mantel Tests of Community Turnover and Phylogenetic Diversity#
###############################################################
#Create Environmental Distance Matrix for Seasonal and Annual Temperature and Precipitation
season_mean<-as.matrix(ecodist::distance(ENV$SEASON_MEANTEMP, method = "euclidean"))
season_max<-as.matrix(ecodist::distance(ENV$SEASON_MAXTEMP, method = "euclidean"))
season_precip<-as.matrix(ecodist::distance(ENV$SEASON_PRECIP, method = "euclidean"))
annual_mean<-as.matrix(ecodist::distance(ENV$YEAR_MEANTEMP, method = "euclidean"))
annual_max<-as.matrix(ecodist::distance(ENV$YEAR_MAXTEMP, method = "euclidean"))
annual_precip<-as.matrix(ecodist::distance(ENV$YEAR_PRECIP, method = "euclidean"))

#Calculate Jaccard Dissimilarity Index for Whole Pool Community
clusSeason_Year<-FULL_COMMUNITY_csv[466:492,9:50]
names<-FULL_COMMUNITY_csv[466:492,1]
rownames(clusSeason_Year)<-names$Clusters
jaccardClus<-as.matrix(vegdist(clusSeason_Year, method="jaccard"))

#Calculate Faith's PD for Whole Pool Community
Faithpd <- ses.pd(clusSeason_Year,SynthTree,null.model="richness",runs=1000)
faithClus<-as.matrix(ecodist::distance(Faithpd$pd.obs.z, method = "euclidean"))
clusSeason_Year<-as.matrix(FULL_COMMUNITY_csv[466:492,9:50])
names<-FULL_COMMUNITY_csv[466:492,1]
rownames(clusSeason_Year)<-names$Clusters
mpd<- ses.mpd(clusSeason_Year,phydist,null.model="taxa.labels",abundance.weighted=FALSE,runs=1000)
MPDClus<-as.matrix(ecodist::distance(mpd$mpd.obs.z, method = "euclidean"))
mntd<- ses.mntd(clusSeason_Year,phydist,null.model="taxa.labels",abundance.weighted=FALSE,runs=1000)
MNTDClus<-as.matrix(ecodist::distance(mntd$mntd.obs.z, method = "euclidean"))


#Perform Mantel Test for Whole Pool Community
s_mean_mantel_J<-mantel.test(season_mean,jaccardClus,nperm=1000)
s_max_mantel_J<-mantel.test(season_max,jaccardClus,nperm=1000)
s_precip_mantel_J<-mantel.test(season_precip,jaccardClus,nperm=1000)
a_mean_mantel_J<-mantel.test(annual_mean,jaccardClus,nperm=1000)
a_max_mantel_J<-mantel.test(annual_max,jaccardClus,nperm=1000)
a_precip_mantel_J<-mantel.test(annual_precip,jaccardClus,nperm=1000)

s_mean_mantel_F<-mantel.test(season_mean,faithClus,nperm=1000)
s_max_mantel_F<-mantel.test(season_max,faithClus,nperm=1000)
s_precip_mantel_F<-mantel.test(season_precip,faithClus,nperm=1000)
a_mean_mantel_F<-mantel.test(annual_mean,faithClus,nperm=1000)
a_max_mantel_F<-mantel.test(annual_max,faithClus,nperm=1000)
a_precip_mantel_F<-mantel.test(annual_precip,faithClus,nperm=1000)

s_mean_mantel_MPD<-mantel.test(season_mean,MPDClus,nperm=1000)
s_max_mantel_MPD<-mantel.test(season_max,MPDClus,nperm=1000)
s_precip_mantel_MPD<-mantel.test(season_precip,MPDClus,nperm=1000)
a_mean_mantel_MPD<-mantel.test(annual_mean,MPDClus,nperm=1000)
a_max_mantel_MPD<-mantel.test(annual_max,MPDClus,nperm=1000)
a_precip_mantel_MPD<-mantel.test(annual_precip,MPDClus,nperm=1000)

s_mean_mantel_MNTD<-mantel.test(season_mean,MNTDClus,nperm=1000)
s_max_mantel_MNTD<-mantel.test(season_max,MNTDClus,nperm=1000)
s_precip_mantel_MNTD<-mantel.test(season_precip,MNTDClus,nperm=1000)
a_mean_mantel_MNTD<-mantel.test(annual_mean,MNTDClus,nperm=1000)
a_max_mantel_MNTD<-mantel.test(annual_max,MNTDClus,nperm=1000)
a_precip_mantel_MNTD<-mantel.test(annual_precip,MNTDClus,nperm=1000)

#######################################################################
#Mantel Tests of Community Turnover and Phylogenetic Diversity by Zone#
#######################################################################

#Calculate Jaccard Dissimilarity Index for each Zone
clusSeason_Year<-FULL_COMMUNITY_csv[376:456,]
names<-FULL_COMMUNITY_csv[376:456,1]
rownames(clusSeason_Year)<-names$Clusters
jaccardClus_B<-as.matrix(vegdist(clusSeason_Year[clusSeason_Year$Zone=="B",9:50], method="jaccard"))
jaccardClus_E<-as.matrix(vegdist(clusSeason_Year[clusSeason_Year$Zone=="E",9:50], method="jaccard"))
jaccardClus_U<-as.matrix(vegdist(clusSeason_Year[clusSeason_Year$Zone=="U",9:50], method="jaccard"))

#Calculate Faith's PD for Whole Pool Community
Faithpd <- ses.pd(clusSeason_Year[,9:50],SynthTree,null.model="richness",runs=1000)
fdf<-Faithpd
fdf$Zone<-FULL_COMMUNITY_csv[376:456,"Zone"]
faithClus_B<-as.matrix(ecodist::distance(fdf[fdf$Zone=="B","pd.obs.z"], method = "euclidean"))
faithClus_E<-as.matrix(ecodist::distance(fdf[fdf$Zone=="E","pd.obs.z"], method = "euclidean"))
faithClus_U<-as.matrix(ecodist::distance(fdf[fdf$Zone=="U","pd.obs.z"], method = "euclidean"))
clusSeason_Year<-as.matrix(FULL_COMMUNITY_csv[376:456,9:50])
names<-FULL_COMMUNITY_csv[376:456,1]
rownames(clusSeason_Year)<-names$Clusters
mpd<- ses.mpd(clusSeason_Year,phydist,null.model="taxa.labels",abundance.weighted=FALSE,runs=1000)
mpd$Zone<-FULL_COMMUNITY_csv[376:456,"Zone"]
mpd<-mpd[!is.na(mpd$mpd.obs),]
mntd<- ses.mntd(clusSeason_Year,phydist,null.model="taxa.labels",abundance.weighted=FALSE,runs=1000)
mntd$Zone<-FULL_COMMUNITY_csv[376:456,"Zone"]
mntd<-mntd[!is.na(mntd$mntd.obs),]
MPDClus_B<-as.matrix(ecodist::distance(mpd[mpd$Zone=="B","mpd.obs.z"], method = "euclidean"))
MNTDClus_B<-as.matrix(ecodist::distance(mntd[mntd$Zone=="B","mntd.obs.z"], method = "euclidean"))
MPDClus_E<-as.matrix(ecodist::distance(mpd[mpd$Zone=="E","mpd.obs.z"], method = "euclidean"))
MNTDClus_E<-as.matrix(ecodist::distance(mntd[mntd$Zone=="E","mntd.obs.z"], method = "euclidean"))
MPDClus_U<-as.matrix(ecodist::distance(mpd[mpd$Zone=="U","mpd.obs.z"], method = "euclidean"))
MNTDClus_U<-as.matrix(ecodist::distance(mntd[mntd$Zone=="U","mntd.obs.z"], method = "euclidean"))


#Perform Mantel Test for Each Zonal Community
s_mean_mantel_B_J<-mantel.test(season_mean,jaccardClus_B,nperm=1000)
s_max_mantel_B_J<-mantel.test(season_max,jaccardClus_B,nperm=1000)
s_precip_mantel_B_J<-mantel.test(season_precip,jaccardClus_B,nperm=1000)
a_mean_mantel_B_J<-mantel.test(annual_mean,jaccardClus_B,nperm=1000)
a_max_mantel_B_J<-mantel.test(annual_max,jaccardClus_B,nperm=1000)
a_precip_mantel_B_J<-mantel.test(annual_precip,jaccardClus_B,nperm=1000)

s_mean_mantel_B_F<-mantel.test(season_mean,faithClus_B,nperm=1000)
s_max_mantel_B_F<-mantel.test(season_max,faithClus_B,nperm=1000)
s_precip_mantel_B_F<-mantel.test(season_precip,faithClus_B,nperm=1000)
a_mean_mantel_B_F<-mantel.test(annual_mean,faithClus_B,nperm=1000)
a_max_mantel_B_F<-mantel.test(annual_max,faithClus_B,nperm=1000)
a_precip_mantel_B_F<-mantel.test(annual_precip,faithClus_B,nperm=1000)

s_mean_mantel_E_J<-mantel.test(season_mean,jaccardClus_E,nperm=1000)
s_max_mantel_E_J<-mantel.test(season_max,jaccardClus_E,nperm=1000)
s_precip_mantel_E_J<-mantel.test(season_precip,jaccardClus_E,nperm=1000)
a_mean_mantel_E_J<-mantel.test(annual_mean,jaccardClus_E,nperm=1000)
a_max_mantel_E_J<-mantel.test(annual_max,jaccardClus_E,nperm=1000)
a_precip_mantel_E_J<-mantel.test(annual_precip,jaccardClus_E,nperm=1000)

s_mean_mantel_E_F<-mantel.test(season_mean,faithClus_E,nperm=1000)
s_max_mantel_E_F<-mantel.test(season_max,faithClus_E,nperm=1000)
s_precip_mantel_E_F<-mantel.test(season_precip,faithClus_E,nperm=1000)
a_mean_mantel_E_F<-mantel.test(annual_mean,faithClus_E,nperm=1000)
a_max_mantel_E_F<-mantel.test(annual_max,faithClus_E,nperm=1000)
a_precip_mantel_E_F<-mantel.test(annual_precip,faithClus_E,nperm=1000)

s_mean_mantel_U_J<-mantel.test(season_mean,jaccardClus_U,nperm=1000)
s_max_mantel_U_J<-mantel.test(season_max,jaccardClus_U,nperm=1000)
s_precip_mantel_U_J<-mantel.test(season_precip,jaccardClus_U,nperm=1000)
a_mean_mantel_U_J<-mantel.test(annual_mean,jaccardClus_U,nperm=1000)
a_max_mantel_U_J<-mantel.test(annual_max,jaccardClus_U,nperm=1000)
a_precip_mantel_U_J<-mantel.test(annual_precip,jaccardClus_U,nperm=1000)

s_mean_mantel_U_F<-mantel.test(season_mean,faithClus_U,nperm=1000)
s_max_mantel_U_F<-mantel.test(season_max,faithClus_U,nperm=1000)
s_precip_mantel_U_F<-mantel.test(season_precip,faithClus_U,nperm=1000)
a_mean_mantel_U_F<-mantel.test(annual_mean,faithClus_U,nperm=1000)
a_max_mantel_U_F<-mantel.test(annual_max,faithClus_U,nperm=1000)
a_precip_mantel_U_F<-mantel.test(annual_precip,faithClus_U,nperm=1000)

s_mean_mantel_B_MPD<-mantel.test(season_mean,MPDClus_B,nperm=1000)
s_max_mantel_B_MPD<-mantel.test(season_max,MPDClus_B,nperm=1000)
s_precip_mantel_B_MPD<-mantel.test(season_precip,MPDClus_B,nperm=1000)
a_mean_mantel_B_MPD<-mantel.test(annual_mean,MPDClus_B,nperm=1000)
a_max_mantel_B_MPD<-mantel.test(annual_max,MPDClus_B,nperm=1000)
a_precip_mantel_B_MPD<-mantel.test(annual_precip,MPDClus_B,nperm=1000)
s_mean_mantel_B_MNTD<-mantel.test(season_mean,MNTDClus_B,nperm=1000)
s_max_mantel_B_MNTD<-mantel.test(season_max,MNTDClus_B,nperm=1000)
s_precip_mantel_B_MNTD<-mantel.test(season_precip,MNTDClus_B,nperm=1000)
a_mean_mantel_B_MNTD<-mantel.test(annual_mean,MNTDClus_B,nperm=1000)
a_max_mantel_B_MNTD<-mantel.test(annual_max,MNTDClus_B,nperm=1000)
a_precip_mantel_B_MNTD<-mantel.test(annual_precip,MNTDClus_B,nperm=1000)

season_mean<-as.matrix(ecodist::distance(ENV[-10,"SEASON_MEANTEMP"], method = "euclidean"))
season_max<-as.matrix(ecodist::distance(ENV[-10,"SEASON_MAXTEMP"], method = "euclidean"))
season_precip<-as.matrix(ecodist::distance(ENV[-10,"SEASON_PRECIP"], method = "euclidean"))
annual_mean<-as.matrix(ecodist::distance(ENV[-10,"YEAR_MEANTEMP"], method = "euclidean"))
annual_max<-as.matrix(ecodist::distance(ENV[-10,"YEAR_MAXTEMP"], method = "euclidean"))
annual_precip<-as.matrix(ecodist::distance(ENV[-10,"YEAR_PRECIP"], method = "euclidean"))
s_mean_mantel_U_MPD<-mantel.test(season_mean,MPDClus_U,nperm=1000)
s_max_mantel_U_MPD<-mantel.test(season_max,MPDClus_U,nperm=1000)
s_precip_mantel_U_MPD<-mantel.test(season_precip,MPDClus_U,nperm=1000)
a_mean_mantel_U_MPD<-mantel.test(annual_mean,MPDClus_U,nperm=1000)
a_max_mantel_U_MPD<-mantel.test(annual_max,MPDClus_U,nperm=1000)
a_precip_mantel_U_MPD<-mantel.test(annual_precip,MPDClus_U,nperm=1000)
s_mean_mantel_U_MNTD<-mantel.test(season_mean,MNTDClus_U,nperm=1000)
s_max_mantel_U_MNTD<-mantel.test(season_max,MNTDClus_U,nperm=1000)
s_precip_mantel_U_MNTD<-mantel.test(season_precip,MNTDClus_U,nperm=1000)
a_mean_mantel_U_MNTD<-mantel.test(annual_mean,MNTDClus_U,nperm=1000)
a_max_mantel_U_MNTD<-mantel.test(annual_max,MNTDClus_U,nperm=1000)
a_precip_mantel_U_MNTD<-mantel.test(annual_precip,MNTDClus_U,nperm=1000)
season_mean<-as.matrix(ecodist::distance(ENV[-1,"SEASON_MEANTEMP"], method = "euclidean"))
season_max<-as.matrix(ecodist::distance(ENV[-1,"SEASON_MAXTEMP"], method = "euclidean"))
season_precip<-as.matrix(ecodist::distance(ENV[-1,"SEASON_PRECIP"], method = "euclidean"))
annual_mean<-as.matrix(ecodist::distance(ENV[-1,"YEAR_MEANTEMP"], method = "euclidean"))
annual_max<-as.matrix(ecodist::distance(ENV[-1,"YEAR_MAXTEMP"], method = "euclidean"))
annual_precip<-as.matrix(ecodist::distance(ENV[-1,"YEAR_PRECIP"], method = "euclidean"))
s_mean_mantel_E_MPD<-mantel.test(season_mean,MPDClus_E,nperm=1000)
s_max_mantel_E_MPD<-mantel.test(season_max,MPDClus_E,nperm=1000)
s_precip_mantel_E_MPD<-mantel.test(season_precip,MPDClus_E,nperm=1000)
a_mean_mantel_E_MPD<-mantel.test(annual_mean,MPDClus_E,nperm=1000)
a_max_mantel_E_MPD<-mantel.test(annual_max,MPDClus_E,nperm=1000)
a_precip_mantel_E_MPD<-mantel.test(annual_precip,MPDClus_E,nperm=1000)
s_mean_mantel_E_MNTD<-mantel.test(season_mean,MNTDClus_E,nperm=1000)
s_max_mantel_E_MNTD<-mantel.test(season_max,MNTDClus_E,nperm=1000)
s_precip_mantel_E_MNTD<-mantel.test(season_precip,MNTDClus_E,nperm=1000)
a_mean_mantel_E_MNTD<-mantel.test(annual_mean,MNTDClus_E,nperm=1000)
a_max_mantel_E_MNTD<-mantel.test(annual_max,MNTDClus_E,nperm=1000)
a_precip_mantel_E_MNTD<-mantel.test(annual_precip,MNTDClus_E,nperm=1000)


##########################################################
#Faith's PD of observation week per pool and year by zone#
##########################################################
week_clus<-FULL_COMMUNITY_csv[133:348,]
names<-FULL_COMMUNITY_csv[133:348,1]
rownames(week_clus)<-names
Zone_list <- c("U","E","B")
cluster_zone.list=list()
cluster.zone.name_list=list()
for (i in Zone_list){
    name <- paste("cluster.groups.",i,sep="")
    cluster_zone.list[[name]]<-week_clus[(week_clus$Zone==i)&!(is.na(week_clus$Week))&!(is.na(week_clus$Year)),]
}
for (q in seq (1,3)){
        name.name <- paste("cluster.group.names.",q,sep="")
        cluster.zone.name_list[[name.name]]<-cluster_zone.list[[q]][,1:5]
        cluster_zone.list[[q]]<-as.matrix(cluster_zone.list[[q]][,9:50])
    }
ses.pd_results <- list()
PD_ZONE <- list()

for (i in Zone_list){
    group <- paste("cluster.groups.",i,sep="")
    ses <- ses.pd(week_clus[week_clus$Zone==i,9:50],SynthTree,null.model="richness",runs=1000,iterations=1000,include.root=TRUE)
    ses.pd_results[[group]]<-ses
}
ses$Clusters <- week_clus[week_clus$Zone==i,1]
for (q in seq (1,3)){
    name_cluster <- paste("cluster.group.names.",q,sep="")
    ses.info<- cbind(Clusters=cluster.zone.name_list[[name_cluster]][,"Clusters"],ses.pd_results[[q]])
    ses.info<- cbind(zone=cluster.zone.name_list[[name_cluster]][,"Zone"],ses.info)
    ses.info<- cbind(year=cluster.zone.name_list[[name_cluster]][,"Year"],ses.info)
    ses.info<- cbind(week=cluster.zone.name_list[[name_cluster]][,"Week"],ses.info)
    PD_ZONE[[name_cluster]]<- ses.info[ses.info$ntaxa>=1,]
    }
PD_shuff_LIST=list()
for (i in seq (1,3)){
    name <- paste("PD_ZONE_",i,sep="")
    PD_shuff_LIST[[name]]<- PD_ZONE[[i]][sample(nrow(PD_ZONE[[i]])),]
}

#Find best fitting polynomial function

degree = 3
K = 10
k_list=list()
regress_list=list()
Breusch_Pagan_list=list()
mse=matrix(data=NA,nrow=K,ncol=degree)
plot_list=list()
for (x in PD_shuff_LIST){
    for(i in 1:K){
    folds <- cut(seq(1,nrow(x)),breaks=K,labels=FALSE)
    index <- which(folds==i,arr.ind=TRUE)
    test <- x[index,]
    train <- x[-index,]

    for(j in 1:degree){
        fittrain = lm(pd.obs.z ~ poly(Week,j),data = train)
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
            regress <- lm(pd.obs.z~poly(Week,k,raw=TRUE),data = y)
            #Perform Breusch-Pagan test to determine if model residuals exhibit heteroscedasticity
            bpmodel <- bptest(regress)
            #Plot the Regression
            plot_list[[name]]<-ggplot(y,aes(x=Week,y=pd.obs.z))+
            geom_point()+
            stat_smooth(method="lm", formula = y ~ poly(x,k,raw=TRUE),linewidth=1)+theme_light()+
            xlab('Week')+
            ylab('Phylogenetic Diversity')+scale_x_continuous(breaks = round(seq(min(y$Week), max(y$Week), by = 1),1))
            name=paste(d,"Polynomial_regression.pdf",sep="_")
            ggsave(
                filename=name,
                plot=last_plot(),
                device=pdf,
                dpi=300,
                limitsize=TRUE
            )
            #Save Regressions as pdf
            }
        regress_list[[d]]<-regress
        Breusch_Pagan_list[[d]]<- bpmodel
    }
}

name=paste(d,"Polynomial_regression.pdf",sep="_")
            ggsave(
                filename=name,
                plot=last_plot(),
                device=pdf,
                dpi=300,
                limitsize=TRUE
#Get Summaries of Each Regression
for (i in regression_list){
    print(summary(i))
}

################################
#Clustering Modoc and Merced####
################################
gosejohan_plus_mercedspecies <- read_csv("~/Documents/Github_Projects/Vernal_Pool_CommunityDiversity/Data/csv_files/gosejohan_plus_mercedspecies.csv")
SynthPhylo_GM<-V.PhyloMaker2::phylo.maker(gosejohan_plus_mercedspecies,tree=GBOTB.extended.TPL,nodes=nodes.info.1.TPL,output.sp.list = TRUE,output.tree = TRUE,scenarios = "S3")
SynthTree_GM<-SynthPhylo_GM$scenario.3
phydist_GM<-cophenetic(SynthTree_GM)
clusterZone_Pool_GM<-as.matrix(FULL_COMMUNITY_csv[c(10,11,12,16,17,18,19,20,21,22,23,24),9:])
names_GM<-FULL_COMMUNITY_csv[c(10,11,12,16,17,18,19,20,21,22,23,24),1]
rownames(clusterZone_Pool_GM)<-names_GM$Clusters
beta_zone.year_GM <- comdist(clusterZone_Pool_GM,phydist_GM)
beta_zone.year.clusters_GM <- hclust(beta_zone.year_GM)

cluster_Pool_GM<-as.matrix(FULL_COMMUNITY_csv[c(4:12),9:])
names_GM<-FULL_COMMUNITY_csv[c(4:12),1]
rownames(cluster_Pool_GM)<-names_GM$Clusters
mpd.result_GM <- ses.mpd(cluster_Pool_GM,phydist_GM,null.model="taxa.labels",abundance.weighted=FALSE,runs=1000)
mntd.result_GM <- ses.mntd(cluster_Pool_GM,phydist_GM,null.model="taxa.labels",abundance.weighted=FALSE,runs=1000)
pd.result_GM <- ses.pd(cluster_Pool_GM,SynthTree,null.model="richness",runs=1000)
