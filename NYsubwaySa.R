library(RgoogleMaps)
library(dplyr)
library(vegan)
library(ade4)
library(e1071)
#library(biomod2)
library(RColorBrewer)
NYCdata <- (read.csv("~/Dropbox/DATA/Data_for_IBS574_metagenome_practical/DataTable5-metaphlan-metadata_v19.csv",stringsAsFactors = FALSE, header = TRUE))[1:4]
strain_SRA <- read.table("~/Dropbox/ARTICLES_BY_TDR/2015-staph-metagenome/STAPH/runs-to-samples.txt", header = TRUE)
colnames(strain_SRA) <- c("Run", "Sample.ID") #to make join easier
NYCdata_SRA <- left_join(NYCdata, strain_SRA, by = "Sample.ID")
Staph_betas <- read.csv("~/Dropbox/ARTICLES_BY_TDR/2015-staph-metagenome/STAPH/Final_Beta_NYC_STAPH.csv",stringsAsFactors = FALSE, header = TRUE)

locations_with_runs <- filter(NYCdata_SRA, !is.na(Run) ) %>%
  select(Latitude, Logitude)

lats <- as.numeric(as.character(locations_with_runs$Latitude))
lons <- as.numeric(as.character(locations_with_runs$Logitude))
cols_staph <- rep(NULL,length(lats))
#all the statiitons which had S. aureus
cols_staph[which(NYCdata_SRA$Run %in% Staph_betas$Run)] <- "red"

#stations reporting CC8_8
CC_8_8_runs <- filter(Staph_betas, CC_8_8 > 0) %>%
  select(Run)
cols_CC_8_8 <- rep(NULL,length(lats))
cols_CC_8_8[which(NYCdata_SRA$Run %in% CC_8_8_runs$Run)] <- "green"

#stations reporting CC_239_239
CC_239_239_runs <- filter(Staph_betas, CC_239_239 > 0) %>%
  select(Run)
cols_CC_239_239 <- rep(NULL,length(lats))
cols_CC_239_239[which(NYCdata_SRA$Run %in% CC_239_239_runs$Run)] <- "green"

#all stations
alllats <- as.numeric(as.character(NYCdata$Latitude)) #every non-numeric value is converted to NA
alllons <- as.numeric(as.character(NYCdata$Logitude))

#centered on NYC
gmap1 <- GetMap(center = c(lat = 40.7127, lon = -74.0059), size = c(640, 640), zoom = 11, GRAYSCALE = TRUE)
#centered on Queens
gmap2 <- GetMap(center = c(lat = 40.7500, lon = -73.8667), size = c(640, 640), zoom = 11, GRAYSCALE = TRUE)
gmap3 <- GetMap(center = c(lat = 40.7500, lon = -73.8667), size = c(640, 640), zoom = 10, GRAYSCALE = TRUE)
## function for all types

plot_CC_types <-function(CC, CCcol = "green") {
  CC_df <- filter(Staph_betas, CC > 0) %>%
    select(Run)
  cols <- rep(NULL,length(lats))
  cols[which(NYCdata_SRA$Run %in% CC_df$Run)] <- CCcol
  plotname <- paste("~/Dropbox/ARTICLES_BY_TDR/2015-staph-metagenome/",CC,"z11.png",collapse = "_")
  print(plotname)
}


##PLOTTING

png("~/Dropbox/ARTICLES_BY_TDR/2015-staph-metagenome/all_stations_z11.png",width=640, height =640, res = 75)
PlotOnStaticMap(gmap2, lat = lats , lon = lons, cex=1.5,pch=20, col = "blue")
dev.off()

png("~/Dropbox/ARTICLES_BY_TDR/2015-staph-metagenome/all_stations_z10.png",width=640, height =640, res = 75)
PlotOnStaticMap(gmap3, lat = lats , lon = lons, cex=1.5,pch=20, col = "blue")
dev.off()

png("~/Dropbox/ARTICLES_BY_TDR/2015-staph-metagenome/all_staph_z11.png",width=640, height =640, res = 75)
PlotOnStaticMap(gmap2, lat = lats , lon = lons, cex=1.5,pch=20, col = cols_staph)
dev.off()

png("~/Dropbox/ARTICLES_BY_TDR/2015-staph-metagenome/all_staph_z10.png",width=640, height =640, res = 75)
PlotOnStaticMap(gmap3, lat = lats , lon = lons, cex=1.5,pch=20, col = cols_staph)
dev.off()

#function to automate plotting
plot_CC_types <-function(CC, CCcol = "red") {
  crows <- which(Staph_betas[[CC]] > 0)
  cruns <-Staph_betas[crows,1]
  
  CC_df <- filter(NYCdata_SRA, Run %in% cruns) %>% select(Latitude,Logitude)
  
  plotname <- paste("~/Dropbox/ARTICLES_BY_TDR/2015-staph-metagenome/",CC,"_z11.png",sep= "", collapse = "")
  png(plotname,width=640, height =640, res = 75)
  PlotOnStaticMap(gmap2, lat = as.numeric(CC_df$Latitude) , lon = as.numeric(CC_df$Logitude), cex=1.5,pch=20, col = CCcol)
  dev.off()
  
  plotname <- paste("~/Dropbox/ARTICLES_BY_TDR/2015-staph-metagenome/",CC,"_z10.png",sep= "", collapse = "")
  png(plotname,width=640, height =640, res = 75)
  PlotOnStaticMap(gmap3, lat = as.numeric(CC_df$Latitude) , lon = as.numeric(CC_df$Logitude), cex=1.5,pch=20, col = CCcol)
  dev.off()
}

#plot all CC types
for (i in colnames(Staph_betas)){
  if(i != "Run" && as.numeric(sum(Staph_betas[[i]])) > 0){
    plot_CC_types(i)
  }  
}

##Mantel test for spatial autocorrelation of Staph aureus genotypes
##First get Bray curtis matrix of dissimilarities
braymat <- vegdist(as.matrix(Staph_betas[-1,-1]))
jacmat <- vegdist(as.matrix(Staph_betas[-1,-1]), method = "jaccard")
##now get geographical distance
geogdist_staph <- filter(NYCdata_SRA, Run %in% Staph_betas$Run) %>% select(Latitude,Logitude) %>% dist
set.seed(45678)
mantel.rtest(geogdist_staph,braymat,nrepet = 9999)
#not significant p = 0.6911
mantel.rtest(geogdist_staph,jacmat,nrepet = 9999)
#not significant p = 0.7129

bintr <- function(mat,cutoff){
  mat[which(mat > cutoff)] <- 1
  mat[which(!(mat > cutoff))] <- 0
  return(mat)
}
#convert staph betas into binary with a cutoff of 0.1
binmat <- bintr(as.matrix(Staph_betas[-1,-1]),0.2)

#calculate hamming distances as an alternative using e1071 function
hammmingmat <- hamming.distance(binmat)
mantel.rtest(geogdist_staph,as.dist(hammmingmat),nrepet = 9999)
#Simulated p-value: 0.2591





genotypes_plot <- function(mat,tit) {
  top_genos <- c("CC_30_36","CC_8_8","CC_45","CC_398_398","CC_133_133","CC_30_30","CC_59_59","CC_15","CC_97","CC_5_5","CC_8_72")
  cS <- colSums(mat)
  cStop <- cS[top_genos]
  cStop[["others"]] <- sum(cS) - sum(cStop)  
  barplot(cStop, main = tit, las = 3, cex.names = 0.8, col = (brewer.pal(12,name = 'Set3')))
}
png("~/Dropbox/ARTICLES_BY_TDR/2015-staph-metagenome/NYC_barchart.png",width=640, height =640, res = 75)
genotypes_plot(Staph_betas[,-1],"NYC All samples")
dev.off()