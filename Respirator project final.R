setwd("E:/OneDrive/BH_Data/Research/CLP/FB2/Analysis/2020")

library(geomorph)
library(shapes)
source("Rfunctions1.txt")
library(ggplot2)
library(Morpho)
library(corrplot)
library(plyr)
library(dplyr)
library(rrcov)
library(Evomorph)
library(kimisc)
library(animation)
#library(exrafont)
#library(Stack)
library(gridExtra)
library(png)
library(raster)
library(sp)
library(rgdal)
library(magick)
library(ecodist)
library(vegan)
library(Momocs)
library(phytools)
library(stringr)
library(grid)
library(ggpubr)
library(ggimage)
library(cowplot)
library(reshape)
library(reshape2)
library(ggrepel)
library(comprehenr)
library(psych)
library(reshape2)
library(pls)
#library(Rpdb)
library(misc3d)
library(plot3D)
library(rgl)
library(forcats)
library(grDevices)
library(ggExtra)
library(svglite)


load("RGL_posing_transforms.rdata")

# Functions

resize.win <- function(Width=6, Height=6)
{
  # works for windows
  dev.off(); # dev.new(width=6, height=6)
  windows(record=TRUE, width=Width, height=Height)
}



#Define multiple plot function (from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/)



#Define datasets



# Start with unstandardized data:




v<-read.csv("Respirator Data.csv", header=T)
v_data<-v[,17:212]


left<-c(40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65)
right<-c(14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39)
pairedLM<-cbind(left,right)

lms<-65
coords<-3


v_class<-v[,1:17]
v_data<-v[,18:212]
v_cov<-cbind(v[,3], v[,c(10,11,12,17)])
colnames(v_cov)<-c("Age","HT","WT","HC","CS")
#v_cov <- read.csv("imputed covariates.csv", header = T)


# Visualization



#Set visualization

mean3D<-mshape(dat_pr3d)
face<-file2mesh("Atlases/atlas.ply")
shade3d(face, color="gray", alpha=1)
landtemplate<-read.table("Atlases/reference65.txt", header=FALSE)
reference<-as.matrix(landtemplate)
plot3d(reference,aspect="iso", type="s",size=1.4,col="red", add=T)

#ref2 <- coords(reference[1,1], reference[1,2], reference[1,3])


# Create landmark labels

places <- seq(from = 1, to = 193, by = 3)
lm_names <- colnames(data_raw)[27:221]
lms <- lm_names[places]
lms <- str_remove(lms, "[_x]")
lms <- str_remove(lms, "[x]")

colnames(v_data) <- lm_names


rgl.close()
mean3D<-mshape(dat_pr3d)
face<-file2mesh("Atlases/atlas.ply")
shade3d(face, color="gray", alpha=1)
landtemplate<-read.table("Atlases/reference65.txt", header=FALSE)
reference<-as.matrix(landtemplate)
plot3d(reference,aspect="iso", type="s",size=1.4,col="red", add=T)

reference <- as.data.frame(reference)

i <- 1
for(i in 1:67)
  {
  #text3d(x = reference2[i,1], y = reference2[i,2], z = reference2[i,3],  cex = 1.2, lms2$lms[i], adj = c(0.5), pos = 2, offset = 1, usePlotmath = is.language(lms2$lms[i])) 
  text3d(x = reference[i,1], y = reference[i,2], z = reference[i,3],  cex = 1.2, lms[i], adj = c(0.5), pos = 2, offset = 1, usePlotmath = is.language(lms[i])) 
  
}

lms_n <- length(lms)
mean_face_lm <- mshape(arrayspecs(v_data, lms_n, coords)) 
mean_face_lm <- as.matrix(mean_face_lm)

# Calculate measurements

measurements <- read.csv("measurements.csv", header = T)
measurements2 <- read.csv("measurements2.csv", header = T)


# Create null vector for measurements

v_data_meas <- data.frame(array(NA, dim = c(nrow(v_data), 0)))
i<-1
for(i in 1:nrow(measurements))
  
{
  meas_x <- measurements$Measurement[i]
  origin_x <- measurements$Origin[i]
  End_x <- measurements$End[i]
  
  o_x_i <- v_data[,paste0(origin_x, "_x")]
  o_y_i <- v_data[,paste0(origin_x, "_y")]
  o_z_i <- v_data[,paste0(origin_x, "_z")]
  
  e_x_i <- v_data[,paste0(End_x, "_x")]
  e_y_i <- v_data[, paste0(End_x, "_y")]
  e_z_i <- v_data[,paste0(End_x, "_z")]
  
  dist_x <- ((((o_x_i-e_x_i)^2)+((o_y_i-e_y_i)^2)+((o_z_i-e_z_i)^2))^(1/2)) * v_cov$CS
  
  
  
}


# Using raw data


# Add pupil landmark
R_P_x <- (data_raw$R03_x + data_raw$R04_x)/2
R_P_y <- (data_raw$R03_y + data_raw$R04_y)/2
R_P_z <- (data_raw$R03_z + data_raw$R04_z)/2

L_P_x <- (data_raw$L03_x + data_raw$L04_x)/2
L_P_y <- (data_raw$L03_y + data_raw$L04_y)/2
L_P_z <- (data_raw$L03_z + data_raw$L04_z)/2

data_raw <- cbind(data_raw, R_P_x, R_P_y, R_P_z, L_P_x, L_P_y, L_P_z)

v_data_meas <- data.frame(array(NA, dim = c(nrow(data_raw), 0)))


raw_lms <- data_raw[, 27:227]
raw_met <- data_raw[,1:26]


i <- 1
for(i in 1:nrow(measurements2))
  
{
  meas_x <- measurements2$Measurement[i]
  origin_x <- measurements2$Origin[i]
  End_x <- measurements2$End[i]
  
  o_x_i <- raw_lms[,paste0(origin_x, "_x")]
  o_y_i <- raw_lms[,paste0(origin_x, "_y")]
  o_z_i <- raw_lms[,paste0(origin_x, "_z")]
  
  e_x_i <- raw_lms[,paste0(End_x, "_x")]
  e_y_i <- raw_lms[, paste0(End_x, "_y")]
  e_z_i <- raw_lms[,paste0(End_x, "_z")]
  
  dist_x <- ((((o_x_i-e_x_i)^2)+((o_y_i-e_y_i)^2)+((o_z_i-e_z_i)^2))^(1/2))
  dist_x <- as.data.frame(dist_x)
  colnames(dist_x) <- meas_x
  
  v_data_meas <- cbind(v_data_meas, dist_x)
  
}

#pupil_dist <- (v_data_meas$Pupil_dist_1 + v_data_meas$Pupil_dist_2)/2
#v_data_meas <- cbind(v_data_meas, pupil_dist)

# Measurement data for analysis

meas_data <- cbind(raw_met, v_data_meas)

# remove outliers

meas_data <- meas_data[which(meas_data$Max_front_breadth > 50),]
meas_data <- meas_data[which(meas_data$Pupil_dist > 30),]

# Make figures showing measurements in 3D


# Set views

rgl.close()
par3d("windowRect"= c(0,0,1000,700))
mean3D<-mshape(dat_pr3d)
face<-file2mesh("Atlases/atlas.ply")
shade3d(face, color="gray", alpha=0.3)

# Add pupil landmark

R_3 <- reference[which(lms == "R03"),]
R_4 <- reference[which(lms == "R04"),]

L_3 <- reference[which(lms == "L03"),]
L_4 <- reference[which(lms == "L04"),]


R_P <- colMeans(rbind(R_3, R_4))
L_P <- colMeans(rbind(L_3, L_4))

reference2 <- rbind(reference, R_P, L_P)

add_P <- c("11", "Pupil_dist", "R_P", "L_P")
measurements2 <- rbind(measurements, add_P)
measurements2 <- measurements2[c(1:4,6:12),]
lms <- as.data.frame(lms)
lms2 <- rbind(lms, "R_P", "L_P")
reference2 <- as.data.frame(reference2)
row.names(reference2) <- lms2$lms

face<-file2mesh("Atlases/atlas.ply")

#model<-par3d()$modelMatrix 
#front<-par3d()$userMatrix
#side<-par3d()$userMatrix
#T3_4_view<-par3d()$userMatrix
#save(front, side, file = "RGL_posing_transforms.rdata")

v_list<-c("front","side", "T3_4_view")


create_views<-function(v_list)
  
{
  n_views<-length(v_list)
  views <- array(NA, dim = c(4,4,n_views))
  
  for(i in 1:n_views)
  {
    views[,,i]<-get(v_list[i])
  }
  
}

create_views(v_list)


# Create directory for measurement images
dir.create("measurement_images")  


i <- 1
z <- 1

for(i in 1:nrow(measurements2))
{
for(z in 1:dim(views)[3])
{


open3d(zoom=0.75, userMatrix = views[,,z])
par3d("windowRect"= c(0,0,1000,700))
shade3d(face, color="grey", alpha = 0.3)
  

measure_name <- measurements2$Measurement[i]
origin_lm_name <- measurements2$Origin[i]
end_lm_name <- measurements2$End[i]
origin_lm_no <- which(lms2 == origin_lm_name)

origin_l <- reference2[origin_lm_name,]
end_l <- reference2[end_lm_name,]

x0 <- origin_l$V1
y0 <- origin_l$V2
z0 <- origin_l$V2

x1 <- end_l$V1
y1 <- end_l$V2
z1 <- end_l$V2

lines3d(rbind(origin_l, end_l), col="blue", lwd=8)

rgl.snapshot(paste0("measurement_images/", measure_name, "_", v_list[z],".png"), top = TRUE )
rgl.close()

}}
        
# Check landmarks

rgl.close()
mean3D<-mean_sym
face<-file2mesh("Atlases/atlas.ply")
shade3d(face, color="gray", alpha=1)
plot3d(reference2,aspect="iso", type="s",size=1.4,col="red", add=T)

i <- 1
for(i in 1:67)
{
  text3d(x = reference2[i,1], y = reference2[i,2], z = reference2[i,3],  cex = 1.2, lms2$lms[i], adj = c(0.5), pos = 2, offset = 1, usePlotmath = is.language(lms2$lms[i])) 
  
}

# Assemble into a panel

#Create page withPCs1-8

#measurements2 <- measurements2[-which(measurements2$Measurement=="Pupil_dist_2"),]

i<-1
z <- 1

image_stack <- c()

for(i in 1:nrow(measurements2))
  {
  meas_i <- measurements2$Measurement[i]
  meas_label <- measurements2$full_meas_name[i]
  view_z<- v_list[z]
  image_row <- c()
  
    for(z in 1:dim(views)[3])
      {
        
    temp_i <- image_read(paste0("measurement_images/", meas_i, "_", v_list[z],".png"))
    temp_i <- image_crop(temp_i, geometry = "800x0+0")
    assign(paste0(meas_i, "_", view_z), temp_i)
    image_row <- append(image_row, (assign(paste0(meas_i, "_", view_z), temp_i)))
    reg_row_i<-image_append(image_scale(image_row, "x200"))
      }
  reg_row_i<-image_border(reg_row_i, "white", "0x20")
  reg_row_i<-image_annotate(reg_row_i, paste0(meas_label), font = "sans", size = 30)
  image_stack <- append(image_stack, (assign(paste0(meas_i, "_row"), reg_row_i)))
  meas_fig<-image_append(image_scale(image_stack), stack=TRUE)
  
  }

  
## Check image
image_browse(meas_fig)
image_write(meas_fig, path = paste0("measurement_images/measurement figure.png"), format = "png")



# Put all measurements on same image
i <- 1
z <- 1

for(i in 1:nrow(measurements2))
{
  for(z in 1:dim(views)[3])
  {
    
    
    open3d(zoom=0.75, userMatrix = views[,,z])
    par3d("windowRect"= c(0,0,1000,700))
    shade3d(face, color="grey", alpha = 0.3)
    
    
    measure_name <- measurements2$Measurement[i]
    origin_lm_name <- measurements2$Origin[i]
    end_lm_name <- measurements2$End[i]
    origin_lm_no <- which(lms2 == origin_lm_name)
    
    origin_l <- reference2[origin_lm_name,]
    end_l <- reference2[end_lm_name,]
    
    x0 <- origin_l$V1
    y0 <- origin_l$V2
    z0 <- origin_l$V2
    
    x1 <- end_l$V1
    y1 <- end_l$V2
    z1 <- end_l$V2
    
    lines3d(rbind(origin_l, end_l), col="blue", lwd=8)
    
    rgl.snapshot(paste0("measurement_images/", measure_name, "_", v_list[z],".png"), top = TRUE )
    rgl.close()
    
  }}

# Check landmarks

rgl.close()
mean3D<-mean_sym
face<-file2mesh("Atlases/atlas.ply")
shade3d(face, color="gray", alpha=1)
plot3d(reference2,aspect="iso", type="s",size=1.4,col="red", add=T)

i <- 1
for(i in 1:67)
{
  text3d(x = reference2[i,1], y = reference2[i,2], z = reference2[i,3],  cex = 1.2, lms2$lms[i], adj = c(0.5), pos = 2, offset = 1, usePlotmath = is.language(lms2$lms[i])) 
  
}


# Creaet a page of analyses for each measurements

# Curves for all measurements for all data

meas_data <- read.csv("meas_data.csv", header = TRUE)
#meas_data$Syndrome_Class[which(meas_data$Syndrome_Class == "")] <- "Syndrome"

meas_data <- meas_data[-which(meas_data$Sex=="U"),]
meas_data <- meas_data[which(meas_data$Age > 0),]
meas_data <- meas_data[-which(meas_data$Race == "Native Hawiian/Pacific Islander"),]
meas_data <- meas_data[-which(meas_data$Race == "American Indian/Alaska Native"),]
meas_data$Race[which(meas_data$Race=="unknown/not reported")] <- "Unknown/Not reported"
meas_data$Race[which(meas_data$Race=="Unknown/not reported")] <- "unknown"


meas_met <- meas_data[, c("Race", "Race2", "Syndrome_Class", "Age", "Sex", "Ethnicity", "HT", "WT")]
#meas_met$Race[which(meas_met$Race=="unknown/not reported")] <- "Unknown/Not reported"
meas_data_meas <- meas_data[, measurements2$Measurement]

meas_data_export <- cbind(meas_met, meas_data_meas)
write.csv(meas_data_export, "Respirator linear distances.csv")

# Create corrplot

meas_cov <- cov(meas_data_meas)
meas_cor <- cor(meas_data_meas)

colnames(meas_cor) <- measurements2$full_meas_name
rownames(meas_cor) <- measurements2$full_meas_name
meas_cor <- as.matrix(meas_cor)
meas_data$Race2[which(meas_data$Race2 == "unknown/not reported")] <- "unknown"

dev.new()
cplot <-corrplot.mixed(meas_cor, lower = "circle", upper = "number", upper.col = "black", number.cex = .7, order = "hclust")
cplot <-corrplot(meas_cor, type = "lower", order = "hclust")

dir.create("respirator_plots")

race_counts<- meas_met %>% count(meas_met$Race)
colours = c("#313695", "#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8", "#FFFFBF", "#FEE090", "#FDAE61", "#F46D43", "#D73027", "#A50026") 





i <- 1

for(i in 1:nrow(measurements2))
{
measurement <- meas_data[,measurements2$Measurement[i]]
meas_dat_i <- data.frame(meas_met, measurement)
measurement_name <- measurements2$full_meas_name[i]

meas_kids_dat <- meas_dat_i[which(meas_dat_i$Age<19),]

# Plots by sex
# adults
dev.off()
dev.new()
plot<-ggplot(meas_dat_i, aes(Age, measurement, colour = factor(Sex)))
plot<-plot+geom_point(alpha=0.3, size=2)+geom_smooth(method='loess')
plot<-plot + scale_color_manual(values=c("#46b4af", "#b44682"),labels=c("Female","Male"))
plot<-plot+theme_bw()
#plot<-plot + geom_density2d(aes(colour = factor(Sex)), size=0.15)
#plot<-plot + stat_density2d(aes(fill = ..level..), geom="polygon", alpha=0.3)
plot<-plot+guides(colour=guide_legend(title="Sex"))
plot<-plot+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot<-plot+ylab(measurement_name)
plot<-plot+labs(title = paste0("A) ", measurement_name, " Sex All Ages"))#+theme(plot.title = element_text(hjust = -0.2, vjust=2.12))
plot_age_sex_all<-plot	
assign(paste0(measurement_name, "_sex_all_ages"), plot_age_sex_all)
#plot_age_sex_all
ggsave(paste0("respirator_plots/", paste0(measurement_name, "_sex_all_ages"), ".pdf"))

# kids

#dev.new()
plot<-ggplot(meas_kids_dat, aes(Age, measurement, colour = factor(Sex)))
plot<-plot+geom_point(alpha=0.3, size=2)+geom_smooth(method='loess')
plot<-plot + scale_color_manual(values=c("#46b4af", "#b44682"),labels=c("Female","Male"))
plot<-plot+theme_bw()
plot<-plot+guides(colour=guide_legend(title="Sex"))
plot<-plot+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot<-plot+ylab(measurement_name)
plot<-plot+labs(title = paste0("B) ", measurement_name, " Sex 18 and younger"))#+theme(plot.title = element_text(hjust = -0.2, vjust=2.12))
plot_age_sex_kids<-plot	
assign(paste0(measurement_name, "_sex_kids"), plot_age_sex_kids)

ggsave(paste0("respirator_plots/", paste0(measurement_name, "_sex_kids"), ".pdf"))


# Plots by race

# adults
#dev.new()
plot<-ggplot(meas_dat_i, aes(Age, measurement, colour = factor(Race)))
plot<-plot+geom_point(alpha=0.2, size=2)+geom_smooth(method='loess')
#plot<-plot + scale_color_manual(values=colours[1:length(unique(Race))])
plot<-plot + scale_color_manual(values=c("blue", "lightblue", "orchid", "slateblue", "aquamarine", "slategrey"))
plot<-plot+theme_bw()
#plot<-plot + geom_density2d(aes(colour = factor(Sex)), size=0.15)
#plot<-plot + stat_density2d(aes(fill = ..level..), geom="polygon", alpha=0.3)
plot<-plot+guides(colour=guide_legend(title="Race")) 
plot<-plot+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot<-plot+ylab(measurement_name)
plot<-plot+labs(title = paste0("C) ", measurement_name, " Race All Ages"))#+theme(plot.title = element_text(hjust = -0.3, vjust=2.12))
plot_age_race_all<-plot	
#plot_age_race_all
assign(paste0(measurement_name, "_race_all_ages"), plot_age_race_all)

ggsave(paste0("respirator_plots/", paste0(measurement_name, "_race_all_ages"), ".pdf"))

# kids
#dev.new()
plot<-ggplot(meas_kids_dat, aes(Age, measurement, colour = factor(Race)))
plot<-plot+geom_point(alpha=0.2, size=2)+geom_smooth(method='loess')
#plot<-plot + scale_color_manual(values=colours[1:length(unique(Race))])
plot<-plot + scale_color_manual(values=c("blue", "lightblue", "orchid", "slateblue", "aquamarine", "slategrey"))
plot<-plot+theme_bw()
#plot<-plot + geom_density2d(aes(colour = factor(Sex)), size=0.15)
#plot<-plot + stat_density2d(aes(fill = ..level..), geom="polygon", alpha=0.3)
plot<-plot+guides(colour=guide_legend(title="Race"))
plot<-plot+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot<-plot+ylab(measurement_name)
plot<-plot+labs(title = paste0("D) ", measurement_name, " Race kids"))#+theme(plot.title = element_text(hjust = -0.3, vjust=2.12))
plot_age_race_kids<-plot	
#plot_age_race_kids
assign(paste0(measurement_name, "_race_kids"), plot_age_race_kids)

ggsave(paste0("respirator_plots/", paste0(measurement_name, "_race_all_kids"), ".pdf"))

#Race2 plots

# adults
#dev.new()
plot<-ggplot(meas_dat_i, aes(Age, measurement, colour = factor(Race2)))
plot<-plot+geom_point(alpha=0.1, size=2)+geom_smooth(method='loess')
#plot<-plot + scale_color_manual(values=colours[1:length(unique(Race))])
plot<-plot + scale_color_manual(values=c("blue", "lightblue", "orchid", "slateblue", "aquamarine", "slategrey", "gray"))
plot<-plot+theme_bw()
#plot<-plot + geom_density2d(aes(colour = factor(Sex)), size=0.15)
#plot<-plot + stat_density2d(aes(fill = ..level..), geom="polygon", alpha=0.3)
plot<-plot+guides(colour=guide_legend(title="Superpopulation"))
plot<-plot+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot<-plot+ylab(measurement_name)
plot<-plot+labs(title = paste0("E) ", measurement_name, " 10k Genomes Superpopulation All Ages"))#+theme(plot.title = element_text(hjust = -0.3, vjust=2.12))
plot_age_race2_all<-plot	
#plot_age_race2_all
assign(paste0(measurement_name, "_race2_all_ages"), plot_age_race2_all)

ggsave(paste0("respirator_plots/", paste0(measurement_name, "_race2_all_ages"), ".pdf"))

# kids
#dev.new()
plot<-ggplot(meas_kids_dat, aes(Age, measurement, colour = factor(Race2)))
plot<-plot+geom_point(alpha=0.2, size=2)+geom_smooth(method='loess')
#plot<-plot + scale_color_manual(values=colours[1:length(unique(Race))])
plot<-plot + scale_color_manual(values=c("blue", "lightblue", "orchid", "slateblue", "aquamarine", "slategrey", "grey"))
plot<-plot+theme_bw()
#plot<-plot + geom_density2d(aes(colour = factor(Sex)), size=0.15)
#plot<-plot + stat_density2d(aes(fill = ..level..), geom="polygon", alpha=0.3)
plot<-plot+guides(colour=guide_legend(title="Superpopulation"))
plot<-plot+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot<-plot+ylab(measurement_name)
plot<-plot+labs(title = paste0("F) ", measurement_name, " 10k Genomes Superpopulation kids"))#+theme(plot.title = element_text(hjust = -0.3, vjust=2.12))
plot_age_race2_kids<-plot	
#plot_age_race2_kids
assign(paste0(measurement_name, "_race2_kids"), plot_age_race2_kids)

ggsave(paste0("respirator_plots/", paste0(measurement_name, "_race2_kids"), ".pdf"))

# Syndrome Class plots

# adults
#dev.new()
plot<-ggplot(meas_dat_i, aes(Age, measurement, colour = factor(Syndrome_Class)))
plot<-plot+geom_point(alpha=0.1, size=2)+geom_smooth(method='loess')
#plot<-plot + scale_color_manual(values=colours[1:length(unique(Race))])
plot<-plot + scale_color_manual(values=c("blue", "slategrey"))
plot<-plot+theme_bw()
#plot<-plot + geom_density2d(aes(colour = factor(Sex)), size=0.15)
#plot<-plot + stat_density2d(aes(fill = ..level..), geom="polygon", alpha=0.3)
plot<-plot+guides(colour=guide_legend(title="Syndrome Class"))
plot<-plot+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot<-plot+ylab(measurement_name)
plot<-plot+labs(title = paste0("G) ", measurement_name, " Syndromic vs Nonsyndr. All Ages"))#+theme(plot.title = element_text(hjust = -0.3, vjust=2.12))
plot_age_syndrome_all<-plot	
#plot_age_race2_all
assign(paste0(measurement_name, "_syndrome_all_ages"), plot_age_syndrome_all)

ggsave(paste0("respirator_plots/", paste0(measurement_name, "_race2_all_ages"), ".pdf"))

# kids
#dev.new()
plot<-ggplot(meas_kids_dat, aes(Age, measurement, colour = factor(Syndrome_Class)))
plot<-plot+geom_point(alpha=0.2, size=2)+geom_smooth(method='loess')
#plot<-plot + scale_color_manual(values=colours[1:length(unique(Race))])
plot<-plot + scale_color_manual(values=c("blue", "slategrey"))
plot<-plot+theme_bw()
#plot<-plot + geom_density2d(aes(colour = factor(Sex)), size=0.15)
#plot<-plot + stat_density2d(aes(fill = ..level..), geom="polygon", alpha=0.3)
plot<-plot+guides(colour=guide_legend(title="Syndrome Class"))
plot<-plot+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot<-plot+ylab(measurement_name)
plot<-plot+labs(title = paste0("H) ", measurement_name, " Syndromic vs Nonsynd. kids"))#+theme(plot.title = element_text(hjust = -0.3, vjust=2.12))
plot_age_syndrome_kids<-plot	
#plot_age_syndrome_kids
assign(paste0(measurement_name, "_syndrome_kids"), plot_age_syndrome_kids)

ggsave(paste0("respirator_plots/", paste0(measurement_name, "_race2_kids"), ".pdf"))

dev.new()

plot1 <-  ggarrange(plot_age_sex_all, plot_age_sex_kids, ncol = 2)
plot2 <-  ggarrange(plot_age_race_all, plot_age_race_kids, ncol = 2)
plot3 <-  ggarrange(plot_age_race2_all, plot_age_race2_kids, ncol = 2)
plot4 <-  ggarrange(plot_age_syndrome_all, plot_age_syndrome_kids, ncol = 2)

plot <- ggarrange(plot1, plot2, plot3, plot4, ncol = 1, nrow = 4)
resize.win(12.5,16)
plot
ggsave(paste0("respirator_plots/", paste0(measurement_name, "_panel"), ".pdf"))


#for (i in dev.list()[2]:dev.list()[length(dev.list())]) {
 # dev.off()}
}







# Perfrom ANOVAS for all 18 measurements


race_coefficients <- as.data.frame(array(NA, dim = c(0,11)))
colnames(race_coefficients) <- c("Measurement",  "Intercept", "Age" , "Age^2", "Age^3", "Sex (Male)",  "Black/African American" ,"More than one race", "Unknown/Not reported", "White", "Syndromic")

race2_coefficients <- as.data.frame(array(NA, dim = c(0,13)))
colnames(race2_coefficients) <- c("Measurement",  "Intercept", "Age" , "Age^2", "Age^3", "Sex (Male)", "AMR","ASA","EAS","EUR","SAS","unknown","Syndromic")



i <- 1

for(i in 1:nrow(measurements2))
{
  measurement <- meas_data[,measurements2$Measurement[i]]
  meas_dat_i <- data.frame(meas_met, measurement)
  measurement_name <- measurements2$full_meas_name[i]
  
  meas_kids_dat <- meas_dat_i[which(meas_dat_i$Age<19),]

  Age_i <- meas_dat_i$Age
  Age2_i <- meas_dat_i$Age^2
  Age3_i <- meas_dat_i$Age^3
  
  Sex_i <- meas_dat_i$Sex
  Race_i <- meas_dat_i$Race
  Race2_i <- meas_dat_i$Race2
  Synd_i <- meas_dat_i$Syndrome_Class
  
  anova_race_i <- lm(measurement ~ Age_i + Age2_i + Age3_i + Sex_i + Race_i + Synd_i)
  coef_race_i <- as.data.frame(t(anova_race_i$coefficients))
  coef_race_i <- cbind(measurement_name, coef_race_i)
  colnames(coef_race_i) <- colnames(race_coefficients)
  
  aov_race_i <- anova(anova_race_i)
  
  race_coefficients <- rbind(race_coefficients, coef_race_i)
  
  
  anova_race2_i <- lm(measurement ~ Age_i + Age2_i + Age3_i + Sex_i + Race2_i + Synd_i)
  coef_race2_i <- as.data.frame(t(anova_race2_i$coefficients))
  coef_race2_i <- as.data.frame(t(anova_race2_i$coefficients))
  coef_race2_i <- cbind(measurement_name, coef_race2_i)
  colnames(coef_race2_i) <- colnames(race2_coefficients)
  aov_race2_i <- anova(anova_race2_i)
    
  race2_coefficients <- rbind(race2_coefficients, coef_race2_i)
}

write.csv(race_coefficients, "measurement race cofficients.csv")
write.csv(race2_coefficients, "measurement race2 cofficients.csv")



## Perform ANOVAs on measurements without syndromic data

# Perfrom ANOVAS for all 18 measurements and only for children


race_coefficients_ns <- as.data.frame(array(NA, dim = c(0,15)))
colnames(race_coefficients_ns) <- c("Measurement",  "Intercept", "Age" , "Age^2", "Age^3", "Sex (Male)",  "Black/African American" ,"More than one race", "Unknown/Not reported", "White", "Age * Black/African American", "Age * More than one race", "Age * Unknown/Not reported", "Age * White", "Age * Sex (Male)")

race2_coefficients_ns <- as.data.frame(array(NA, dim = c(0,15)))
colnames(race2_coefficients_ns) <- c("Measurement",  "Intercept", "Age" , "Age^2", "Age^3", "Sex (Male)", "AMR","ASA","EAS","EUR","Age*AMR", "Age*ASA", "Age * EAS", "Age*EUR", "Age *Sex (Male)")

race_aov <- as.data.frame(array(NA, dim = c(0,7)))
colnames(race_aov) <- c("Measurement", "Factor", "Df", "Sumsq" , "meansq", "F_value", "Prob")

race2_aov <- as.data.frame(array(NA, dim = c(0,7)))
colnames(race2_aov) <- c("Measurement", "Factor", "Df", "Sumsq" , "meansq", "F_value", "Prob")

meas_data_ns <- meas_data[which(meas_data$Syndrome_Class=="Control"),]
meas_data_ns <- meas_data_ns[which(meas_data_ns$Age<19),]
meas_data_ns <- meas_data_ns[which(meas_data_ns$Mand_height> 25 & meas_data_ns$Age >4),]
meas_data_ns <- meas_data_ns[which(meas_data_ns$Max_front_breadth> 80 & meas_data_ns$Age >4),]
meas_data_ns <- meas_data_ns[which(meas_data_ns$Nose_protrusion> 8 & meas_data_ns$Age >4),]
meas_data_ns <- meas_data_ns[which(meas_data_ns$Race2 !="unknown"),]

# To calclulate variance components, use the realimpo package as this addresses colinearity

library(relaimpo)


i <- 1

for(i in 1:nrow(measurements2))
{
  measurement <- meas_data_ns[,measurements2$Measurement[i]]
  meas_dat_i <- data.frame(meas_data_ns, measurement)
  measurement_name <- measurements2$full_meas_name[i]
  measurement_name_R <- measurements2$Measurement[i]
  
  #meas_kids_dat <- meas_dat_i[which(meas_dat_i$Age<19),]
  
  Age_i <- meas_dat_i$Age
  Age2_i <- meas_dat_i$Age^2
  Age3_i <- meas_dat_i$Age^3
  
  Sex_i <- meas_dat_i$Sex
  Race_i <- meas_dat_i$Race
  Race2_i <- meas_dat_i$Race2
  Synd_i <- meas_dat_i$Syndrome_Class
  
  anova_race_i <- lm(measurement ~ Age_i + Age2_i + Age3_i + Sex_i + Race_i + Race_i * Age_i + Sex_i * Age_i)
  anova_race_ib <- lm(measurement ~ Age_i + Race_i + factor(Sex_i))
  coef_race_i <- as.data.frame(t(anova_race_i$coefficients))
  coef_race_i <- cbind(measurement_name, coef_race_i)
  colnames(coef_race_i) <- colnames(race_coefficients_ns)
  label_i <- array(measurement_name, dim=c(8,1))
  
  aov_race_i <- as.data.frame(anova(anova_race_i))
  factor_i <- row.names(aov_race_i)
  aov_race_i <- cbind(label_i, factor_i, aov_race_i)
  race_aov_b <- cbind(Fact_2, aov_race_i)
  colnames(race_aov_b)[2:8] <- colnames(race_aov)
  
  colnames(aov_race_i) <- colnames(race_aov)
  
  race_aov <- rbind(race_aov, aov_race_i) 
  
  Fact_2 <- c("Age", "Age^2", "Age^3", "Sex", "Race", "Age * Race", "Age * Sex", "Residuals")
  
  var_comp_i <- calc.relimp(anova_race_ib, type = c("lmg"), rela = TRUE)
  var_factors <- c("Race", "Age", "Sex")
  var_comp_t <- data.frame(var_factors, var_comp_i@lmg)
  colnames(var_comp_t) <- c("Factor", "Relative_Variance")
  
  # Calculate variance components with relaimpo
  
  perc_expl <- paste0("Tot. Expl.: ", round((var_comp_i@R2*100), digits = 2), "%")
  
  # create small barplot
  dev.new(width=3.5, height=2, noRStudioGD = TRUE, units = "inch")
  ms_plot_i <- ggplot(race_aov_b, aes(x = reorder(Fact_2, -meansq), y = meansq))
  ms_plot_i <- ms_plot_i + geom_col(fill = "darkblue")
  ms_plot_i <- ms_plot_i + labs(x = "Factor", y = "Mean square") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ms_plot_i<- ms_plot_i+labs(title = paste0("C) ", measurement_name, " variance explained"), hjust = -2)
  ms_plot_i <- ms_plot_i + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  ms_plot_i
  
  ggsave(paste0("respirator_plots/", paste0(measurement_name, "race_mean_square"), ".pdf"))
  assign(paste0(measurement_name_R, "race_mean_square"),  ms_plot_i)
  dev.off()
  
  # create small barplot for relative variance
  dev.new(width=2, height=2.5, noRStudioGD = TRUE, units = "inch")
  ms_plot_v <- ggplot(var_comp_t, aes(x = reorder(Factor, -Relative_Variance), y = Relative_Variance))
  ms_plot_v <- ms_plot_v + geom_col(fill = "darkblue")
  ms_plot_v <- ms_plot_v + labs(x = "Factor", y = "Relative Variance") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ms_plot_v<- ms_plot_v+labs(title = paste0("C) ", measurement_name, " variance explained"), hjust = -2)
  ms_plot_v <- ms_plot_v + annotate("text", x = 2.0, y = 1.0, label =  perc_expl, size = 3)
  ms_plot_v <- ms_plot_v + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
  ms_plot_v

  ggsave(paste0("respirator_plots/", paste0(measurement_name, "race_rel_var"), ".pdf"))
  assign(paste0(measurement_name_R, "race_rel_var"),  ms_plot_v)
  dev.off()
  
  race_coefficients_ns <- rbind(race_coefficients_ns, coef_race_i)
  
  
  anova_race2_i <- lm(measurement ~ Age_i + Age2_i + Age3_i + Sex_i + Race2_i+ Race2_i * Age_i + Sex_i * Age_i)
  anova_race2_ib <- lm(measurement ~ Age_i + Race2_i + factor(Sex_i))
  coef_race2_i <- as.data.frame(t(anova_race2_i$coefficients))
  coef_race2_i <- as.data.frame(t(anova_race2_i$coefficients))
  coef_race2_i <- cbind(measurement_name, coef_race2_i)
  colnames(coef_race2_i) <- colnames(race2_coefficients_ns)
  aov_race2_i <- anova(anova_race2_i)
  
  aov_race2_i <- as.data.frame(anova(anova_race2_i))
  factor_i <- row.names(aov_race2_i)
  aov_race2_i <- cbind(label_i, factor_i, aov_race2_i)
  colnames(aov_race2_i) <- colnames(race2_aov)
  
  race2_aov <- rbind(race2_aov, aov_race2_i)
  race2_aov_b <- cbind(Fact_2, aov_race2_i)
  colnames(race2_aov_b)[2:8] <- colnames(race2_aov)
  
  var_comp_i <- calc.relimp(anova_race2_ib, type = c("lmg"), rela = TRUE)
  var_factors <- c("Race", "Age", "Sex")
  var_comp_t <- data.frame(var_factors, var_comp_i@lmg)
  colnames(var_comp_t) <- c("Factor", "Relative_Variance")
  
  # create small barplot
  dev.new(width=3.5, height=2, noRStudioGD = TRUE, units = "inch")
  ms_plot_i <- ggplot(aov_race2_i, aes(x = reorder(Fact_2, -meansq), y = meansq))
  ms_plot_i <- ms_plot_i + geom_col(fill = "darkblue")
  ms_plot_i <- ms_plot_i + labs(x = "Factor", y = "Mean square") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ms_plot_i <- ms_plot_i + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  ms_plot_i
  
  assign(paste0(measurement_name_R, "race2_mean_square"),  ms_plot_i)

  ggsave(paste0("respirator_plots/", paste0(measurement_name, "race2_mean_square"), ".pdf"))
  dev.off()
  
  perc_expl <- paste0("Tot. Expl.: ", round((var_comp_i@R2*100), digits = 2), "%")
  
  # create small barplot for relative variance
  dev.new(width=2, height=2.5, noRStudioGD = TRUE, units = "inch")
  ms_plot_v <- ggplot(var_comp_t, aes(x = reorder(Factor, -Relative_Variance), y = Relative_Variance))
  ms_plot_v <- ms_plot_v + geom_col(fill = "darkblue")
  ms_plot_v <- ms_plot_v + labs(x = "Factor", y = "Relative Variance") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ms_plot_v<- ms_plot_v+labs(title = paste0("C) ", measurement_name, " variance explained"), hjust = -2)
  ms_plot_v <- ms_plot_v + annotate("text", x = 2.0, y = 1.0, label =  perc_expl, size = 3)
  ms_plot_v <- ms_plot_v + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
   ms_plot_v
  
  
  ggsave(paste0("respirator_plots/", paste0(measurement_name, "race2_rel_var"), ".pdf"))
  assign(paste0(measurement_name_R, "race2_rel_var"),  ms_plot_v)
  dev.off()
  
  race2_coefficients_ns <- rbind(race2_coefficients_ns, coef_race2_i)
}

write.csv(race_coefficients_ns, "measurement race cofficients_ns.csv")
write.csv(race2_coefficients_ns, "measurement race2 cofficients_ns.csv")


write.csv(race_aov, "measurement race_aov.csv")
write.csv(race2_aov, "measurement race2_aov.csv")


# Create new scatterplots for kids only
library("RColorBrewer")


race_cols <- brewer.pal(5, "Spectral")

i <- 1

for(i in 1:nrow(measurements2))
{
  measurement <- meas_data_ns[,measurements2$Measurement[i]]
  meas_dat_i <- data.frame(meas_data_ns, measurement)
  measurement_name <- measurements2$full_meas_name[i]
  measurement_name_R <- measurements2$Measurement[i]
  
  #meas_kids_dat <- meas_dat_i[which(meas_dat_i$Age<19),]
  
  Age_i <- meas_dat_i$Age
  Age2_i <- meas_dat_i$Age^2
  Age3_i <- meas_dat_i$Age^3
  
  Sex_i <- meas_dat_i$Sex
  Race_i <- meas_dat_i$Race
  Race2_i <- meas_dat_i$Race2
  Synd_i <- meas_dat_i$Syndrome_Class

  # Sex plot
  dev.new(width=3, height=3, noRStudioGD = TRUE, units = "inch")
  plot<-ggplot(meas_dat_i, aes(Age, measurement, colour = factor(Sex)))
  plot<-plot+geom_point(alpha=0.8, size=1.5)+geom_smooth(method='loess')
  plot<-plot + scale_color_manual(values=c("#FDAE6B", "#54278F"),labels=c("Female","Male"))
  plot<-plot+theme_bw()
  plot<-plot+guides(colour=guide_legend(title="Sex"))
  plot<-plot+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  plot<-plot+ylab(paste0(measurement_name, " (mm)")) + xlab("Age (years)")
  plot<-plot+labs(title = paste0("A) ", measurement_name, " by Sex"))#+theme(plot.title = element_text(hjust = -0.2, vjust=2.12))
  plot_age_sex_kids<-plot	+ theme(legend.position = "none")
  assign(paste0(measurement_name_R, "_sex_kids"), plot_age_sex_kids) 
  plot_age_sex_kids
  ggsave(paste0("respirator_plots/", paste0(measurement_name, "_sex_kids"), ".pdf"))
  dev.off()
  
  # Race2 plot
  dev.new(width=5, height=3, noRStudioGD = TRUE, units = "inch")
  plot<-ggplot(meas_dat_i, aes(Age, measurement, colour = factor(Race2)))
  plot<-plot+geom_point(alpha=0.4, size=1.5)+geom_smooth(method='loess')
  plot<-plot + scale_color_manual(values=c(race_cols))
  plot<-plot+theme_bw()
  plot<-plot+guides(colour=guide_legend(title="Race (1000 Genomes Superpopulation)"))
  plot<-plot+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  plot<-plot+ylab(paste0(measurement_name, " (mm)")) + xlab("Age (years)")
  plot<-plot+labs(title = paste0("B) ", measurement_name, " by Race"))#+theme(plot.title = element_text(hjust = -0.2, vjust=2.12))
  plot_age_race2_kids<-plot	+ theme(legend.position = "none")
  assign(paste0(measurement_name_R, "_race_kids"), plot_age_race2_kids)
  plot_age_race2_kids
  ggsave(paste0("respirator_plots/", paste0(measurement_name, "_race_kids"), ".pdf"))
  dev.off()
  
}

# Make figures into rows. need to figure out the do.call or assign thing

library(Biobase)
library(grobblR)

i <- 1

for(i in 1:nrow(measurements2))
{
  measurement <- meas_data_ns[,measurements2$Measurement[i]]
  meas_dat_i <- data.frame(meas_data_ns, measurement)
  measurement_name <- measurements2$full_meas_name[i]
  measurement_name_R <- measurements2$Measurement[i]
  
  #plot_1 <- openPDF(paste0("respirator_plots/", measurement_name, "_sex_kids.pdf"))
  
  #plot_a <- rlang::syms(paste0(measurement_name_R, "_sex_kids"))[[1]]
  
  plot_a <- str2expression(paste0(measurement_name_R, "_sex_kids"))
  plot_b <- str2expression(paste0(measurement_name_R, "_race_kids"))
  plot_c <- str2expression(paste0(measurement_name_R, "race2_rel_var"))                    
                           
  #eval(str2expression(paste0(measurement_name_R, "_sex_kids"))) 
  dev.new(width=7.5, height=3, noRStudioGD = TRUE, units = "inch")
  plot1 <-  ggarrange(eval(plot_a), eval(plot_b), eval(plot_c), ncol = 3)
  plot1
  
  assign(paste0(measurement_name_R, "_panel"), plot1)
  ggsave(paste0("respirator_plots/", paste0(measurement_name, "_panel"), ".pdf"))
  dev.off()
}

# Stack panels

 plot_names <- str2expression(paste0(measurements2$Measurement, "_panel"))
  dev.new() 
  plot1 <-  ggarrange(eval(plot_names[1]), eval(plot_names[2]),eval(plot_names[3]),eval(plot_names[4]),eval(plot_names[5]),eval(plot_names[6]),eval(plot_names[7]),eval(plot_names[8]),eval(plot_names[9]),eval(plot_names[10]),eval(plot_names[11]),eval(plot_names[12]),eval(plot_names[13]),eval(plot_names[14]),eval(plot_names[15]),eval(plot_names[16]),eval(plot_names[17]),eval(plot_names[18]), ncol = 2, nrow = 9)
  plot1
  ggsave("respirator_plots/measurement_panel.pdf")
  
# Age distributions
#dev.new(2,5)
dev.new()
plot<- ggplot(meas_data, aes(Age))
plot<-plot + geom_density(aes(fill=factor(Sex)), alpha=0.5)+scale_fill_manual(values=	c("#0F084B", "#A0D2E7", "#3D60A7"))
#plot<-plot+xlim(c(-5,15))
plot <- plot + theme_bw()
plot<-plot + labs(title=c("Age distributions by sex"), x="Age (years)", fill="Group", size=1)
plot

dev.new()
plot<- ggplot(meas_data, aes(Age))
plot<-plot + geom_density(aes(fill=factor(Sex)), alpha=0.5)+scale_fill_manual(values=	c("#0F084B", "#A0D2E7", "#3D60A7"))
#plot<-plot+xlim(c(-5,15))
plot <- plot + theme_bw()
plot<-plot + labs(title=c("Age distributions by sex"), x="Age (years)", fill="Group", size=1)
plot


# barplots

Sex_counts<- meas_met %>% count(meas_met$Sex)
colnames(Sex_counts) <- c("Sex", "N")
race_counts<- meas_met %>% count(meas_met$Race)
colnames(race_counts) <- c("Race", "N")
race2_counts<- meas_met %>% count(meas_met$Race2)
colnames(race2_counts) <- c("Superpopulation", "N")
Synd_counts<- meas_met %>% count(meas_met$Syndrome_Class)
colnames(Synd_counts) <- c("Syndrome_Class", "N")


dev.new()
plot<-ggplot(Sex_counts, aes(Sex, N))
#plot<-plot+geom_bar(stat="identity", fill="#4682B4", "#A0D2E7")
plot<-plot+geom_bar(stat="identity")
plot<-plot+xlab("Sex") + ylab("N")+theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot<-plot+ggtitle("Sample Size by Sex")
plot
plot_sex<-plot


dev.new()
plot<-ggplot(race_counts, aes(Race, N))
#plot<-plot+geom_bar(stat="identity", fill="#4682B4", "#A0D2E7")
plot<-plot+geom_bar(stat="identity")
plot<-plot+xlab("Race") + ylab("N")+theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot<-plot+ggtitle("Sample Size by Race")
plot
plot_race<-plot



dev.new()
plot<-ggplot(race2_counts, aes(Superpopulation, N))
#plot<-plot+geom_bar(stat="identity", fill="#4682B4", "#A0D2E7")
plot<-plot+geom_bar(stat="identity")
plot<-plot+xlab("Superpopulation") + ylab("N")+theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot<-plot+ggtitle("Sample Size by 10K Superpopulation")
plot
plot_race2<-plot


dev.new()
plot<-ggplot(Synd_counts, aes(Syndrome_Class, N))
#plot<-plot+geom_bar(stat="identity", fill="#4682B4", "#A0D2E7")
plot<-plot+geom_bar(stat="identity")
plot<-plot+xlab("Syndrome Status") + ylab("N")+theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot<-plot+ggtitle("Sample Size by Syndrome Status")
plot
plot_synd<-plot

dev.new()
plot <- ggarrange(plot_sex, plot_race, plot_race2, plot_synd, ncol = 2, nrow = 2)
resize.win(5,7)
plot



# ANOVA for all measurements

# stack measurements

meas_long <- as.data.frame(array(NA, dim = c(0, 10)))
long_col_names <- c("Individual", "Measurement", "Age", "Sex", "Race2", "Syndrome_R", "Syndrome_Class", "Height","Weight","Value" )
colnames(meas_long) <- long_col_names

i <-1
for(i in 1:nrow(measurements2))
{
  Value <- meas_data[,measurements2$Measurement[i]]
  Individual <- meas_data$Individual
  Age <- meas_data$Age
  Sex <- meas_data$Sex
  Race2 <- meas_data$Race2
  Syndrome_R <- meas_data$Syndrome_R
  Syndrome_Class <- meas_data$Syndrome_Class
  Height <- meas_data$HT
  Weight <- meas_data$WT
  Measurement <- measurements2$Measurement[i]
  
  chunk <- data.frame(Individual, Measurement, Age, Sex, Race2, Syndrome_R, Syndrome_Class, Height, Weight, Value)
  meas_long <- rbind(meas_long, chunk)
}

rm(Value, Individual, Age, Sex, Race2, Syndrome_R, Syndrome_Class, Height,Weight, Measurement)
  
# Create sample of nonsyndromic kids
meas_long_ns <- meas_long[which(meas_long$Syndrome_Class == "Control"),]
meas_long_ns <- meas_long_ns[which(meas_long_ns$Age <19),]

all_anova <- lm(meas_long_ns$Value ~ meas_long_ns$Measurement + meas_long_ns$Age + meas_long_ns$Sex + meas_long_ns$Race2 + meas_long_ns$Measurement * meas_long_ns$Age + meas_long_ns$Measurement * meas_long_ns$Sex + meas_long_ns$Race2)


# Perform Procrustes ANOVA 



#dat_ns <- v_data[which(v_cov$Age < 19 & v_class$Ctr_Syndrome == "Control" & v_class$Race2 != "unknown/not reported"),]
#v_cov_ns <- v_cov[which(v_cov$Age < 19 & v_class$Ctr_Syndrome == "Control"& v_class$Race2 != "unknown/not reported"),]
#v_class_ns <- v_class[which(v_cov$Age < 19 & v_class$Ctr_Syndrome == "Control"& v_class$Race2 != "unknown/not reported"),]

# Get data from Race module. This is to make sure that we are using the data corrected for camera and method using the Harry Comparison module

dat_kids <- read.csv("dat_kids.csv", header = T)
dat_kids <- dat_kids[,-1]
dat_kids_3d <- arrayspecs(dat_kids, 65, 3)
cov_kids <- read.csv("cov_kids.csv", header = T)
class_kids <- read.csv("class_kids.csv", header = T)

dat_kids <- dat_kids[-which(class_kids$Race2 =="unknown/not reported"), ]
cov_kids <- cov_kids[-which(class_kids$Race2 =="unknown/not reported"), ]
class_kids <- class_kids[-which(class_kids$Race2 =="unknown/not reported"), ]

#dat_ns_3d <- arrayspecs(dat_ns, 65, 3)
#dat_ns_3d <- procGPA(dat_ns_3d)
#dat_ns_2d <- two.d.array(dat_ns_3d$rotated)

kids_mshape <- mshape(dat_kids_3d)

outliers <- plotOutliers(dat_kids_3d)


# Calculate Procrustes Distance from mean
dists <- array(NA, dim = c(nrow(dat_kids), 1))
dists <- as.data.frame(dists)
colnames(dists) <- "Procdist"
  

for(i in 1:nrow(dat_kids))
{
  dist_i <- procdist(dat_kids_3d[,,i], kids_mshape)
  dists[i,] <- dist_i
}

# remove outliers

dat_kids <- dat_kids[which(dists$Procdist <0.09),]
cov_kids <- cov_kids[which(dists$Procdist <0.09),]
class_kids <- class_kids[which(dists$Procdist <0.09),]
dat_kids_3d <- arrayspecs(dat_kids, 65, 3)

#Age2_ns <- v_cov_ns$Age^2
#Age3_ns <- v_cov_ns$Age^3

Age2_ns <- cov_kids$Age^2
Age3_ns <- cov_kids$Age^3


#pr_gm_ns<-geomorph.data.frame(shape=dat_ns_3d$rotated, Age=v_cov_ns$Age, Age_2=Age2_ns,Age_3=Age3_ns, Sex=v_class_ns$Sex, Race2= v_class_ns$Race2, CS=v_cov_ns$CS, height = v_cov_ns$HT, weight = v_cov_ns$WT)
pr_gm_ns<-geomorph.data.frame(shape=dat_kids_3d, Age=cov_kids$Age, Age_2=Age2_ns,Age_3=Age3_ns, Sex=class_kids$Sex, Race2 = class_kids$Race2, CS=cov_kids$CS, height = cov_kids$HT, weight = cov_kids$WT)

ac_pr_ns <- procD.lm(shape ~ Age + Age_2 + Age_3 +Sex + Race2 + Age * Sex + Age * Race2, data=pr_gm_ns, iter=999)


# Create graph of change in Sex and Race-related variance with age


Age_groups_min <- seq(0:17)-1
Age_groups_max <- Age_groups_min + 2
Age_groups <- as.data.frame(cbind(Age_groups_min, Age_groups_max))
colnames(Age_groups) <- c("min", "max")


sex_race_results <- as.data.frame(array(NA, dim = c(0, 7)))
colnames(sex_race_results) <- c("age_label", "age_val", "df", "Race_rsq", "Race_prob", "Sex_rsq", "Sex_prob")

i <-1
for (i in 1:nrow(Age_groups))
{
  
  
  age_group_i <- dat_kids[which(cov_kids$Age > (Age_groups$min[i]-1) & cov_kids$Age < (Age_groups$max[i]+1)),]
  age_cov_i <- cov_kids[which(cov_kids$Age > (Age_groups$min[i]-1) & cov_kids$Age < (Age_groups$max[i]+1)),]
  age_class_i <- class_kids[which(cov_kids$Age > (Age_groups$min[i]-1) & cov_kids$Age < (Age_groups$max[i]+1)),]
  age_label <- paste0(Age_groups$min[i], "_", Age_groups$max[i])
  age_val <- mean(Age_groups$min[i], Age_groups$max[i])
  
  dat_i_3d <- arrayspecs(age_group_i, 65, 3)
  
  Age2_ns_i <- age_cov_i$Age^2
  Age3_ns_i <- age_cov_i$Age^3
  
  Race_i <- age_class_i$Race2
  Sex_i <- age_class_i$Sex
  
  pr_gm_ns_i<-geomorph.data.frame(shape=  dat_i_3d, Age = age_cov_i$Age, Age_2=Age2_ns_i ,Age_3=Age3_ns_i, Sex=Sex_i, Race2= Race_i)
  
  ac_pr_ns_i <- procD.lm(shape ~ Age + Sex + Race2, data=pr_gm_ns_i, iter=999)
  aov_table_i <- as.data.frame(ac_pr_ns_i$aov.table)
  term_i <- rownames(aov_table_i)
  aov_table_i <- cbind(term_i, aov_table_i)
  colnames(aov_table_i) <- c("term", "DF", "SS", "MS", "Rsq", "F", "Z", "Prob")
  Race_rsq_i <- aov_table_i$Rsq[which(aov_table_i$term == "Race2")]
  Race_prob_i <- aov_table_i$Prob[which(aov_table_i$term == "Race2")]
  Sex_rsq_i <- aov_table_i$Rsq[which(aov_table_i$term == "Sex")]
  Sex_prob_i <- aov_table_i$Prob[which(aov_table_i$term == "Sex")]
  df_i <- aov_table_i$DF[which(aov_table_i$term == "Total")]
  
  result_i <- data.frame(age_label, age_val, df_i, Race_rsq_i, Race_prob_i, Sex_rsq_i, Sex_prob_i)
  colnames(result_i) <- c("age_label", "age_val", "df", "Race_rsq", "Race_prob", "Sex_rsq", "Sex_prob")
  sex_race_results <-rbind(sex_race_results, result_i)
  
}


# Plot these results
# Create dataframe for ggplot

Race_res <- sex_race_results[,c(1,2,3,4,5),]
race_label <- array("Race", dim = c(nrow(Race_res), 1))
Race_res <- cbind(race_label, Race_res)
colnames(Race_res) <- c("Variable", "Age_group", "Age", "df", "Rsq", "prob")

Sex_res <- sex_race_results[,c(1,2,3,6,7),]
sex_label <- array("Sex", dim = c(nrow(Sex_res), 1))
Sex_res <- cbind(sex_label, Sex_res)
colnames(Sex_res) <- c("Variable", "Age_group", "Age", "df", "Rsq", "prob")

SR_results <- rbind(Race_res, Sex_res)

# Filter out nonsignificant
SR_results <- SR_results[-which(SR_results$prob > 0.05),]

plot <- ggplot(SR_results, aes( x = age_val, y = ))


dev.new()
plot<-ggplot(SR_results, aes(Age, Rsq, colour = factor(Variable)))
plot<-plot+geom_point(alpha=0.3, size=2)+geom_smooth(method='loess') + geom_text_repel(label= SR_results$df, size=2, max.overlaps = 5)
plot<-plot + scale_color_manual(values=c("royalblue2", "khaki3"),labels=c("Race","Sex"))
plot<-plot+theme_bw()
plot<-plot+guides(colour=guide_legend(title="Variable"))
plot<-plot+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot<-plot+ylab("Variance Explained (R2)")
plot<-plot+labs(title = "Facial Shape Variance due to Race and Sex by Age") #+theme(plot.title = element_text(hjust = -0.2, vjust=2.12))
plot

ggsave("respirator_plots/variance explained by sex and race by age.pdf")

library(ggExtra)

# Create scatter/density plots for sex race and age

CVA_race <- CVA(dat_kids, class_kids$Race2, iter = 999)

sex_race_pca <- prcomp(dat_kids)
pc_scores <- as.data.frame(sex_race_pca$x)


scores <- as.data.frame(CVA_race$CVscores)
colnames(scores) <- c("CV1", "CV2","CV3", "CV4")
Race <- class_kids$Race2
Sex <- class_kids$Sex

#CV1 vs CV2
 
dev.new(width=6, height=5, noRStudioGD = TRUE, units = "inch")
plot<-ggplot(scores, aes(CV1,CV2, colour= factor(Race)))
plot<-plot+geom_point(alpha=0.4, size = class_kids$Age/10)
plot<-plot + scale_color_manual(values=c("lightblue","blue", "gray45", "darkred", "darkgreen"))
plot<-plot + geom_density2d(aes(colour = factor(Race)), size=0.15)
#plot<-plot + stat_density2d(aes(fill = ..level..), geom="polygon", alpha=0.1)+ylim(-0.08,0.08)
plot<- plot + xlim(-8,5) + ylim(-4, 7)
plot<-plot+guides(colour=guide_legend(title="Genotype Group"))
plot<-plot+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot<-plot+labs(title = "CVA for facial shape by Race (1000 Genomes Superpopulation)")+theme(plot.title = element_text(hjust = -0.15, vjust=2.12))
#plot<- ggMarginal(plot, scores, CV1, CV2, type = c("density"), margins = c("both","x", "y"), size = 3, groupColour = TRUE, groupFill = TRUE)
#plot<- ggMarginal(plot, scores, CV2, type = c("density"), margins = c("y"), size = 3, groupColour = TRUE, groupFill = TRUE)

plot


ggsave("respirator_plots/Race CV1_2 y margin.svg", width = 6, height = 5)

#CV3 vs CV4

dev.new(width=6, height=5, noRStudioGD = TRUE, units = "inch")
plot<-ggplot(scores, aes(CV3,CV4, colour= factor(Race)))
plot<-plot+geom_point(alpha=0.4, size = class_kids$Age/10)
plot<-plot + scale_color_manual(values=c("lightblue","blue", "gray45", "darkred", "darkgreen"))
plot<-plot + geom_density2d(aes(colour = factor(Race)), size=0.15)
#plot<-plot + stat_density2d(aes(fill = ..level..), geom="polygon", alpha=0.1)+ylim(-0.08,0.08)
plot<-plot+guides(colour=guide_legend(title="Genotype Group"))
plot<- plot + xlim(-7,7) + ylim(-7, 7)
plot<-plot+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot<-plot+labs(title = "CVA for facial shape by Race (1000 Genomes Superpopulation)")+theme(plot.title = element_text(hjust = -0.15, vjust=2.12))
#plot<- ggMarginal(plot, scores, CV3, CV4, type = c("density"), margins = c("both", "x", "y"), size = 3, groupColour = TRUE, groupFill = TRUE)
#plot<- ggMarginal(plot, scores, CV4, type = c("density"), margins = c("y"), size = 3, groupColour = TRUE, groupFill = TRUE)

plot


ggsave("respirator_plots/Race CV3_4 y margin.svg", width = 6, height = 5)


#PC plots


dev.new(width=6, height=5, noRStudioGD = TRUE, units = "inch")
plot<-ggplot(pc_scores, aes(PC1,PC2, colour= factor(Race)))
plot<-plot+geom_point(alpha=0.4, size = class_kids$Age/10)
plot<-plot + scale_color_manual(values=c("lightblue","blue", "gray45", "darkred", "darkgreen"))
plot<-plot + geom_density2d(aes(colour = factor(Race)), size = 0.015)
plot<-plot + stat_density2d(aes(fill = ..level..), geom="polygon", alpha=0.1)+ylim(-0.08,0.08)
plot<-plot+guides(colour=guide_legend(title="Race"))
#plot<- plot + xlim(-10,10) + ylim(-0.10, 0.10)
plot<-plot+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot<-plot+labs(title = "PCA for facial shape by Race (1000 Genomes Superpopulation)")+theme(plot.title = element_text(hjust = -0.15, vjust=2.12))
#plot<- ggMarginal(plot, scores, PC1, PC2, type = c("density"), margins = c("both", "x", "y"), size = 3, groupColour = TRUE, groupFill = TRUE)
plot<- ggMarginal(plot, scores, PC1, type = c("density"), margins = c("x"), size = 3, groupColour = TRUE, groupFill = TRUE)

plot


ggsave("respirator_plots/Race x margin.svg", width = 6, height = 5)


dev.new(width=5, height=5, noRStudioGD = TRUE, units = "inch")
plot<-ggplot(pc_scores, aes(PC3,PC4, colour= factor(Race)))
plot<-plot+geom_point(alpha=0.4, size = class_kids$Age/10)
plot<-plot + scale_color_manual(values=c("lightblue","blue", "gray45", "darkred", "darkgreen"))
plot<-plot + geom_density2d(aes(colour = factor(Race)), size = 0.015)
plot<-plot + stat_density2d(aes(fill = ..level..), geom="polygon", alpha=0.1)+ylim(-0.08,0.08)
plot<-plot+guides(colour=guide_legend(title="Race"))
plot<- plot + xlim(-0.05,0.05) + ylim(-0.05, 0.05)
plot<-plot+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot<-plot+labs(title = "PCA for facial shape by Race (1000 Genomes Superpopulation)")+theme(plot.title = element_text(hjust = -0.15, vjust=2.12))
#plot<- ggMarginal(plot, scores, PC3, PC4, type = c("density"), margins = c("both", "x", "y"), size = 3, groupColour = TRUE, groupFill = TRUE)
plot

ggsave("respirator_plots/Race PC3_4.svg", width = 6, height = 5)


dev.new(width=5, height=5, noRStudioGD = TRUE, units = "inch")
plot<-ggplot(pc_scores, aes(PC5,PC6, colour= factor(Race)))
plot<-plot+geom_point(alpha=0.4, size = class_kids$Age/10)
plot<-plot + scale_color_manual(values=c("lightblue","blue", "gray45", "darkred", "darkgreen"),labels=c("EUR" ,"AMR", "SAS", "AFR", "EAS"))
plot<-plot + geom_density2d(aes(colour = factor(Race)), size = 0.015)
plot<-plot + stat_density2d(aes(fill = ..level..), geom="polygon", alpha=0.1)+ylim(-0.08,0.08)
plot<-plot+guides(colour=guide_legend(title="Race"))
plot<- plot + xlim(-0.03,0.03) + ylim(-0.03, 0.03)
plot<-plot+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot<-plot+labs(title = "PCA for facial shape by Race (1000 Genomes Superpopulation)")+theme(plot.title = element_text(hjust = -0.15, vjust=2.12))
plot<- ggMarginal(plot, scores, PC5, PC6, type = c("density"), margins = c("both", "x", "y"), size = 3, groupColour = TRUE, groupFill = TRUE)
plot

ggsave("respirator_plots/Race PC5_6.pdf")


#.  by sex

dev.new(width=5, height=5, noRStudioGD = TRUE, units = "inch")
plot<-ggplot(pc_scores, aes(PC1,PC2, colour= factor(Sex)))
plot<-plot+geom_point(alpha=0.4, size = class_kids$Age/10)
plot<-plot + scale_color_manual(values=c("gray45","darkred", "white"),labels=c("F" ,"M", "NA"))
plot<-plot + geom_density2d(aes(colour = factor(Sex)), size = 0.015)
plot<-plot + stat_density2d(aes(fill = ..level..), geom="polygon", alpha=0.1)+ylim(-0.08,0.08)
plot<-plot+guides(colour=guide_legend(title="Sex"))
#plot<- plot + xlim(-10,10) + ylim(-0.10, 0.10)
plot<-plot+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot<-plot+labs(title = "PCA for facial shape by Sex")+theme(plot.title = element_text(hjust = -0.15, vjust=2.12))
plot<- ggMarginal(plot, scores, PC1, PC2, type = c("density"), margins = c("both", "x", "y"), size = 3, groupColour = TRUE, groupFill = TRUE)
plot


ggsave("respirator_plots/Sex PC1_2.pdf")


# Histogram of eigenvalues

sum_eig <- sum(sex_race_pca$sdev)
PC1_per <- sex_race_pca$sdev[1]/sum_eig
PC2_per <- sex_race_pca$sdev[2]/sum_eig


#ac_pr3 <- procD.lm(shape ~ Age + Age_2 + Age_3 +Sex+Race, data=pr_gm, iter=1000)
#ac_pr_a<-plotAllometry(ac_pr3, v_cov$imp_Age, logsz=FALSE)


#age_CAC<-data.frame(ac_pr_a$RegScore,v_cov$imp_Age)
#colnames(age_CAC)<-c("CAC","Age")



plot<-ggplot(age_CAC, aes(Age, CAC, colour = factor(Sex)))
plot<-plot+geom_point(alpha=0.3, size=2)+geom_smooth(method='loess')
plot<-plot + scale_color_manual(values=c("blue", "lightblue", "orchid", "slateblue", "aquamarine", "slategrey"),labels=c("EUR" ,"AMR", "SAS", "AFR", "EAS"))
plot<-plot+theme_bw()
plot<-plot + geom_density2d(aes(colour = factor(v_class$Sex)), size=0.15)
plot<-plot + stat_density2d(aes(fill = ..level..), geom="polygon", alpha=0.1)
dev.new(width=8, height=5)
plot<-plot+guides(colour=guide_legend(title="Sex"))
plot<-plot+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot<-plot+ylab("Regresssion Score (CAC)")
plot<-plot+labs(title = "CVA for facial shape by Race (1000 Genomes Superpopulation)")+theme(plot.title = element_text(hjust = -0.15, vjust=2.12))
plot_race_pca<-plot	
plot_race_pca









#Explore PCA of data. This is script from the syndrome analysis.

g_mean<-mshape(dat_kids_3d)

data_2d<-dat_kids

#Generate morphs for first 30PCs

pr_dat<-data_2d
pca<-prcomp(pr_dat)


pc_names<-colnames(pca$rotation)
scale<-1
mean_pc<-g_mean

i<-1




#BEGINNING LOOP. This loop automatically positions faces in front and side views.
for (i in 1:30)
{
  clear3d()
  
  pc_rot<-pca$rotation[,i]
  pc_rot<-rbind(pc_rot,pc_rot,pc_rot)
  pc_rot<-mshape(arrayspecs(pc_rot, 65,3))
  
  PC_min_x<-min(pca$x[,i])	
  PC_max_x<-max(pca$x[,i])
  PC<-pc_names[i]
  
  pc_min<-g_mean+(PC_min_x*pc_rot)
  pc_max<-g_mean+(PC_max_x*pc_rot)
  
  
  c_sup_m<-pPsup(reference, pc_min)
  pc_min_face<-tps3d(face,reference,c_sup_m$Mp2)
  
  c_sup_M<-pPsup(reference, pc_max)
  pc_max_face<-tps3d(face,reference,c_sup_M$Mp2)
  
  mean_face<-tps3d(face,reference,g_mean)
  
  clear3d()
  rgl.close()
  dev.new(width=1, height=3)
  
  #Front views
  open3d(zoom=0.75, userMatrix = front)
  par3d("windowRect"= c(0,100,600,800))
  shade3d(pc_min_face, color="grey")
  
  rgl.snapshot(paste0("respirator_plots/pc_morphs/",PC,"_","min_front",".png"), top = TRUE )
  
  clear3d()
  rgl.close()
  
  #Front views
  open3d(zoom=0.75, userMatrix = front)
  par3d("windowRect"= c(0,100,600,800))
  shade3d(pc_max_face, color="grey")
  
  rgl.snapshot(paste0("respirator_plots/pc_morphs/",PC,"_","max_front",".png"), top = TRUE )
  
  clear3d()
  rgl.close()
  
  #Side views
  open3d(zoom=0.75, userMatrix = side)
  par3d("windowRect"= c(0,100,600,800))
  shade3d(pc_min_face, color="grey")
  
  rgl.snapshot(paste0("respirator_plots/pc_morphs/",PC,"_","min_side",".png"), top = TRUE )
  
  clear3d()
  rgl.close()
  
  #Side views
  open3d(zoom=0.75, userMatrix = side)
  par3d("windowRect"= c(0,100,600,800))
  shade3d(pc_max_face, color="grey")
  
  rgl.snapshot(paste0("respirator_plots/pc_morphs/",PC,"_","max_side",".png"), top = TRUE )
  
  clear3d()
  rgl.close()
  
  #dev.new(width=1.8, height=4)
  
  #Heatmaps
  
  open3d(zoom=0.75, userMatrix = front)
  par3d("windowRect"= c(0,100,600,800))
  mD<-meshDist(mean_face, pc_max_face, rampcolors = c("blue","green","red", "dark red"))
  
  rgl.snapshot(paste0("respirator_plots/pc_morphs/heat_",PC,"_","front",".png"), top = TRUE )
  dev.off()
  
  
  clear3d()
  rgl.close()
  
  open3d(zoom=0.75, userMatrix = side)
  par3d("windowRect"= c(0,100,600,800))
  mD<-meshDist(mean_face, pc_max_face, rampcolors = c("blue","green","red", "dark red"))
  
  rgl.snapshot(paste0("respirator_plots/pc_morphs/heat_",PC,"_","side",".png"), top = TRUE )
  
  dev.copy(pdf,paste0("respirator_plots/pc_morphs/heat_",PC,"_scale_"), width=1.7, height=4)
  dev.off()
  
  
  clear3d()
  rgl.close()
  
}
#END LOOP



image_write(syn_img, path = paste0("morphs/",synd_i,"_","panel",".png"), format = "png")	
}
#End loop


#Create page withPCs1-8
i<-1
PC_min_img_f<- image_read(paste0("res_pc_morphs/PC",i,"_min_front.png"))
PC_max_img_f<- image_read(paste0("res_pc_morphs/PC",i,"_max_front.png"))
PC_min_img_s<- image_read(paste0("res_pc_morphs/PC",i,"_min_side.png"))
PC_max_img_s<- image_read(paste0("res_pc_morphs/PC",i,"_max_side.png"))
PC_HT_img_f<- image_read(paste0("res_pc_morphs/heat_PC",i,"_front.png"))
PC_HT_img_s<- image_read(paste0("res_pc_morphs/heat_PC",i,"_side.png"))

PC_row <- c(PC_min_img_f,PC_max_img_f,PC_min_img_s,PC_max_img_s,PC_HT_img_f,PC_HT_img_s)
imgs<-image_append(image_scale(PC_row, "x200"))
imgs<-image_border(imgs, "white", "0x15")
imgs<-image_annotate(imgs, paste0("PC",i), font = "sans", size = 30)


for (i in 2:8)
{
  PC_min_img_f<- image_read(paste0("res_pc_morphs/PC",i,"_min_front.png"))
  PC_max_img_f<- image_read(paste0("res_pc_morphs/PC",i,"_max_front.png"))
  PC_min_img_s<- image_read(paste0("res_pc_morphs/PC",i,"_min_side.png"))
  PC_max_img_s<- image_read(paste0("res_pc_morphs/PC",i,"_max_side.png"))
  PC_HT_img_f<- image_read(paste0("res_pc_morphs/heat_PC",i,"_front.png"))
  PC_HT_img_s<- image_read(paste0("res_pc_morphs/heat_PC",i,"_side.png"))
  
  
  PC_row <- c(PC_min_img_f,PC_max_img_f,PC_min_img_s,PC_max_img_s,PC_HT_img_f,PC_HT_img_s)
  img2<-image_append(image_scale(PC_row, "x200"))
  img2<-image_border(img2, "white", "0x15")
  img2<-image_annotate(img2, paste0("PC",i), font = "sans", size = 30)
  
  imgs<-c(imgs,img2)
  
  imgs<-image_append(image_scale(imgs), stack=TRUE)
}

image_browse(imgs)
image_write(imgs, path = "res_pc_morphs/PC1_8.png", format = "png")

#Create page withPCs 9-16
i<-9
PC_min_img_f<- image_read(paste0("res_pc_morphs/PC",i,"_min_front.png"))
PC_max_img_f<- image_read(paste0("res_pc_morphs/PC",i,"_max_front.png"))
PC_min_img_s<- image_read(paste0("res_pc_morphs/PC",i,"_min_side.png"))
PC_max_img_s<- image_read(paste0("res_pc_morphs/PC",i,"_max_side.png"))
PC_HT_img_f<- image_read(paste0("res_pc_morphs/heat_PC",i,"_front.png"))
PC_HT_img_s<- image_read(paste0("res_pc_morphs/heat_PC",i,"_side.png"))

PC_row <- c(PC_min_img_f,PC_max_img_f,PC_min_img_s,PC_max_img_s,PC_HT_img_f,PC_HT_img_s)
imgs<-image_append(image_scale(PC_row, "x200"))
img2<-image_border(img2, "white", "0x15")
img2<-image_annotate(img2, paste0("PC",i), font = "sans", size = 30)

for (i in 10:16)
{
  PC_min_img_f<- image_read(paste0("res_pc_morphs/PC",i,"_min_front.png"))
  PC_max_img_f<- image_read(paste0("res_pc_morphs/PC",i,"_max_front.png"))
  PC_min_img_s<- image_read(paste0("res_pc_morphs/PC",i,"_min_side.png"))
  PC_max_img_s<- image_read(paste0("res_pc_morphs/PC",i,"_max_side.png"))
  PC_HT_img_f<- image_read(paste0("res_pc_morphs/heat_PC",i,"_front.png"))
  PC_HT_img_s<- image_read(paste0("res_pc_morphs/heat_PC",i,"_side.png"))
  
  
  PC_row <- c(PC_min_img_f,PC_max_img_f,PC_min_img_s,PC_max_img_s,PC_HT_img_f,PC_HT_img_s)
  img2<-image_append(image_scale(PC_row, "x200"))
  img2<-image_border(img2, "white", "0x15")
  img2<-image_annotate(img2, paste0("PC",i), font = "sans", size = 30)
  imgs<-c(imgs,img2)
  
  imgs<-image_append(image_scale(imgs), stack=TRUE)
}

image_browse(imgs)
image_write(imgs, path = "res_pc_morphs/PC9_16.png", format = "png")	

#Create page withPCs 17-24
i<-17
PC_min_img_f<- image_read(paste0("res_pc_morphs/PC",i,"_min_front.png"))
PC_max_img_f<- image_read(paste0("res_pc_morphs/PC",i,"_max_front.png"))
PC_min_img_s<- image_read(paste0("res_pc_morphs/PC",i,"_min_side.png"))
PC_max_img_s<- image_read(paste0("res_pc_morphs/PC",i,"_max_side.png"))
PC_HT_img_f<- image_read(paste0("res_pc_morphs/heat_PC",i,"_front.png"))
PC_HT_img_s<- image_read(paste0("res_pc_morphs/heat_PC",i,"_side.png"))

PC_row <- c(PC_min_img_f,PC_max_img_f,PC_min_img_s,PC_max_img_s,PC_HT_img_f,PC_HT_img_s)
imgs<-image_append(image_scale(PC_row, "x200"))
img2<-image_border(img2, "white", "0x15")
img2<-image_annotate(img2, paste0("PC",i), font = "sans", size = 30)


for (i in 18:24)
{
  PC_min_img_f<- image_read(paste0("res_pc_morphs/PC",i,"_min_front.png"))
  PC_max_img_f<- image_read(paste0("res_pc_morphs/PC",i,"_max_front.png"))
  PC_min_img_s<- image_read(paste0("res_pc_morphs/PC",i,"_min_side.png"))
  PC_max_img_s<- image_read(paste0("res_pc_morphs/PC",i,"_max_side.png"))
  PC_HT_img_f<- image_read(paste0("res_pc_morphs/heat_PC",i,"_front.png"))
  PC_HT_img_s<- image_read(paste0("res_pc_morphs/heat_PC",i,"_side.png"))
  
  
  PC_row <- c(PC_min_img_f,PC_max_img_f,PC_min_img_s,PC_max_img_s,PC_HT_img_f,PC_HT_img_s)
  img2<-image_append(image_scale(PC_row, "x200"))
  img2<-image_border(img2, "white", "0x15")
  img2<-image_annotate(img2, paste0("PC",i), font = "sans", size = 30)
  
  
  imgs<-c(imgs,img2)
  
  imgs<-image_append(image_scale(imgs), stack=TRUE)
}

image_browse(imgs)
image_write(imgs, path = "res_pc_morphs/PC17_24.png", format = "png")	




