# Timbrell, L. et al. (2022). Testing interobserver error under a collaborative research framework for studying shape lithic variability

# Data available to download at: https://github.com/lucytimbrell/error_analysis_lithics/ 

# To reproduce the analysis, download the data as a zip file, extract it to a folder, and set the
# working directory below to the folder where you extracted the data.


# Clear R environment 
rm(list=ls())

####  Set working directory and load packages ####

setwd("...")

if(!require("tidyverse")) install.packages('Momocs', repos='http://cran.us.r-project.org')  
if(!require("ggpubr")) install.packages('tidyverse', repos='http://cran.us.r-project.org')
if(!require("rstatic")) install.packages('rio', repos='http://cran.us.r-project.org')
if(!require("psych")) install.packages('ggplot2', repos='http://cran.us.r-project.org')
if(!require("Momocs")) install.packages('ggpubr', repos='http://cran.us.r-project.org')
if(!require("rio")) install.packages('patchwork', repos='http://cran.us.r-project.org')
if(!require("patchwork")) install.packages('MASS', repos='http://cran.us.r-project.org')

library(tidyverse)
library(ggpubr)
library(rstatix)
library(psych)
library(Momocs)
library(rio)
library(patchwork)

####  Metric analysis ####

# Load meta-data
error_data <- read.csv("error_data.csv")
colnames(error_data) <- c("ID", "Museum", "Artefact", "Length", "Width", "Thickness", "Photo_ID")

error_data2 <- as_tibble(error_data)
error_data2$Tool <- factor(error_data2$Artefact)

# Summary statistics 

length_sum <- error_data2 %>%
  group_by(Artefact) %>%
  get_summary_stats(Length, type = "mean_sd")

length_sum$CV <- rep(0,6)

for(i in 1:nrow(length_sum)){
  cv <- (length_sum[i,5]/length_sum[i,4])*100
  length_sum[i,6] <- cv
}

length_sum$Variance <- length_sum$sd^2

width_sum <- error_data2 %>%
  group_by(Artefact) %>%
  get_summary_stats(Width, type = "mean_sd")

width_sum$CV <- rep(0,6)

for(i in 1:nrow(width_sum)){
  cv <- (width_sum[i,5]/width_sum[i,4])*100
  width_sum[i,6] <- cv
}

width_sum$Variance <- width_sum$sd^2

thickness_sum <- error_data2 %>%
  group_by(Artefact) %>%
  get_summary_stats(Thickness, type = "mean_sd")

thickness_sum$CV <- rep(0,6)

for(i in 1:nrow(thickness_sum)){
  cv <- (thickness_sum[i,5]/thickness_sum[i,4])*100
  thickness_sum[i,6] <- cv
}

thickness_sum$Variance <- thickness_sum$sd^2

# Boxplots
length_bxp <- ggboxplot(error_data2,
                        x = "Tool",
                        y = "Length", 
                        fill = c("#44AA99","#AA4499", "#6699CC", "#999933", "#888888", "#000000"), 
                        bxp.errorbar = TRUE)
length_bxp

width_bxp <- ggboxplot(error_data2,
                       x = "Tool",
                       y = "Width", 
                       fill = c("#44AA99","#AA4499", "#6699CC", "#999933", "#888888", "#000000"), 
                       bxp.errorbar = TRUE)
width_bxp

thickness_bxp <- ggboxplot(error_data2,
          x = "Tool",
          y = "Thickness", 
          fill = c("#44AA99","#AA4499", "#6699CC", "#999933", "#888888", "#000000"), 
          bxp.errorbar = TRUE)

thickness_bxp

length_bxp/width_bxp/thickness_bxp +
  plot_annotation(tag_levels = "a",
                  theme = theme(plot.title = element_text(size = 20))) 


# ANOVA and Tukey HSD
res <- aov(Length~factor(Artefact), error_data2)
summary(res)
TukeyHSD(res)

res1 <- aov(Width ~ factor(Artefact), error_data2)
summary(res1)
TukeyHSD(res1)

res2 <- aov(Thickness ~ factor(Artefact), error_data2)
summary(res2)
TukeyHSD(res2)

# Intra-class correlation coefficient 

data2 <- read.csv("rearranged_error_data.csv") # rearranged data

rownames(data2 ) <- data2$X # Assign first column as row name rather than as a variable
data2  <- data2[,-1] # remove first column 

ICC(data2) # we use ICC3

# R calculation

output <- matrix(0, nrow = ncol(data2), ncol = ncol(data2 ))
for(i in 1:ncol(data2 )){
  for(m in 1:ncol(data2)){
    a <- cbind(data2[,i],data2[,m])
    b <- var(c(a[,1], a[,2]), use = "all.obs") # sample variance of the whole of matrix 'a'
    b <- b[1]
    c <- (a[,1])-(a[,2])
    d <- var(c)
    output[i,m] = b/(b+d)
  }
}

rownames(output) <- colnames(data2)
colnames(output) <- colnames(data2)
write.csv(output, "R_scores_metric.csv")

#### FUNCTION FOR CALCULATING TEAM TEM ######

calculate_team_tem <- function(n, k, m) {
  
  if(class(m) != "data.frame") {
    stop(
      paste(
        strwrap("Measurements should be supplied as a data frame.
                Please check and try again.",
                width = 80),
        "\n",
        collapse = "\n"
      )
    )
  }
  
  if(ncol(m) != k | nrow(m) != n) {
    stop(
      paste(
        strwrap("Measurements data frame should have rows equal to number of
                subjects and columns equal to number of observers.
                Please check and try again.",
                width = 80),
        "\n",
        collapse = "\n"
      )
    )
  }
  
  ## Function to square a value
  fun_squared <- function(x) { x ^ 2 }
  
  ## Square of all measurements and sum per subject measurement squares
  firstPart <- rowSums(sapply(m, FUN = fun_squared))
  
  ## Sum per subject measurements and then squared divided by number of
  ## observers
  secondPart <- (rowSums(m) ^ 2) / k
  
  ## Calculate TEM
  tem <- sqrt(sum(firstPart - secondPart) / (n * (k - 1)))
  
  ## Return output
  return(tem)
}

###### FUNCTION FOR CALCULATING RELATIVE TEM  ######
calculate_relative_tem <- function(tem, mean_value) {
  ## Calculate relative TEM
  relative_tem <- (tem / mean_value) * 100
  
  ## Return output
  return(relative_tem)
}

# TEM and RTEM 
TEM <- calculate_team_tem(nrow(data2), ncol(data2), data2)
RTEM <- calculate_relative_tem(tem = TEM, mean_value = as.matrix(data2))

## COMPARATIVE STUDY ##

# Load data
comparative_data <- read.csv("comparative_data.csv")
colnames(comparative_data) <- c("ID", "Museum", "Artefact", "Length", "Width", "Thickness", "Photo_ID")

comparative_data2 <- as_tibble(comparative_data)


# Summary statistics
length_sum <- comparative_data2 %>%
  group_by(Artefact) %>%
  get_summary_stats(Length, type = "mean_sd")

length_sum$CV <- rep(0,6)

for(i in 1:nrow(length_sum)){
  cv <- (length_sum[i,5]/length_sum[i,4])*100
  length_sum[i,6] <- cv
}

width_sum <- comparative_data2 %>%
  group_by(Artefact) %>%
  get_summary_stats(Width, type = "mean_sd")

width_sum$CV <- rep(0,6)

for(i in 1:nrow(width_sum)){
  cv <- (width_sum[i,5]/width_sum[i,4])*100
  width_sum[i,6] <- cv
}


thickness_sum <- comparative_data2 %>%
  group_by(Artefact) %>%
  get_summary_stats(Thickness, type = "mean_sd")

thickness_sum$CV <- rep(0,6)

for(i in 1:nrow(thickness_sum)){
  cv <- (thickness_sum[i,5]/thickness_sum[i,4])*100
  thickness_sum[i,6] <- cv
}


# t-test and f-test
res_length <- matrix(0, nrow = 6, ncol = 2)

for(i in 1:nrow(res_length)){
  a <- subset(error_data2, Artefact == i)
  a1 <- a$Length
  b <- subset(comparative_data2, Artefact == i)
  b1 <- b$Length
  t <- t.test(a1,b1)
  f <- var.test(a1,b1)
  p1 <- t$p.value
  p2 <- f$p.value
  res_length[i,1:2] <- c(p1, p2) 
}

colnames(res_length) <- c("T", "F")
rownames(res_length) <- seq(1,6,1)


res_width <- matrix(0, nrow = 6, ncol = 2)

for(i in 1:nrow(res_width)){
  a <- subset(error_data2, Artefact == i)
  a1 <- a$Width
  b <- subset(comparative_data2, Artefact == i)
  b1 <- b$Width
  t <- t.test(a1,b1)
  f <- var.test(a1,b1)
  p1 <- t$p.value
  p2 <- f$p.value
  res_width[i,1:2] <- c(p1, p2) 
}

colnames(res_width) <- c("T", "F")
rownames(res_width) <- seq(1,6,1)

res_thickness <- matrix(0, nrow = 6, ncol = 2)

for(i in 1:nrow(res_thickness)){
  a <- subset(error_data2, Artefact == i)
  a1 <- a$Thickness
  b <- subset(comparative_data2, Artefact == i)
  b1 <- b$Thickness
  t <- t.test(a1,b1)
  f <- var.test(a1,b1)
  p1 <- t$p.value
  p2 <- f$p.value
  res_thickness[i,1:2] <- c(p1, p2) 
}

colnames(res_thickness) <- c("T", "F")
rownames(res_thickness) <- seq(1,6,1)

final_res <- cbind(res_length, res_width, res_thickness) 
write.csv(final_res, "T_F_results_metric.csv")


# R calculation
compare_data2 <- read.csv("rearranged_comparative_data.csv") # rearranged data

rownames(compare_data2 ) <- compare_data2$X # Assign first column as row name rather than as a variable
compare_data2  <- compare_data2[,-1] # remove first column 

ICC(compare_data2) # we use ICC3

# R calculation

output <- matrix(0, nrow = ncol(compare_data2), ncol = ncol(compare_data2 ))
for(i in 1:ncol(compare_data2 )){
  for(m in 1:ncol(compare_data2)){
    a <- cbind(compare_data2[,i],compare_data2[,m])
    b <- var(c(a[,1], a[,2]), use = "all.obs") # sample variance of the whole of matrix 'a'
    b <- b[1]
    c <- (a[,1])-(a[,2])
    d <- var(c)
    output[i,m] = b/(b+d)
  }
}

rownames(output) <- colnames(compare_data2)
colnames(output) <- colnames(compare_data2)

range(output)
write.csv(output, "R_compare_scores_metric.csv")

##### GMM analysis ####

# Load data 
tpslines <- import_tps("error_data.tps")

saveRDS(tpslines, file =  "ER_tpslines.rds")
saveRDS(error_data, file =  "ER_database.rds")

# Create ouline data and normalisation
tpsdata <- import("ER_tpslines.rds") 
database <- import("ER_database.rds")   

shape <- Out(tpsdata$coo, fac = database)
names(shape) <- database$Photo_ID

panel(shape, main = "Outline data") # Visualization of points in their original orientation

shapenorm <- shape %>% 
  coo_centre() %>% 
  coo_scale() %>% 
  coo_alignxax()%>% 
  coo_rotate(pi/2) %>%
  coo_slidedirection("right") %>% 
  coo_close() 

panel(shapenorm, main = "Normalised outline data") # Visualization of points in their original orientation
stack(shapenorm)

saveRDS(shapenorm, file = "Normalised_outlinesER.rds") # save normalised and transformed landmarks

# RUN FROM HERE USING RDS FILE - Reload normalised data and display according to tool type

ER_outlines <- import("Normalised_outlinesER.rds")
database <- import("ER_database.rds")  

ER_outlines$fac$Artefact <- as.factor(ER_outlines$fac$Artefact)
panel(ER_outlines, main = "Error data", fac = "Artefact")

## Calculate mean Artefact of points representing outlines
meanpoints <- matrix(0, nrow = length(ER_outlines), ncol = 1)
for(i in 1:length(ER_outlines)){
  artefact <- ER_outlines[i]
  nlandmarks <- length(unlist(artefact))/2
  meanpoints[i,] <- nlandmarks
}

mean(meanpoints)

# EFA
calibrate_harmonicpower_efourier(ER_outlines, nb.h = 20, plot = FALSE) # 8 harmonics

efashape <- efourier(ER_outlines, nb.h = 8, norm = FALSE)

coo_oscillo(ER_outlines[1], "efourier") # figure demonstrating EFA

# PCA

pcashape <- PCA(efashape) 


# Scree and PCA contrib plots

plot.new()
gg <- PCcontrib(pcashape, nax = 1:3, plot = FALSE)
gg$gg + 
  geom_polygon(fill="gray", col="black") 

p1 <- scree_plot(pcashape, nax =1:10) # PC1-3 gives over 95% of cum variance in the data, PC1 = 79% of variance
p1 <- p1  + theme_minimal()
p1

# PCA scatterplots
tidy.table <- cbind(as.tibble(pcashape$x[,1:4]), as.tibble(database))
tidy.table$Tool <- as.factor(tidy.table$Artefact)

write.csv(tidy.table, "PC_scores.csv")

a <- ggplot(tidy.table, aes(PC1, PC2, colour = Tool)) + 
  geom_point(size = 3) +
  stat_ellipse() +
  scale_colour_manual(values = c("#44AA99","#AA4499", "#6699CC", "#999933", "#888888", "#000000")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray")+
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgray")+
  theme_pubr()+
  theme(legend.position = "none")

b <- ggplot(tidy.table, aes(PC1, PC3, colour = Tool)) + 
  geom_point(size = 3) +
  stat_ellipse() +
  scale_colour_manual(values = c("#44AA99","#AA4499", "#6699CC", "#999933", "#888888", "#000000")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray")+
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgray")+
  theme_pubr()

c <- ggplot(tidy.table, aes(PC2, PC3, colour = Tool)) + 
  geom_point(size = 3) +
  stat_ellipse() +
  scale_colour_manual(values = c("#44AA99","#AA4499", "#6699CC", "#999933", "#888888", "#000000")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray")+
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgray")+
  theme_pubr()+
  theme(legend.position = "none")


d <- ggboxplot(tidy.table, "Tool", "PC1", fill = c("#44AA99","#AA4499", "#6699CC", "#999933", "#888888", "#000000"),
               bxp.errorbar = TRUE)+ theme(legend.position = "none")

e <- ggboxplot(tidy.table, "Tool", "PC2", fill = c("#44AA99","#AA4499", "#6699CC", "#999933", "#888888", "#000000"),
               bxp.errorbar = TRUE)+  theme(legend.position = "none")

f<- ggboxplot(tidy.table, "Tool", "PC3", fill = c("#44AA99","#AA4499", "#6699CC", "#999933", "#888888", "#000000"),
              bxp.errorbar = TRUE) +  theme(legend.position = "none")

(a|b|c)/(d|e|f) + 
  plot_annotation(theme = theme(plot.title = element_text(size = 20))) 

# Summary statistics

PC1_sum <- tidy.table %>%
  group_by(Artefact) %>%
  get_summary_stats(PC1, type = "mean_sd")

PC1_sum$Variance <- PC1_sum$sd^2


PC2_sum <- tidy.table %>%
  group_by(Artefact) %>%
  get_summary_stats(PC2, type = "mean_sd")

PC2_sum$Variance <- PC2_sum$sd^2

PC3_sum <- tidy.table %>%
  group_by(Artefact) %>%
  get_summary_stats(PC3, type = "mean_sd")

PC3_sum$Variance <- PC3_sum$sd^2

# ICC

rearranged_pca_data <- read.csv("rearranged_pca_data.csv")
rownames(rearranged_pca_data) <- rearranged_pca_data$X # Assign first column as row name rather than as a variable
rearranged_pca_data <- rearranged_pca_data[,-1] # remove first column 


ICC(rearranged_pca_data)

# R calculation

output <- matrix(0, nrow = ncol(rearranged_pca_data), ncol = ncol(rearranged_pca_data))
for(i in 1:ncol(rearranged_pca_data)){
  for(m in 1:ncol(rearranged_pca_data)){
    a <- cbind(rearranged_pca_data[,i],rearranged_pca_data[,m])
    b <- var(c(a[,1], a[,2]), y = NULL, use = "all") # sample variance of the whole of matrix 'a'
    c <- (a[,1])-(a[,2])
    d <- var(c)
    output[i,m] = b/(b+d)
  }
}

rownames(output) <- colnames(rearranged_pca_data)
colnames(output) <- colnames(rearranged_pca_data)
write.csv(output, "R_scores_PCA.csv")

## LDA

dashape90 <- LDA(pcashape, ~Museum, prior = c(6,6,6,6,6,6)/nrow(tidy.table), retain = 0.90, cv = TRUE) # 95% cum var PC scores
dashape90$CV.correct
dashape90$CV.ce

## Tukey HSD

res <- aov(PC1~Artefact, tidy.table)
res1 <- aov(PC2~Artefact, tidy.table)
res2 <- aov(PC3~Artefact, tidy.table)

TukeyHSD(res)
TukeyHSD(res1)
TukeyHSD(res2)


## COMPARATIVE ## 

# Load data
tpslines <- import_tps("comp_error_data.tps")
database <- read.csv("comparative_data.csv")
saveRDS(tpslines, file =  "CER_tpslines.rds")
saveRDS(database, file =  "CER_database.rds")

# Create ouline object and normalisation
tpsdata <- import("CER_tpslines.rds") 
database2 <- import("CER_database.rds")   

shape <- Out(tpsdata$coo, fac = database2)
names(shape) <- database2$Photo_ID

panel(shape, main = "Outline data") # Visualization of points in their original orientation

shapenorm <- shape %>% 
  coo_centre() %>% 
  coo_scale() %>% 
  coo_alignxax()%>% 
  coo_rotate(pi/2) %>%
  coo_slidedirection("right") %>% 
  coo_close() 

panel(shapenorm, main = "Normalised outline data") # Visualization of points in their original orientation
stack(shapenorm)

saveRDS(shapenorm, file = "Normalised_outlinesCER.rds") # save normalised and transformed landmarks

# RUN FROM HERE USING RDS FILE - Reload normalised data and display according to tool type

CER_outlines <- import("Normalised_outlinesCER.rds")
CER_database <- import("CER_database.rds")  

CER_outlines$fac$Artefact <- as.factor(CER_outlines$fac$Artefact)
panel(CER_outlines, main = "Normalised outline data", fac = "Artefact") # Visualization of points in their original orientation

# Combine with multiple observer data

combined <- Momocs:: combine(CER_outlines, ER_outlines)
panel(combined, main = "Outline data") # Visualization of points in their original orientation

database3 <- rbind(CER_database, error_data)

# EFA
calibrate_harmonicpower_efourier(combined, nb.h = 20, plot = FALSE) # 8 harmonics

efashape2 <- efourier(combined, nb.h = 8, norm = FALSE)

# PCA  

pcashape2 <- PCA(efashape2) 

# Scree and PCA contrib plots 

plot.new()
gg <- PCcontrib(pcashape2, nax = 1:3, plot = FALSE)
gg$gg + 
  geom_polygon(fill="gray", col="black") 

p1 <- scree_plot(pcashape2, nax =1:10) # PC1-3 gives over 95% of cum variance in the data, PC1 = 79% of variance
p1 <- p1  + theme_minimal()
p1


# PCA plots
tidy.table <- cbind(as.tibble(pcashape2$x[,1:3]), as.tibble(database3))
tidy.table$Tool <- as.factor(tidy.table$Artefact)

CodeID <- c(rep("Single", 36), rep("Multiple", 36))
tidy.table$Observer <- CodeID
tidy.table$Observer<- factor(tidy.table$Observer)

write.csv(tidy.table, "comparative_PC_scores.csv")

a <- ggplot(tidy.table, aes(PC1, PC2, colour = Observer, pch = Tool)) + 
  geom_point(size = 3) +
  scale_colour_manual(values = c("#44AA99","#AA4499")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray")+
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgray")+
  theme_pubr() +
  theme(legend.position = "none")

b <- ggplot(tidy.table, aes(PC1, PC3, colour = Observer, pch = Tool)) + 
  geom_point(size = 3) +
  scale_colour_manual(values = c("#44AA99","#AA4499")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray")+
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgray")+
  theme_pubr()

c <- ggplot(tidy.table, aes(PC2, PC3, colour = Observer, pch = Tool)) + 
  geom_point(size = 3) +
  scale_colour_manual(values = c("#44AA99","#AA4499")) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray")+
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgray")+
  theme_pubr()+
  theme(legend.position = "none")

d <- ggboxplot(tidy.table, "Tool", "PC1", fill = 'Observer', palette = c("#44AA99","#AA4499"),
                     bxp.errorbar = TRUE)+ theme(legend.position = "none")

e <- ggboxplot(tidy.table, "Tool", "PC2", fill  = 'Observer', palette = c("#44AA99","#AA4499"),
                     bxp.errorbar = TRUE)+  theme(legend.position = "none")

f<- ggboxplot(tidy.table, "Tool", "PC3", fill  = 'Observer', palette =  c("#44AA99","#AA4499"),
                     bxp.errorbar = TRUE) +  theme(legend.position = "none")

(a|b|c)/(d|e|f) + 
  plot_annotation(theme = theme(plot.title = element_text(size = 20))) 

# Summary statistics

comp_PC1_sum <- tidy.table %>%
  group_by(Artefact,Observer) %>%
  get_summary_stats(PC1, type = "mean_sd")

comp_PC1_sum$Variance <- comp_PC1_sum$sd^2

comp_PC2_sum <- tidy.table%>%
  group_by(Artefact,Observer) %>%
  get_summary_stats(PC2, type = "mean_sd")

comp_PC2_sum$Variance <- comp_PC2_sum$sd^2

comp_PC3_sum <- tidy.table %>%
  group_by(Artefact,Observer) %>%
  get_summary_stats(PC3, type = "mean_sd")

comp_PC3_sum$Variance <- comp_PC3_sum$sd^2


PC3_m_bxp <- ggboxplot(subset(tidy.table, Observer == "Multiple"),
                       x = "Tool",
                       y = "PC3", 
                       fill = c("#44AA99","#AA4499", "#6699CC", "#999933", "#888888", "#000000"), 
                       bxp.errorbar = TRUE)

PC1_s_bxp <- ggboxplot(subset(tidy.table, Observer == "Single"),
                       x = "Tool",
                       y = "PC1", 
                       fill = c("#44AA99","#AA4499", "#6699CC", "#999933", "#888888", "#000000"), 
                       bxp.errorbar = TRUE)

PC2_s_bxp <- ggboxplot(subset(tidy.table, Observer == "Single"),
                       x = "Tool",
                       y = "PC2", 
                       fill = c("#44AA99","#AA4499", "#6699CC", "#999933", "#888888", "#000000"), 
                       bxp.errorbar = TRUE)

PC3_s_bxp <- ggboxplot(subset(tidy.table, Observer == "Single"),
                       x = "Tool",
                       y = "PC3", 
                       fill = c("#44AA99","#AA4499", "#6699CC", "#999933", "#888888", "#000000"), 
                       bxp.errorbar = TRUE)



(PC1_m_bxp| PC1_s_bxp) /width_bxp/thickness_bxp +
  plot_annotation(tag_levels = "a",
                  theme = theme(plot.title = element_text(size = 20))) 


# t-test and f-test
res_PC1 <- matrix(0, nrow = 6, ncol = 2)

for(i in 1:nrow(res_PC1)){
  a <- subset(tidy.table, Artefact == i & Observer == "Single")
  a1 <- a$PC1
  b <- subset(tidy.table, Artefact == i & Observer == "Multiple")
  b1 <- b$PC1
  t <- t.test(a1,b1)
  f <- var.test(a1,b1)
  p1 <- t$p.value
  p2 <- f$p.value
  res_PC1[i,1:2] <- c(p1, p2) 
}

colnames(res_PC1) <- c("T", "F")
rownames(res_PC1) <- seq(1,6,1)


res_PC2 <- matrix(0, nrow = 6, ncol = 2)

for(i in 1:nrow(res_PC2)){
  a <- subset(tidy.table, Artefact == i & Observer == "Single")
  a1 <- a$PC2
  b <- subset(tidy.table, Artefact == i & Observer == "Multiple")
  b1 <- b$PC2
  t <- t.test(a1,b1)
  f <- var.test(a1,b1)
  p1 <- t$p.value
  p2 <- f$p.value
  res_PC2[i,1:2] <- c(p1, p2) 
}

colnames(res_PC2) <- c("T", "F")
rownames(res_PC2) <- seq(1,6,1)

res_PC3 <- matrix(0, nrow = 6, ncol = 2)

for(i in 1:nrow(res_PC3)){
  a <- subset(tidy.table, Artefact == i & Observer == "Single")
  a1 <- a$PC3
  b <- subset(tidy.table, Artefact == i & Observer == "Multiple")
  b1 <- b$PC3
  t <- t.test(a1,b1)
  f <- var.test(a1,b1)
  p1 <- t$p.value
  p2 <- f$p.value
  res_PC3[i,1:2] <- c(p1, p2) 
}

colnames(res_PC3) <- c("T", "F")
rownames(res_PC3) <- seq(1,6,1)

final_res <- cbind(res_PC1, res_PC2, res_PC3) 
write.csv(final_res, "T_F_results_PCA.csv")



# R calculation for single observer

rearranged_pca_data <- read.csv("rearranged_comparative_PC_scores.csv")
rownames(rearranged_pca_data) <- rearranged_pca_data$X # Assign first column as row name rather than as a variable
rearranged_pca_data <- rearranged_pca_data[,-1] # remove first column 


ICC(rearranged_pca_data)

output <- matrix(0, nrow = ncol(rearranged_pca_data), ncol = ncol(rearranged_pca_data))
for(i in 1:ncol(rearranged_pca_data)){
  for(m in 1:ncol(rearranged_pca_data)){
    a <- cbind(rearranged_pca_data[,i],rearranged_pca_data[,m])
    b <- var(c(a[,1], a[,2]), y = NULL, use = "all") # sample variance of the whole of matrix 'a'
    c <- (a[,1])-(a[,2])
    d <- var(c)
    output[i,m] = b/(b+d)
  }
}

rownames(output) <- colnames(rearranged_pca_data)
colnames(output) <- colnames(rearranged_pca_data)
write.csv(output, "comparative_R_scores_PCA.csv")


