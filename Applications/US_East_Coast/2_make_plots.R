rm(list = ls())
library(perla)
library(ggplot2)
library(reshape2)
library(sp)
library(latex2exp)

saving_directory <- "/mnt/callisto/Sottosanti/perla/Applications/US_East_Coast/Results"


# DIC3 --------------------------------------------------------------

gender <- c("males","females")
infoCr <- data.frame(matrix(NA,0,4))
colnames(infoCr) <- c("DIC","Cluster","Penalisation","Gender")
for(i in gender){
  dir_name <- paste(saving_directory,i,sep="/")
  for(j in 1:length(dir(dir_name))){
    load(paste(dir_name,"/",dir(dir_name)[j],sep=""))
    infoCr <- rbind(infoCr,
                    data.frame(DIC = results$DIC3, Cluster = dim(results$Z)[2], Penalisation = paste(results$mean.penalty, collapse = "_"), Gender = i))
  }
}
infoCr$Penalisation <- factor(infoCr$Penalisation)
levels(infoCr$Penalisation) <- c("(c)","(c,d)","(cd)","(d)","(d,cd)")
infoCr$Penalisation <- factor(infoCr$Penalisation, levels=c("(c)","(d)","(c,d)","(cd)","(d,cd)"))
infoCr$Gender <- factor(infoCr$Gender)
levels(infoCr$Gender) <- c("Females","Males")

v <- 16
ggplot(aes(x=factor(Cluster), y=DIC, group=Penalisation),data=infoCr[which(infoCr$Gender=="Males"),])+
  geom_point(aes(shape=Penalisation, color=Penalisation), size=2)+
  geom_line(aes(color=Penalisation))+
  theme_bw()+
  xlab("Number of clusters")+ylab(expression(~DIC[3]))+
  theme(axis.text = element_text(size=v),
        axis.title = element_text(size=v),
        legend.title = element_text(size=v),
        legend.text = element_text(size=v),
        strip.text = element_text(size=v),
        legend.key.size = unit(1, 'cm'))+
  scale_color_manual(values=c("#CC79A7","#009E73","#E69F00","#0072B2","#F0E442"))
ggsave(file = "~/perla/Applications/US_East_Coast/Graphs/males/DIC.pdf", width = 12, height = 9)

v <- 16
ggplot(aes(x=factor(Cluster), y=DIC, group=Penalisation),data=infoCr[which(infoCr$Gender=="Females"),])+
  geom_point(aes(shape=Penalisation, color=Penalisation), size=2)+
  geom_line(aes(color=Penalisation))+
  theme_bw()+
  xlab("Number of clusters")+ylab(expression(~DIC[3]))+
  theme(axis.text = element_text(size=v),
        axis.title = element_text(size=v),
        legend.title = element_text(size=v),
        legend.text = element_text(size=v),
        strip.text = element_text(size=v),
        legend.key.size = unit(1, 'cm'))+
  scale_color_manual(values=c("#CC79A7","#009E73","#E69F00","#0072B2","#F0E442"))
ggsave(file = "~/perla/Applications/US_East_Coast/Graphs/females/DIC.pdf", width = 12, height = 9)



# MAP OF CLUSTERS ---------------------------------------------------------

# ---MALES
rm(list = ls())
saving_directory <- "/mnt/callisto/Sottosanti/perla/Applications/US_East_Coast/Results"
colors <- c("#0072B2","#E69F00","#009E73","#CC79A7")
gender <- "males"
select_pen <- "c"
select_K <- 3
load(paste(saving_directory,"/",gender,"/Results_",select_K,"_",select_pen,".RData",sep=""))
map <- readRDS(paste("~/perla/Applications/US_East_Coast/Data/",gender,".RDS",sep=""))
map@data$cluster <-  relab$relabelling$clusters[2,]
colore <- rep(NA, nrow(map@data))
for(i in 1:length(colors)) colore[map@data$cluster == i] <- colors[i]
map@data$colore <- colore

pdf(file = paste("~/perla/Applications/US_East_Coast/Graphs/",gender,"/Map_",select_K,"_",select_pen,".pdf",sep=""), width = 12, height = 9)
plot(map, col=map@data$colore)
if(gender == "males") title("Males")
if(gender == "females") title("Females")
dev.off()


# ---FEMALES
rm(list = ls())
saving_directory <- "/mnt/callisto/Sottosanti/perla/Applications/US_East_Coast/Results"
colors <- c("#0072B2","#E69F00","#009E73","#CC79A7")
gender <- "females"
select_pen <- "d_cd"
select_K <- 4
load(paste(saving_directory,"/",gender,"/Results_",select_K,"_",select_pen,".RData",sep=""))
map <- readRDS(paste("~/perla/Applications/US_East_Coast/Data/",gender,".RDS",sep=""))
map@data$cluster <-  relab$relabelling$clusters[2,]
colore <- rep(NA, nrow(map@data))
for(i in 1:length(colors)) colore[map@data$cluster == i] <- colors[i]
map@data$colore <- colore

pdf(file = paste("~/perla/Applications/US_East_Coast/Graphs/",gender,"/Map_",select_K,"_",select_pen,".pdf",sep=""), width = 12, height = 9)
plot(map, col=map@data$colore)
if(gender == "males") title("Males")
if(gender == "females") title("Females")
dev.off()


# DENSITY PLOT ------------------------------------------------------------

# ---MALES
rm(list = ls())
saving_directory <- "/mnt/callisto/Sottosanti/perla/Applications/US_East_Coast/Results"
colors <- c("#0072B2","#E69F00","#009E73","#CC79A7")
gender <- "males"
select_pen <- "c"
select_K <- 3
load(paste(saving_directory,"/",gender,"/Results_",select_K,"_",select_pen,".RData",sep=""))
diseases <- c("Circulatory","Respiratory","Tumours")
dd <- data.frame(matrix(NA,0,3))
names(dd) <- c("value","Cluster","Disease")
for(k in 1:select_K){
  for(j in 1:length(diseases)){
    dd <- rbind(dd,
                data.frame(value = exp(relab$results$ECR$M[,k,j]),
                        Cluster = rep(k, length(relab$results$ECR$M[,k,j])),
                        Disease = rep(diseases[j],length(relab$results$ECR$M[,k,j])))
    )
  }
}
dd$Cluster <- factor(dd$Cluster)
dd$Disease <- factor(dd$Disease)

v <- 16
ggplot(aes(x=value, fill=Cluster),data=dd)+geom_density(linewidth=0.7, alpha=0.7)+theme_bw()+
  ylab("Density")+xlab(expression(exp(mu[kj])))+facet_wrap(~Disease)+
  theme(axis.text = element_text(size=v),
        axis.title = element_text(size=v),
        legend.title = element_text(size=v),
        legend.text = element_text(size=v),
        strip.text = element_text(size=v),
        legend.key.size = unit(1, 'cm'))+
  scale_fill_manual(values=colors)
ggsave(file = paste("~/perla/Applications/US_East_Coast/Graphs/",gender,"/Density_",select_K,"_",select_pen,".pdf",sep=""), width = 12, height = 9)

# ---FEMALES
rm(list = ls())
saving_directory <- "/mnt/callisto/Sottosanti/perla/Applications/US_East_Coast/Results"
colors <- c("#0072B2","#E69F00","#009E73","#CC79A7")
gender <- "females"
select_pen <- "d_cd"
select_K <- 4
load(paste(saving_directory,"/",gender,"/Results_",select_K,"_",select_pen,".RData",sep=""))
diseases <- c("Circulatory","Respiratory","Tumours")
dd <- data.frame(matrix(NA,0,3))
names(dd) <- c("value","Cluster","Disease")
for(k in 1:select_K){
  for(j in 1:length(diseases)){
    dd <- rbind(dd,
                data.frame(value = exp(relab$results$PRA$M[,k,j]),
                           Cluster = rep(k, length(relab$results$PRA$M[,k,j])),
                           Disease = rep(diseases[j],length(relab$results$PRA$M[,k,j])))
    )
  }
}
dd$Cluster <- factor(dd$Cluster)
dd$Disease <- factor(dd$Disease)

v <- 16
ggplot(aes(x=value, fill=Cluster),data=dd)+geom_density(linewidth=0.7, alpha=0.7)+theme_bw()+
  ylab("Density")+xlab(expression(exp(mu[kj])))+facet_wrap(~Disease)+
  theme(axis.text = element_text(size=v),
        axis.title = element_text(size=v),
        legend.title = element_text(size=v),
        legend.text = element_text(size=v),
        strip.text = element_text(size=v),
        legend.key.size = unit(1, 'cm'))+
  scale_fill_manual(values=colors)
ggsave(file = paste("~/perla/Applications/US_East_Coast/Graphs/",gender,"/Density_",select_K,"_",select_pen,".pdf",sep=""), width = 12, height = 9)




# POSTERIOR SHRINKAGE FACTORS ---------------------------------------------

# ---FEMALES
rm(list = ls())
saving_directory <- "/mnt/callisto/Sottosanti/perla/Applications/US_East_Coast/Results"
colors <- c("#0072B2","#E69F00","#009E73","#CC79A7")
gender <- "females"
select_pen <- "d_cd"
select_K <- 4
load(paste(saving_directory,"/",gender,"/Results_",select_K,"_",select_pen,".RData",sep=""))

# the parametrisation for the shrinkage factors is the inverse of the one reported in the article.
# For instance, to access the posterior distribution of the disease shrinkage factors, one should look at 1/results$shrinkage.parameters$Zeta.d

sum_zeta_d <- rowSums(1/results$shrinkage.parameters$Zeta.d)
zeta_d_norm <- results$shrinkage.parameters$Zeta.d
for(i in 1:nrow(zeta_d_norm)){
  zeta_d_norm[i,] <- (1/results$shrinkage.parameters$Zeta.d[i,])/sum_zeta_d[i]
}

d <- ncol(zeta_d_norm)
zetaF <- matrix(apply(zeta_d_norm,2,median), select_K, d, byrow=T)
colnames(zetaF) <- c("Circulatory\n(j = 1)","Respiratory\n(j = 2)","Tumours\n(j = 3)")
rownames(zetaF) <- 1:select_K
longDataB <- melt(zetaF)
v <- 16
Shrinkage_plots <- list()
Shrinkage_plots[[1]] <- ggplot(longDataB, aes(x = Var2, y = Var1)) +
  geom_tile(aes(fill = value), color = NA) +  # use geom_tile
  scale_fill_continuous(low="orange", high="blue4",limits=c(0,0.6), breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5))+
  labs(x="Cause of death", y="Cluster") +
  theme_bw() + guides(fill=guide_legend(title=TeX(r'($\zeta_j/\sum_l{\zeta_l}$)') ))+
  theme(axis.text = element_text(size=v),
        axis.title = element_text(size=v),
        legend.title = element_text(size=v),
        plot.title = element_text(size=v, hjust = 0.5),
        legend.text = element_text(size=v),
        strip.text = element_text(size=v),
        legend.key.size = unit(1, 'cm'))
Shrinkage_plots[[1]]
ggsave(file = paste("~/perla/Applications/US_East_Coast/Graphs/",gender,"/ShrinkFactorD_",select_K,"_",select_pen,".pdf",sep=""), width = 12, height = 9)


C <- matrix(apply(1/results$shrinkage.parameters$Zeta.cd,2,median), select_K, 3, byrow=F)
colnames(C) <- c("Circulatory\n(j = 1)","Respiratory\n(j = 2)","Tumours\n(j = 3)")
rownames(C) <- 1:select_K
longDataC <- melt(C)
longDataC$Var1 <- paste("k =",longDataC$Var1)
Shrinkage_plots[[2]] <- ggplot(longDataC, aes(x = Var2, y = Var1)) +
  geom_tile(aes(fill = value), color = NA) +  # use geom_tile
  scale_fill_continuous(low = "orange", high = "blue4",
                        limits = c(0, 13),
                        breaks = c(0, 2.5, 5, 7.5, 10, 12.5)) +
  labs(x = "Cause of death", y = "Cluster") +
  theme_bw() +
  guides(fill = guide_legend(title = TeX(r'($\gamma_{kj}$)'))) +
  theme(
    axis.text = element_text(size = v),
    axis.title = element_text(size = v),
    legend.title = element_text(size = v),
    plot.title = element_text(size = v, hjust = 0.5),
    legend.text = element_text(size = v),
    strip.text = element_text(size = v),
    legend.key.size = unit(1, 'cm')
  )
Shrinkage_plots[[2]]
ggsave(file = paste("~/perla/Applications/US_East_Coast/Graphs/",gender,"/ShrinkFactorCD_",select_K,"_",select_pen,".pdf",sep=""), width = 12, height = 9)

(Shrinkage_plots[[1]]|Shrinkage_plots[[2]])+
  plot_annotation(title = "Females", theme = theme(title = element_text(size = 16),
                                                   plot.title = element_text(hjust = .5)))
ggsave(file = paste("~/perla/Applications/US_East_Coast/Graphs/",gender,"/ShrinkFactors_",select_K,"_",select_pen,".pdf",sep=""),
       width = 24, height = 9)






# Probability maps with cities --------------------------------------------

library(sf)
library(coda)
library(MCMCglmm)
library(tidyr)

rm(list = ls())
saving_directory <- "/mnt/callisto/Sottosanti/perla/Applications/US_East_Coast/Results"
colors <- c("#0072B2","#E69F00","#009E73","#CC79A7")
gender <- "females"
select_pen <- "d_cd"
select_K <- 4
load(paste(saving_directory,"/",gender,"/Results_",select_K,"_",select_pen,".RData",sep=""))
map <- readRDS(paste("~/perla/Applications/US_East_Coast/Data/",gender,".RDS",sep=""))
map@data <- as.data.frame(t(apply(relab$results$ECR$Responsabilities, c(2,3), function(x) posterior.mode(mcmc(x)))))
colore <- rep(NA, nrow(map@data))
spdf_sf <- st_as_sf(map)

# Pivot longer to plot multiple variables
map_long <- spdf_sf %>%
  pivot_longer(
    cols = setdiff(names(spdf_sf), attr(spdf_sf, "sf_column")),
    names_to = "Cluster",
    values_to = "value"
  )


library(ggrepel)  # For better label placement

# Major East Coast cities with coordinates (lat, long)
east_coast_cities <- maps::us.cities[,names(maps::us.cities) %in% c("name","long","lat")]

east_coast_cities <- east_coast_cities[east_coast_cities$long > -80 & east_coast_cities$lat > 34,]
east_coast_cities <- st_as_sf(east_coast_cities, coords = c("long", "lat"), crs = 4326)  # Convert to sf object



# Create faceted plot
ggplot(map_long[map_long$Cluster == "V2",]) +
  geom_sf(aes(fill = value)) + # Remove borders with color=NA
  scico::scale_fill_scico(
                           palette = "navia")+
  facet_wrap(~Cluster, ncol = 4) +       # Arrange in grid
  theme_bw() +
  labs()+
  geom_sf(data = east_coast_cities, color = "red", size = 1) +
  # Add city labels with repelling
  geom_label_repel(
    data = east_coast_cities,
    aes(label = name, geometry = geometry),
    stat = "sf_coordinates",
    size = 2,
    min.segment.length = 1,
    segment.color = "gray50"
  )



