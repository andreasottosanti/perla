library(perla)
library(ggplot2)
library(reshape2)
library(sp)


# DIC3 --------------------------------------------------------------

gender <- c("males","females")
infoCr <- data.frame(matrix(NA,0,4))
colnames(infoCr) <- c("DIC","Cluster","Penalisation","Gender")
for(i in gender){
  dir_name <- paste("~/perla/Applications/US_East_Coast/Results/",i,sep="")
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
colors <- c("#0072B2","#E69F00","#009E73","#CC79A7")
gender <- "males"
select_pen <- "c"
select_K <- 3
load(paste("~/perla/Applications/US_East_Coast/Results/",gender,"/Results_",select_K,"_",select_pen,".RData",sep=""))
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
colors <- c("#0072B2","#E69F00","#009E73","#CC79A7")
gender <- "females"
select_pen <- "d_cd"
select_K <- 4
load(paste("~/perla/Applications/US_East_Coast/Results/",gender,"/Results_",select_K,"_",select_pen,".RData",sep=""))
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
colors <- c("#0072B2","#E69F00","#009E73","#CC79A7")
gender <- "males"
select_pen <- "c"
select_K <- 3
load(paste("~/perla/Applications/US_East_Coast/Results/",gender,"/Results_",select_K,"_",select_pen,".RData",sep=""))
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
colors <- c("#0072B2","#E69F00","#009E73","#CC79A7")
gender <- "females"
select_pen <- "d_cd"
select_K <- 4
load(paste("~/perla/Applications/US_East_Coast/Results/",gender,"/Results_",select_K,"_",select_pen,".RData",sep=""))
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




# POSTERIOR SHRINKAGE FACTORS ---------------------------------------------

# ---FEMALES
rm(list = ls())
colors <- c("#0072B2","#E69F00","#009E73","#CC79A7")
gender <- "females"
select_pen <- "d_cd"
select_K <- 4
load(paste("~/perla/Applications/US_East_Coast/Results/",gender,"/Results_",select_K,"_",select_pen,".RData",sep=""))


sum_zeta_d <- rowSums(1/results$shrinkage.parameters$Zeta.d)
zeta_d_norm <- results$shrinkage.parameters$Zeta.d
for(i in 1:nrow(zeta_d_norm)){
  zeta_d_norm[i,] <- (1/results$shrinkage.parameters$Zeta.d[i,])/sum_zeta_d[i]
}




################################ GGPLOT #############################################################################


#females
d <- ncol(zeta_d_norm)
zetaF <- matrix(apply(zeta_d_norm,2,median), select_K, d, byrow=T)
zetaFgg <- matrix(NA, 12, 3)
zetaFgg[,1] <- zetaF
zetaFgg[,2] <- rep(1:4,3)
zetaFgg[,3] <- c(rep(1,4),rep(2,4),rep(3,4))
zetaFgg <- as.data.frame(zetaFgg)
colnames(zetaFgg) <- c("Value","Cluster","Disease")
zetaFgg
B <- zetaF
colnames(B) <- c("Circulatory","Respiratory","Tumours")#c("Disease 1","Disease 2 ","Disease 3")
rownames(B) <- 1:4#c("Cluster 1","Cluster 2","Cluster 3","Cluster 4")
longDataB<-melt(B)
pdf("C:/Users/enric/OneDrive/Desktop/ArticoloFuturo/Boxplots/Plot_GG_shrink_USAzetaF202.pdf", width = 12, height = 9)
ggplot(longDataB, aes(x = Var2, y = Var1)) +
  geom_raster(aes(fill=value)) +
  #  scale_fill_gradient(low="orange", high="blue4")+
  scale_fill_continuous(low="orange", high="blue4",limits=c(0,0.5), breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5))+
  labs(x="Cause of death", y="Cluster", title="Females") +
  theme_bw() + guides(fill=guide_legend(title=TeX(r'($\zeta_j/\sum_l{\zeta_l}$)') ))+
  theme(axis.text = element_text(size=v),
        axis.title = element_text(size=v),
        legend.title = element_text(size=v),
        plot.title = element_text(size=v, hjust = 0.5),
        legend.text = element_text(size=v),
        strip.text = element_text(size=v),
        legend.key.size = unit(1, 'cm'))
dev.off()


gammaF <- matrix(apply(try,2,median),4,3,byrow=T)
C <- gammaF
colnames(C) <- c("Circulatory","Respiratory","Tumours")#c("Disease 1","Disease 2 ","Disease 3")
rownames(C) <- 1:4#c("Cluster 1","Cluster 2","Cluster 3","Cluster 4")
longDataC<-melt(C)
pdf("C:/Users/enric/OneDrive/Desktop/ArticoloFuturo/Boxplots/Plot_GG_shrink_USAgammaF202.pdf", width = 12, height = 9)
ggplot(longDataC, aes(x = Var2, y = Var1)) +
  geom_raster(aes(fill=value)) +
  #  scale_fill_gradient(low="orange", high="blue4")+
  scale_fill_continuous(low="orange", high="blue4",limits=c(0,13), breaks = c(0, 2.5, 5, 7.5, 10, 12.5))+
  labs(x="Cause of death", y="Cluster", title="Females") +
  theme_bw() + guides(fill=guide_legend(title=TeX(r'($\gamma_{kj}$)') )) +
  theme(axis.text = element_text(size=v),
        axis.title = element_text(size=v),
        legend.title = element_text(size=v),
        plot.title = element_text(size=v, hjust = 0.5),
        legend.text = element_text(size=v),
        strip.text = element_text(size=v),
        legend.key.size = unit(1, 'cm'))
dev.off()


#males
zetaM <- matrix(apply(zeta_j_t_tilde,2,median),3,3,byrow=T)
A <- zetaM
colnames(A) <- c("Circulatory","Respiratory","Tumours")#c("Disease 1","Disease 2 ","Disease 3")
rownames(A) <- 1:3#c("Cluster 1","Cluster 2","Cluster 3","Cluster 4")
longDataA<-melt(A)
pdf("C:/Users/enric/OneDrive/Desktop/ArticoloFuturo/Boxplots/Plot_GG_shrink_USAzetaM202.pdf", width = 12, height = 9)
ggplot(longDataA, aes(x = Var2, y = Var1)) +
  geom_raster(aes(fill=value)) +
  #  scale_fill_gradient(low="orange", high="blue4")+
  scale_fill_continuous(low="orange", high="blue4",limits=c(0,1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1))+
  labs(x="Cause of death", y="Cluster", title="Males") +
  theme_bw() + guides(fill=guide_legend(title=TeX(r'($\zeta_j/\sum_l{\zeta_l^}$)') ))+
  theme(axis.text = element_text(size=v),
        axis.title = element_text(size=v),
        legend.title = element_text(size=v),
        plot.title = element_text(size=v, hjust = 0.5),
        legend.text = element_text(size=v),
        strip.text = element_text(size=v),
        legend.key.size = unit(1, 'cm'))
dev.off()
