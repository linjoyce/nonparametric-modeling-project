library(meanShiftR)
library(RColorBrewer)
library(pheatmap)
library(Compositional)
library(MASS)
library(dplyr)
library(randomForest)
library(proxy)
library(parallel)
library(fpc)
library(aricode)

##########################
#      Input Data        #
##########################
setwd("/Users/joycelin/Documents/UW/Course/STAT527/Project")
set.seed(1108)
dt = read.csv("data/project-data.csv", header = F)
#sdt = dt[sample(nrow(dt)/10), ]


##########################
#      Visualization     #
##########################
jpeg(file="images/lineplot.jpeg", width=14, height=8, units="cm", res=300, pointsize=6)
par(mar = c(5.1, 4.5, 2.1, 2.1))
plot(0,0, xlim=c(1, ncol(dt)), ylim=c(min(dt), max(dt)), xlab = "Dimensions", ylab="Value",
     type="n", cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
for(i in 1:nrow(dt)){
  cat(i, "\t")
  par(new=T)
  plot(1:ncol(dt), dt[i,], ylim=c(min(dt), max(dt)), lwd = 0.5,
       xaxt="n", yaxt="n", xlab = "", ylab="", type="l", col=i)
}
dev.off()

jpeg(file="images/histogram.jpeg", width=12, height=8, units="cm", res=300, pointsize=6)
par(mar=c(2,2,1,1), mfrow = c(8,8))
for(i in 1:ncol(dt)){
  cat(i, "\t")
  hist(dt[,i], main = "")
}
dev.off()

##########################
#       Mean-Shift       #
##########################
# Bandwidth Selection
# Tuning of the bandwidth h of the kernel using the maximum likelihood cross validation.
h = mkde.tune(as.matrix(dt), c(0.1, 3));h

# Number of Neighbors Selection by WCSS
WCSS = NULL
for (nNbr in 1:30) {
  cat(nNbr, "\t")
  classification <- meanShift(as.matrix(dt), nNeighbors = nNbr,
                              bandwidth = rep(h$hopt, ncol(dt)))
  grp = classification$assignment
  
  group.data = cbind(grp, dt)
  group.mean  = as.data.frame(group.data %>% group_by(grp) %>% summarise_at(vars(colnames(dt)), list(name = mean)))[,-1]
  colnames(group.mean) = colnames(dt)
  
  wcss = 0
  for (i in 1:nrow(group.mean)) {
    wcss = wcss + sum(dist(dt[which(grp == i),], group.mean[i,]))
  }
  
  WCSS[nNbr] = wcss
}

jpeg(file="images/MS_wcss.jpeg", width=14, height=8, units="cm", res=300, pointsize=6)
par(mar = c(5.1, 4.5, 2.1, 2.1))
plot(1:30, WCSS, type = "b", xlab = "Number of Neighbors",
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
points(21, WCSS[21], pch = 16, col = "red")
dev.off()

# Mean-shift Classification
nNbr = 21
classification <- meanShift(as.matrix(dt), nNeighbors = nNbr,
                            bandwidth = rep(h$hopt, ncol(dt)))
grp = classification$assignment
length(unique(classification$assignment))

write.table(grp, file = "alg1_out.txt", sep = "\n", row.names = F, col.names = F)

# Heatmap
gdt = dt[order(grp), 1:10]

colnames(gdt) = 1:ncol(gdt)
brk = 50
color = colorRampPalette(c("red","white","blue"))(brk)
myBreaks <- c(seq(min(gdt), 0, length.out=ceiling(brk/2) + 1), 
              seq(max(gdt)/brk, max(gdt), length.out=floor(brk/2)))
cls_label = rep("", nrow(dt))
cls_label[c(0, cumsum(table(grp))[-nNbr]) + round((cumsum(table(grp)) - c(1, cumsum(table(grp))[-nNbr]+1))/2)] =
  paste("C", 1:nNbr, sep = "")
pheatmap(gdt, cluster_cols = F, cluster_row = F, color = color, breaks=myBreaks, fontsize=13.5,
         labels_row = cls_label, main="Heatmap", gaps_row = cumsum(table(grp)))


# Evaluate Stability
set.seed(1108)
sc.dt = dt[which(grp == which.min(table(grp))), ]
nosie_sigma = apply(sc.dt, 2, function(x) max(dist(x))*0.001)
noise.dt = matrix(NA, nrow(dt), ncol(dt))
for (i in 1:ncol(dt)) {
  noise.dt[,i] = rnorm(nrow(dt), sd = nosie_sigma[i])
}

perturb.dt = dt + noise.dt
write.table(perturb.dt, "data/pertubdata.csv", sep = ",", row.names = FALSE, col.names = FALSE)

jpeg(file="images/perturbdt_histogram.jpeg", width=12, height=8, units="cm", res=300, pointsize=6)
par(mar=c(2,2,1,1), mfrow = c(8,8))
for(i in 1:ncol(dt)){
  cat(i, "\t")
  hist(perturb.dt[,i], main = "")
}
dev.off()

perturb.h = mkde.tune(as.matrix(perturb.dt), c(0.1, 3));perturb.h

# Number of Neighbors Selection by WCSS
WCSS = NULL
for (nNbr in 1:30) {
  cat(nNbr, "\t")
  classification <- meanShift(as.matrix(perturb.dt), nNeighbors = nNbr,
                              bandwidth = rep(perturb.h$hopt, ncol(perturb.dt)))
  perturb.grp = classification$assignment
  
  group.data = cbind(perturb.grp, perturb.dt)
  group.mean  = as.data.frame(group.data %>% group_by(perturb.grp) %>% summarise_at(vars(colnames(perturb.dt)), list(name = mean)))[,-1]
  colnames(group.mean) = colnames(perturb.dt)
  
  wcss = 0
  for (i in 1:nrow(group.mean)) {
    wcss = wcss + sum(dist(perturb.dt[which(perturb.grp == i),], group.mean[i,]))
  }
  
  WCSS[nNbr] = wcss
}

jpeg(file="images/MS_wcss_perturb.jpeg", width=14, height=8, units="cm", res=300, pointsize=6)
par(mar = c(5.1, 4.5, 2.1, 2.1))
plot(1:30, WCSS, type = "b", xlab = "Number of Neighbors")
points(21, WCSS[21], pch = 16, col = "red")
dev.off()

# Mean-shift Classification
nNbr = 21
classification <- meanShift(as.matrix(perturb.dt), nNeighbors = nNbr,
                            bandwidth = rep(perturb.h$hopt, ncol(perturb.dt)))
perturb.grp = classification$assignment
length(unique(classification$assignment))

d.ME = 1-sum(diag(table(grp, perturb.grp)))/nrow(dt);d.ME


# Heatmap
perturb.gdt = perturb.dt[order(perturb.grp), 1:10]

colnames(perturb.gdt) = 1:ncol(perturb.gdt)
brk = 50
color = colorRampPalette(c("red","white","blue"))(brk)
myBreaks <- c(seq(min(perturb.gdt), 0, length.out=ceiling(brk/2) + 1), 
              seq(max(perturb.gdt)/brk, max(perturb.gdt), length.out=floor(brk/2)))
cls_label = rep("", nrow(perturb.dt))
cls_label[c(0, cumsum(table(perturb.grp))[-nNbr]) + round((cumsum(table(perturb.grp)) - c(1, cumsum(table(perturb.grp))[-nNbr]+1))/2)] =
  paste("C", 1:nNbr, sep = "")
pheatmap(perturb.gdt, cluster_cols = F, cluster_row = F, color = color, breaks=myBreaks,
         labels_row = cls_label, main="Heatmap", gaps_row = cumsum(table(perturb.grp)))



##########################
#      Random Forest     #
##########################
# https://nishanthu.github.io/articles/ClusteringUsingRandomForest.html
# ISOMAP Data Reduction
reconner = read.csv("data/reconerror.csv")

jpeg(file="images/reconerror_screeplot.jpeg", width=14, height=8, units="cm", res=300, pointsize=6)
plot(reconner$nNbr, reconner$reconerror, type = "b", cex.lab=1.5, cex.axis=1.5, cex.main=1.5,
     xlab = "Number of Neighbors", ylab = "Reconstruction Error")
points(reconner$nNbr[3], reconner$reconerror[3], pch = 16, col = "red")
dev.off()

# Random Forest
# https://ithelp.ithome.com.tw/articles/10303882
# Number of Tree Selection by WCSS
isopy.data = read.csv("data/isomap_data_dim3.csv")[,-1]
ntr = seq(500, 3000, 100)
lWCSS = mclapply(1:length(ntr),
                 function(j){
                   #set.seed(5813)
                   set.seed(1108)
                   clust.urf <- randomForest(isopy.data, proximity = TRUE, oob.prox = TRUE, ntree = ntr[j])
                   hclust.rf <- hclust(as.dist(1-clust.urf$proximity), method = "mcquitty")
                   grp = cutree(hclust.rf, k=15)
                   
                   group.data = cbind(grp, isopy.data)
                   group.mean  = as.data.frame(group.data %>%group_by(grp) %>%
                                                 summarise_at(vars(colnames(isopy.data)), list(name = mean)))[,-1]
                   colnames(group.mean) = colnames(isopy.data)
                   
                   wcss = 0
                   for (i in 1:nrow(group.mean)) {
                     wcss = wcss + sum(dist(isopy.data[which(grp == i),], group.mean[i,]))
                   }
                   return(list(ntr[j], wcss))
                 }, mc.cores = detectCores(), mc.allow.recursive = TRUE)

plot(sapply(lWCSS, "[[", 1), sapply(lWCSS, "[[", 2), type = "b",
     xlab = "Number of Trees", ylab = "WCSS")

write.csv(data.frame(ntree = sapply(lWCSS, "[[", 1), WCSS = sapply(lWCSS, "[[", 2)),
          "ntree_selection_1108.csv", row.names = F)

ntree_selection_1108 = read.csv("data/ntree_selection_1108.csv")
jpeg(file="images/RF_ntree_wcss.jpeg", width=14, height=8, units="cm", res=300, pointsize=6)
par(mar = c(5.1, 4.5, 2.1, 2.1))
plot(ntree_selection_1108$ntree, ntree_selection_1108$WCSS, type = "b",
     xlab = "Number of Trees", ylab = "WCSS", cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
points(ntree_selection_1108$ntree[which.min(ntree_selection_1108$WCSS)],
       min(ntree_selection_1108$WCSS), pch = 16, col = "red")
dev.off()


# ntree = 1000 (set.seed(1108))
# ntree = 1200 (set.seed(5813))

# k Selection by WCSS
set.seed(1108)
clust.urf <- randomForest(isopy.data, proximity = TRUE, oob.prox = TRUE, ntree = 1000)
hclust.rf <- hclust(as.dist(1-clust.urf$proximity), method = "mcquitty")
importance(clust.urf)

WCSS = NULL
for (k in 1:30) {
  cat(k, "\t")
  grp = cutree(hclust.rf, k=k)
  
  group.data = cbind(grp, isopy.data)
  group.mean  = as.data.frame(group.data %>% group_by(grp) %>% summarise_at(vars(colnames(isopy.data)), list(name = mean)))[,-1]
  colnames(group.mean) = colnames(isopy.data)
  
  wcss = 0
  for (i in 1:nrow(group.mean)) {
    wcss = wcss + sum(dist(isopy.data[which(grp == i),], group.mean[i,]))
  }
  
  WCSS[k] = wcss
}

jpeg(file="images/RF_k_wcss.jpeg", width=14, height=8, units="cm", res=300, pointsize=6)
par(mar = c(5.1, 4.5, 2.1, 2.1))
plot(1:30, WCSS, type = "b", xlab = "Number of Neighbors",
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
points(14, WCSS[14], pch = 16, col = "red")
dev.off()



#dim = 10; k = 16
#dim = 4; k = 11
#dim = 3; k = 13
#k = 17 ("mcquitty")
#k = 19 ("average")
#k = 11 (set.seed(1108))
k = 14
grp = cutree(hclust.rf, k=k)
write.table(grp, file = "alg2_out.txt", sep = "\n", row.names = F, col.names = F)
grp = unlist(as.vector(read.csv("alg2_out.txt", header = F)))

write.csv(grp, "isoclass_dim3.csv", row.names = F)

# Heatmap
gdt = dt[order(grp), 1:10]

colnames(gdt) = 1:ncol(gdt)
brk = 50
color = colorRampPalette(c("red","white","blue"))(brk)
myBreaks <- c(seq(min(gdt), 0, length.out=ceiling(brk/2) + 1), 
              seq(max(gdt)/brk, max(gdt), length.out=floor(brk/2)))
cls_label = rep("", nrow(dt))
cls_label[c(0, cumsum(table(grp))[-k]) + round((cumsum(table(grp)) - c(1, cumsum(table(grp))[-k]+1))/2)] =
  paste("C", 1:k, sep = "")
pheatmap(gdt, cluster_cols = F, cluster_row = F, color = color, breaks=myBreaks, fontsize=13.5,
         labels_row = cls_label, main="Heatmap", gaps_row = cumsum(table(grp)))



# Raw data
set.seed(1108)
clust.urf <- randomForest(dt, proximity = TRUE, oob.prox = TRUE, ntree = 1000)
hclust.rf <- hclust(as.dist(1-clust.urf$proximity), method = "mcquitty")
rf.cluster = cutree(hclust.rf, k=10)
table(rf.cluster)

# Number of Group Selection by WCSS
WCSS = NULL
for (k in 1:30) {
  cat(k, "\t")
  grp = cutree(hclust.rf, k=k)
  
  group.data = cbind(grp, dt)
  group.mean  = as.data.frame(group.data %>% group_by(grp) %>% summarise_at(vars(colnames(dt)), list(name = mean)))[,-1]
  colnames(group.mean) = colnames(dt)
  
  wcss = 0
  for (i in 1:nrow(group.mean)) {
    wcss = wcss + sum(dist(dt[which(grp == i),], group.mean[i,]))
  }
  
  WCSS[k] = wcss
}

plot(1:30, WCSS, type = "b")


# Random Forest Classification
k = 16
grp = cutree(hclust.rf, k=k)

# Heatmap
gdt = dt[order(grp), 1:10]

colnames(gdt) = 1:ncol(gdt)
brk = 50
color = colorRampPalette(c("red","white","blue"))(brk)
myBreaks <- c(seq(min(gdt), 0, length.out=ceiling(brk/2) + 1), 
              seq(max(gdt)/brk, max(gdt), length.out=floor(brk/2)))
cls_label = rep("", nrow(dt))
cls_label[c(0, cumsum(table(grp))[-k]) + round((cumsum(table(grp)) - c(1, cumsum(table(grp))[-k]+1))/2)] =
  paste("C", 1:k, sep = "")
pheatmap(gdt, cluster_cols = F, cluster_row = F, color = color, breaks=myBreaks, fontsize=13.5,
         labels_row = cls_label, main="Heatmap", gaps_row = cumsum(table(grp)))


# Evaluate clustering
set.seed(1108)

sc.isodt = isopy.data[which(grp == which.min(table(grp))), ]
nosie_sigma = apply(sc.isodt, 2, function(x) max(dist(x))*0.001)
noise.isodt = matrix(NA, nrow(isopy.data), ncol(isopy.data))
for (i in 1:ncol(isopy.data)) {
  noise.isodt[,i] = rnorm(nrow(isopy.data), sd = nosie_sigma[i])
}

perturb.isodt = isopy.data + noise.isodt
perturb.isodt = read.csv("data/isomap_perturbeddata_dim3.csv")[,-1]
perturb.isodt = isopy.data

set.seed(1108)
clust.urf <- randomForest(perturb.isodt, proximity = TRUE, oob.prox = TRUE, ntree = 1000)
hclust.rf <- hclust(as.dist(1-clust.urf$proximity), method = "mcquitty")

WCSS = NULL
for (k in 1:30) {
  cat(k, "\t")
  perturb.grp = cutree(hclust.rf, k=k)
  
  group.data = cbind(perturb.grp, perturb.isodt)
  group.mean  = as.data.frame(group.data %>% group_by(perturb.grp) %>% summarise_at(vars(colnames(perturb.isodt)), list(name = mean)))[,-1]
  colnames(group.mean) = colnames(perturb.isodt)
  
  wcss = 0
  for (i in 1:nrow(group.mean)) {
    wcss = wcss + sum(dist(perturb.isodt[which(perturb.grp == i),], group.mean[i,]))
  }
  
  WCSS[k] = wcss
}

plot(1:30, WCSS, type = "b", xlab = "Number of Neighbors")


# Random Forest Classification
k = 17
perturb.grp = cutree(hclust.rf, k=k)

# Heatmap
perturb.gdt = dt[order(perturb.grp), 1:10]

colnames(perturb.gdt) = 1:ncol(perturb.gdt)
brk = 50
color = colorRampPalette(c("red","white","blue"))(brk)
myBreaks <- c(seq(min(perturb.gdt), 0, length.out=ceiling(brk/2) + 1), 
              seq(max(perturb.gdt)/brk, max(perturb.gdt), length.out=floor(brk/2)))
cls_label = rep("", nrow(dt))
cls_label[c(0, cumsum(table(perturb.grp))[-k]) + round((cumsum(table(perturb.grp)) - c(1, cumsum(table(perturb.grp))[-k]+1))/2)] =
  paste("C", 1:k, sep = "")
pheatmap(perturb.gdt, cluster_cols = F, cluster_row = F, color = color, breaks=myBreaks, fontsize=13.5,
         labels_row = cls_label, main="Heatmap", gaps_row = cumsum(table(perturb.grp)))


M = as.matrix(table(grp, perturb.grp))
apply(M, 2, which.max)

library(fpc)

NMI(grp, perturb.grp, variant = "sum")
cont_table <- table(grp, perturb.grp)
vi.dist(grp, perturb.grp)


