setwd("~/Study/research/scisoform/")
allpvalues2 <- readRDS("transset.rds")
allpvalues3 <- readRDS("transset_depth.rds")
pvalues <- readRDS("pvalues.rds")
pvalues <- t(matrix(unlist(pvalues), nrow = 2))
pvalues0 <- readRDS("pvalues0.rds")
pvalues0 <- t(matrix(unlist(pvalues0), nrow = 2))
pvalues.binom <- readRDS("pvalues_binom.rds")
pvalues.binom <- t(matrix(unlist(pvalues.binom), nrow = 2))
pvalues.binom.depth <- readRDS("pvalues_binom_depth.rds")
pvalues.binom.depth <- t(matrix(unlist(pvalues.binom.depth), nrow = 3))



gene.depth <- read.table("SingleCell_CellType_Gene_Abundance.txt", header = T)
gene.name <- gene.depth$annot_gene_name
gene.depth <- gene.depth[,-1]

for (j in 1:length(allpvalues3)){
  if (is.null(allpvalues3[[j]])) {
    allpvalues3[[j]] <- c(1,1,1)
  }
  if (allpvalues3[[j]] == c(1,1)) {
    allpvalues3[[j]] <- c(1,1,1)
  }
}
p <- t(matrix(unlist(allpvalues3), nrow = 3))
p3 <- p[,3]

p1 <- pvalues[,2]
p0 <- pvalues0[,2]
p.binom <- pvalues.binom[,2]
p.binom.depth <- pvalues.binom.depth[,3]
p.binom.depth[is.na(p.binom.depth)] <- 1
index3 <- which(p.adjust(p3, method = "BH") < 0.05)
index1 <- which(p.adjust(p1, method = "BH") < 0.05)
index0 <- which(p.adjust(p0, method = "BH") < 0.05)
index.binom <- which(p.adjust(p.binom, method = "BH") < 0.05)
index.binom.depth <- which(p.adjust(p.binom.depth, method = "BH") < 0.05)

intersect(order(p.binom)[1:100], order(p.binom.depth)[1:100])

length(intersect(index.binom, index.binom.depth))

num.nonzero <- sapply(1:length(index1), function(i){
  length(unique(celltype[which(gene.iso[index1[i],]!=0)]))
  })
saveRDS(num.nonzero, "index1_cell_non0.rds")
index1.non0 <- index1[num.nonzero>1]
gene.iso.interesting <- gene.iso[index3, ]
rownames(gene.iso.interesting) <- transcript[index3]
saveRDS(gene.iso.interesting, "transset_interesting_depth.rds")
index_top3 <- order(p3[index3])[1:20]

pvalue.table <- data.frame(transcript = transcript, pvalue = p3)
write.csv(pvalue.table, file = "pvalues_nbreg.csv", row.names = F)
a <- read.csv("pvalues_nbreg.csv")

covariates <- names(summary(model)$coefficients[,4])
sum(allpvalues < 0.05)

index2 <- which(p.adjust(allpvalues2[,2], method = "BH") < 0.05)
gene.iso.interesting <- gene.iso[index2, ]
saveRDS(gene.iso.interesting, "transset_interesting.rds")
allpvalues2 <- allpvalues2[index2,]
index_top <- order(allpvalues[,2])[1:20]
index_top2 <- order(allpvalues2[,2])[1:20]
index_top1 <- order(p1)[1:20]
index_top0 <- order(p0)[1:20]
index_top.binom <- order(p.binom)[1:20]


index.top1.non0 <-  intersect( index_top1, index1.non0)
allpvalues2 <- readRDS(file = "transset1.rds")
for (i in 2:129){
  allpvalues3 <- readRDS(paste0("transset",i,".rds"))
  
}
sapply(1:129, function(i){
  nrow(readRDS(paste0("transset",i,".rds")))
})
index <- order(allpvalues)[1:20]
trans_intresting <- index%/%18 +1
table(transcript[index%/%18 +1])
table(covariants[index%%18])
ct_interesting <- strsplit(covariants[index%%18],"celltype")


library(ggplot2)
for (i in 1:5){
  data <- data.frame(expression = as.numeric(gene.iso[index_top[i],]), celltype = celltype)
  fig <- ggplot(data, aes(x=celltype, y=expression)) + 
    geom_boxplot() + ggtitle(paste0("Plot of ", i, "th transcript") )
  print(fig)
}

for (i in 1:5){
  data <- data.frame(expression = as.numeric(gene.iso.interesting[index_top2[i],]), celltype = celltype)
  fig <- ggplot(data, aes(x=celltype, y=expression)) + theme(axis.text.x = element_text(angle = 90))+
    geom_boxplot() + ggtitle(paste0("Transcript ", rownames(gene.iso.interesting)[index_top2[i]]) )
  print(fig)
}

for (i in 1:20){
  data <- data.frame(expression = as.numeric(gene.iso[i,]), celltype = celltype)
  fig <- ggplot(data, aes(x=celltype, y=expression)) + 
    geom_boxplot() + ggtitle(paste0("Plot of ", i, "th transcript") )
  print(fig)
}

for (i in 1:5){
  data <- data.frame(expression = as.numeric(gene.iso.interesting[index_top3[i],]), celltype = celltype)
  fig <- ggplot(data, aes(x=celltype, y=expression)) + theme(axis.text.x = element_text(angle = 90))+
    geom_boxplot() + ggtitle(paste0("Transcript ", rownames(gene.iso.interesting)[index_top3[i]]) )
  print(fig)
}

for (i in 1:5){
  data <- data.frame(expression = as.numeric(gene.iso[index_top1[i],]), celltype = celltype)
  fig <- ggplot(data, aes(x=celltype, y=expression)) + theme(axis.text.x = element_text(angle = 90))+
    geom_boxplot() + ggtitle(paste0("Transcript ", rownames(gene.iso.interesting)[index_top3[i]]) )
  print(fig)
}
for (i in 1:5){
  j <- index_top.binom[i]
  y <- as.numeric(gene.iso[j,])
  genename <- trans.gene.name[j]
  y1 <- as.numeric(gene.depth[match(genename, gene.name),])
  y <- y/y1
  data <- data.frame(expression = y, celltype = celltype)
  fig <- ggplot(data, aes(x=celltype, y=expression)) + theme(axis.text.x = element_text(angle = 90))+
    geom_boxplot() + ggtitle(paste0("Transcript ", trans.id[index_top.binom[i]]) )
  print(fig)
}
for (i in 1:5){
  j <- index_top1[i]
  y <- as.numeric(gene.iso[j,])
  genename <- trans.gene.name[j]
  y1 <- as.numeric(gene.depth[match(genename, gene.name),])
  y <- y/y1
  data <- data.frame(proportion = y, celltype = celltype)
  fig <- ggplot(data, aes(x=celltype, y=proportion)) + theme(axis.text.x = element_text(angle = 90))+
    geom_boxplot() + ggtitle(paste0("Transcript ", rownames(trans.id)[index_top1[i]]))
  print(fig)
}
#hierachical anova
#nb regression
data <- data.frame(expression = as.numeric(gene.iso[i,]), celltype = celltype)
fig <- ggplot(data, aes(x=celltype, y=expression)) + 
  geom_boxplot() + ggtitle(paste0("Plot of ", i, "th transcript") )
print(fig)

model <- glm(cbind(y,y1-y) ~ donar + celltype, family = quasibinomial)
model0 <- glm(cbind(y,y1-y) ~ donar, family = quasibinomial)
pchisq(deviance(model0) - deviance(model), df.residual(model0) - df.residual(model), lower.tail = F)
pvalues[j,]
model <- glm(cbind(y,y1-y) ~ donar + celltype, family = binomial)
model0 <- glm(cbind(y,y1-y) ~ donar, family = binomial)
pchisq(deviance(model0) - deviance(model), df.residual(model0) - df.residual(model), lower.tail = F)
pvalues.binom[j,]

pvalues.new <- data.frame(trans.id, p.binom, p.binom.depth)
write.csv(pvalues.new, file = "pvalues_binom.csv", row.names = F)
pvalues.new <- read.csv("pvalues_binom.csv")
pvalues.binom.pair <- readRDS('pvalues_binom_3pairs.rds')
pvalues.binom.pair <- t(matrix(unlist(pvalues.binom.pair), nrow = 3))

p1 <- pvalues.binom.pair[,1]
length(which(p.adjust(p1, method = "BH") < 0.05))
p2 <- p1[p1<=0.99]
length(which(p.adjust(p2, method = "BH") < 0.05))

pvalues.nbinom.pair <- readRDS('pvalues_negbinom_3pairs.rds')
pvalues.nbinom.pair <- t(matrix(unlist(pvalues.nbinom.pair), nrow = 3))
#vRG vs oRG
index_top <- order(pvalues.binom.pair[,1])[1:20]
for (i in 1:5){
  j <- index_top[i]
  y <- as.numeric(gene.iso[j,])
  genename <- trans.gene.name[j]
  y1 <- as.numeric(gene.depth[match(genename, gene.name),])
  y <- y/y1
  data <- data.frame(proportion = y, celltype = celltype)
  fig <- ggplot(data, aes(x=celltype, y=proportion)) + theme(axis.text.x = element_text(angle = 90))+
    geom_boxplot() + ggtitle(paste0("Transcript ", trans.id[index_top[i]]))
  print(fig)
}
#vRG vs IP
index_top <- order(pvalues.binom.pair[,2])[1:20]
for (i in 1:5){
  j <- index_top[i]
  y <- as.numeric(gene.iso[j,])
  genename <- trans.gene.name[j]
  y1 <- as.numeric(gene.depth[match(genename, gene.name),])
  y <- y/y1
  data <- data.frame(proportion = y, celltype = celltype)
  fig <- ggplot(data, aes(x=celltype, y=proportion)) + theme(axis.text.x = element_text(angle = 90))+
    geom_boxplot() + ggtitle(paste0("Transcript ", trans.id[index_top[i]]))
  print(fig)
}
#oRG vs IP
index_top <- order(pvalues.binom.pair[,3])[1:20]
for (i in 1:5){
  j <- index_top[i]
  y <- as.numeric(gene.iso[j,])
  genename <- trans.gene.name[j]
  y1 <- as.numeric(gene.depth[match(genename, gene.name),])
  y <- y/y1
  data <- data.frame(proportion = y, celltype = celltype)
  fig <- ggplot(data, aes(x=celltype, y=proportion)) + theme(axis.text.x = element_text(angle = 90))+
    geom_boxplot() + ggtitle(paste0("Transcript ", trans.id[index_top[i]]))
  print(fig)
}
#vRG vs oRG
index_top <- order(pvalues.nbinom.pair[,1])[1:20]
for (i in 1:5){
  j <- index_top[i]
  y <- as.numeric(gene.iso[j,])
  data <- data.frame(expression = y, celltype = celltype)
  fig <- ggplot(data, aes(x=celltype, y=expression)) + theme(axis.text.x = element_text(angle = 90))+
    geom_boxplot() + ggtitle(paste0("Transcript ", trans.id[index_top[i]]))
  print(fig)
}
#vRG vs IP
index_top <- order(pvalues.nbinom.pair[,2])[1:20]
for (i in 1:5){
  j <- index_top[i]
  y <- as.numeric(gene.iso[j,])
  data <- data.frame(expression = y, celltype = celltype)
  fig <- ggplot(data, aes(x=celltype, y=expression)) + theme(axis.text.x = element_text(angle = 90))+
    geom_boxplot() + ggtitle(paste0("Transcript ", trans.id[index_top[i]]))
  print(fig)
}
#oRG vs IP
index_top <- order(pvalues.nbinom.pair[,3])[1:20]
for (i in 1:5){
  j <- index_top[i]
  y <- as.numeric(gene.iso[j,])
  data <- data.frame(expression = y, celltype = celltype)
  fig <- ggplot(data, aes(x=celltype, y=expression)) + theme(axis.text.x = element_text(angle = 90))+
    geom_boxplot() + ggtitle(paste0("Transcript ", trans.id[index_top[i]]))
  print(fig)
}
pvalues.pair <- cbind(pvalues.binom.pair, pvalues.nbinom.pair)
colnames(pvalues.pair) <- c("vRG_oRG_prop", "vRG_IP_prop", "oRG_IP_prop", "vRG_oRG_expr", "vRG_IP_expr", "oRG_IP_expr")
rownames(pvalues.pair) <- trans.id
write.csv(pvalues.pair, file = "pvalues_3pairs.csv", row.names = T)


pvalues.multi <- readRDS('pvalues_multinom.rds')
pvalues.multi <- unlist(pvalues.multi)
hist(pvalues.multi)
length(which(p.adjust(pvalues.multi, method = "BH") < 0.01))
pvalues.multi[1:20]

####### get gene index of each transcript #####

transindex <- mclapply(1:nrow(gene.iso), function(j){
  y <- as.numeric(gene.iso[j,])
  genename <- trans.gene.name[j]
  y1 <- match(genename, gene.name)
  return(y1)
 }, mc.cores = ncores)
transindex <- unlist(transindex)
saveRDS(transindex,'trans_geneindex.rds')
transindex <- readRDS('trans_geneindex.rds')

transindex[order(p.binom)[1:20]]
pvalues.multi[transindex[order(p.binom)[1:20]]]
p3528 <- p.binom[which(transindex==3528)]
pchisq(-2*sum(log(p3528)), 2*length(p3528), lower.tail = F)
     
 
plot(-log(pvalues.multi[transindex]), -log(p.binom))
pvalues.multi[order(pvalues.multi)[1:20]]

readRDS('pvalues_multinom.rds')
write.csv(pvalues.multi, file = "pvalues_gene.csv", row.names = F)
