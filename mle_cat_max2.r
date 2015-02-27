######
##  ML methods for ancestral nodes inference using paralog numbers

library("ape")
tree<-read.tree("all_prot_with_outgroups_nodes_labeled.tree")
genomes<-tree$tip.label
#source("/n/projects/lka/LUCAcog103/Sept2010/Optimization Codes/aceOptim_nonUniformPrior.r")
source("aceOptim_nonequalpriors.r")

 ############# parameter models ##############
 
data = read.csv("all_prot_with_outgroups.csv", h=T, row.names=1)
outer = c("PseAer","SheOne","AerHyd","VibCho","VibFis","PhoPro", "HaeInf","PasMul","ManSuc","HaeDuc")
tree = drop.tip(tree, outer)
cogs = data[,1:34]
cogs = cogs[, !(colnames(cogs) %in% outer)]
ID <- as.vector(row.names(cogs))
cogs = cogs[apply(cogs,1, function(x) sum(x != 0) > 1),]

cogs[cogs >1] = 2
ncogs<-dim(cogs)[1]
 
###model 1
 modelx<-cbind(c(0,2,4), c(1,0, 5), c(3,5,0))
 colnames(modelx)<-c("state0", "state1", "state2")
 rownames(modelx)<-c("state0", "state1", "state2")

 modelx
#       state0 state1 state2
#state0      0      1      3
#state1      2      0      5
#state2      4      5      0




out_M1<-list()
length(out_M1)<-ncogs
for (i in 1:ncogs){
 print(i)
 cogi<-data.frame(colnames(cogs), t(cogs[i, ]))
 colnames(cogi)<-c("genome", "state")
 # remove the X heading the two genomes that have a number as a name
 cogi$genome =  sub("X", "", cogi$genome)
 tips<-data.frame(seq(1:length(tree$tip.label)), tree$tip.label)
 colnames(tips)<-c("id", "label")

 tmp<-merge(tips, cogi, by.x="label", by.y="genome")
 tmp<-tmp[order(tmp$id),]
 x<-tmp$state

 nstates<-length(unique(x))
 modeli<-modelx[1:nstates, 1:nstates]
 out_M1[[i]]<-aceOptim_nonUniformPrior(x, phy=tree,  model=modeli)
 out_M1[[i]]$ID <- ID[i]
}
save(out_M1, file="/n/projects/lka/LUCAcog103/Sept2010/Results/out_M1")


###model 2
 modelx<-cbind(c(0,2,4), c(1,0, 6), c(3,5,0))
 colnames(modelx)<-c("state0", "state1", "state2")
 rownames(modelx)<-c("state0", "state1", "state2")

 modelx
#       state0 state1 state2
#state0      0      1      3
#state1      2      0      5
#state2      4      6      0

out_M2<-list()
length(out_M2)<-ncogs
for (i in 1:ncogs){
 print(i)
 cogi<-data.frame(colnames(cogs)[-1], t(cogs[i, -1]))
 colnames(cogi)<-c("genome", "state")

 tips<-data.frame(seq(1:103), tree$tip.label)
 colnames(tips)<-c("id", "label")

 tmp<-merge(tips, cogi, by.x="label", by.y="genome")
 tmp<-tmp[order(tmp$id),]
 x<-tmp$state

 nstates<-length(unique(x))
 modeli<-modelx[1:nstates, 1:nstates]
 out_M2[[i]]<-aceOptim_nonUniformPrior(x, phy=tree,  model=modeli)
 out_M2[[i]]$ID <- ID[i]
}

save(out_M2, file="/n/projects/lka/LUCAcog103/Sept2010/Results/out_M2")
