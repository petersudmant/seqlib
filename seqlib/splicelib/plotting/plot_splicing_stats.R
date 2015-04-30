
t=read.table("/Users/psudmant/Documents/science/mounts/c3ddb01/scratch/users/psudmant/indexes/full_micro_exon_index/MISO_indexes/mus_musculus/stats/stats.df",header=T,sep="\t")

names(t)
g=ggplot(t)
g+geom_bar(aes(x=stype,y=n,fill=stype),stat='identity')+theme_bw()+scale_fill_brewer(palette="Set1")+coord_flip()+theme(legend.position="none")+scale_y_continuous(breaks=c(0,10000,20000))+scale_x_discrete("")+ggtitle("mouse")


t=read.table("/Users/psudmant/Documents/science/mounts/c3ddb01/scratch/users/psudmant/indexes/full_micro_exon_index/MISO_indexes/homo_sapiens/stats/stats.df",header=T,sep="\t")
g=ggplot(t)
g+geom_bar(aes(x=stype,y=n,fill=stype),stat='identity')+theme_bw()+scale_fill_brewer(palette="Set1")+coord_flip()+theme(legend.position="none")+scale_x_discrete("")+scale_y_continuous(breaks=c(0,25000,50000,75000))+ggtitle("human")


