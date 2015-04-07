library(ggplot2)

#t=read.table("/Users/psudmant/Documents/science/mounts/c3ddb01/home/psudmant/projects/ancient_exons/timecourse_expression/RbFox_levels/output/RbFOX_expression.df",header=T,sep="\t")

args <- commandArgs(trailingOnly=TRUE)
input.file <- args[1]
output.file <- args[2]


t=read.table(input.file,header=T,sep="\t")

xlabels=unique(t$alt_name[order(t$group)])
xbreaks=seq(1,length(xlabels))

pdf(output.file,width=2.5,height=3.5)
g=ggplot(t)
#+geom_line(aes(x=group,y=fpkm,color=GeneName))+geom_point(aes(x=group,y=fpkm,group=group,color=GeneName),size=1.5)
g+stat_summary(aes(x=group,y=fpkm,color=GeneName),fun.y=mean,geom="line")+stat_summary(fun.data="mean_cl_boot",aes(x=group,y=fpkm,color=GeneName))+theme_bw()+scale_color_brewer(palette="Set1")+theme_bw()+scale_color_brewer(palette="Set1")+theme(legend.title=element_blank(),legend.key=element_blank(),legend.position="bottom",legend.direction='horizontal',axis.text.x=element_text(angle=45,hjust=1,vjust=1))+guides(color=guide_legend(nrow=2))+scale_x_continuous("",breaks=xbreaks,labels=xlabels)+scale_y_continuous("FPKM")

dev.off()
