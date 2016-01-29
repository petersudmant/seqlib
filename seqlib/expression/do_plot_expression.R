library(R.utils)
library(ggplot2)
library(MASS)
library(scales)

#t=read.table("/Users/psudmant/Documents/science/mounts/c3ddb01/home/psudmant/projects/trap_development/analysis/expression/D1_D2_genes_of_interest/output/striatum/Drd.df",header=T,sep="\t")

args = cmdArgs()
input.file = cmdArg("fn_input")
output.file = cmdArg("fn_output")
timecourse = cmdArg("timecourse",TRUE)
log = cmdArg("log",FALSE)
add1 = cmdArg("add1",FALSE)

print(input.file)
print(output.file)

t=read.table(input.file,header=T,sep="\t")
if (add1){
    t$fpkm = t$fpkm+1
}

xlabels=unique(t$alt_name[order(t$group)])
xbreaks=unique(t$group[order(t$group)])
#xbreaks=seq(1,length(xlabels))

pdf(output.file,width=2.5,height=3)
g=ggplot(t)
if (timecourse){
#+geom_line(aes(x=group,y=fpkm,color=GeneName))+geom_point(aes(x=group,y=fpkm,group=group,color=GeneName),size=1.5)
    g=g+stat_summary(aes(x=group,y=fpkm,color=GeneName),fun.y=mean,geom="line")+stat_summary(fun.data="mean_cl_normal",aes(x=group,y=fpkm,color=GeneName))+scale_color_brewer(palette="Set1")+theme_classic()+scale_color_brewer(palette="Set1")+theme(legend.title=element_blank(),legend.key=element_blank(),legend.position="bottom",legend.direction='horizontal',axis.text.x=element_text(angle=45,hjust=1,vjust=1))+guides(color=guide_legend(nrow=2))+scale_x_continuous("",breaks=xbreaks,labels=xlabels)+scale_y_continuous("FPKM")

}else{
	g=g+stat_summary(aes(x=group,y=fpkm,fill=GeneName),fun.y=mean,geom='bar',position='dodge')+stat_summary(geom="errorbar",width=.5,fun.data="mean_cl_normal",position=position_dodge(width=.93),aes(x=group,y=fpkm,color=GeneName,group=GeneName))+theme_classic()+scale_fill_brewer(palette="Set1")+scale_color_brewer(palette="Set1")+theme(legend.title=element_blank(),legend.key=element_blank(),legend.position="bottom",legend.direction='horizontal',axis.text.x=element_text(angle=45,hjust=1,vjust=1))+guides(color=guide_legend(nrow=2),fill=guide_legend(nrow=2))+scale_x_continuous("",breaks=xbreaks,labels=xlabels)+scale_y_continuous("FPKM")
	
    #g+stat_summary(fun.data="mean_cl_boot",aes(x=group,y=fpkm,color=GeneName,group=GeneName),position='dodge')+theme_bw()+scale_color_brewer(palette="Set1")+theme(legend.title=element_blank(),legend.key=element_blank(),legend.position="bottom",legend.direction='horizontal',axis.text.x=element_text(angle=45,hjust=1,vjust=1))+guides(color=guide_legend(nrow=2))+scale_x_continuous("",breaks=xbreaks,labels=xlabels)+scale_y_continuous("FPKM")
}

if (log){
    g=g+scale_y_log10("FPMK", breaks = trans_breaks("log10", function(x) 10^x),
                      labels = trans_format("log10", math_format(10^.x)))
}
print(g)
dev.off()
