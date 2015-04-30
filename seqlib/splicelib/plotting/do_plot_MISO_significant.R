library(R.utils)
library(ggplot2)
library(MASS)
library(scales)

#input.file2="/Users/psudmant/Documents/science/mounts/c3ddb01/home/psudmant/projects/trap_development/analysis/splicing/D1_D2/output/D1_vs_D2.delta_psi.df"
#input.file1="/Users/psudmant/Documents/science/mounts/c3ddb01/home/psudmant/projects/trap_development/analysis/splicing/D1_D2/output/D1_vs_D2.df"

args = cmdArgs()
input.file1 = cmdArg("fn_input_1")
output.file1 = cmdArg("fn_output_1")
input.file2 = cmdArg("fn_input_2")
output.file2 = cmdArg("fn_output_2")
output.file3 = cmdArg("fn_output_3")
title=cmdArg("title")

t1=read.table(input.file1,header=TRUE,sep='\t')
t2=read.table(input.file2,header=TRUE,sep='\t')

level_order = rev(unique(t1$comparison))
t1$comparison = factor(t1$comparison,levels=level_order)
t2$comparison = factor(t2$comparison,levels=level_order)

pdf(output.file1,width=5,height=4)
g=ggplot(t1)
g+geom_bar(aes(x=stype,y=significant/total,fill=comparison),stat='identity',position='dodge')+geom_text(aes(x=stype,y=significant/total,label=significant,group=comparison),fontface=2,stat='identity',position=position_dodge(.9),color='white',size=3.5,hjust=1)+theme_bw()+coord_flip()+scale_fill_brewer(palette="Set1")+scale_y_continuous(bquote(paste(hat(p)," BF>10 ",.(title))))+scale_x_discrete("")
dev.off()

s1=t2$sample1[1]
s2=t2$sample2[1]
pdf(output.file2,width=4,height=4)
g=ggplot(t2)
g+geom_boxplot(aes(x=stype,y=delta_psi,color=comparison),scale='width',position=position_dodge(width=1))+geom_point(aes(x=stype,y=delta_psi,color=comparison,fill=comparison),size=.8,alpha=.1,position=position_jitterdodge(dodge.width=1,jitter.height=0))+theme_bw()+scale_color_brewer(palette="Set1")+coord_flip()+scale_y_continuous(bquote(paste(Delta, Psi, " ",.(title))))+scale_x_discrete("")
dev.off()

pdf(output.file3,width=4,height=4)
g=ggplot(t2)
g+geom_boxplot(aes(x=stype,y=delta_psi,color=comparison),notch=TRUE,scale='width',position=position_dodge(width=1))+geom_point(aes(x=stype,y=delta_psi,color=comparison,fill=comparison),size=.8,alpha=.1,position=position_jitterdodge(dodge.width=1,jitter.height=0))+theme_bw()+scale_color_brewer(palette="Set1")+coord_flip()+scale_y_continuous(bquote(paste(Delta, Psi, " ",.(title))))+scale_x_discrete("")
dev.off()

#g=ggplot(t)
#g+geom_bar(aes(x=stype,y=significant,fill=comparison),stat='identity',position='dodge')+theme_bw()+coord_flip()+scale_fill_brewer(palette="Set1")
