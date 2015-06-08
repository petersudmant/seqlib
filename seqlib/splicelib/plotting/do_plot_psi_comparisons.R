library(R.utils)
library(ggplot2)
library(MASS)
library(scales)

args = cmdArgs()
input.file1 = cmdArg("fn_input_1")
output.file1 = cmdArg("fn_output_1")
input.file2 = cmdArg("fn_input_2")
output.file2 = cmdArg("fn_output_2")
comp1 = cmdArg("comp1")
comp2 = cmdArg("comp2")

title=cmdArg("title")

t1=read.table(input.file1,header=TRUE,sep='\t')
t2=read.table(input.file2,header=TRUE,sep='\t')

#t=read.table("/Users/psudmant/Documents/science/mounts/c3ddb01/home/psudmant/projects/trap_development/analysis/splicing/D1_D2/output/transcript_comparisons/E16_vs_PN42.transcript.df",header=T,sep="\t")

pdf(output.file1,width=6,height=2)
g=ggplot(t1)
g+geom_point(aes(x=D1,y=D2),alpha=.7,size=1)+geom_smooth(aes(x=D1,y=D2),method="lm")+theme_classic(base_size=6)+scale_color_brewer(palette="Set1")+facet_grid(~stype)+scale_x_continuous(bquote(paste(.(comp1)," ", Delta, Psi, " ",.(title))),lim=c(-1,1))+scale_y_continuous(bquote(paste(.(comp2)," ", Delta, Psi, " ",.(title))),lim=c(-1,1))
dev.off()

pdf(output.file2,width=4,height=4)
g=ggplot(t2)
g+geom_point(aes(x=D1,y=D2),alpha=.7,size=1)+geom_smooth(aes(x=D1,y=D2),method="lm")+theme_classic(base_size=6)+scale_color_brewer(palette="Set1")+facet_grid(stype_1~stype_2)+scale_x_continuous(bquote(paste(.(comp1)," ", Delta, Psi, " ",.(title))),lim=c(-1,1))+scale_y_continuous(bquote(paste(.(comp2)," ", Delta, Psi, " ",.(title))),lim=c(-1,1))
dev.off()

#AFE_ALE = (t$stype_2=="AFE")&(t$stype_1=="ALE")
#ALE_AFE = (t$stype_1=="AFE")&(t$stype_2=="ALE")

#summary(lm(t$D2[AFE_ALE]~t$D1[AFE_ALE]))
#summary(lm(t$D2[ALE_AFE]~t$D1[ALE_AFE]))
#t=read.table("/Users/psudmant/Documents/science/mounts/c3ddb01/home/psudmant/projects/trap_development/analysis/splicing/D1_D2/output/transcript_comparisons/E16_vs_PN42.event.df",header=T,sep="\t")
#t=read.table("/Users/psudmant/Documents/science/mounts/c3ddb01/home/psudmant/projects/trap_development/analysis/splicing/D1_D2/output/transcript_comparisons/PN42_vs_2yr.event.df",header=T,sep="\t")


