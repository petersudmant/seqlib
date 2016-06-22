library(ggplot2)
library(dplyr)
library(scales)
library(tidyr)
library(R.utils)
library(gridExtra)

args = cmdArgs()
fn_cvg = cmdArg("fn_cvg")
fn_gene = cmdArg("fn_gene_model")
fn_output = cmdArg("fn_output")

t=read.table(fn_cvg, header=TRUE,sep="\t")
t$age = t$time
t_gene_model = read.table(fn_gene, header=TRUE, sep="\t")

g1=ggplot(t)
g1=g1+geom_area(aes(x=t_pos,y=log((cvg)+1)/log(10),fill=age,group=sample))+
    theme_classic(base_size=8)+         
    scale_fill_gradient2("age",midpoint=30,low=muted("red"),high=muted('blue'))+
    theme(axis.text.x=element_blank())+
    theme(axis.ticks.x=element_blank())+
    theme(axis.line.x=element_blank())+
    guides(fill=guide_colorbar(barwidth=.3,barheight=3,title.size=6))+
    scale_x_continuous("")+scale_y_continuous("log(cvg)")+
    theme(legend.position=c(.1,.7))+
    theme(plot.margin = unit(c(0,.2,.5,0),"cm"))

# g1=ggplot(filter(t_plot))
# g1+geom_area(aes(x=t_pos,y=log((cvg/CDS_mu)+1)/log(10),fill=age,group=sample))+
    # theme_classic()+
    # scale_fill_gradient2(midpoint=30,low=muted("red"),high=muted('blue'))
    
g2=ggplot(t_gene_model)+
    geom_polygon(aes(x=x, y=y, group = id), fill="steelblue")+
    theme_classic()+
    theme(axis.text.x=element_text(size=6), axis.line.x=element_blank())+
    scale_x_continuous("")+scale_y_continuous("")+
    theme(axis.line.y=element_blank())+
    theme(axis.text.y=element_blank())+
    theme(axis.ticks.y=element_blank())+
    theme(plot.margin = unit(c(-1,.2,0,0),"cm"))+
    theme(panel.background=element_blank())
gA = ggplotGrob(g1)
gB = ggplotGrob(g2)
gB$widths = gA$widths

pdf(fn_output, width=2.5, height=1.8)
grid.arrange(gA, gB, ncol=1, heights=c(.95,.05))
dev.off()    
    
    

    
