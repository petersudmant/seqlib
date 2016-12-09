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
midpoint = cmdArg("midpoint",30)

t=read.table(fn_cvg, header=TRUE,sep="\t")
t_gene_model = read.table(fn_gene, header=TRUE, sep="\t")

t$age = t$time

levs = t %>% select(sample,age) %>%
            group_by(sample,age) %>%
            unique() %>%
            arrange(age)
    
t$sample = factor(t$sample,levels=levs$sample)

g1=ggplot(t)
g1=g1+geom_area(aes(x=t_pos,
                    y=log((cvg)+1)/log(10),
                    fill=age,
                    group=sample),
                position=position_stack(reverse=FALSE))+
    theme_classic(base_size=8)+
    scale_fill_gradient2("age",midpoint=30,low="red",high=muted('blue'))+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
              axis.line.x=element_blank())+
    guides(fill=guide_colorbar(barwidth=.3,barheight=3,title.size=6))+
    scale_x_continuous("")+
    scale_y_continuous(expression(log[10](cvg)))+
    #theme(legend.position=c(.1,.7))+
    theme(plot.margin = unit(c(0,.2,.5,0),"cm"))

    
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

pdf(fn_output, width=2.7, height=1.8)
grid.arrange(gA, gB, ncol=1, heights=c(.95,.05))
dev.off()    
    
    

    
