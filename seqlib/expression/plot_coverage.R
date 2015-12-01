library(R.utils)
library(ggplot2)
library(MASS)
library(scales)
library(gridExtra)
library(gtable)
library(RColorBrewer)
library(dplyr)

args = cmdArgs()
fn_cvg = cmdArg("fn_cvg")
fn_gene_model = cmdArg("fn_gene_model")
fn_output = cmdArg("fn_output")
gene = cmdArg("gene")

height = cmdArg("height",-1)
width = cmdArg("width",4)
scale = cmdArg("scale",.5)

t_cvg = read.table(fn_cvg, header=TRUE, sep="\t")

samples = levels(t_cvg$sample)
if (height==-1){
    n = length(samples) 
    height = n * scale+.3 
}

t_gene_model = read.table(fn_gene_model, header=TRUE, sep="\t")
xmin = min(t_gene_model$x)
xmax = max(t_gene_model$x)

annot = t_cvg %>% 
            group_by(sample) %>%
            summarize(y=max(cvg))

annot$x = xmin 

theme_no_x_axis = theme(axis.text.x=element_blank(),
                        axis.line.x=element_blank(),
                        axis.ticks.x=element_blank())

theme_no_y_axis = theme(axis.text.y=element_blank(),
                        axis.line.y=element_blank(),
                        axis.ticks.y=element_blank())

theme_no_legend = theme(legend.position="none")
theme_no_facet_box = theme(strip.background=element_blank(),
                           strip.text=element_blank())

colors = brewer.pal(7, "Set1")
colors = c(brewer.pal(11, "Spectral"), brewer.pal(12, "Paired"))
colors = c(colors,colors,colors) 

g1=ggplot(t_cvg)+
    geom_area(aes(x=pos,y=cvg,fill=celltype))+
    facet_grid(sample~., scales="free_y")+
    scale_x_continuous("",lim=c(xmin,xmax))+scale_y_continuous("")+
    scale_fill_manual(values=colors)+
    geom_text(data=annot, aes(x=x, y=y+5, label=sample),size=1.5,hjust=0) +
    theme_classic()+
    theme(axis.text.y=element_text(size=5))+
    theme_no_x_axis +
    theme_no_legend +
    theme_no_facet_box

g2=ggplot(t_gene_model)+
    geom_polygon(aes(x=x, y=y, group = id), fill="steelblue")+
    theme_classic()+
    theme(axis.text.x=element_text(size=5), axis.line.x=element_blank())+
    scale_x_continuous("",lim=c(xmin,xmax))+scale_y_continuous("")+
    theme(plot.margin = unit(c(-1,0,0,0),"cm"))+
    annotate("text", x=xmin, y=0, label=gene, hjust=0, size=1.5)+
    theme_no_y_axis

gA = ggplotGrob(g1)
gB = ggplotGrob(g2)
gB$widths = gA$widths

pdf(fn_output, width=width, height=height)
grid.arrange(gA, gB, ncol=1, heights=c(.95,.05))
dev.off()

