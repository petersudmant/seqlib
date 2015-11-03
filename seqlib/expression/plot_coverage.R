library(R.utils)
library(ggplot2)
library(MASS)
library(scales)
library(gridExtra)
library(gtable)

args = cmdArgs()
fn_cvg = cmdArg("fn_cvg")
fn_gene_model = cmdArg("fn_gene_model")
fn_output = cmdArg("fn_output")

t_cvg = read.table(fn_cvg, header=TRUE, sep="\t")

theme_no_x_axis = theme(axis.text.x=element_blank(),
                        axis.line.x=element_blank(),
                        axis.ticks.x=element_blank())

theme_no_y_axis = theme(axis.text.y=element_blank(),
                        axis.line.y=element_blank(),
                        axis.ticks.y=element_blank())
g1=ggplot(t_cvg)+
    geom_area(aes(x=pos,y=cvg),fill="steelblue")+
    theme_classic()+
    facet_grid(sample~.,scales="free_y")+
    theme(axis.text.y=element_text(size=5))+
    scale_x_continuous("")+scale_y_continuous("")+
    theme_no_x_axis  



t_gene_model = read.table(fn_gene_model, header=TRUE, sep="\t")
g2=ggplot(t_gene_model)+
    geom_polygon(aes(x=x, y=y, group = id), fill="steelblue")+
    theme_classic()+
    theme(axis.text.x=element_text(size=5), axis.line.x=element_blank())+
    scale_x_continuous("")+scale_y_continuous("")+
    theme(plot.margin = unit(c(-1,0,0,0),"cm"))+
    theme_no_y_axis

print("1")
gA = ggplotGrob(g1)
print("2")
gB = ggplotGrob(g2)

gB$widths = gA$widths

print("3")
pdf(fn_output, width=4, height=6)
grid.arrange(gA, gB, ncol=1, heights=c(.95,.05))
dev.off()

