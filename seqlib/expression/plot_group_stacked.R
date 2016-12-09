library(R.utils)
library(ggplot2)
library(dplyr)
library(broom)
library(scales)
library(RColorBrewer)
library(grid)
library(gridExtra)

#fn = "/Users/psudmant/Documents/science/mounts/c3ddb01/home/psudmant/projects/trap_development/analysis/3p_UTR_frags/cvg_assess/plot_cvg/mouse/output/cvg_plots/DRD_all/Cdpf1.cvg.df.gz"
#fn_gene="/Users/psudmant/Documents/science/mounts/c3ddb01/home/psudmant/projects/trap_development/analysis/3p_UTR_frags/cvg_assess/plot_cvg/mouse/output/cvg_plots/DRD_2yr/Cdpf1.model.df.gz"
#fn_output = "~/Documents/science/burgelab/projects/trap_development/paper/figure_assembly/meta_plots/plot_D1_D2_examples/Cdpf1.pdf"

args = cmdArgs()
fn_cvg = cmdArg("fn_cvg")
fn_gene = cmdArg("fn_gene_model")
fn_output = cmdArg("fn_output")
gene_name = cmdArg("gene_name")

t=read.table(fn_cvg,header=T,sep="\t")
t_gene_model = read.table(fn_gene, header=TRUE, sep="\t")

t$age = t$time
t = t %>% filter(group %in% c("PN42","2yr"))
t$group = factor(t$group)

t$t_celltype = interaction(t$celltype, t$group)
t$t_celltype = factor(t$t_celltype, levels=c("Drd1a.PN42",
											 "Drd2.PN42",									
										     "Drd1a.2yr",
											 "Drd2.2yr"))
											 											
labels = data.frame(t_celltype=c("Drd1a.PN42","Drd2.2yr","Drd2.PN42","Drd1a.2yr"),
					label=c("D1 PN42", "D2 2yr","D2 PN42","D1 2yr"))


Reds=brewer.pal(6,"Reds")[3:6]
Blues=brewer.pal(6,"Blues")[3:6]
Greens=brewer.pal(6,"Greens")[3:6]
Purples=brewer.pal(6,"Purples")[3:6]

R1=brewer.pal(9,"Reds")[2:5]
R2=brewer.pal(9,"Reds")[6:9]
B1=brewer.pal(9,"Blues")[2:5]
B2=brewer.pal(9,"Blues")[6:9]

colors=c(Reds,Blues)
colors=c(Reds,Blues,Purples,Greens)

#2yD1 2yrD2 Pn42D1 Pn42D2
C1=brewer.pal(9,"Greys")[2:5]
greys=brewer.pal(9,"Greys")[3:6]

YlOr = brewer.pal(9,"YlOrRd")[2:5]
OrRd = brewer.pal(9,"YlOrRd")[6:9]

Pu = brewer.pal(9,"PuBu")[2:5]
Bu = brewer.pal(9,"PuBu")[6:9]

D2_PN42 = Pu
D2_PN42 = Purples
D2_2yr = Bu
D1_PN42 = YlOr
D1_2yr = OrRd

colors=c(D1_2yr, D2_2yr, D1_PN42, D2_PN42)

STOP_CODON = t_gene_model %>% filter(type=='UTR_3p') %>% select(x) %>% min()

g1=ggplot(t)
g1=g1+geom_area(aes(x=t_pos,y=log((cvg)+1)/log(10),fill=sample,group=sample),color="white")+
    theme_bw(base_size=10)+
    theme(axis.text.x=element_blank())+
    theme(axis.ticks.x=element_blank())+
    theme(axis.line.x=element_blank())+
    theme()+
    facet_grid(t_celltype~.)+
    scale_x_continuous("")+
    scale_y_continuous("log(cvg)")+
    scale_fill_manual(values=colors)+
    theme(strip.background=element_blank())+
    geom_vline(xintercept=STOP_CODON,linetype=4)+
    theme(legend.position="none")+
    theme(panel.grid=element_blank())+
    ggtitle(gene_name)
    
g2=ggplot(t_gene_model)+
    geom_polygon(aes(x=x, y=y, group = id), fill="steelblue")+
    theme_classic()+
    theme(axis.text.x=element_text(size=8), axis.line.x=element_blank())+
    scale_x_continuous("")+scale_y_continuous("")+
    theme(axis.line.y=element_blank())+
    theme(axis.text.y=element_blank())+
    theme(axis.ticks.y=element_blank())+
    theme(plot.margin = unit(c(-1,.2,0,0),"cm"))+
    theme(panel.background=element_blank())
    
gA = ggplotGrob(g1)
gB = ggplotGrob(g2)
gB$widths = gA$widths

gC <- grid.rect(gp=gpar(col="white"))
gC.widths=gA$widths

pdf(fn_output, width=2.5, height=4.5)
grid.arrange(gA, gC, gB, ncol=1, heights=c(.94,.03,.03))
dev.off()
