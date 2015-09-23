library(grid)
library(ggplot2)
library(wesanderson)
library(R.utils)
library(RColorBrewer)

args = cmdArgs()

fn_input  = cmdArg("fn_input")
fn_output = cmdArg("fn_output")
fn_output_2 = cmdArg("fn_output_2")
x = cmdArg("x")

t=read.table(fn_input,header=T,sep="\t")

#t$dataset = factor(t$dataset, c("Brawand_2011", "Merkin_2012", "ENCODE_set1", "ENCODE_set2"))
