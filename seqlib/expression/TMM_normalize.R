library(R.utils)
library(edgeR)

args = cmdArgs()
fn_input = cmdArg("fn_input")
fn_output = cmdArg("fn_output")
counts_start_col = cmdArg("exp_start_col",default=6)

counts_t = read.table(fn_input,header=T,sep="\t")
counts_end_col = dim(counts_t)[2]

y=DGEList(counts=counts_t[,counts_start_col:counts_end_col])
rownames(y)=counts_t$gene_name
y=calcNormFactors(y)
adj=t(t(y$counts)/(y$samples$norm.factors * y$samples$lib.size))

write.table(adj, file=fn_output, sep="\t", col.names=NA)
