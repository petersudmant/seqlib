import sys


if __name__=="__main__":

    for l in sys.stdin:
        if l[0] == "#":
            continue
        l=l.rstrip()
        sl = l.split()
        if "geneID" in sl[-1]:
            d={e.split("=")[0]:e.split("=")[1] for e in sl[-1].split(";")}
            d["gene_id"] = d['geneID']
            sl[-1] = ";".join(["%s=%s"%(k,v) for k,v in d.iteritems()])
        print("\t".join(sl))     

