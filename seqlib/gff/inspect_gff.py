import BCBio.GFF as GFF
import argparse
import pdb
import time
import pprint

if __name__=="__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--fn_gff")
    o = parser.parse_args()

    GFF_examiner = GFF.GFFExaminer()
    t = time.time()
    pprint.pprint(GFF_examiner.parent_child_map(open(o.fn_gff)))
    print "----------------------------"
    pprint.pprint(GFF_examiner.available_limits(open(o.fn_gff)))
    #print time.time()-t
