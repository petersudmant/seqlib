import argparse
import primer3
from fastahack import FastaHack
import os
import pdb


class PrimerDesign(object):
    
    def __init__(self, **kwargs):
        
        mispriming_species = kwargs.get("mispriming_lib", "rodent")
        target_size_range = kwargs.get('size_range', [150-250])
        self.target_size_range = "-".join(["%s"%s for s in target_size_range])

        cwd = os.path.dirname(os.path.abspath(__file__))
        mispriming_libs = {"human":"%s/mispriming_libs/HUMAN_AND_SIMPLE.fa"%cwd, 
                           "mouse":"%s/mispriming_libs/RODENT_AND_SIMPLE.fa"%cwd,
                           "rodent":"%s/mispriming_libs/RODENT_AND_SIMPLE.fa"%cwd}
        
        assert mispriming_species in mispriming_libs
        
        self.load_mispriming_lib(mispriming_libs[mispriming_species]) 
        
        self.return_vals = ["PRIMER_LEFT_{n}",
                            "PRIMER_LEFT_{n}_SEQUENCE",
                            "PRIMER_LEFT_{n}_GC_PERCENT",
                            "PRIMER_LEFT_{n}_TM",
                            "PRIMER_RIGHT_{n}",
                            "PRIMER_RIGHT_{n}_SEQUENCE",
                            "PRIMER_RIGHT_{n}_GC_PERCENT",
                            "PRIMER_RIGHT_{n}_TM",
                            "PRIMER_PAIR_{n}_PRODUCT_SIZE"]
        

    def load_mispriming_lib(self, fn_mispriming_fa):
        fa = FastaHack(fn_mispriming_fa)
        self.mispriming_dict = {}
        for contig in fa.names:
            self.mispriming_dict[contig] = fa.get_sequence(contig)[:]

    def get_primers(self, **kwargs):
        """
        addditional_args that could be added
            #SEQUENCE_TARGET 37,89
            #SEQUENCE_EXCLUDED_REGION 37,89
        """
        
        seq = kwargs.get("seq")
        contig = kwargs.get("contig")
        self.start = kwargs.get("start")
        self.end = kwargs.get("end")
        self.strand = kwargs.get("strand")
        additional_args = kwargs.get("additional_args", {})
        name_prefix = kwargs.get("name_prefix", "")
        nested = kwargs.get("nested", False)
        n_max = kwargs.get("n_max", 1)

        ret_val_lambdas = {}
        ret_val_lambdas["PRIMER_PAIR_{n}_PRODUCT_SIZE"] = kwargs.get("PRIMER_PRODUCT_SIZE", lambda x: int(x))
        
        if self.strand:
            length = self.end-self.start
            ret_val_lambdas["PRIMER_LEFT_{n}"] = kwargs.get("PRIMER_LEFT", lambda x: int(x[0])+self.start)
            ret_val_lambdas["PRIMER_RIGHT_{n}"] = kwargs.get("PRIMER_RIGHT", lambda x: self.end-(length-int(x[0]))+1)
        else:
            length = len(seq)
            ret_val_lambdas["PRIMER_LEFT_{n}"] = kwargs.get("PRIMER_LEFT", lambda x: self.end-int(x[0]))
            ret_val_lambdas["PRIMER_RIGHT_{n}"] = kwargs.get("PRIMER_RIGHT", lambda x: self.start+(length-int(x[0]))-1)
        
        seq_args = {"SEQUENCE_ID":"",
                    "SEQUENCE_TEMPLATE":seq, 
                    "PRIMER_PRODUCT_SIZE_RANGE": self.target_size_range}
        
        seq_args.update(additional_args)
        p3_ret = primer3.designPrimers(seq_args, self.mispriming_dict)
        #try:
        #except:
        #    return []
    
        if int(p3_ret['PRIMER_PAIR_NUM_RETURNED']) == 0:
            return []
        
        ret_dicts = []

        for i_primer_ret in range(n_max):
            ret_dict = {}
            ret_dicts.append(ret_dict)
            p_left_key = 'PRIMER_LEFT_{n}'.format(n=i_primer_ret)
            p_right_key = 'PRIMER_RIGHT_{n}'.format(n=i_primer_ret)

            seg_s=ret_val_lambdas['PRIMER_LEFT_{n}'](p3_ret[p_left_key])
            seg_e=ret_val_lambdas['PRIMER_RIGHT_{n}'](p3_ret[p_right_key])
            seg_s,seg_e = sorted([seg_s,seg_e])
            
            for ret_val in self.return_vals:
                p3_key = ret_val.format(n=i_primer_ret)
                key = p3_key.replace("_{n}".format(n=i_primer_ret),"")
                if ret_val in ret_val_lambdas:
                    ret_dict[key] = ret_val_lambdas[ret_val](p3_ret[p3_key])
                else:
                    ret_dict[key] = p3_ret[p3_key]
        
            ret_dict['contig'] = contig
            ret_dict['strand'] = self.strand
            ret_dict['name'] = "{prefix}".format(prefix=name_prefix)
            ret_dict['locus'] = "{contig}_{s}_{e}".format(contig=contig,
                                                          s=seg_s,
                                                          e=seg_e)
            if nested:
                add_args = dict(additional_args)
                excl1 = (0, int(p3_ret[p_left_key][0]+p3_ret[p_left_key][1]))
                excl2 = (int(p3_ret[p_right_key][0]),int(len(seq)-p3_ret[p_right_key][0]))
                add_args.update({"SEQUENCE_EXCLUDED_REGION":[excl1,excl2]})
                nested = self.get_primers(seq=seq,
                                          contig=contig,
                                          start=self.start,
                                          end=self.end,
                                          name_prefix="%s_nested"%name_prefix,
                                          strand=self.strand,
                                          nested=False,
                                          additional_args=add_args,
                                          PRIMER_LEFT=ret_val_lambdas["PRIMER_LEFT_{n}"],
                                          PRIMER_RIGHT=ret_val_lambdas["PRIMER_RIGHT_{n}"])
                ret_dicts.extend(nested)

        return ret_dicts

if __name__=="__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", required=True)
    o = parser.parse_args()
    
