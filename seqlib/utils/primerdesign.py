import argparse
import primer3
from fastahack import FastaHack
import os
import pdb


class PrimerDesign(object):
    
    def __init__(self, **kwargs):
        
        mispriming_species = kwargs.get("mispriming_lib", "rodent")
        
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
        self.col_order = ['name',
                          'locus',
                          'strand',
                          'contig',
                          'PRIMER_LEFT_SEQUENCE',
                          'PRIMER_RIGHT_SEQUENCE',
                          'PRIMER_LEFT_GC_PERCENT',
                          'PRIMER_RIGHT_GC_PERCENT',
                          'PRIMER_LEFT_TM',
                          'PRIMER_RIGHT_TM']

    def load_mispriming_lib(self, fn_mispriming_fa):
        fa = FastaHack(fn_mispriming_fa)
        self.mispriming_dict = {}
        for contig in fa.names:
            self.mispriming_dict[contig] = fa.get_sequence(contig)[:]

    def get_primers(self, **kwargs):
        """
        addditional_args that could be added
            #SEQUENCE_TARGET [37,1] ###NOTE< secondd number is len
            #SEQUENCE_EXCLUDED_REGION 37,89
        """
        
        seq = kwargs.get("seq")
        contig = kwargs.get("contig", "-")
        start = kwargs.get("start", 0)
        end = kwargs.get("end", 0)
        strand = kwargs.get("strand", "+")
        name_prefix = kwargs.get("name_prefix", "")
        nested = kwargs.get("nested", False)
        n_max = kwargs.get("n_max", 1)
        
        target_size_range = kwargs.get("target_size_range", [150,250])
        target_size_range = [target_size_range]
        
        _seq_args = kwargs.get("seq_args", {})
        _global_args = kwargs.get("global_args", {})
        
        ret_val_lambdas = {}
        ret_val_lambdas["PRIMER_PAIR_{n}_PRODUCT_SIZE"] = kwargs.get("PRIMER_PRODUCT_SIZE", lambda x: int(x))
        
        if strand:
            length = end-start
            ret_val_lambdas["PRIMER_LEFT_{n}"] = kwargs.get("PRIMER_LEFT", lambda x: int(x[0])+start)
            ret_val_lambdas["PRIMER_RIGHT_{n}"] = kwargs.get("PRIMER_RIGHT", lambda x: end-(length-int(x[0]))+1)
        else:
            length = len(seq)
            ret_val_lambdas["PRIMER_LEFT_{n}"] = kwargs.get("PRIMER_LEFT", lambda x: end-int(x[0]))
            ret_val_lambdas["PRIMER_RIGHT_{n}"] = kwargs.get("PRIMER_RIGHT", lambda x: start+(length-int(x[0]))-1)
        
        seq_args = {"SEQUENCE_ID":"ID_%s"%name_prefix,
                    "SEQUENCE_TEMPLATE":seq}
        
        seq_args.update(_seq_args)

        global_args = {
            'PRIMER_OPT_SIZE': 20,
            'PRIMER_PICK_INTERNAL_OLIGO': 1,
            'PRIMER_INTERNAL_MAX_SELF_END': 8,
            'PRIMER_MIN_SIZE': 18,
            'PRIMER_MAX_SIZE': 25,
            'PRIMER_OPT_TM': 60.0,
            'PRIMER_MIN_TM': 57.0,
            'PRIMER_MAX_TM': 63.0,
            'PRIMER_MIN_GC': 20.0,
            'PRIMER_MAX_GC': 80.0,
            'PRIMER_MAX_POLY_X': 100,
            'PRIMER_INTERNAL_MAX_POLY_X': 100,
            'PRIMER_SALT_MONOVALENT': 50.0,
            'PRIMER_DNA_CONC': 50.0,
            'PRIMER_MAX_NS_ACCEPTED': 0,
            'PRIMER_MAX_SELF_ANY': 12,
            'PRIMER_MAX_SELF_END': 8,
            'PRIMER_PAIR_MAX_COMPL_ANY': 12,
            'PRIMER_PAIR_MAX_COMPL_END': 8,
            'PRIMER_PRODUCT_SIZE_RANGE': target_size_range
        }
        global_args.update(_global_args)

        if len(seq)<100:
            return []
        p3_ret = primer3.designPrimers(seq_args, 
                                       global_args=global_args, 
                                       misprime_lib=self.mispriming_dict)
         
        if int(p3_ret['PRIMER_PAIR_NUM_RETURNED']) == 0:
            return []
        
        ret_dicts = []
        i_count = 0
        
        for i_primer_ret in range(n_max):
            ret_dict = {}
            ret_dicts.append(ret_dict)
            p_left_key = 'PRIMER_LEFT_{n}'.format(n=i_primer_ret)
            p_right_key = 'PRIMER_RIGHT_{n}'.format(n=i_primer_ret)
            
            if not p_left_key in p3_ret:
                #break if fewer than requested found
                break

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
            ret_dict['strand'] = strand
            ret_dict['name'] = "{prefix}".format(prefix=name_prefix)
            ret_dict['locus'] = "{contig}_{s}_{e}".format(contig=contig,
                                                          s=seg_s,
                                                          e=seg_e)
        
            if nested:
                add_args = dict(additional_args)
                excl1 = (0, int(p3_ret[p_left_key][0]+p3_ret[p_left_key][1]))
                excl2 = (int(p3_ret[p_right_key][0]),int(len(seq)-p3_ret[p_right_key][0]))
                t_range  = [max(50, ret_dict['PRIMER_PAIR_PRODUCT_SIZE']-100), target_size_range[1]]
                add_args.update({"SEQUENCE_EXCLUDED_REGION":[excl1, excl2],
                                 "target_size_range":t_range})
                                 
                nested = self.get_primers(seq=seq,
                                          contig=contig,
                                          start=start,
                                          end=end,
                                          name_prefix="%s_nested"%name_prefix,
                                          strand=strand,
                                          nested=False,
                                          additional_args=add_args,
                                          PRIMER_LEFT=ret_val_lambdas["PRIMER_LEFT_{n}"],
                                          PRIMER_RIGHT=ret_val_lambdas["PRIMER_RIGHT_{n}"])
                ret_dicts.extend(nested)
                i_count += 1

        return ret_dicts

if __name__=="__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", required=True)
    o = parser.parse_args()
    
