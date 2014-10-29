import argparse
import xmltodict
import os
import json
import pdb

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--fn_xml",default="SraExperimentPackage.xml")
    o = parser.parse_args()

    xml = open(o.fn_xml).read()
    xml_dict = xmltodict.parse(xml)

    cwd = os.getcwd()
    STUDY_ACC=cwd.split("/")[-1].split("_")[0] 

    j_out = {"study_accession":STUDY_ACC}
    j_out["symlink_pattern"] = ["accession","disease state","tissue"]
    j_out["samples"] = {}

    for d in xml_dict["EXPERIMENT_PACKAGE_SET"]["EXPERIMENT_PACKAGE"]:
        """
        [u'EXPERIMENT', u'SUBMISSION', u'STUDY', u'SAMPLE', u'RUN_SET']
        """
        STUDY_accession = d["STUDY"]["@accession"]
        assert STUDY_accession==STUDY_ACC
        
        sample_accession = d["SAMPLE"]["@accession"]
        sample_alias = d["SAMPLE"]["@alias"]
        sample_title = d["SAMPLE"]["TITLE"]
        print(type(d["RUN_SET"]["RUN"])), len(d["RUN_SET"]["RUN"])
        #iterate over the 
        run_accessions = []
        if type(d["RUN_SET"]["RUN"]) == list:
            for st in d["RUN_SET"]["RUN"]:
                run_accessions.append(st["@accession"])
        else:
            run_accessions = [d["RUN_SET"]["RUN"]["@accession"]]

        j_out["samples"][sample_accession] = {"sample_accession":sample_accession,
                                              "sample_alias":sample_alias,
                                              "sample_title":sample_title,
                                              "run_accessions":run_accessions
                                              }
        for attribute in d["SAMPLE"]["SAMPLE_ATTRIBUTES"]["SAMPLE_ATTRIBUTE"]:
            attr = attribute["TAG"]
            val = attribute["VALUE"]
            j_out["samples"][sample_accession][attr] = val
    
    FOUT = open("%s/config.json"%cwd,'w')
    FOUT.write(json.dumps(j_out, indent=4, separators=(",", ": ")))
    FOUT.close()

    
    

