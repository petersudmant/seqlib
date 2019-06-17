import argparse
import xmltodict
import os
import json
import pdb


def make_json_config(xml_dict, naming_pattern, fn_add_to_config):
    
    add = {}
    if fn_add_to_config:
        add = json.load(open(fn_add_to_config))

    
    cwd = os.getcwd()
    STUDY_ACC=cwd.split("/")[-1].split("_")[0] 

    j_out = {"study_accession":STUDY_ACC}
    j_out["naming_pattern"] = naming_pattern
    #j_out["samples"] = {}
    j_out["experiments"] = {}

    for d in xml_dict["EXPERIMENT_PACKAGE_SET"]["EXPERIMENT_PACKAGE"]:
        """
        [u'EXPERIMENT', u'SUBMISSION', u'STUDY', u'SAMPLE', u'RUN_SET']
        """
        STUDY_accession = d["STUDY"]["@accession"]
        #print STUDY_accession, STUDY_ACC
        #assert STUDY_accession==STUDY_ACC
        
        experiment_accession = d['EXPERIMENT']["@accession"]
        sample_accession = d["SAMPLE"]["@accession"]
        sample_alias = d["SAMPLE"]["@alias"]
        sample_title = d["SAMPLE"]["TITLE"]
        #iterate over the 
        #pdb.set_trace()
        run_accessions = []
        print(sample_accession, experiment_accession, sample_alias, sample_title)
        if type(d["RUN_SET"]["RUN"]) == list:
            for st in d["RUN_SET"]["RUN"]:
                run_accessions.append(st["@accession"])
        else:
            run_accessions = [d["RUN_SET"]["RUN"]["@accession"]]

        #j_out["samples"][sample_accession] = {"sample_accession":sample_accession,
        #                                      "sample_alias":sample_alias,
        #                                      "sample_title":sample_title,
        #                                      "experiment_accession":experiment_accession,
        #                                      "run_accessions":run_accessions
        #                                      }
        j_out["experiments"][experiment_accession] = {"sample_accession":sample_accession,
                                              "sample_alias":sample_alias,
                                              "experiment_accession":experiment_accession,
                                              "sample_title":sample_title,
                                              "run_accessions":run_accessions
                                              }
        if "SAMPLE_ATTRIBUTES" in d["SAMPLE"]:
            for attribute in d["SAMPLE"]["SAMPLE_ATTRIBUTES"]["SAMPLE_ATTRIBUTE"]:
                attr = attribute["TAG"]
                val = attribute["VALUE"]
                j_out["experiments"][experiment_accession][attr] = val
                #j_out["samples"][sample_accession][attr] = val

        if "SAMPLE_NAME" in d["SAMPLE"]:
            for attribute_pair in d["SAMPLE"]["SAMPLE_NAME"].items():
                attr, val = attribute_pair
                #print attr, val
                #j_out["samples"][sample_accession][attr] = val
                j_out["experiments"][experiment_accession][attr] = val

        if "samples" in add.keys():
            for key, value in add["samples"][sample_accession].items():
                #j_out["samples"][sample_accession][key] = value
                j_out["experiments"][experiment_accession][attr] = val

    FOUT = open("%s/config.json"%cwd,'w')
    FOUT.write(json.dumps(j_out, indent=4, separators=(",", ": ")))
    FOUT.close()

   

def manual_parse(xml_dict, manually_currate = []):

    cwd = os.getcwd()
    STUDY_ACC=cwd.split("/")[-1].split("_")[0] 

    j_out = {"study_accession":STUDY_ACC}
    #j_out["naming_pattern"] = ["sample_accession", "tissue"]
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
        
        j_out["samples"][sample_accession] = {"sample_accession":sample_accession}

        if "SAMPLE_ATTRIBUTES" in manually_currate:
            #print "OK"
            #print d["SAMPLE"]["DESCRIPTION"]
            #attr = raw_input("attribute: ")
            #val = raw_input("value: ")
            #j_out["samples"][sample_accession][attr] = val
            j_out["samples"][sample_accession]["description"] = d["SAMPLE"]["DESCRIPTION"]
             
        else:
            for attribute in d["SAMPLE"]["SAMPLE_ATTRIBUTES"]["SAMPLE_ATTRIBUTE"]:
                attr = attribute["TAG"]
                val = attribute["VALUE"]
                j_out["samples"][sample_accession][attr] = val
    F = open("manual_add.json",'w') 
    F.write(json.dumps(j_out, indent=4, separators=(",", ": ")))
    F.close()

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--fn_xml",default="SraExperimentPackage.xml")
    parser.add_argument("--manual_curration", default=None)
    parser.add_argument("--add_to_config", default=None)
    parser.add_argument("--naming_pattern", default="sample_accession:tissue")
    o = parser.parse_args()

    xml = open(o.fn_xml).read()
    xml_dict = xmltodict.parse(xml)

    if o.manual_curration:
        manual_parse(xml_dict, o.manual_curration.split(":"))
    else:
        make_json_config(xml_dict, o.naming_pattern.split(":"), o.add_to_config)
    

