import argparse
import xmltodict
import os
import json

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--fn_xml",default="SraExperimentPackage.xml")
    o = parser.parse_args()

    xml = open(o.fn_xml).read()
    xml_dict = xmltodict.parse(xml)

    cwd = os.getcwd()
    STUDY=cwd.split("/")[-1].split("_")[0] 

    j_out = {"STUDY":STUDY}
    j_out["samples"] = {}

    for d in xml_dict["EXPERIMENT_PACKAGE_SET"]["EXPERIMENT_PACKAGE"]:
        """
        [u'EXPERIMENT', u'SUBMISSION', u'STUDY', u'SAMPLE', u'RUN_SET']
        """
        STUDY_accession = d["STUDY"]["@accession"]
        assert STUDY_accession==STUDY

        sample_accession = d["SAMPLE"]["@accession"]
        sample_alias = d["SAMPLE"]["@alias"]
        sample_title = d["SAMPLE"]["TITLE"]

        j_out["samples"][sample_accession] = {"accession":sample_accession,
                                              "alias":sample_alias,
                                              "title":sample_title}
        for attribute in d["SAMPLE"]["SAMPLE_ATTRIBUTES"]["SAMPLE_ATTRIBUTE"]:
            attr = attribute["TAG"]
            val = attribute["VALUE"]
            j_out["samples"][sample_accession][attr] = val
    
    FOUT = open("%s/config.json"%cwd,'w')
    FOUT.write(json.dumps(j_out, indent=4, separators=(",", ": ")))
    FOUT.close()

    
    

