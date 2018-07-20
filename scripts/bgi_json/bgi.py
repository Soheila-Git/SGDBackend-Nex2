from src.models import DBSession, Base, Colleague, ColleagueLocus, Dbentity, Goextension, Locusdbentity, LocusUrl, LocusAlias, Dnasequenceannotation, So, Locussummary, Phenotypeannotation, PhenotypeannotationCond, Phenotype, Goannotation, Go, Goslimannotation, Goslim, Apo, Straindbentity, Strainsummary, Reservedname, GoAlias, Goannotation, Referencedbentity, Referencedocument, Referenceauthor, ReferenceAlias, Chebi
from sqlalchemy import create_engine, and_, inspect
import os
import json
import re
import time
import sys
from random import randint
#from pycallgraph import PyCallGraph
#from pycallgraph.output import GraphvizOutput
from datetime import datetime
from threading import Thread
import concurrent.futures

engine = create_engine(os.environ['NEX2_URI'], pool_recycle=3600)
DBSession.configure(bind=engine)
Base.metadata.bind = engine


def get_goext():
    ext_data = DBSession.query(Goextension).all()
    if ext_data:
        obj = {}
        for item in ext_data:
            if item.ro.roid == 'NTR:0000002' or item.ro.roid == 'RO:0002092':
                if item.annotation_id in obj:
                    obj[item.annotation_id].append(item.annotation_id)
                else:
                    obj[item.annotation_id] = []
        return obj
    return None


def get_goannotation_data():
    flag_ship = {}
    data = []
    metadata = {}
    main_obj = {}
    xcontainer = {
        'GO:0000749': 'GO:0000749',
        'GO:0006974': 'GO:0006974',
        'GO:0006995': 'GO:0006995',
        'GO:0009408': 'GO:0009408',
        'GO:0033554': 'GO:0033554',
        'GO:0034198': 'GO:0034198',
        'GO:0034599': 'GO:0034599',
        'GO:0034605': 'GO:0034605',
        'GO:0034614': 'GO:0034614',
        'GO:0042149': 'GO:0042149',
        'GO:0070301': 'GO:0070301',
        'GO:0071311': 'GO:0071311',
        'GO:0071333': 'GO:0071333',
        'GO:0071361': 'GO:0071361',
        'GO:0071400': 'GO:0071400',
        'GO:0071456': 'GO:0071456',
        'GO:0071470': 'GO:0071470',
        'GO:0071474': 'GO:0071474',
        'GO:0097185': 'GO:0097185'
    }
    ext_temp = get_goext()
    go_annot_data = DBSession.query(Goannotation).filter(
        Goannotation.annotation_type != 'computational',
        Goannotation.eco_id != 236340).all()
    for item in go_annot_data:

        obj = {}
        ro_flag = ext_temp.get(item.annotation_id, None)
        go_flag = xcontainer.get(item.go.goid, None)
        if ro_flag is None and go_flag is None and item.go.go_namespace == "cellular component":
            obj = {
                "geneId":
                    "SGD:" + str(item.dbentity.sgdid),
                "evidence": {},
                "whenExpressedStage":
                    "N/A",
                "assay":
                    "MMO:0000642",
                "dateAssigned":
                    item.date_created.strftime("%Y-%m-%dT%H:%m:%S-00:00"),
                "wildtypeExpressionTermIdentifiers": {
                    "anatomicalStructureTermId": "N/A"
                }
            }
            if item.reference.pmid:
                obj["evidence"]["pubMedId"] = "PMID:" + str(item.reference.pmid)
            else:
                obj["evidence"][
                    "modPublicationId"] = "SGD:" + str(item.reference.sgdid)
            if obj not in data:
                data.append(obj)

    if data:
        main_obj['data'] = data
        main_obj['metaData'] = {
            "dataProvider": [{
                "crossReference": {
                    "id": "SGD"
                },
                "type": "curated"
            }],
            "dateProduced":
                datetime.utcnow().strftime("%Y-%m-%dT%H:%m:%S-00:00")
        }

        return main_obj
    return None


def get_expression_data():
    with concurrent.futures.ProcessPoolExecutor(max_workers=128) as executor:
        go_data = get_goannotation_data()
        if go_data:
            fileStr = './scripts/bgi_json/data_dump/SGD.1.0.0.4_expression_' + str(randint(0, 1000)) + '.json'
            with open(fileStr, 'w+') as res_file:
                res_file.write(json.dumps(go_data))


# populate text file with sgdis to be used to retrieve panther data
def get_sgdids_for_panther():
    new_data = Locusdbentity.get_s288c_genes()
    temp = []
    for loc in new_data:
        temp.append(loc.sgdid)
    result = json.dumps(temp, ensure_ascii=False)
    with open('./scripts/bgi_json/data_dump/sgd_ids_for_panther.txt',
              'w+') as res_file:
        res_file.write(
            result.replace('"', '').replace('[', '').replace(']', ''))


# pair pantherIds with corresponding sgdids
def get_panther_sgdids():
    data_dict = {}
    with open('./scripts/bgi_json/data_dump/panther/panther_search_results.json') as json_data_file:
        json_data = json.load(json_data_file)
        for item in json_data:
            temp_str = ','.join(map(str, item))
            reg_pattern = r'(SGD=S\d+)|(PTHR\d+)'
            reg_result = sorted(list(set(re.findall(reg_pattern, temp_str))))
            if(len(reg_result) > 1):
                item_str1 = ''.join(reg_result[0])
                item_str2 = ''.join(reg_result[1]).split("=")
                data_dict[item_str2[1]] = item_str1
            elif(len(reg_result) == 1):
                item_str1 = ''.join(reg_result[0]).split("=")
                data_dict[item[1]] = None
            else:
                continue

        return data_dict


# pupulate json file with basic gene infromation(bgi)
def get_bgi_data(soFlag=False):
    combined_list = combine_panther_locus_list(get_panther_sgdids(), Locusdbentity.get_s288c_genes())
    print("computing " + str(len(combined_list)) + " genes")
    result = []
    if(len(combined_list) > 0):
        with concurrent.futures.ProcessPoolExecutor(max_workers=128) as executor:
            for item_key in combined_list:
                obj = {
                    "crossReferences":
                        [],
                    "primaryId":
                        "",
                    "symbol":
                        "",
                    "genomeLocations": [{
                        "startPosition": 0,
                        "chromosome": "",
                        "assembly": "R64-2-1",
                        "endPosition": 0,
                        "strand": ""
                    }],
                    "soTermId":
                        "",
                    "taxonId":
                        "NCBITaxon:559292",
                    "synonyms":
                        [],
                    "geneSynopsis":
                        ""
                }
                item = combined_list[item_key]["locus_obj"]
                temp_itm = ["gene"]
                temp_itm.append("gene/references")
                temp_itm.append("homepage")
                if(item.has_expression):
                    temp_itm.append("gene/expression")
                    temp_itm.append("gene/spell")
                if(item.has_interaction):
                    temp_itm.append("gene/interaction")

                obj["crossReferences"].append({"id": "SGD:"+item.sgdid, "pages": temp_itm})
                item_panther = combined_list[item_key]["panther_id"]
                locus_alias_data = DBSession.query(LocusAlias).filter(
                    LocusAlias.locus_id == item.dbentity_id).all()

                if(len(locus_alias_data) > 0):
                    dna_seq_annotation_obj = DBSession.query(
                        Dnasequenceannotation).filter(
                            Dnasequenceannotation.dbentity_id == item.dbentity_id,
                            Dnasequenceannotation.taxonomy_id == 274901,
                            Dnasequenceannotation.dna_type == "GENOMIC").all()

                    if(len(dna_seq_annotation_obj) > 0):
                        strnd = ""
                        if dna_seq_annotation_obj[0].strand == "0":
                            strnd = "."
                        else:
                            strnd = dna_seq_annotation_obj[0].strand
                        chromosome = dna_seq_annotation_obj[0].contig.display_name.split(" ")
                        obj["genomeLocations"][0]["startPosition"] = dna_seq_annotation_obj[0].start_index
                        obj["genomeLocations"][0]["endPosition"] = dna_seq_annotation_obj[0].end_index
                        obj["genomeLocations"][0]["strand"] = strnd
                        obj["genomeLocations"][0]["startPosition"] = dna_seq_annotation_obj[0].start_index
                        obj["genomeLocations"][0]["chromosome"] = "chr"+chromosome[1]
                        if dna_seq_annotation_obj[0].so.so_id == 263757:
                            obj["soTermId"] = "SO:0001217"
                        else:
                            obj["soTermId"] =  dna_seq_annotation_obj[0].so.soid
                    mod_locus_alias_data = get_locus_alias_data(locus_alias_data, item.dbentity_id, item)

                    for mod_item in mod_locus_alias_data:
                        mod_value = mod_locus_alias_data.get(mod_item)
                        if (type(mod_value) is list):
                            if (mod_locus_alias_data.get("aliases") is not None):
                                obj["synonyms"] = mod_locus_alias_data.get(
                                    "aliases")

                        else:
                            if (mod_value.get("secondaryIds") is not None):
                                temp_sec_item = mod_value.get("secondaryIds")
                                if(len(temp_sec_item) > 0):
                                    if(item.name_description is not None):
                                        obj["name"] = item.name_description
                                    if(len(temp_sec_item) > 1):
                                        obj["secondaryIds"] = [str(x) for x in temp_sec_item]
                                    else:
                                        if(len(temp_sec_item) == 1):
                                            obj["secondaryIds"] = [str(temp_sec_item[0])]
                            if (mod_value.get("crossReferences") is not None):
                                temp_cross_item = mod_value.get("crossReferences")
                                if(len(temp_cross_item) > 1):
                                    for x_ref in temp_cross_item:
                                        obj["crossReferences"].append({"id": str(x_ref), "pages": ["generic_cross_reference"]})
                                else:
                                    if(len(temp_cross_item) == 1):
                                        obj["crossReferences"].append({"id": str(temp_cross_item[0]), "pages": ["generic_cross_reference"]})
                                        #obj["crossReferences"] = [str(temp_cross_item[0])]
                    if(item_panther is not None):
                        obj["crossReferences"].append({"id": "PANTHER:" + item_panther})
                        #obj["crossReferences"].append("PANTHER:" + item_panther)
                        obj["primaryId"] = "SGD:" + item.sgdid
                        item = combined_list[item_key]["locus_obj"]
                        obj["geneSynopsis"] = item.description
                        obj["symbol"] = item.gene_name if item.gene_name is not None else item.systematic_name
                        obj["synonyms"].append(item.systematic_name)
                        result.append(obj)

                    else:
                        obj["primaryId"] = "SGD:" + item.sgdid
                        item = combined_list[item_key]["locus_obj"]
                        obj["geneSynopsis"] = item.description
                        obj["symbol"] = item.gene_name if item.gene_name is not None else item.systematic_name
                        obj["synonyms"].append(item.systematic_name)
                        result.append(obj)
            if(len(result) > 0):
                output_obj = {
                    "data": result,
                    "metaData": {
                        "dataProvider": [
                            {
                                "crossReference": {
                                    "id": "SGD"
                                }
                            }
                        ],
                        "dateProduced": datetime.utcnow().strftime("%Y-%m-%dT%H:%m:%S-00:00")
                    }
                }
                fileStr = './scripts/bgi_json/data_dump/SGD.1.0.0.4_basicGeneInformation_' + str(randint(0, 1000)) + '.json'
                with open(fileStr, 'w+') as res_file:
                    res_file.write(json.dumps(output_obj))

# create single gene list containing genes with pantherids and genes without pantherids
def combine_panther_locus_list(panther_list, locus_list):
    combined_list = {}
    if(len(panther_list) > 0 and len(locus_list) > 0):
        for item in locus_list:
            obj = {
                "panther_id": "",
                "locus_obj": ""
            }
            if(panther_list.get(item.sgdid) is not None):
                obj["panther_id"] = panther_list.get(item.sgdid)
                obj["locus_obj"] = item
                combined_list[item.dbentity_id] = obj
            else:
                obj["panther_id"] = None
                obj["locus_obj"] = item
                combined_list[item.dbentity_id] = obj
    return combined_list


# helper method to get locus_alias data
def get_locus_alias_data(locus_alias_list, id, item_obj):
    data_container = {}
    aliases = []
    aliases_types = ["Uniform", "Non-uniform"]
    aliases_types_other = ["SGDID Secondary", "UniProtKB ID", "Gene ID"]
    obj = {
        "secondaryIds": [],
         "crossReferences": []
        }
    flag = False
    for item in locus_alias_list:
        if(item_obj):
            obj["crossReferences"].append
        if(item.alias_type in aliases_types):
            aliases.append(item.display_name)
        if(item.alias_type == "SGDID Secondary"):
            obj["secondaryIds"].append(item.source.display_name+":" + item.display_name)
            flag = True
        if (item.alias_type == "UniProtKB ID"):
            obj["crossReferences"].append("UniProtKB:" + item.display_name)
            flag = True
        if (item.alias_type == "Gene ID" and item.source.display_name == 'NCBI'):
            obj["crossReferences"].append("NCBI_Gene:" + item.display_name)
            flag = True

    if(flag):
        data_container[id] = obj
    data_container["aliases"] = aliases
    return data_container


def get_phenotype_data():
    _data = DBSession.query(Phenotypeannotation).all()
    result = []
    print("computing " + str(len(_data)) + " phenotypes")
    with concurrent.futures.ProcessPoolExecutor(max_workers=4) as executor:
        for item in _data:
            obj = {
                "objectId": "",
                "phenotypeTermIdentifiers": [],
                "phenotypeStatement": "",
                "dateAssigned": "",
                "evidence": {}
            }
            if item.reference.pmid:
                obj["evidence"]["pubMedId"] = "PMID:" + str(item.reference.pmid)
            else:
                obj["evidence"]["pubModId"] = "SGD:" + str(item.reference.sgdid)
            if item.phenotype.qualifier:
                pString = item.phenotype.qualifier.display_name
                obj["phenotypeTermIdentifiers"].append({
                "termId": str(item.phenotype.qualifier.apoid),
                "termOrder": 1
                })
                if item.phenotype.observable:
                    pString = pString + " " +item.phenotype.observable.display_name
                    obj["phenotypeTermIdentifiers"].append({
                    "termId": str(item.phenotype.observable.apoid),
                    "termOrder": 2
                    })

            else:
                if item.phenotype.observable:
                    pString = item.phenotype.observable.display_name
                    obj["phenotypeTermIdentifiers"].append({
                    "termId": str(item.phenotype.observable.apoid),
                    "termOrder": 1
                    })
            obj["objectId"] = "SGD:" + str(item.dbentity.sgdid)
            obj["phenotypeStatement"] = pString
            obj["dateAssigned"] = item.date_created.strftime(
                "%Y-%m-%dT%H:%m:%S-00:00")
            result.append(obj)

        if len(result) > 0:
            output_obj = {
                "data": result,
                "metaData": {
                    "dataProvider": [{
                        "crossReference": {
                            "id": "SGD",
                            "pages": ["homepage"]
                        },
                        "type": "curated"
                    }],
                    "dateProduced":
                        datetime.utcnow().strftime("%Y-%m-%dT%H:%m:%S-00:00")
                }
            }
            fileStr = './scripts/bgi_json/data_dump/SGD.1.0.0.4_2_phenotype.json'
            with open(fileStr, 'w+') as res_file:
                res_file.write(json.dumps(output_obj))



# entry point
if __name__ == '__main__':
    '''print "--------------start computing data--------------"
    start_time = time.time()
    get_bgi_data()
    time_taken = "time taken: " + ("--- %s seconds ---" % (time.time() - start_time))
    print "------------------ bgi time taken: " + time_taken + " --------------------"
    with open('./scripts/bgi_json/data_dump/log_time_bgi.txt', 'w+') as res_file:
        time_taken = "time taken: " + ("--- %s seconds ---" %
                                       (time.time() - start_time))
        res_file.write(time_taken)

    ### phenotype ###
    second_start_time = time.time()'''
    get_phenotype_data()
    '''second_time_taken = "time taken: " + ("--- %s seconds ---" %
                                   (time.time() - second_start_time))
    print "------------------ phenotype time taken: " + second_time_taken + " --------------------"
    with open('./scripts/bgi_json/data_dump/log_time_pheno.txt', 'w+') as res_file_2:
        second_time_taken = "time taken: " + ("--- %s seconds ---" %
                                              (time.time() - second_start_time))
        res_file_2.write(second_time_taken)

    ### expression ###
    third_start_time = time.time()
    get_expression_data()
    third_time_taken = "time taken: " + ("--- %s seconds ---" %
                                          (time.time() - third_start_time))
    print "------------------ expression time taken: " + third_time_taken + " --------------------"
    with open('./scripts/bgi_json/data_dump/log_time_expression.txt',
              'w+') as res_file_2:
        third_start_time = "time taken: " + ("--- %s seconds ---" %
                                             (time.time() - third_start_time))
        res_file_2.write(third_start_time)'''
