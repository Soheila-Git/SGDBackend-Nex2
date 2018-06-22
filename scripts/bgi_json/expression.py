from src.models import DBSession, Base, Goextension, Dbentity, Locusdbentity, LocusUrl, LocusAlias, Dnasequenceannotation, So, Locussummary, Goannotation, Go, Eco, Goslimannotation,Goslim, GoAlias, Goannotation, Referencedbentity
from sqlalchemy import create_engine, and_, inspect, func
from sqlalchemy.sql import label
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
        Goannotation.annotation_type != 'computational', Goannotation.eco_id != 236340).all()
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
                    "MMO:0000658",
                "dateAssigned":
                    item.date_created.strftime("%Y-%m-%dT%H:%m:%S-00:00"),
                "wildtypeExpressionTermIdentifiers": {"anatomicalStructureTermId": "N/A"}
            }
            if item.reference.pmid:
                obj["evidence"]["pubMedId"] = "PMID:" + str(item.reference.pmid)
            else:
                obj["evidence"]["modPublicationId"] = "SGD:" + str(item.reference.sgdid)
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
            "dateProduced": datetime.utcnow().strftime("%Y-%m-%dT%H:%m:%S-00:00")
        }

        return main_obj
    return None


def get_expression_data():
    with concurrent.futures.ProcessPoolExecutor(max_workers=8) as executor:
        go_data = get_goannotation_data()
        if go_data:
            with open('expression_2.json', 'w') as file:
                file.write(json.dumps(go_data))


if __name__ == '__main__':
    get_expression_data()
