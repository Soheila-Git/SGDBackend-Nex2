import urllib
import logging
import os
from datetime import datetime
import sys
reload(sys)  # Reload does the trick!
from src.models import Locusdbentity, Referencedbentity, Source, Taxonomy, \
                       Phenotype, Geninteractionannotation, Apo
from scripts.loading.database_session import get_session

__author__ = 'sweng66'

logging.basicConfig()
log = logging.getLogger()
log.setLevel(logging.INFO)

CREATED_BY = os.environ['DEFAULT_USER']

TAXON_ID = "TAX:4932"

def load_data(infile, logfile):

    nex_session = get_session()

    taxonomy =  nex_session.query(Taxonomy).filter_by(taxid=TAXON_ID).one_or_none()
    taxonomy_id = taxonomy.taxonomy_id
    
    src = nex_session.query(Source).filter_by(display_name='BioGRID').one_or_none()
    source_id = src.source_id
    
    sgdid_to_dbentity_id = dict([(x.sgdid, x.id) for x in nex_session.query(Locusdbentity).all()])
    pmid_to_reference_id = dict([(x.pmid, x.id) for x in nex_session.query(Referencedbentity).all()])
    phenotype_to_phenotype_id = dict([(x.display_name, x.id) for x in nex_session.query(Phenotype).all()])
    key_to_interaction = dict([((x.dbentity1_id, x.dbentity2_id, x.bait_hit, x.biogrid_experimental_system, x.reference_id), x) for x in nex_session.query(Geninteractionannotation).all()])
    apo_to_id = dict([(x.display_name, x.id) for x in nex_session.query(Apo).all()])

    apo_id = apo_to_id['unspecified']
    
    genetic_type_to_phenotype = get_genetic_type_to_pheno_mapping()
    old_to_new_phenotype = get_old_to_new_pheno_mapping()

    fw = open(logfile, "w")    
    f = open(infile)

    genetic_type = {}
    header = True
    for line in f:
        row = line.strip().split('\t')
        if len(row) < 12 or genetic_type_to_phenotype.get(row[4]) == None:
            continue
        if len(row) == 12:
            if header:
                header = False
            else:
                oldPheno = genetic_type_to_phenotype.get(row[4])
                newPheno = old_to_new_phenotype.get(oldPheno)

                reference_id = pmid_to_reference_id.get(int(row[6]))
                if reference_id is None:
                    print "The PMID: ", row[6], " is not in the REFERENCEDBENTITY table."
                    continue

                annotation_type = 'high-throughput' if 'HTP' in row[9] else 'manually curated'

                biogrid_experimental_system = row[4]

                dbentity1_id = sgdid_to_dbentity_id.get(row[2])
                dbentity2_id = sgdid_to_dbentity_id.get(row[3])

                if dbentity1_id is None:
                    print "The SGDID: ", row[2], " is not in the LOCUSDBENTITY table."
                    continue
                if dbentity2_id is None:
                    print "The SGDID: ", row[3], " is not in the LOCUSDBENTITY table."
                    continue

                id1 = int(row[2].lstrip("S").lstrip("0"))
                id2 = int(row[3].lstrip("S").lstrip("0"))
                bait_hit = 'bait_hit'
                if id1 >= id2:
                    (dbentity1_id, dbentity2_id) = (dbentity2_id, dbentity1_id)
                    bait_hit = 'Hit-Bait'

                description = None
                if row[10] != '-':
                    description = row[10]

                phenotype_id = None
                if newPheno is not None:
                    phenotype_id = phenotype_to_phenotype_id.get(newPheno)
                    if phenotype_id is None:
                        print "The phenotype: ", newPheno, " is not in the Phenotype table."
                    else:
                        phenotype_id = phenotype_id
                        mutant_id = apo_id
                else:
                    print "The phenotype: ", oldPheno, " can not be mapped to a new one."

                key = (dbentity1_id, dbentity2_id, bait_hit, biogrid_experimental_system, reference_id)
                
                if key in key_to_interaction:
                    x = key_to_interaction[key]
                    update_interaction(nex_session, fw, key, x, annotation_type, 
                                       phenotype_id, mutant_id, description)
                else:    
                    insert_interaction(nex_session, fw, source_id, taxonomy_id, reference_id,
                                       biogrid_experimental_system, annotation_type, dbentity1_id,
                                       dbentity2_id, bait_hit, phenotype_id, mutant_id, description)
        
                nex_session.rollback()
                # nex_session.commit() 
    f.close()

    nex_session.rollback()
    # nex_session.commit()

    fw.close()
    f.close()


def update_interaction(nex_session, fw, key, x, annotation_type, phenotype_id, mutant_id, description):
    
    update = 0
    if annotation_type is not None and annotation_type != x.annotation_type:
        x.annotation_type = annotation_type
        update = 1
    if phenotype_id is not None and phenotype_id != x.phenotype_id:
        x.phenotype_id = phenotype_id
        update = 1
    if mutant_id is not None and mutant_id != x.mutant_id:
        x.mutant_id = mutant_id
        update = 1
    if description is not None and description != x.description:
        x.description = description
        update = 1

    if update == 1:
        nex_session.add(x)
        fw.write("Update Geninteractionannotation row for key = " + str(key) + "\n") 
        

def insert_interaction(nex_session, fw, source_id, taxonomy_id, reference_id, biogrid_experimental_system, annotation_type, dbentity1_id, dbentity2_id, bait_hit, phenotype_id, mutant_id, description):

    x = Geninteractionannotation(source_id = source_id,
                                 taxonomy_id = taxonomy_id,
                                 reference_id = reference_id,
                                 biogrid_experimental_system = biogrid_experimental_system,
                                 annotation_type = annotation_type,
                                 dbentity1_id = dbentity1_id,
                                 dbentity2_id = dbentity2_id,
                                 bait_hit = bait_hit,
                                 phenotype_id = phenotype_id,
                                 mutant_id = mutant_id,
                                 description = description,
                                 created_by = CREATED_BY)
    nex_session.add(x)

    fw.write("Insert Geninteractionannotation row for dbentity1_id = " + str(dbentity1_id) + ", dbentity2_id = " + str(dbentity2_id) + ", reference_id = " + str(reference_id) + ", bait_hit = " + bait_hit + ", biogrid_experimental_system = " + biogrid_experimental_system + "\n")
           
    
def get_genetic_type_to_pheno_mapping():

    return {'Dosage Lethality'               : 'inviable',
            'Dosage Rescue'                  : 'wildtype',
            'Dosage Growth Defect'           : 'slow growth',
            'Epistatic MiniArray Profile'    : 'Not available',
            'Negative Genetic'               : 'Not available',
            'Positive Genetic'               : 'Not available',
            'Phenotypic Enhancement'         : 'Not available',
            'Phenotypic Suppression'         : 'Not available',
            'Synthetic Growth Defect'        : 'slow growth',
            'Synthetic Haploinsufficiency'   : 'Not available',
            'Synthetic Lethality'            : 'inviable',
            'Synthetic Rescue'               : 'wildtype'}


def get_old_to_new_pheno_mapping():

    return { 'inviable'    : 'inviable',
             'slow growth' : 'vegetative growth: decreased',
             'wildtype'    : 'vegetative growth: normal' }


if __name__ == '__main__':

    # https://www.thebiogrid.org/downloads/datasets/SGD.tab.txt

    url_path = 'https://www.thebiogrid.org/downloads/datasets/'
    infile = 'SGD.tab.txt'

    urllib.urlretrieve(url_path + infile, infile)
            
    logfile = "scripts/loading/biogrid/logs/geninteractionannotation.log"
    
    load_data(infile, logfile)
                                                  
