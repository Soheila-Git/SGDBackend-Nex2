from src.models import DBSession, Base, Colleague, ColleagueLocus, Locusdbentity, LocusAlias, Dnasequenceannotation, So, Locussummary, Phenotypeannotation, PhenotypeannotationCond, Phenotype, Goannotation, Go, Goslimannotation, Goslim, Apo, Straindbentity, Strainsummary, Reservedname, GoAlias, Goannotation, Referencedbentity, Referencedocument, Referenceauthor, ReferenceAlias, Chebi, Proteindomain, Contig, Dataset, Keyword

from sqlalchemy import create_engine, and_
import os
import redis

engine = create_engine(os.environ['NEX2_URI'], pool_recycle=3600)
DBSession.configure(bind=engine)
Base.metadata.bind = engine

disambiguation = redis.Redis()

ignoring = []

def table_set(key, value, prefix):
    key = str("/" + prefix + "/" + str(key)).upper()

    redis_value = disambiguation.get(key)
    
    if redis_value and redis_value != str(value):
        ignoring.append(key)
    else:
        disambiguation.set(key, value)

def load_locus():
    print("Loading genes into Redis...")

    genes = DBSession.query(Locusdbentity).all()

    aliases = DBSession.query(LocusAlias.locus_id, LocusAlias.display_name).filter(LocusAlias.alias_type.in_(['Uniform', 'Non-uniform'])).all()
    ids_to_aliases = {}
    for alias in aliases:
        if alias.locus_id in ids_to_aliases:
            ids_to_aliases[alias.locus_id].append(alias.display_name)
        else:
            ids_to_aliases[alias.locus_id] = [alias.display_name]

    # table_set will ignore in case of collisions
    # so by indexing each field separately guarantees priority
    for gene in genes:
        table_set(str(gene.dbentity_id), gene.dbentity_id, "locus")

    for gene in genes:
        table_set(str(gene.display_name.upper()), gene.dbentity_id, "locus")

    for gene in genes:
        table_set(str(gene.systematic_name.upper()), gene.dbentity_id, "locus")

    for gene in genes:
        table_set(str(gene.sgdid.upper()), gene.dbentity_id, "locus")

    for gene in genes:
        for alias in ids_to_aliases.get(gene.dbentity_id, []):
            table_set(str(alias).upper(), gene.dbentity_id, "locus")

def load_reserved_names():
    print("Loading reserve names into Redis...")

    reserved_names = DBSession.query(Reservedname).all()

    for name in reserved_names:
        table_set(str(name.reservedname_id), name.reservedname_id, "reservedname")
        table_set(str(name.format_name).upper(), name.reservedname_id, "reservedname")

def load_references():
    print("Loading references into Redis...")

    references = DBSession.query(Referencedbentity).all()

    for ref in references:
        table_set(ref.sgdid, ref.dbentity_id, "reference")
        table_set(ref.pmid, ref.dbentity_id, "reference")

def load_author():
    print("Loading authors into Redis...")

    authors = DBSession.query(Referenceauthor).all()

    ignoring = []
    
    for author in authors:
        format_name = author.obj_url.split("/")[2]

        try:
            table_set(format_name.upper(), author.referenceauthor_id, "author")
            table_set(str(author.referenceauthor_id), author.referenceauthor_id, "author")
        except UnicodeEncodeError:
            ignoring.append(format_name)

    print "Ignoring for special characters: " + ",".join(ignoring)

def load_chemical():
    print("Loading chemicals into Redis...")

    chemicals = DBSession.query(Chebi).all()

    ignoring = []
    
    for chemical in chemicals:
        try:
            table_set(str(chemical.display_name.replace(" ", "_")).upper(), chemical.chebi_id, "chebi")
            table_set(chemical.format_name.upper(), chemical.chebi_id, "chebi")
            table_set(chemical.chebi_id, chemical.chebi_id, "chebi")
        except UnicodeEncodeError:
            ignoring.append(chemical.display_name)

    print "Ignoring for special characters: " + ",".join(ignoring)

def load_phenotype():
    print("Loading phenotypes into Redis...")

    phenotypes = DBSession.query(Phenotype).all()

    for phenotype in phenotypes:
        table_set(phenotype.format_name.upper(), phenotype.phenotype_id, "phenotype")
        table_set(phenotype.phenotype_id, phenotype.phenotype_id, "phenotype")

def load_observables():
    print("Loading observables into Redis...")

    apos = DBSession.query(Apo).all()

    for apo in apos:
        table_set(apo.apo_id, apo.apo_id, "apo")
        table_set(apo.format_name.upper(), apo.apo_id, "apo")
        table_set(apo.display_name.replace(" ", "_").upper(), apo.apo_id, "apo")

def load_go():
    print("Loading go into Redis...")

    gos = DBSession.query(Go).all()

    for go in gos:
        table_set(go.format_name.upper(), go.go_id, "go")
        table_set(go.go_id, go.go_id, "go")
        table_set(go.display_name.replace(" ", "_").upper(), go.go_id, "go")

def load_protein_domain():
    print("Loading protein domains into Redis...")

    pds = DBSession.query(Proteindomain).all()

    for pd in pds:
        table_set(pd.proteindomain_id, pd.proteindomain_id, "proteindomain")
        table_set(pd.format_name.upper(), pd.proteindomain_id, "proteindomain")

def load_contigs():
    print("Loading contigs into Redis...")

    contigs = DBSession.query(Contig).all()

    for contig in contigs:
        table_set(contig.format_name.upper(), contig.contig_id, "contig")
        table_set(contig.contig_id, contig.contig_id, "contig")

def load_dataset():
    print("Loading datasets into Redis...")

    datasets = DBSession.query(Dataset).all()

    for dataset in datasets:
        table_set(dataset.format_name.upper(), dataset.dataset_id, "dataset")
        table_set(dataset.dataset_id, dataset.dataset_id, "dataset")

def load_keyword():
    print("Loading Keywords into Redis...")

    keywords = DBSession.query(Keyword).all()

    for keyword in keywords:
        table_set(keyword.keyword_id, keyword.keyword_id, "keyword")
        table_set(keyword.format_name.upper(), keyword.keyword_id, "keyword")

def load_strains():
    print("Loading strains into Redis...")

    strains = DBSession.query(Straindbentity).all()

    for strain in strains:
        table_set(strain.dbentity_id, strain.dbentity_id, "strain")
        table_set(strain.sgdid, strain.dbentity_id, "strain")
        table_set(strain.display_name.replace(" ", "_"), strain.dbentity_id, "strain")
        

load_locus()
load_references()
load_reserved_names()
load_author()
load_chemical()
load_phenotype()
load_observables()
load_go()
load_protein_domain()
load_contigs()
load_dataset()
load_keyword()
load_strains()