-- Generated by Ora2Pg, the Oracle database Schema converter, version 17.4
-- Copyright 2000-2016 Gilles DAROLD. All rights reserved.
-- DATASOURCE: dbi:Oracle:host=sgd-nex2-db.stanford.edu;sid=SGD

SET client_encoding TO 'UTF8';

\set ON_ERROR_STOP ON

-- Set priviledge on TABLE ALLELE
ALTER TABLE allele OWNER TO nex;
GRANT ALL ON  allele TO nex;
REVOKE ALL ON allele FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON allele TO curator;
GRANT SELECT ON allele TO PUBLIC;

-- Set priviledge on TABLE APO
ALTER TABLE apo OWNER TO nex;
GRANT ALL ON  apo TO nex;
REVOKE ALL ON apo FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON apo TO curator;
GRANT SELECT ON apo TO PUBLIC;

-- Set priviledge on TABLE APO_ALIAS
ALTER TABLE apo_alias OWNER TO nex;
GRANT ALL ON  apo_alias TO nex;
REVOKE ALL ON apo_alias FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON apo_alias TO curator;
GRANT SELECT ON apo_alias TO PUBLIC;

-- Set priviledge on TABLE APO_RELATION
ALTER TABLE apo_relation OWNER TO nex;
GRANT ALL ON  apo_relation TO nex;
REVOKE ALL ON apo_relation FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON apo_relation TO curator;
GRANT SELECT ON apo_relation TO PUBLIC;

-- Set priviledge on TABLE APO_URL
ALTER TABLE apo_url OWNER TO nex;
GRANT ALL ON  apo_url TO nex;
REVOKE ALL ON apo_url FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON apo_url TO curator;
GRANT SELECT ON apo_url TO PUBLIC;

-- Set priviledge on TABLE ARCH_CONTIG
ALTER TABLE arch_contig OWNER TO nex;
GRANT ALL ON  arch_contig TO nex;
REVOKE ALL ON arch_contig FROM PUBLIC;
GRANT SELECT,INSERT ON arch_contig TO curator;
GRANT SELECT ON arch_contig TO PUBLIC;

-- Set priviledge on TABLE ARCH_CONTIGCHANGE
ALTER TABLE arch_contigchange OWNER TO nex;
GRANT ALL ON  arch_contigchange TO nex;
REVOKE ALL ON arch_contigchange FROM PUBLIC;
GRANT SELECT,INSERT ON arch_contigchange TO curator;
GRANT SELECT ON arch_contigchange TO PUBLIC;

-- Set priviledge on TABLE ARCH_DNASEQUENCEANNOTATION
ALTER TABLE arch_dnasequenceannotation OWNER TO nex;
GRANT ALL ON  arch_dnasequenceannotation TO nex;
REVOKE ALL ON arch_dnasequenceannotation FROM PUBLIC;
GRANT SELECT,INSERT ON arch_dnasequenceannotation TO curator;
GRANT SELECT ON arch_dnasequenceannotation TO PUBLIC;

-- Set priviledge on TABLE ARCH_DNASUBSEQUENCE
ALTER TABLE arch_dnasubsequence OWNER TO nex;
GRANT ALL ON  arch_dnasubsequence TO nex;
REVOKE ALL ON arch_dnasubsequence FROM PUBLIC;
GRANT SELECT,INSERT ON arch_dnasubsequence TO curator;
GRANT SELECT ON arch_dnasubsequence TO PUBLIC;

-- Set priviledge on TABLE ARCH_LITERATUREANNOTATION
ALTER TABLE arch_literatureannotation OWNER TO nex;
GRANT ALL ON  arch_literatureannotation TO nex;
REVOKE ALL ON arch_literatureannotation FROM PUBLIC;
GRANT SELECT,INSERT ON arch_literatureannotation TO curator;
GRANT SELECT ON arch_literatureannotation TO PUBLIC;

-- Set priviledge on TABLE ARCH_LOCUSCHANGE
ALTER TABLE arch_locuschange OWNER TO nex;
GRANT ALL ON  arch_locuschange TO nex;
REVOKE ALL ON arch_locuschange FROM PUBLIC;
GRANT SELECT,INSERT ON arch_locuschange TO curator;
GRANT SELECT ON arch_locuschange TO PUBLIC;

-- Set priviledge on TABLE ARCH_PROTEINSEQUENCEANNOTATION
ALTER TABLE arch_proteinsequenceannotation OWNER TO nex;
GRANT ALL ON  arch_proteinsequenceannotation TO nex;
REVOKE ALL ON arch_proteinsequenceannotation FROM PUBLIC;
GRANT SELECT,INSERT ON arch_proteinsequenceannotation TO curator;
GRANT SELECT ON arch_proteinsequenceannotation TO PUBLIC;

-- Set priviledge on TABLE AUTHORRESPONSE
ALTER TABLE authorresponse OWNER TO nex;
GRANT ALL ON  authorresponse TO nex;
REVOKE ALL ON authorresponse FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON authorresponse TO curator;
GRANT SELECT ON authorresponse TO PUBLIC;

-- Set priviledge on TABLE BINDINGMOTIFANNOTATION
ALTER TABLE bindingmotifannotation OWNER TO nex;
GRANT ALL ON  bindingmotifannotation TO nex;
REVOKE ALL ON bindingmotifannotation FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON bindingmotifannotation TO curator;
GRANT SELECT ON bindingmotifannotation TO PUBLIC;

-- Set priviledge on TABLE BOOK
ALTER TABLE book OWNER TO nex;
GRANT ALL ON  book TO nex;
REVOKE ALL ON book FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON book TO curator;
GRANT SELECT ON book TO PUBLIC;

-- Set priviledge on TABLE CHEBI
ALTER TABLE chebi OWNER TO nex;
GRANT ALL ON  chebi TO nex;
REVOKE ALL ON chebi FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON chebi TO curator;
GRANT SELECT ON chebi TO PUBLIC;

-- Set priviledge on TABLE CHEBI_ALIAS
ALTER TABLE chebi_alias OWNER TO nex;
GRANT ALL ON  chebi_alias TO nex;
REVOKE ALL ON chebi_alias FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON chebi_alias TO curator;
GRANT SELECT ON chebi_alias TO PUBLIC;

-- Set priviledge on TABLE CHEBI_RELATION
ALTER TABLE chebi_relation OWNER TO nex;
GRANT ALL ON  chebi_relation TO nex;
REVOKE ALL ON chebi_relation FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON chebi_relation TO curator;
GRANT SELECT ON chebi_relation TO PUBLIC;

-- Set priviledge on TABLE CHEBI_URL
ALTER TABLE chebi_url OWNER TO nex;
GRANT ALL ON  chebi_url TO nex;
REVOKE ALL ON chebi_url FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON chebi_url TO curator;
GRANT SELECT ON chebi_url TO PUBLIC;

-- Set priviledge on TABLE COLLEAGUE
ALTER TABLE colleague OWNER TO nex;
GRANT ALL ON  colleague TO nex;
REVOKE ALL ON colleague FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON colleague TO curator;
GRANT SELECT ON colleague TO PUBLIC;

-- Set priviledge on TABLE COLLEAGUETRIAGE
ALTER TABLE colleaguetriage OWNER TO nex;
GRANT ALL ON  colleaguetriage TO nex;
REVOKE ALL ON colleaguetriage FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON colleaguetriage TO curator;
GRANT SELECT ON colleaguetriage TO PUBLIC;

-- Set priviledge on TABLE COLLEAGUE_RELATION
ALTER TABLE colleague_relation OWNER TO nex;
GRANT ALL ON  colleague_relation TO nex;
REVOKE ALL ON colleague_relation FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON colleague_relation TO curator;
GRANT SELECT ON colleague_relation TO PUBLIC;

-- Set priviledge on TABLE COLLEAGUE_KEYWORD
ALTER TABLE colleague_keyword OWNER TO nex;
GRANT ALL ON  colleague_keyword TO nex;
REVOKE ALL ON colleague_keyword FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON colleague_keyword TO curator;
GRANT SELECT ON colleague_keyword TO PUBLIC;

-- Set priviledge on TABLE COLLEAGUE_LOCUS
ALTER TABLE colleague_locus OWNER TO nex;
GRANT ALL ON  colleague_locus TO nex;
REVOKE ALL ON colleague_locus FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON colleague_locus TO curator;
GRANT SELECT ON colleague_locus TO PUBLIC;

-- Set priviledge on TABLE COLLEAGUE_REFERENCE
ALTER TABLE colleague_reference OWNER TO nex;
GRANT ALL ON  colleague_reference TO nex;
REVOKE ALL ON colleague_reference FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON colleague_reference TO curator;
GRANT SELECT ON colleague_reference TO PUBLIC;

-- Set priviledge on TABLE COLLEAGUE_URL
ALTER TABLE colleague_url OWNER TO nex;
GRANT ALL ON  colleague_url TO nex;
REVOKE ALL ON colleague_url FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON colleague_url TO curator;
GRANT SELECT ON colleague_url TO PUBLIC;

-- Set priviledge on TABLE CONTIG
ALTER TABLE contig OWNER TO nex;
GRANT ALL ON  contig TO nex;
REVOKE ALL ON contig FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON contig TO curator;
GRANT SELECT ON contig TO PUBLIC;

-- Set priviledge on TABLE CONTIGNOTEANNOTATION
ALTER TABLE contignoteannotation OWNER TO nex;
GRANT ALL ON  contignoteannotation TO nex;
REVOKE ALL ON contignoteannotation FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON contignoteannotation TO curator;
GRANT SELECT ON contignoteannotation TO PUBLIC;

-- Set priviledge on TABLE CONTIG_URL
ALTER TABLE contig_url OWNER TO nex;
GRANT ALL ON  contig_url TO nex;
REVOKE ALL ON contig_url FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON contig_url TO curator;
GRANT SELECT ON contig_url TO PUBLIC;

-- Set priviledge on TABLE CURATION
ALTER TABLE curation OWNER TO nex;
GRANT ALL ON  curation TO nex;
REVOKE ALL ON curation FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON curation TO curator;
GRANT SELECT ON curation TO PUBLIC;

-- Set priviledge on TABLE DATASET
ALTER TABLE dataset OWNER TO nex;
GRANT ALL ON  dataset TO nex;
REVOKE ALL ON dataset FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON dataset TO curator;
GRANT SELECT ON dataset TO PUBLIC;

-- Set priviledge on TABLE DATASETLAB
ALTER TABLE datasetlab OWNER TO nex;
GRANT ALL ON  datasetlab TO nex;
REVOKE ALL ON datasetlab FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON datasetlab TO curator;
GRANT SELECT ON datasetlab TO PUBLIC;

-- Set priviledge on TABLE DATASETSAMPLE
ALTER TABLE datasetsample OWNER TO nex;
GRANT ALL ON  datasetsample TO nex;
REVOKE ALL ON datasetsample FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON datasetsample TO curator;
GRANT SELECT ON datasetsample TO PUBLIC;

-- Set priviledge on TABLE DATASETTRACK
ALTER TABLE datasettrack OWNER TO nex;
GRANT ALL ON  datasettrack TO nex;
REVOKE ALL ON datasettrack FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON datasettrack TO curator;
GRANT SELECT ON datasettrack TO PUBLIC;

-- Set priviledge on TABLE DATASET_FILE
ALTER TABLE dataset_file OWNER TO nex;
GRANT ALL ON  dataset_file TO nex;
REVOKE ALL ON dataset_file FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON dataset_file TO curator;
GRANT SELECT ON dataset_file TO PUBLIC;

-- Set priviledge on TABLE DATASET_KEYWORD
ALTER TABLE dataset_keyword OWNER TO nex;
GRANT ALL ON  dataset_keyword TO nex;
REVOKE ALL ON dataset_keyword FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON dataset_keyword TO curator;
GRANT SELECT ON dataset_keyword TO PUBLIC;

-- Set priviledge on TABLE DATASET_REFERENCE
ALTER TABLE dataset_reference OWNER TO nex;
GRANT ALL ON  dataset_reference TO nex;
REVOKE ALL ON dataset_reference FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON dataset_reference TO curator;
GRANT SELECT ON dataset_reference TO PUBLIC;

-- Set priviledge on TABLE DATASET_URL
ALTER TABLE dataset_url OWNER TO nex;
GRANT ALL ON  dataset_url TO nex;
REVOKE ALL ON dataset_url FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON dataset_url TO curator;
GRANT SELECT ON dataset_url TO PUBLIC;

-- Set priviledge on TABLE DBENTITY
ALTER TABLE dbentity OWNER TO nex;
GRANT ALL ON  dbentity TO nex;
REVOKE ALL ON dbentity FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON dbentity TO curator;
GRANT SELECT ON dbentity TO PUBLIC;

-- Set priviledge on TABLE DBUSER
ALTER TABLE dbuser OWNER TO nex;
GRANT ALL ON  dbuser TO nex;
REVOKE ALL ON dbuser FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE ON dbuser TO curator;
GRANT SELECT ON dbuser TO PUBLIC;

-- Set priviledge on TABLE DELETELOG
ALTER TABLE deletelog OWNER TO nex;
GRANT ALL ON  deletelog TO nex;
REVOKE ALL ON deletelog FROM PUBLIC;
GRANT SELECT,INSERT ON deletelog TO curator;
GRANT SELECT ON deletelog TO PUBLIC;

-- Set priviledge on TABLE DNASEQUENCEANNOTATION
ALTER TABLE dnasequenceannotation OWNER TO nex;
GRANT ALL ON  dnasequenceannotation TO nex;
REVOKE ALL ON dnasequenceannotation FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON dnasequenceannotation TO curator;
GRANT SELECT ON dnasequenceannotation TO PUBLIC;

-- Set priviledge on TABLE DNASUBSEQUENCE
ALTER TABLE dnasubsequence OWNER TO nex;
GRANT ALL ON  dnasubsequence TO nex;
REVOKE ALL ON dnasubsequence FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON dnasubsequence TO curator;
GRANT SELECT ON dnasubsequence TO PUBLIC;

-- Set priviledge on TABLE DISEASE
ALTER TABLE disease OWNER TO nex;
GRANT ALL ON  disease TO nex;
REVOKE ALL ON disease FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON disease TO curator;
GRANT SELECT ON disease TO PUBLIC;

-- Set priviledge on TABLE DISEASEANNOTATION
ALTER TAble diseaseannotation OWNER TO nex;
GRANT ALL ON  diseaseannotation TO nex;
REVOKE ALL ON diseaseannotation FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON diseaseannotation TO curator;
GRANT SELECT ON diseaseannotation TO PUBLIC;

-- Set priviledge on TABLE DISEASESUBSET
ALTER TABLE diseasesubset OWNER TO nex;
GRANT ALL ON  diseasesubset TO nex;
REVOKE ALL ON diseasesubset FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON diseasesubset TO curator;
GRANT SELECT ON diseasesubset TO PUBLIC;

-- Set priviledge on TABLE DISEASESUBSETANNOTATION
ALTER TABLE diseasesubsetannotation OWNER TO nex;
GRANT ALL ON  diseasesubsetannotation TO nex;
REVOKE ALL ON diseasesubsetannotation FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON diseasesubsetannotation TO curator;
GRANT SELECT ON diseasesubsetannotation TO PUBLIC;

-- Set priviledge on TABLE DISEASESUPPORTINGEVIDENCE
ALTER TABLE diseasesupportingevidence OWNER TO nex;
GRANT ALL ON  diseasesupportingevidence TO nex;
REVOKE ALL ON diseasesupportingevidence FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON diseasesupportingevidence TO curator;
GRANT SELECT ON diseasesupportingevidence TO PUBLIC;

-- Set priviledge on TABLE DISEASE_ALIAS
ALTER TABLE disease_alias OWNER TO nex;
GRANT ALL ON  disease_alias TO nex;
REVOKE ALL ON disease_alias FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON disease_alias TO curator;
GRANT SELECT ON disease_alias TO PUBLIC;

-- Set priviledge on TABLE DISEASE_RELATION
ALTER TABLE disease_relation OWNER TO nex;
GRANT ALL ON  disease_relation TO nex;
REVOKE ALL ON disease_relation FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON disease_relation TO curator;
GRANT SELECT ON disease_relation TO PUBLIC;

-- Set priviledge on TABLE DISEASE_URL
ALTER TABLE disease_url OWNER TO nex;
GRANT ALL ON  disease_url TO nex;
REVOKE ALL ON disease_url FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON disease_url TO curator;
GRANT SELECT ON disease_url TO PUBLIC;

-- Set priviledge on TABLE EC
ALTER TABLE ec OWNER TO nex;
GRANT ALL ON  ec TO nex;
REVOKE ALL ON ec FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON ec TO curator;
GRANT SELECT ON ec TO PUBLIC;

-- Set priviledge on TABLE ECO
ALTER TABLE eco OWNER TO nex;
GRANT ALL ON  eco TO nex;
REVOKE ALL ON eco FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON eco TO curator;
GRANT SELECT ON eco TO PUBLIC;

-- Set priviledge on TABLE ECO_ALIAS
ALTER TABLE eco_alias OWNER TO nex;
GRANT ALL ON  eco_alias TO nex;
REVOKE ALL ON eco_alias FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON eco_alias TO curator;
GRANT SELECT ON eco_alias TO PUBLIC;

-- Set priviledge on TABLE ECO_RELATION
ALTER TABLE eco_relation OWNER TO nex;
GRANT ALL ON  eco_relation TO nex;
REVOKE ALL ON eco_relation FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON eco_relation TO curator;
GRANT SELECT ON eco_relation TO PUBLIC;

-- Set priviledge on TABLE ECO_URL
ALTER TABLE eco_url OWNER TO nex;
GRANT ALL ON  eco_url TO nex;
REVOKE ALL ON eco_url FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON eco_url TO curator;
GRANT SELECT ON eco_url TO PUBLIC;

-- Set priviledge on TABLE EC_ALIAS
ALTER TABLE ec_alias OWNER TO nex;
GRANT ALL ON  ec_alias TO nex;
REVOKE ALL ON ec_alias FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON ec_alias TO curator;
GRANT SELECT ON ec_alias TO PUBLIC;

-- Set priviledge on TABLE EC_URL
ALTER TABLE ec_url OWNER TO nex;
GRANT ALL ON  ec_url TO nex;
REVOKE ALL ON ec_url FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON ec_url TO curator;
GRANT SELECT ON ec_url TO PUBLIC;

-- Set priviledge on TABLE EDAM
ALTER TABLE edam OWNER TO nex;
GRANT ALL ON  edam TO nex;
REVOKE ALL ON edam FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON edam TO curator;
GRANT SELECT ON edam TO PUBLIC;

-- Set priviledge on TABLE EDAM_ALIAS
ALTER TABLE edam_alias OWNER TO nex;
GRANT ALL ON  edam_alias TO nex;
REVOKE ALL ON edam_alias FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON edam_alias TO curator;
GRANT SELECT ON edam_alias TO PUBLIC;

-- Set priviledge on TABLE EDAM_RELATION
ALTER TABLE edam_relation OWNER TO nex;
GRANT ALL ON  edam_relation TO nex;
REVOKE ALL ON edam_relation FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON edam_relation TO curator;
GRANT SELECT ON edam_relation TO PUBLIC;

-- Set priviledge on TABLE EDAM_URL
ALTER TABLE edam_url OWNER TO nex;
GRANT ALL ON  edam_url TO nex;
REVOKE ALL ON edam_url FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON edam_url TO curator;
GRANT SELECT ON edam_url TO PUBLIC;

-- Set priviledge on TABLE ENZYMEANNOTATION
ALTER TABLE enzymeannotation OWNER TO nex;
GRANT ALL ON  enzymeannotation TO nex;
REVOKE ALL ON enzymeannotation FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON enzymeannotation TO curator;
GRANT SELECT ON enzymeannotation TO PUBLIC;

-- Set priviledge on TABLE EXPRESSIONANNOTATION
ALTER TABLE expressionannotation OWNER TO nex;
GRANT ALL ON  expressionannotation TO nex;
REVOKE ALL ON expressionannotation FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON expressionannotation TO curator;
GRANT SELECT ON expressionannotation TO PUBLIC;

-- Set priviledge on TABLE FILEDBENTITY
ALTER TABLE filedbentity OWNER TO nex;
GRANT ALL ON  filedbentity TO nex;
REVOKE ALL ON filedbentity FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON filedbentity TO curator;
GRANT SELECT ON filedbentity TO PUBLIC;

-- Set priviledge on TABLE FILEPATH
ALTER TABLE filepath OWNER TO nex;
GRANT ALL ON  filepath TO nex;
REVOKE ALL ON filepath FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON filepath TO curator;
GRANT SELECT ON filepath TO PUBLIC;

-- Set priviledge on TABLE FILE_KEYWORD
ALTER TABLE file_keyword OWNER TO nex;
GRANT ALL ON  file_keyword TO nex;
REVOKE ALL ON file_keyword FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON file_keyword TO curator;
GRANT SELECT ON file_keyword TO PUBLIC;

-- Set priviledge on TABLE GENINTERACTIONANNOTATION
ALTER TABLE geninteractionannotation OWNER TO nex;
GRANT ALL ON  geninteractionannotation TO nex;
REVOKE ALL ON geninteractionannotation FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON geninteractionannotation TO curator;
GRANT SELECT ON geninteractionannotation TO PUBLIC;

-- Set priviledge on TABLE GENOMERELEASE
ALTER TABLE genomerelease OWNER TO nex;
GRANT ALL ON  genomerelease TO nex;
REVOKE ALL ON genomerelease FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON genomerelease TO curator;
GRANT SELECT ON genomerelease TO PUBLIC;

-- Set priviledge on TABLE GO
ALTER TABLE go OWNER TO nex;
GRANT ALL ON  go TO nex;
REVOKE ALL ON go FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON go TO curator;
GRANT SELECT ON go TO PUBLIC;

-- Set priviledge on TABLE GOANNOTATION
ALTER TABLE goannotation OWNER TO nex;
GRANT ALL ON  goannotation TO nex;
REVOKE ALL ON goannotation FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON goannotation TO curator;
GRANT SELECT ON goannotation TO PUBLIC;

-- Set priviledge on TABLE GOEXTENSION
ALTER TABLE goextension OWNER TO nex;
GRANT ALL ON  goextension TO nex;
REVOKE ALL ON goextension FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON goextension TO curator;
GRANT SELECT ON goextension TO PUBLIC;

-- Set priviledge on TABLE GOSLIM
ALTER TABLE goslim OWNER TO nex;
GRANT ALL ON  goslim TO nex;
REVOKE ALL ON goslim FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON goslim TO curator;
GRANT SELECT ON goslim TO PUBLIC;

-- Set priviledge on TABLE GOSLIMANNOTATION
ALTER TABLE goslimannotation OWNER TO nex;
GRANT ALL ON  goslimannotation TO nex;
REVOKE ALL ON goslimannotation FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON goslimannotation TO curator;
GRANT SELECT ON goslimannotation TO PUBLIC;

-- Set priviledge on TABLE GOSUPPORTINGEVIDENCE
ALTER TABLE gosupportingevidence OWNER TO nex;
GRANT ALL ON  gosupportingevidence TO nex;
REVOKE ALL ON gosupportingevidence FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON gosupportingevidence TO curator;
GRANT SELECT ON gosupportingevidence TO PUBLIC;

-- Set priviledge on TABLE GO_ALIAS
ALTER TABLE go_alias OWNER TO nex;
GRANT ALL ON  go_alias TO nex;
REVOKE ALL ON go_alias FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON go_alias TO curator;
GRANT SELECT ON go_alias TO PUBLIC;

-- Set priviledge on TABLE GO_RELATION
ALTER TABLE go_relation OWNER TO nex;
GRANT ALL ON  go_relation TO nex;
REVOKE ALL ON go_relation FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON go_relation TO curator;
GRANT SELECT ON go_relation TO PUBLIC;

-- Set priviledge on TABLE GO_URL
ALTER TABLE go_url OWNER TO nex;
GRANT ALL ON  go_url TO nex;
REVOKE ALL ON go_url FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON go_url TO curator;
GRANT SELECT ON go_url TO PUBLIC;

-- Set priviledge on TABLE JOURNAL
ALTER TABLE journal OWNER TO nex;
GRANT ALL ON  journal TO nex;
REVOKE ALL ON journal FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON journal TO curator;
GRANT SELECT ON journal TO PUBLIC;

-- Set priviledge on TABLE KEYWORD
ALTER TABLE keyword OWNER TO nex;
GRANT ALL ON  keyword TO nex;
REVOKE ALL ON keyword FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON keyword TO curator;
GRANT SELECT ON keyword TO PUBLIC;

-- Set priviledge on TABLE LITERATUREANNOTATION
ALTER TABLE literatureannotation OWNER TO nex;
GRANT ALL ON  literatureannotation TO nex;
REVOKE ALL ON literatureannotation FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON literatureannotation TO curator;
GRANT SELECT ON literatureannotation TO PUBLIC;

-- Set priviledge on TABLE LOCUSDBENTITY
ALTER TABLE locusdbentity OWNER TO nex;
GRANT ALL ON  locusdbentity TO nex;
REVOKE ALL ON locusdbentity FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON locusdbentity TO curator;
GRANT SELECT ON locusdbentity TO PUBLIC;

-- Set priviledge on TABLE LOCUSNOTEANNOTATION
ALTER TABLE locusnoteannotation OWNER TO nex;
GRANT ALL ON  locusnoteannotation TO nex;
REVOKE ALL ON locusnoteannotation FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON locusnoteannotation TO curator;
GRANT SELECT ON locusnoteannotation TO PUBLIC;

-- Set priviledge on TABLE LOCUS_ALIAS
ALTER TABLE locus_alias OWNER TO nex;
GRANT ALL ON  locus_alias TO nex;
REVOKE ALL ON locus_alias FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON locus_alias TO curator;
GRANT SELECT ON locus_alias TO PUBLIC;

-- Set priviledge on TABLE LOCUS_RELATION
ALTER TABLE locus_relation OWNER TO nex;
GRANT ALL ON  locus_relation TO nex;
REVOKE ALL ON locus_relation FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON locus_relation TO curator;
GRANT SELECT ON locus_relation TO PUBLIC;

-- Set priviledge on TABLE LOCUSSUMMARY
ALTER TABLE locussummary OWNER TO nex;
GRANT ALL ON  locussummary TO nex;
REVOKE ALL ON locussummary FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON locussummary TO curator;
GRANT SELECT ON locussummary TO PUBLIC;

-- Set priviledge on TABLE LOCUSSUMMARY_REFERENCE
ALTER TABLE locussummary_reference OWNER TO nex;
GRANT ALL ON  locussummary_reference TO nex;
REVOKE ALL ON locussummary_reference FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON locussummary_reference TO curator;
GRANT SELECT ON locussummary_reference TO PUBLIC;

-- Set priviledge on TABLE LOCUS_URL
ALTER TABLE locus_url OWNER TO nex;
GRANT ALL ON  locus_url TO nex;
REVOKE ALL ON locus_url FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON locus_url TO curator;
GRANT SELECT ON locus_url TO PUBLIC;

-- Set priviledge on TABLE OBI
ALTER TABLE obi OWNER TO nex;
GRANT ALL ON  obi TO nex;
REVOKE ALL ON obi FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON obi TO curator;
GRANT SELECT ON obi TO PUBLIC;

-- Set priviledge on TABLE OBI_RELATION
ALTER TABLE obi_relation OWNER TO nex;
GRANT ALL ON  obi_relation TO nex;
REVOKE ALL ON obi_relation FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON obi_relation TO curator;
GRANT SELECT ON obi_relation TO PUBLIC;

-- Set priviledge on TABLE OBI_URL
ALTER TABLE obi_url OWNER TO nex;
GRANT ALL ON  obi_url TO nex;
REVOKE ALL ON obi_url FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON obi_url TO curator;
GRANT SELECT ON obi_url TO PUBLIC;

-- Set priviledge on TABLE PATHWAYANNOTATION
ALTER TABLE pathwayannotation OWNER TO nex;
GRANT ALL ON  pathwayannotation TO nex;
REVOKE ALL ON pathwayannotation FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON pathwayannotation TO curator;
GRANT SELECT ON pathwayannotation TO PUBLIC;

-- Set priviledge on TABLE PATHWAYDBENTITY
ALTER TABLE pathwaydbentity OWNER TO nex;
GRANT ALL ON  pathwaydbentity TO nex;
REVOKE ALL ON pathwaydbentity FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON pathwaydbentity TO curator;
GRANT SELECT ON pathwaydbentity TO PUBLIC;

-- Set priviledge on TABLE PATHWAY_ALIAS
ALTER TABLE pathway_alias OWNER TO nex;
GRANT ALL ON  pathway_alias TO nex;
REVOKE ALL ON pathway_alias FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON pathway_alias TO curator;
GRANT SELECT ON pathway_alias TO PUBLIC;

-- Set priviledge on TABLE PATHWAYSUMMARY
ALTER TABLE pathwaysummary OWNER TO nex;
GRANT ALL ON  pathwaysummary TO nex;
REVOKE ALL ON pathwaysummary FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON pathwaysummary TO curator;
GRANT SELECT ON pathwaysummary TO PUBLIC;

-- Set priviledge on TABLE PATHWAYSUMMARY_REFERENCE
ALTER TABLE pathwaysummary_reference OWNER TO nex;
GRANT ALL ON  pathwaysummary_reference TO nex;
REVOKE ALL ON pathwaysummary_reference FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON pathwaysummary_reference TO curator;
GRANT SELECT ON pathwaysummary_reference TO PUBLIC;

-- Set priviledge on TABLE PATHWAY_URL
ALTER TABLE pathway_url OWNER TO nex;
GRANT ALL ON  pathway_url TO nex;
REVOKE ALL ON pathway_url FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON pathway_url TO curator;
GRANT SELECT ON pathway_url TO PUBLIC;

-- Set priviledge on TABLE PHENOTYPE
ALTER TABLE phenotype OWNER TO nex;
GRANT ALL ON  phenotype TO nex;
REVOKE ALL ON phenotype FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON phenotype TO curator;
GRANT SELECT ON phenotype TO PUBLIC;

-- Set priviledge on TABLE PHENOTYPEANNOTATION
ALTER TABLE phenotypeannotation OWNER TO nex;
GRANT ALL ON  phenotypeannotation TO nex;
REVOKE ALL ON phenotypeannotation FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON phenotypeannotation TO curator;
GRANT SELECT ON phenotypeannotation TO PUBLIC;

-- Set priviledge on TABLE PHENOTYPEANNOTATION_COND
ALTER TABLE phenotypeannotation_cond OWNER TO nex;
GRANT ALL ON  phenotypeannotation_cond TO nex;
REVOKE ALL ON phenotypeannotation_cond FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON phenotypeannotation_cond TO curator;
GRANT SELECT ON phenotypeannotation_cond TO PUBLIC;

-- Set priviledge on TABLE PHYSINTERACTIONANNOTATION
ALTER TABLE physinteractionannotation OWNER TO nex;
GRANT ALL ON  physinteractionannotation TO nex;
REVOKE ALL ON physinteractionannotation FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON physinteractionannotation TO curator;
GRANT SELECT ON physinteractionannotation TO PUBLIC;

-- Set priviledge on TABLE POSTTRANSLATIONANNOTATION
ALTER TABLE posttranslationannotation OWNER TO nex;
GRANT ALL ON  posttranslationannotation TO nex;
REVOKE ALL ON posttranslationannotation FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON posttranslationannotation TO curator;
GRANT SELECT ON posttranslationannotation TO PUBLIC;

-- Set priviledge on TABLE PROTEINDOMAIN
ALTER TABLE proteindomain OWNER TO nex;
GRANT ALL ON  proteindomain TO nex;
REVOKE ALL ON proteindomain FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON proteindomain TO curator;
GRANT SELECT ON proteindomain TO PUBLIC;

-- Set priviledge on TABLE PROTEINDOMAINANNOTATION
ALTER TABLE proteindomainannotation OWNER TO nex;
GRANT ALL ON  proteindomainannotation TO nex;
REVOKE ALL ON proteindomainannotation FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON proteindomainannotation TO curator;
GRANT SELECT ON proteindomainannotation TO PUBLIC;

-- Set priviledge on TABLE PROTEINDOMAIN_URL
ALTER TABLE proteindomain_url OWNER TO nex;
GRANT ALL ON  proteindomain_url TO nex;
REVOKE ALL ON proteindomain_url FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON proteindomain_url TO curator;
GRANT SELECT ON proteindomain_url TO PUBLIC;

-- Set priviledge on TABLE PROTEINEXPTANNOTATION
ALTER TABLE proteinexptannotation OWNER TO nex;
GRANT ALL ON  proteinexptannotation TO nex;
REVOKE ALL ON proteinexptannotation FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON proteinexptannotation TO curator;
GRANT SELECT ON proteinexptannotation TO PUBLIC;

-- Set priviledge on TABLE PROTEINEXPTANNOTATION_COND
ALTER TABLE proteinexptannotation_cond OWNER TO nex;
GRANT ALL ON  proteinexptannotation_cond TO nex;
REVOKE ALL ON proteinexptannotation_cond FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON proteinexptannotation_cond TO curator;
GRANT SELECT ON proteinexptannotation_cond TO PUBLIC;

-- Set priviledge on TABLE PROTEINSEQUENCEANNOTATION
ALTER TABLE proteinsequenceannotation OWNER TO nex;
GRANT ALL ON  proteinsequenceannotation TO nex;
REVOKE ALL ON proteinsequenceannotation FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON proteinsequenceannotation TO curator;
GRANT SELECT ON proteinsequenceannotation TO PUBLIC;

-- Set priviledge on TABLE PROTEINSEQUENCE_DETAIL
ALTER TABLE proteinsequence_detail OWNER TO nex;
GRANT ALL ON  proteinsequence_detail TO nex;
REVOKE ALL ON proteinsequence_detail FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON proteinsequence_detail TO curator;
GRANT SELECT ON proteinsequence_detail TO PUBLIC;

-- Set priviledge on TABLE PSIMOD
ALTER TABLE psimod OWNER TO nex;
GRANT ALL ON  psimod TO nex;
REVOKE ALL ON psimod FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON psimod TO curator;
GRANT SELECT ON psimod TO PUBLIC;

-- Set priviledge on TABLE PSIMOD_RELATION
ALTER TABLE psimod_relation OWNER TO nex;
GRANT ALL ON  psimod_relation TO nex;
REVOKE ALL ON psimod_relation FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON psimod_relation TO curator;
GRANT SELECT ON psimod_relation TO PUBLIC;

-- Set priviledge on TABLE PSIMOD_URL
ALTER TABLE psimod_url OWNER TO nex;
GRANT ALL ON  psimod_url TO nex;
REVOKE ALL ON psimod_url FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON psimod_url TO curator;
GRANT SELECT ON psimod_url TO PUBLIC;

-- Set priviledge on TABLE REFERENCEDBENTITY
ALTER TABLE referencedbentity OWNER TO nex;
GRANT ALL ON  referencedbentity TO nex;
REVOKE ALL ON referencedbentity FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON referencedbentity TO curator;
GRANT SELECT ON referencedbentity TO PUBLIC;

-- Set priviledge on TABLE REFERENCETRIAGE
ALTER TABLE referencetriage OWNER TO nex;
GRANT ALL ON  referencetriage TO nex;
REVOKE ALL ON referencetriage FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON referencetriage TO curator;
GRANT SELECT ON referencetriage TO PUBLIC;

-- Set priviledge on TABLE REFERENCE_ALIAS
ALTER TABLE reference_alias OWNER TO nex;
GRANT ALL ON  reference_alias TO nex;
REVOKE ALL ON reference_alias FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON reference_alias TO curator;
GRANT SELECT ON reference_alias TO PUBLIC;

-- Set priviledge on TABLE REFERENCEAUTHOR
ALTER TABLE referenceauthor OWNER TO nex;
GRANT ALL ON  referenceauthor TO nex;
REVOKE ALL ON referenceauthor FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON referenceauthor TO curator;
GRANT SELECT ON referenceauthor TO PUBLIC;

-- Set priviledge on TABLE REFERENCE_RELATION
ALTER TABLE reference_relation OWNER TO nex;
GRANT ALL ON  reference_relation TO nex;
REVOKE ALL ON reference_relation FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON reference_relation TO curator;
GRANT SELECT ON reference_relation TO PUBLIC;

-- Set priviledge on TABLE REFERENCEDELETED
ALTER TABLE referencedeleted OWNER TO nex;
GRANT ALL ON  referencedeleted TO nex;
REVOKE ALL ON referencedeleted FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON referencedeleted TO curator;
GRANT SELECT ON referencedeleted TO PUBLIC;

-- Set priviledge on TABLE REFERENCEDOCUMENT
ALTER TABLE referencedocument OWNER TO nex;
GRANT ALL ON  referencedocument TO nex;
REVOKE ALL ON referencedocument FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON referencedocument TO curator;
GRANT SELECT ON referencedocument TO PUBLIC;

-- Set priviledge on TABLE REFERENCE_FILE
ALTER TABLE reference_file OWNER TO nex;
GRANT ALL ON  reference_file TO nex;
REVOKE ALL ON reference_file FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON reference_file TO curator;
GRANT SELECT ON reference_file TO PUBLIC;

-- Set priviledge on TABLE REFERENCETYPE
ALTER TABLE referencetype OWNER TO nex;
GRANT ALL ON  referencetype TO nex;
REVOKE ALL ON referencetype FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON referencetype TO curator;
GRANT SELECT ON referencetype TO PUBLIC;

-- Set priviledge on TABLE REFERENCEUNLINK
ALTER TABLE referenceunlink OWNER TO nex;
GRANT ALL ON  referenceunlink TO nex;
REVOKE ALL ON referenceunlink FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON referenceunlink TO curator;
GRANT SELECT ON referenceunlink TO PUBLIC;

-- Set priviledge on TABLE REFERENCE_URL
ALTER TABLE reference_url OWNER TO nex;
GRANT ALL ON  reference_url TO nex;
REVOKE ALL ON reference_url FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON reference_url TO curator;
GRANT SELECT ON reference_url TO PUBLIC;

-- Set priviledge on TABLE REGULATIONANNOTATION  
ALTER TABLE regulationannotation OWNER TO nex;
GRANT ALL ON  regulationannotation TO nex;
REVOKE ALL ON regulationannotation FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON regulationannotation TO curator;
GRANT SELECT ON regulationannotation TO PUBLIC;

-- Set priviledge on TABLE REPORTER
ALTER TABLE reporter OWNER TO nex;
GRANT ALL ON  reporter TO nex;
REVOKE ALL ON reporter FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON reporter TO curator;
GRANT SELECT ON reporter TO PUBLIC;

-- Set priviledge on TABLE RESERVEDNAME
ALTER TABLE reservedname OWNER TO nex;
GRANT ALL ON  reservedname TO nex;
REVOKE ALL ON reservedname FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON reservedname TO curator;
GRANT SELECT ON reservedname TO PUBLIC;

-- Set priviledge on TABLE RO
ALTER TABLE ro OWNER TO nex;
GRANT ALL ON  ro TO nex;
REVOKE ALL ON ro FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON ro TO curator;
GRANT SELECT ON ro TO PUBLIC;

-- Set priviledge on TABLE RO_RELATION
ALTER TABLE ro_relation OWNER TO nex;
GRANT ALL ON  ro_relation TO nex;
REVOKE ALL ON ro_relation FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON ro_relation TO curator;
GRANT SELECT ON ro_relation TO PUBLIC;

-- Set priviledge on TABLE RO_URL
ALTER TABLE ro_url OWNER TO nex;
GRANT ALL ON  ro_url TO nex;
REVOKE ALL ON ro_url FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON ro_url TO curator;
GRANT SELECT ON ro_url TO PUBLIC;

-- Set priviledge on TABLE SGDID
ALTER TABLE sgdid OWNER TO nex;
GRANT ALL ON  sgdid TO nex;
REVOKE ALL ON sgdid FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE ON sgdid TO curator;
GRANT SELECT ON sgdid TO PUBLIC;

-- Set priviledge on TABLE SO
ALTER TABLE so OWNER TO nex;
GRANT ALL ON  so TO nex;
REVOKE ALL ON so FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON so TO curator;
GRANT SELECT ON so TO PUBLIC;

-- Set priviledge on TABLE SOURCE
ALTER TABLE source OWNER TO nex;
GRANT ALL ON  source TO nex;
REVOKE ALL ON source FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON source TO curator;
GRANT SELECT ON source TO PUBLIC;

-- Set priviledge on TABLE SO_ALIAS
ALTER TABLE so_alias OWNER TO nex;
GRANT ALL ON  so_alias TO nex;
REVOKE ALL ON so_alias FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON so_alias TO curator;
GRANT SELECT ON so_alias TO PUBLIC;

-- Set priviledge on TABLE SO_RELATION
ALTER TABLE so_relation OWNER TO nex;
GRANT ALL ON  so_relation TO nex;
REVOKE ALL ON so_relation FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON so_relation TO curator;
GRANT SELECT ON so_relation TO PUBLIC;

-- Set priviledge on TABLE SO_URL
ALTER TABLE so_url OWNER TO nex;
GRANT ALL ON  so_url TO nex;
REVOKE ALL ON so_url FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON so_url TO curator;
GRANT SELECT ON so_url TO PUBLIC;

-- Set priviledge on TABLE STRAINDBENTITY
ALTER TABLE straindbentity OWNER TO nex;
GRANT ALL ON  straindbentity TO nex;
REVOKE ALL ON straindbentity FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON straindbentity TO curator;
GRANT SELECT ON straindbentity TO PUBLIC;

-- Set priviledge on TABLE STRAINSUMMARY
ALTER TABLE strainsummary OWNER TO nex;
GRANT ALL ON  strainsummary TO nex;
REVOKE ALL ON strainsummary FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON strainsummary TO curator;
GRANT SELECT ON strainsummary TO PUBLIC;

-- Set priviledge on TABLE STRAINSUMMARY_REFERENCE
ALTER TABLE strainsummary_reference OWNER TO nex;
GRANT ALL ON  strainsummary_reference TO nex;
REVOKE ALL ON strainsummary_reference FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON strainsummary_reference TO curator;
GRANT SELECT ON strainsummary_reference TO PUBLIC;

-- Set priviledge on TABLE STRAIN_URL
ALTER TABLE strain_url OWNER TO nex;
GRANT ALL ON  strain_url TO nex;
REVOKE ALL ON strain_url FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON strain_url TO curator;
GRANT SELECT ON strain_url TO PUBLIC;

-- Set priviledge on TABLE TAXONOMY
ALTER TABLE taxonomy OWNER TO nex;
GRANT ALL ON  taxonomy TO nex;
REVOKE ALL ON taxonomy FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON taxonomy TO curator;
GRANT SELECT ON taxonomy TO PUBLIC;

-- Set priviledge on TABLE TAXONOMY_ALIAS
ALTER TABLE taxonomy_alias OWNER TO nex;
GRANT ALL ON  taxonomy_alias TO nex;
REVOKE ALL ON taxonomy_alias FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON taxonomy_alias TO curator;
GRANT SELECT ON taxonomy_alias TO PUBLIC;

-- Set priviledge on TABLE TAXONOMY_RELATION
ALTER TABLE taxonomy_relation OWNER TO nex;
GRANT ALL ON  taxonomy_relation TO nex;
REVOKE ALL ON taxonomy_relation FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON taxonomy_relation TO curator;
GRANT SELECT ON taxonomy_relation TO PUBLIC;

-- Set priviledge on TABLE TAXONOMY_URL
ALTER TABLE taxonomy_url OWNER TO nex;
GRANT ALL ON  taxonomy_url TO nex;
REVOKE ALL ON taxonomy_url FROM PUBLIC;
GRANT SELECT,INSERT,UPDATE,DELETE ON taxonomy_url TO curator;
GRANT SELECT ON taxonomy_url TO PUBLIC;

-- Set priviledge on TABLE UPDATELOG
ALTER TABLE updatelog OWNER TO nex;
GRANT ALL ON  updatelog TO nex;
REVOKE ALL ON updatelog FROM PUBLIC;
GRANT SELECT,INSERT ON updatelog TO curator;
GRANT SELECT ON updatelog TO PUBLIC;

-- Set priviledge on SEQUENCES
grant select on all sequences in schema nex to curator;

-- Set priviledge on FUNCTIONS
grant execute on all functions in schema nex to curator;
