from src.models import Locusdbentity
from scripts.loading.database_session import get_session

__author__ = 'sweng66'


oldfile = "scripts/dumping/curation/data/SGD_1.0.0.0_gff3.gff3.old"

def dump_data():
 
    nex_session = get_session()

    name_to_desc = dict([(x.systematic_name, x.name_description) for x in nex_session.query(Locusdbentity).all()])

    f = open(oldfile)
    
    for line in f:
        line = line.rstrip()
        # line = line.replace("gene=", "Name=").replace("Note=", "description=")
        line = line.replace("Note=", "description=")            
        field = line.split("\t")
        if len(field) >= 9 and field[8].startswith("ID="):
            ## no gene name without alias example:
            # ID=YAL064W-B;Name=YAL064W-B;Ontology_term=GO:0003674,GO:0005783,GO:0008150,GO:0016021;Note=Fungal-specific%20protein%20of%20unknown%20function;display=Fungal-specific%20protein%20of%20unknown%20function;dbxref=SGD:S000002141;orf_classification=Uncharacterized;curie=SGD:S000002141

            ## no gene name with alias field 
            # ID=YGR161C-D;Name=YGR161C-D;Alias=gag-pol%20fusion%20protein;Ontology_term=GO:0000943,GO:0003723,GO:0003887,GO:0003964,GO:0004540,GO:0005634,GO:0008233,GO:0032197;Note=Retrotransposon%20TYA%20Gag%20and%20TYB%20Pol%20genes%3B%20transcribed%2Ftranslated%20as%20one%20unit%3B%20polyprotein%20is%20processed%20to%20make%20a%20nucleocapsid-like%20protein%20%28Gag%29%2C%20reverse%20transcriptase%20%28RT%29%2C%20protease%20%28PR%29%2C%20and%20integrase%20%28IN%29%3B%20similar%20to%20retroviral%20genes;display=Retrotransposon%20TYA%20Gag%20and%20TYB%20Pol%20genes;dbxref=SGD:S000007368;curie=SGD:S000007368
            ## with a gene name example
            # ID=YFL039C;Name=YFL039C;gene=ACT1;Alias=ACT1,ABY1,END7,actin;Ontology_term=GO:0000001,GO:0000011,GO:0000132,GO:0000142,GO:0000812,GO:0000916,GO:0001300,GO:0004402,GO:0005200,GO:0005884,GO:0006281,GO:0006887,GO:0006897,GO:0007119,GO:0009306,GO:0016573,GO:0030010,GO:0030050,GO:0030476,GO:0030479,GO:0031011,GO:0031505,GO:0032432,GO:0034599,GO:0035267;Note=Actin%3B%20structural%20protein%20involved%20in%20cell%20polarization%2C%20endocytosis%2C%20and%20other%20cytoskeletal%20functions;display=Actin;dbxref=SGD:S000001855;orf_classification=Verified;curie=SGD:S000001855

            # ID=YFL039C;Name=YFL039C;gene=ACT1;full_name=ACTin;Alias=ACT1,ABY1,END7,actin;Ontology_term=GO:0000001,GO:0000011,GO:0000132,GO:0000142,GO:0000812,GO:0000916,GO:0001300,GO:0004402,GO:0005200,GO:0005884,GO:0006281,GO:0006887,GO:0006897,GO:0007119,GO:0009306,GO:0016573,GO:0030010,GO:0030050,GO:0030476,GO:0030479,GO:0031011,GO:0031505,GO:0032432,GO:0034599,GO:0035267;description=Actin%3B%20structural%20protein%20involved%20in%20cell%20polarization%2C%20endocytosis%2C%20and%20other%20cytoskeletal%20functions;display=Actin;dbxref=SGD:S000001855;orf_classification=Verified;curie=SGD:S000001855

            pieces = field[8].split(";")
            ID = pieces[0].split("=")[1]
            if "gene=" in line:
                name = pieces[1].split("=")[1]
                if "Alias=" in line:
                    field[8] = pieces[0] + ";" + pieces[2].replace("gene=", "Name=") + ";" + "Alias=" + name + "," + pieces[3].split("=")[1] + ";" + ";".join(pieces[4:])
                else:
                    field[8] = pieces[0] + ";" + pieces[2].replace("gene=", "Name=") + ";" + "Alias=" + name + ";" + ";".join(pieces[3:])
                
                if name_to_desc.get(ID) is not None:
                    
                    name_desc = name_to_desc[ID].replace(" ", "%20")
                    pieces = field[8].split(";")
                    field[8] = ";".join(pieces[0:3]) + ";full_name=" + name_desc + ";" + ";".join(pieces[3:])
              
            else:
                if name_to_desc.get(ID) is not None:
                    name_desc = name_to_desc[ID].replace(" ", "%20")
                    if "Alias=" in line:
                        field[8] = ";".join(pieces[0:3]) + ";full_name=" + name_desc + ";" + ";".join(pieces[3:])
                    else:
                        field[8] = ";".join(pieces[0:2]) + ";full_name=" + name_desc + ";" + ";".join(pieces[2:])

            line = "\t".join(field)
                
        print line
            
    nex_session.close()


if __name__ == '__main__':
    
    dump_data()
