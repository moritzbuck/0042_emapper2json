from tqdm import tqdm
from ete3 import NCBITaxa
import sqlite3

ncbi = NCBITaxa()

header = ["query_name" , "ortholog" , "evalue" , "score" , "taxo" , "name" , "GO-terms" , "EC" , "KO" , "Pathway" , "Module" , "Reaction" , "rclass" , "BRITE" , "KEGG_TC" , "CAZy " , "BiGG" , "tax_scope" , "OG" , "deprecated" , "COG_category" , "description"]

aa2eggnog = dict()
with open("data/Erken.emapper.annotations") as handle:
     for l in tqdm(handle):
         lls = l[:-1].split("\t")
         aas_id = lls[0]
         aa2eggnog[aas_id] = dict(zip(header, lls))


for k,v in tqdm(aa2eggnog.items()):
         aa2eggnog[k]['OG'] = { vv.split("@")[1] : vv.split("@")[0] for vv in v['OG'].split(",")}

all_tax_ids = {vv for v in aa2eggnog.values() for vv in v['OG'].keys()}

exepctions = {'85005' : '2037',
              '206350' : '32003',
              '1212' : '1213',
              '35237' : '51368',
              '29258' : '35342',
              '358033' : '59732',
              '34383' : '1593277',
              '713636' : '32003'
            }


tax2rank = { tax : len(ncbi.get_lineage(tax) if tax not in exepctions else ncbi.get_lineage(exepctions[tax])) for tax in tqdm(all_tax_ids)}
for k,v in tqdm(aa2eggnog.items()):
         aa2eggnog[k]['rankOG'] = { tax2rank[kk] : vv for kk,vv in v['OG'].items()}

all_OGs = {vv for v in aa2eggnog.values() for vv in v['OG'].values()}
all_OG_tuples = {tuple([vv[1] for vv in sorted(v['rankOG'].items(), key = lambda x: x[0])]) for v in aa2eggnog.values()}
all_OG_tuples = [list(OG) for OG in all_OG_tuples]


og2ko = {og : dict() for og in all_OGs}
og2count = {og : 0 for og in all_OGs}

for aa, vv in tqdm(aa2eggnog.items()):
    for og in vv['OG'].values():
        ko = vv['KO']
        if ko != '':
            for kkoo in ko.split(","):
                    if ko not in og2ko[og]:
                        og2ko[og][kkoo] = 1
                    else :
                        og2ko[og][kkoo] += 1
        og2count[og] +=1

og2freqs = { k: {kk : vv/og2count[k] for kk, vv in d.items()} for k, d in og2ko.items()}

conn = sqlite3.connect('example.db')
conn = sqlite3.connect('/home/moritz/kadath/miniconda3/envs/eggnogmapper/lib/python2.7/site-packages/data/eggnog.db') 

tt = conn.cursor()

# FIND ALL TABLE NAMES
#
#all_tables_req = """SELECT name FROM
#     (SELECT * FROM sqlite_master UNION ALL
#      SELECT * FROM sqlite_temp_master)
#  WHERE type='table'
#  ORDER BY name"""
#tt.execute(all_tables_req)
#tt.execute("PRAGMA table_info(OG);")
