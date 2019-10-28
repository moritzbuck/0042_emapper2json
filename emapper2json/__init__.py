import os
import shutil
import sys
from os.path import join as pjoin

_N_HEADLINES_ = 4
_HEADER_ = ['query_name',
 'seed_eggNOG_ortholog',
 'seed_ortholog_evalue',
 'seed_ortholog_score',
 'best_tax_level',
 'Preferred_name',
 'GOs',
 'EC',
 'KEGG_ko',
 'KEGG_Pathway',
 'KEGG_Module',
 'KEGG_Reaction',
 'KEGG_rclass',
 'BRITE',
 'KEGG_TC',
 'CAZy',
 'BiGG_Reaction',
 'tax_scope',
 'eggNOG_OGs',
 'bestOG',
 'COG_functional_cat',
 'notes'
]

_COG_CATS_ = {
"A": "RNA processing and modification",
"B": "Chromatin Structure and dynamics",
"C": "Energy production and conversion",
"D": "Cell cycle control and mitosis",
"E": "Amino Acid metabolis and transport",
"F": "Nucleotide metabolism and transport",
"G": "Carbohydrate metabolism and transport",
"H": "Coenzyme metabolis",
"I": "Lipid metabolism",
"J": "Translation",
"K": "Transcription",
"L": "Replication and repair",
"M": "Cell wall/membrane/envelop biogenesis",
"N": "Cell motility",
"O": "Post-translational modification, protein turnover, chaperone functions",
"P": "Inorganic ion transport and metabolism",
"Q": "Secondary Structure",
"T": "Signal Transduction",
"U": "Intracellular trafficing and secretion",
"V": "Defense mechanisms",
"W": "Extracellular structures",
"Y": "Nuclear structure",
"Z": "Cytoskeleton",
"R": "General Functional Prediction only",
"S": "Function Unknown"
}

def parse_file(file):

    with open( file ) as handle:
        for i in range(_N_HEADLINES_):
            handle.readline()
        #head = handle.readline()[1:-1].split("\t") + [ 'notes']
        big_dict = { l.split("\t")[0] : {h : vv for h, vv in zip(_HEADER_[1:], l[:-1].split("\t")[1:])} for l in handle}
    return big_dict

def compose_annot(annot_set):
    out_data =  {}
    if type(annot_set) == dict:
        out_data['genes'] = list(annot_set.keys())
        out_data['genomes'] = {"_".join(g.split("_")[:-1]) for g in out_data['genes']}
        annot_set = list(annot_set.values())
    else :
        out_data['genes'] = None
        out_data['genomes'] = None

    comp_dat = {k : [dd  for d in annot_set for dd in d[k].split(",") if dd != ''] for k in annot_set[0]}
    del comp_dat['seed_ortholog_evalue']
    del comp_dat['seed_ortholog_score']
    del comp_dat['bestOG']
    comp_dat = {k : {dd : v.count(dd)/len(annot_set) for dd in set(v)} for k, v in comp_dat.items()}
    for k, v in comp_dat.items():
        if len(v) == 0:
            comp_dat[k] = None
    out_data.update(comp_dat)
    return out_data



def get_blocks(keggs, def_line, andline = True):
    def_line = def_line.split() if andline else def_line.split(",")
    i = 0
    blocks = []
    fwd = None
    rev = None
    for v in def_line:
        blocks += [i]
        if v.startswith("("):
            fwd = v.count("(")
            rev = v.count(")")
        elif fwd:
            fwd += v.count("(")
            rev += v.count(")")
        if fwd == rev:
            i += 1
            fwd = None
            rev = None
    operation = " " if andline else ","
    big_operation = lambda x: sum(x)/len(x) if andline else max(x)
    and_block = [operation.join([b for i,b in enumerate(def_line) if blocks[i] == block]) for block in set(blocks)]
    and_block = [b[1:-1] if b.startswith("(") else b for b in and_block]
    return big_operation([get_blocks(keggs, b, not andline) if ("," in b or " " in b) else int((b in keggs) if not b.startswith("-") else True) for b in and_block])
