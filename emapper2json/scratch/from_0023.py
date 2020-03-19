header = ["query_name" , "ortholog" , "evalue" , "score" , "taxo" , "name" , "GO-terms" , "EC" , "KO" , "Pathway" , "Module" , "Reaction" , "rclass" , "BRITE" , "KEGG_TC" , "CAZy " , "BiGG" , "tax_scope
" , "OG" , "deprecated" , "COG_category" , "description"]


bin2og = dict()
with open("data/output_file.emapper.annotations") as handle:
    for l in tqdm(handle):
        lls = l[:-1].split("\t")
        bin_id = "_".join(lls[0].split("_")[:-1])
        if bin_id not in bin2og:
            bin2og[bin_id] = set()
        bin2og[bin_id].add(lls[18].split(",")[0].split("@")[0])


og2md = dict() 
with open("data/output_file.emapper.annotations") as handle:
    for l in tqdm(handle):
        lls = l[:-1].split("\t")
        og_id = lls[18].split(",")[0].split("@")[0]
        if og_id not in og2md:
            og2md[og_id] = { a : b for a,b in zip(header[5:],lls[5:])}


inter_size = {(k,l) : len(v.intersection(w))/len(v) for k,v in tqdm(bin2og.items()) for l,w  in bin2og.items() if len(w) > 20 and len(v.intersection(w))/len(v) > 0.5}
