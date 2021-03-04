#! /usr/bin/env python3
# author : Picolo Floriane

from Bio import Phylo
from math import *

order = ["Yeast", "Metazoa", "Bilateria", "Chordata", "Vertebrata", "Clupeocephala", "Sarcopterygii", "Tetrapoda", "Amniota", "Mammalia"]
d = {"Yeast": [], "Metazoa": [], "Bilateria": [], "Chordata": [], "Vertebrata": [], "Clupeocephala": [], "Sarcopterygii": [], "Tetrapoda": [], "Amniota": [], "Mammalia": []}

list_genes = []
with open("list_genes.txt") as infile:
    for row in infile: 
        list_genes.append(row.strip())


# list_genes = ["ENSG00000171862"]

t = Phylo.read("vertebrates_species-tree_Ensembl.nh", "newick")
terminals = t.get_terminals()


for specie in terminals:
    path = str(t.get_path(specie))
    if "Mammalia" in path: 
        d["Mammalia"].append(str(specie))
    elif "Amniota" in path:
        d["Amniota"].append(str(specie))
    elif "Tetrapoda" in path:
        d["Tetrapoda"].append(str(specie))
    elif "Sarcopterygii" in path:
        d["Sarcopterygii"].append(str(specie))
    elif "Clupeocephala" in path:
        d["Clupeocephala"].append(str(specie))
    elif "Vertebrata" in path: 
        d["Vertebrata"].append(str(specie))
    elif "Chordata" in path:
        d["Chordata"].append(str(specie))
    elif "Bilateria" in path: 
        d["Bilateria"].append(str(specie))
    elif "Metazoa" in path: 
        d["Metazoa"].append(str(specie))
    else: 
        d["Yeast"].append(str(specie))

print(d)

trees = Phylo.parse("/home/fpicolo/Téléchargements/protein_trees.nhx", "newick")
for tree in trees: # pour chaque arbre du fichier arbres

# tree = Phylo.read("ENSG00000171862.nh", "newick")
    terminals = tree.get_terminals()
    if "Homo.sapiens" in str(terminals): # si on retrouve un gène humain dans l'arbre
        all_human_id = list(tree.find_elements({"comment":".*Homo.sapiens.*"}))
        for geneid in all_human_id:
            if str(geneid) in list_genes:
                dico = {"Yeast": [], "Metazoa": [], "Bilateria": [], "Chordata": [], "Vertebrata": [], "Clupeocephala": [], "Sarcopterygii": [], "Tetrapoda": [], "Amniota": [], "Mammalia": []}
                for genecode in terminals: # pour toutes les espèces de l'arbre du gène de notre liste
                    path = str(tree.get_path(genecode)).split("),") # on récupère les chemins de orthologue de notre gène
                    specie = path[-1].split("=")[2][:-2].replace(".","_")
                    if specie in str(d["Mammalia"]) or "Mammalia" in str(path):
                        dico["Mammalia"].append(specie) # .append(genecode) pour récupérer l'id du gène de l'espèce
                    elif specie in str(d["Amniota"]) or "Amniota" in str(path):
                        dico["Amniota"].append(specie)
                    elif specie in str(d["Tetrapoda"]) or "Tetrapoda" in str(path):
                        dico["Tetrapoda"].append(specie)
                    elif specie in str(d["Sarcopterygii"]) or "Sarcopterygii" in str(path):
                        dico["Sarcopterygii"].append(specie)
                    elif specie in str(d["Clupeocephala"]) or "Clupeocephala" in str(path):
                        dico["Clupeocephala"].append(specie)
                    elif specie in str(d["Vertebrata"]) or "Vertebrata" in str(path): 
                        dico["Vertebrata"].append(specie)
                    elif specie in str(d["Chordata"]) or "Chordata" in str(path):
                        dico["Chordata"].append(specie)
                    elif specie in str(d["Bilateria"]) or "Bilateria" in str(path): 
                        dico["Bilateria"].append(specie)
                    elif specie in str(d["Metazoa"]) or "Metazoa" in str(path): 
                        dico["Metazoa"].append(specie)
                    else: # le reste c'est "forcément" des levures
                        dico["Yeast"].append(specie)

                for clade in order:
                    dico[clade] = list(set(dico[clade])) # on supprime les paralogues
                    if len(dico[clade]) >= ceil((80*len(d[clade]))/100) and len(dico[clade]) > 0: # si le gène est présent dans 
                        print(str(geneid), clade)
                        break # on s'arrête à ce clade

                list_genes.remove(str(geneid)) # on supprime le gène de notre liste

print(list_genes)
    


