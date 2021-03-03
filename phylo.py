#! /usr/bin/env python3
# author : Picolo Floriane

from Bio import Phylo
from math import *

order = ["Yeast", "Metazoa", "Bilateria", "Chordata", "Vertebrata", "Clupeocephala", "Sarcopterygii", "Tetrapoda", "Amniota", "Mammalia"]
yeast = []; metazoa = []; bilateria = []; chordata = []; vertebrata = []; clupeocephala = [] #teleosts
sarcopterygii = []; tetrapoda = []; amniota = []; mammalia = []

d = {"Yeast": yeast, "Metazoa": metazoa, "Bilateria": bilateria, "Chordata": chordata, "Vertebrata": vertebrata, "Clupeocephala": clupeocephala, "Sarcopterygii": sarcopterygii, "Tetrapoda": tetrapoda, "Amniota": amniota, "Mammalia": mammalia}


list_genes = ["ENSG00000171862"]

# list_genes = ["ENSG00000118260","ENSG00000198055","ENSG00000171862","ENSG00000118515","ENSG00000123159","ENSG00000063046","ENSG00000198873","ENSG00000104205","ENSG00000288602","ENSG00000106211","ENSG00000134308","ENSG00000137486","ENSG00000157500","ENSG00000173020","ENSG00000127955","ENSG00000114867","ENSG00000164742","ENSG00000139318","ENSG00000185345","ENSG00000079337","ENSG00000100077","ENSG00000087460","ENSG00000137841","ENSG00000158828"]

t = Phylo.read("vertebrates_species-tree_Ensembl.nh", "newick")
terminals = t.get_terminals()


for specie in terminals:
    path = str(t.get_path(specie))
    if "Mammalia" in path: 
        mammalia.append(str(specie))
    elif "Amniota" in path:
        amniota.append(str(specie))        
    elif "Tetrapoda" in path:
        tetrapoda.append(str(specie))
    elif "Sarcopterygii" in path:
        sarcopterygii.append(str(specie))
    elif "Clupeocephala" in path:
        clupeocephala.append(str(specie))
    elif "Vertebrata" in path: 
        vertebrata.append(str(specie))
    elif "Chordata" in path:
        chordata.append(str(specie))
    elif "Bilateria" in path: 
        bilateria.append(str(specie))
    elif "Metazoa" in path: 
        metazoa.append(str(specie))
    else: 
        yeast.append(str(specie))



# trees = Phylo.parse("/home/fpicolo/Téléchargements/protein_trees.nhx", "newick")
# for tree in trees: # pour chaque arbre du fichier arbres
tree = Phylo.read("ENSG00000171862.nh", "newick")
terminals = tree.get_terminals()
# print(terminals)
if "Homo.sapiens" in str(terminals): # si on retrouve un gène humain dans l'arbre
    path_human = list(tree.get_path({"comment": ".*Homo.sapiens.*"}))  
    print(path_human)
    try:
        B1 = []; B2 = []; B3 = []; B4 = []; B5 = []; B6 = []; B7 = []; B8 = []; B9 = []; B10 = []
        dico = {"Yeast": B1, "Metazoa": B2, "Bilateria": B3, "Chordata": B4, "Vertebrata": B5, "Clupeocephala": B6, "Sarcopterygii": B7, "Tetrapoda": B8, "Amniota": B9, "Mammalia": B10}

        if str(path_human[-1]) in list_genes: # si l'id du gène humain se trouve dans notre liste 
            for genecode in terminals: # pour toutes les espèces de l'arbre du gène de notre liste
                path = str(tree.get_path(genecode)).split("),") # on récupère les chemins de orthologue de notre gène
                specie = path[-1].split("=")[2][:-2].replace(".","_")
                
                if specie in str(mammalia) or "Mammalia" in str(path):
                    B10.append(specie) #B10.append(str(gencode)) pour récupérer l'id du gène de l'espèce
                elif specie in str(amniota) or "Amniota" in str(path):
                    B9.append(specie)
                elif specie in str(tetrapoda) or "Tetrapoda" in str(path):
                    B8.append(specie)
                elif specie in str(sarcopterygii) or "Sarcopterygii" in str(path):
                    B7.append(specie)
                elif specie in str(clupeocephala) or "Clupeocephala" in str(path):
                    B6.append(specie)
                elif specie in str(vertebrata) or "Vertebrata" in str(path): 
                    B5.append(specie)
                elif specie in str(chordata) or "Chordata" in str(path):
                    B4.append(specie)
                elif specie in str(bilateria) or "Bilateria" in str(path): 
                    B3.append(specie)
                elif specie in str(metazoa) or "Metazoa" in str(path): 
                    B2.append(specie)
                else: # le reste c'est "forcément" des levures
                    B1.append(specie)
            
            for clade in order:
                dico[clade] = list(set(dico[clade])) # on supprime les paralogues
                if len(dico[clade]) >= ceil((80*len(d[clade]))/100) and len(dico[clade]) > 0: # si le gène est présent dans 
                    print(clade)
                    break # on s'arrête à ce clade
            list_genes.remove(str(path_human[-1])) # on supprime le gène de notre liste
            # break       
    except IndexError: # contrer les arbres vides
        # print(tree)
        pass
print(list_genes)
    




