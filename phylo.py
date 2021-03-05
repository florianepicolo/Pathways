#! /usr/bin/env python3
# author : Picolo Floriane

from Bio import Phylo
from math import *

def genefile_to_genelist(genefile):
    """
        retrieves a file containing a gene list in Ensembl format () and outputs a list from it
    """
    gene_list = list()
    with open(genefile) as infile:
        for gene in infile: 
            gene_list.append(gene.strip())
    return gene_list

def explore_perfect_tree(perfecttree_file):
    """
        the perfect tree corresponds to the tree including all the species available on Ensembl
    """
    dico_perfecttree = {"Yeast": [], "Metazoa": [], "Bilateria": [], "Chordata": [], "Vertebrata": [], "Clupeocephala": [], "Sarcopterygii": [], "Tetrapoda": [], "Amniota": [], "Mammalia": []}
    perfecttree = Phylo.read(perfecttree_file, "newick")
    all_species = perfecttree.get_terminals() # returns all species of tree

    for specie in all_species:
        path = str(perfecttree.get_path(specie)) # phylogenetic path of specie
        if "Mammalia" in path: 
            dico_perfecttree["Mammalia"].append(str(specie))
        elif "Amniota" in path:
            dico_perfecttree["Amniota"].append(str(specie))
        elif "Tetrapoda" in path:
            dico_perfecttree["Tetrapoda"].append(str(specie))
        elif "Sarcopterygii" in path:
            dico_perfecttree["Sarcopterygii"].append(str(specie))
        elif "Clupeocephala" in path:
            dico_perfecttree["Clupeocephala"].append(str(specie))
        elif "Vertebrata" in path: 
            dico_perfecttree["Vertebrata"].append(str(specie))
        elif "Chordata" in path:
            dico_perfecttree["Chordata"].append(str(specie))
        elif "Bilateria" in path: 
            dico_perfecttree["Bilateria"].append(str(specie))
        elif "Metazoa" in path: 
            dico_perfecttree["Metazoa"].append(str(specie))
        else: 
            dico_perfecttree["Yeast"].append(str(specie))
    return dico_perfecttree


def explore_trees(treesfile, gene_list, orderclade, selected_specie, dico_perfecttree):
    trees = Phylo.parse(treesfile, "newick")
    for tree in trees: 
    # tree = Phylo.read("ENSG00000171862.nh", "newick")
        all_lastclade = tree.get_terminals() # returns all last clade of tree
        if selected_specie in str(all_lastclade): # if our human gene is found in the tree
            all_genesid = list(tree.find_elements({"comment":".*" + selected_specie + ".*"}))
            for geneid in all_genesid:
                if str(geneid) in gene_list:
                    dico_tree = {"Yeast": [], "Metazoa": [], "Bilateria": [], "Chordata": [], "Vertebrata": [], "Clupeocephala": [], "Sarcopterygii": [], "Tetrapoda": [], "Amniota": [], "Mammalia": []}
                    for genecode in all_lastclade:
                        path = str(tree.get_path(genecode)).split("),") # we retrieve the orthologic pathways of our gene
                        specie = path[-1].split("=")[2][:-2].replace(".","_")
                        if specie in str(dico_perfecttree["Mammalia"]) or "Mammalia" in str(path):
                            dico_tree["Mammalia"].append(specie) # .append(genecode) to retrieve the id of the gene of the species
                        elif specie in str(dico_perfecttree["Amniota"]) or "Amniota" in str(path):
                            dico_tree["Amniota"].append(specie)
                        elif specie in str(dico_perfecttree["Tetrapoda"]) or "Tetrapoda" in str(path):
                            dico_tree["Tetrapoda"].append(specie)
                        elif specie in str(dico_perfecttree["Sarcopterygii"]) or "Sarcopterygii" in str(path):
                            dico_tree["Sarcopterygii"].append(specie)
                        elif specie in str(dico_perfecttree["Clupeocephala"]) or "Clupeocephala" in str(path):
                            dico_tree["Clupeocephala"].append(specie)
                        elif specie in str(dico_perfecttree["Vertebrata"]) or "Vertebrata" in str(path): 
                            dico_tree["Vertebrata"].append(specie)
                        elif specie in str(dico_perfecttree["Chordata"]) or "Chordata" in str(path):
                            dico_tree["Chordata"].append(specie)
                        elif specie in str(dico_perfecttree["Bilateria"]) or "Bilateria" in str(path): 
                            dico_tree["Bilateria"].append(specie)
                        elif specie in str(dico_perfecttree["Metazoa"]) or "Metazoa" in str(path): 
                            dico_tree["Metazoa"].append(specie)
                        elif specie in str(dico_perfecttree["Yeast"]) or "Yeast" in str(path):
                            dico_tree["Yeast"].append(specie)
                        else:
                            print(path)

                    for clade in orderclade:
                        dico_tree[clade] = list(set(dico_tree[clade])) # deletion of paralogues
                        if len(dico_tree[clade]) >= ceil((80*len(dico_perfecttree[clade]))/100) and len(dico_tree[clade]) > 0: # if the gene is present in sufficient species of clade
                            # print(str(geneid), clade)
                            break # we have the most ancestral clade

                    gene_list.remove(str(geneid)) # deletion of the gene from our initial list

    print(gene_list) # list of genes without trees
        


### main
genefile = "list_genes.txt"
perfecttreefile = "vertebrates_species-tree_Ensembl.nh"
treesfile = "/home/fpicolo/Téléchargements/protein_trees.nhx"
selectedspecie = "Homo.sapiens"

# genelist = ["ENSG00000171862"]
genelist = genefile_to_genelist(genefile)
dico = explore_perfect_tree(perfecttreefile)

orderclade = ["Yeast", "Metazoa", "Bilateria", "Chordata", "Vertebrata", "Clupeocephala", "Sarcopterygii", "Tetrapoda", "Amniota", "Mammalia"]


explore_trees(treesfile, genelist, orderclade, selectedspecie, dico)
