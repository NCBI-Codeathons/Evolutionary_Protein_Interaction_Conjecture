#!/bin/python
#mirror_tree_test.py

#python mirror_tree_test.py 0.01 0.9

import matplotlib.pyplot as plt
import coevo2 as ce
import itertools as it
import numpy as np
import copy
import time
import sys

## Change the Gene names 
names = ['glgA','glgC', 'cydA', 'cydB']
algPath = 'TestSet/eggNOG_aligns/Slice_09/'
prots = ce.prots_from_scratch(names,path2alg=algPath)
ps = ce.ProtSet(prots,names)

phylo_names = ['aspS','ffh','lepA','pgk','recN','rho','rpoA','ruvB','tig','uvrB']
phylo_prots = ce.prots_from_scratch(phylo_names,path2alg='TestSet/eggNOG_aligns/phylogenes/')


phylo2 = ce.PhyloSet(phylo_prots)
phylo2.set_indexer(thresh=7)
for pt in phylo2.prots:  # temporary fix for duplicated locus ids in the same msa
	pt.msa = pt.msa[~pt.msa.index.duplicated(keep='first')]
phylo2.set_sim_mat()

mt_mat2 = ce.mirror_tree(ps,phylo2, float(sys.argv[1]))
filename = "MirrorTree_Delta_"+sys.argv[1]+"_Similarity_"+sys.argv[2]+".npy"
np.save(filename,[mt_mat2])
