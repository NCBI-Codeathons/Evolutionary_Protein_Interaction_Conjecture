#!/bin/python
#sca_test.py

import matplotlib.pyplot as plt
import coevo2 as ce
import itertools as it
import numpy as np
import copy
import time

reload(ce)

names = ['glgA', 'glgC', 'cydA', 'cydB']
algPath = 'TestSet/eggNOG_aligns/slice_0.9/'
prots = ce.prots_from_scratch(names,path2alg=algPath)
ps = ce.ProtSet(prots,names)

phylo_names = ['aspS','ffh','lepA','pgk','recN','rho','rpoA','ruvB','tig','uvrB']
phylo_prots = ce.prots_from_scratch(phylo_names,path2alg='TestSet/eggNOG_aligns/phylogenes/')

phylo2 = ce.PhyloSet(phylo_prots)
phylo2.set_indexer(thresh=7)
for pt in phylo2.prots:  # temporary fix for duplicated locus ids in the same msa
    pt.msa = pt.msa[~pt.msa.index.duplicated(keep='first')]
phylo2.set_sim_mat()

protsmats,pairmats,pairrandmats,sca_score,sca_score2 = ce.sca(ps,phylo2,delta=0.0001)

for pt,sca in it.izip(ps.prots,protsmats): pt.sca_mat = sca
for pair,sca_cat in it.izip(ps.pairs,pairmats): pair.sca_mat = sca_cat
np.save('GettingStartedSCACalcs.npy',ps)

print(sca_score2)
