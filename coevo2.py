"""
"""

import numpy as np
from scipy.stats import gaussian_kde
import matplotlib.pyplot as plt
import pandas as pd
import itertools as it
import cPickle as pickle
import glob
import random
from multiprocessing import Pool

class Prot(object):
	'''
	A class for Protein Families. Contains their sequence alignments as well corresponding locus and taxonomic identifiers.

		 **Attributes:**
			-  `name`  =  E. coli family name
			-  `msa` =  Pandas DataFrame of the multiple sequence alignment (numerical format) indexed by locus id
			-  `indexer` =  Dictionary mapping from taxid (keys) to loci (strings in a list)
			-  `speciesmap` =  Pandas DataFrame mapping from locus (index) to taxid (values)
			-  `Npos`  =  Number of positions in the alignment

		 **Properties:**
			-  `aln`  =  Outputs the full multiple sequence alignment (numeric format) and a list of corresponding taxids

		 **Example:**
		   >>> pt = Prot(name, msa, indexer, speciesmap)
	'''
	def __init__(self, name, msa, indexer, speciesmap):
		self.name = name
		self.msa = msa # pandas DataFrame where the data is the numeric MSA and the indexes are gene IDs
		self.indexer = indexer # mapping from taxonomic IDs to gene IDs in dictionary form
		self.speciesmap = speciesmap # mapping from gene IDs to taxonomic IDs in DataFrame form, 1:1 with no duplicate indices
		self.Npos = msa.shape[1]

	def filter_loci(self, loci):
		'''
		Filter the multiple sequence alignment by a list of loci. Returns a numeric MSA in an np array.

			 **Arguments:**
			   -  Bound to an instance of Prot
			   -  `loci`  =  A list or 1D array of locus ids for the corresponding protein

			 **Example:**
			  >>> msa_filt = pt.filter_loci(loci)
		'''
		return self.msa.filter(items=loci, axis=0).values

	def filter_aln(self, includes=None):
		'''
		Filter the multiple sequence alignment by a list of taxids. Returns a numeric MSA in an np array, returns a list of taxids corresponding to rows of the MSA. Note that some taxids can be repeated due to paralogs, so the matrix axis will not necessarily match the input of taxids exactly.

			 **Arguments:**
			   -  Bound to an instance of Prot

			 **Keyword Arguments:**
			   -  `includes`  =  (optional) List or array of taxids

			 **Example:**
			   >>> msa_filt, tids = pt.filter_aln(species)
		'''
		if includes is None: includes = self.indexer.keys()
		tids, complete_loci = [],[]
		for (taxid,loci) in it.ifilter(lambda (taxid,loci): taxid in includes, self.indexer.iteritems()):
			tids+=[taxid]*len(loci)
			complete_loci+=loci
		return self.filter_loci(complete_loci), tids

	@ property
	def aln(self):
		'''
		Returns the numeric MSA and taxids for all loci in self.indexer

			 **Example:**
			   >>> msa, tids = pt.filter_aln(species)
		'''
		tids, complete_loci = [],[]
		for (taxid,loci) in self.indexer.iteritems():
			tids+=[taxid]*len(loci)
			complete_loci+=loci
		return self.filter_loci(complete_loci), tids

class ProtPair(object):
	'''
	A wrapper that handles pairs of protein families. Contains a reference to their contents and a bound function for pairing alignments.

		 **Attributes:**
			-  `names`  =  E. coli family names
			-  `prots` =  A list of two Prot Class objects which constitute the pair
			-  `ij` =  Index of this pair in a matrix representation of the ProtSet
			-  `overlap` =  List of overlapping taxids for the two proteins

		 **Example:**
		   >>> pair = ProtPair(pt1, pt2)
	'''
	def __init__(self, prot1, prot2, ij=None):
		self.names = [prot1.name, prot2.name]
		self.prots = [prot1, prot2]
		self.ij = ij # position of that pair in the protein-protein coevolution matrix
		self.overlap = list(set(prot1.indexer.keys()).intersection(prot2.indexer.keys())) # overlap of taxids

	def filter_alns(self, includes=None):
		'''
		Filter the multiple sequence alignment of each protein by a list of taxids. Returns a pair of row matched numeric MSAs in np arrays, returns a list of taxids corresponding to rows of the MSAs. Note that some taxids can be repeated due to paralogs, so the matrix axis will not necessarily match the input of taxids exactly.

			 **Arguments:**
			   -  Bound to an instance of ProtPair

			 **Keyword Arguments:**
			   -  `includes`  =  (optional) List or array of taxids

			 **Example:**
			   >>> msa1, msa2, tids = pt.filter_aln(species)
		'''
		p1,p2 = self.prots
		indexer1, indexer2 = p1.indexer, p2.indexer
		if not includes: overlap = self.overlap
		else: overlap = list(set(self.overlap).intersection(includes))
		loci1, loci2, taxids = match_loci(indexer1, indexer2, overlap)
		return p1.filter_loci(loci1), p2.filter_loci(loci2), taxids

class ProtSet(object):
	'''
	A class for collections of protein families. Contains Prot and ProtPair Class objects. Interaction analyses act on ProtSet.

		 **Attributes:**
			-  `prots`  =  List of protein families (Prot Class)
			-  `names` =  E. coli family names
			-  `pairs` =  List of protein pairs (ProtPair Class)
			-  `ijs` = Matrix indices corresponding to the pair list

		 **Example:**
		   >>> ps = ProtSet(prots)
	'''
	@staticmethod
	def from_file(path):
		'''
		Create a ProtSet object from a collection of Prots saved as a Pickle

			 **Arguments:**
			   -  `path`  =  Path to a Pickle produced by the create_db function

			 **Example:**
			   >>> ps = ProtSet.from_file(test_proteins.db)
		'''
		prots, names, extras = unpack_file(path)
		ps = ProtSet(prots)
		for key in extras: # this will add annotation matrices and what have you
			setattr(ps, key, extras[key])
		return ps

	def __init__(self, prots, names=None):
		self.prots = prots # list of Prot objects
		if names is not None: self.names = names
		else: self.names = [ pt.name for pt in prots ]
		self.pairs, self.ijs = _ppair_list(prots) # creates a list of all possible pairs and corresponding ij matrix indices

	def get_prot(self, name):
		'''
		Returns the Prot with the provided name

			 **Arguments:**
			   -  Bound to an instance of ProtSet
			   -  `name`  =  Name of the desired Prot

			 **Example:**
			   >>> pt = ps.get_prot('trpB')
		'''
		return next(pt for pt in self.prots if name==pt.name)

	def get_pair(self, name1, name2):
		'''
		Returns the ProtPair with the provided names

			 **Arguments:**
			   -  Bound to an instance of ProtSet
			   -  `name1`  =  Name of one of the Prots
			   -  `name2`  =  Name of one of the Prots

			 **Example:**
			   >>> pair = ps.get_pair('trpA','trpB')
		'''
		return next(pair for pair in self.pairs if name1 in pair.names and name2 in pair.names)

class PhyloSet(object):
	'''
	A collection of Prots to be used as a phylogenetic model for mirror-tree

		 **Attributes:**
			-  `prots`  =  List of protein families (Prot Class)
			-  `names` =  E. coli family names

		 **Example:**
		   >>> phylo1 = PhyloSet(prots)
	'''
	@staticmethod
	def from_file(path):
		'''
		Create a PhyloSet object from a collection of Prots saved as a Pickle

			 **Arguments:**
			   -  `path`  =  Path to a Pickle produced by the create_db function

			 **Example:**
			   >>> phylo1 = PhyloSet.from_file(test_proteins.db)
		'''
		prots, names, extras = unpack_file(path)
		phylo = PhyloSet(prots)
		for key in extras: # this will add annotation matrices and what have you
			setattr(phylo, key, extras[key])
		return phylo

	def __init__(self, prots, names=None):
		self.prots = prots
		if names is not None: self.names = names
		else: self.names = [ pt.name for pt in prots ]

	def info(self, tids, delta):
		'''
		Obtains the relevant phylogenetic information for a given list of taxonomic ids

			 **Arguments:**
			   -  Bound to an instance of ProtPair
			   -  `tids`  =  A list of taxonomic IDs

			 **Example:**
			   >>> phylo_vec, specw, pwsw, Meff = phylo.info(tids)
		'''
		phylo_sub = self.get_submat(tids)
		specw = specw_calc(phylo_sub,delta)
		pwsw = triu_flatten(np.outer(specw, specw))
		Meff = pwsw.sum()
		phylo_vec = triu_flatten(phylo_sub)
		phylo_vec = phylo_vec/np.sqrt(wdot(phylo_vec,phylo_vec,pwsw))
		return phylo_vec, specw, pwsw, Meff

	def get_submat(self, tids):
		'''
		Subsample the phylogenetic simlarity matrix based on the provided taxonomic ids. Can handle duplicate ids for paralogs.

			 **Arguments:**
			   -  Bound to an instance of ProtPair
			   -  `tids`  =  A list of taxonomic IDs

			 **Example:**
			   >>> phylo_sub = phylo.get_submat(tids)
		'''
		ix = self.indexer.filter(items=tids,axis=0).values
		return self.sim_mat[ix,ix.T]

	def set_indexer(self, thresh=5):
		'''
		Determine the set of species for which phylogenic similarity will be modeled.

		 **Arguments:**
		   -  Bound to an instance of PhyloSet

		 **Keyword Arguments:**
		   -  `thresh`  =  Minimum number of PhyloSet proteins in a species for it to be included

		 **Example:**
		   >>> phylo.set_indexer(thresh=7)
		'''
		pt = self.prots[0]
		seq_counts = { key: 1 for key in pt.indexer }
		for pt in self.prots[1:]: # I wish I knew how to do this in a flat way
			for key in pt.indexer:
				if key in seq_counts: seq_counts[key]+=1
				else: seq_counts[key]=1
		countsdf = pd.DataFrame.from_dict(seq_counts,orient='index')
		keep = countsdf.index[countsdf.values.flat>=thresh]
		print 'Keeping: '+str(len(keep))+' species'
		self.indexer = pd.DataFrame(range(len(keep)),index=keep)

	def set_sim_mat(self):
		'''
		Calculate phylogenetic similarity over all included species. Called after set_indexer

		 **Arguments:**
		   -  Bound to an instance of PhyloSet

		 **Example:**
		   >>> phylo.set_sim_mat()
		'''
		Nprot,Mtot = len(self.prots),len(self.indexer.index)
		phylo_sim_mat = np.zeros((Mtot,Mtot,Nprot))
		for j,pt in enumerate(self.prots):
			bin_msa = []
			matches = []
			for key in it.ifilter(lambda key: key in pt.indexer.keys(), self.indexer.index):
				temp = binarize_msa(pt.filter_loci(pt.indexer[key]),Naa=21)
				bin_msa.append(temp.mean(axis=0)) # average sequence across paralogs
				matches.append(key)
			bin_msa = np.vstack(bin_msa)
			sim_mat = bin_msa.dot(bin_msa.T)/float(pt.Npos)
			ix = self.indexer.filter(items=matches,axis=0).values
			phylo_sim_mat[ix,ix.T,j] = sim_mat
		phylo_sim_mat = np.divide(phylo_sim_mat.sum(axis=2),(phylo_sim_mat>0).sum(axis=2))
		np.fill_diagonal(phylo_sim_mat,1.)
		self.sim_mat = phylo_sim_mat

	@property
	def index(self):
		'''
		Returns the list of modeled taxonomic ids

			 **Example:**
			   >>> tids = phylo.index
		'''
		return list(self.indexer.index)

def match_loci(indexer1, indexer2, overlap):
	'''
	This function produces two row matched lists of loci based on taxids. A 1:1 matching is assumed among paralogs, which are paired based on the proximity of their gene locus.

		 **Arguments:**
		   -  `indexer1`  =  An indexer dict from a Prot object. Keys are taxids and the elements are lists of the corresponding loci for the protein family.
		   -  `indexer2`  =  A second indexer dict to be matched against the first.
		   -  `overlap`  =  A collection of taxids to filter by. Usually provided by the phylogenetic model.

		 **Example:**
		   >>> loci1, loci2, taxids = match_loci(indexer1, indexer2, overlap)
	'''
	taxids, loci1, loci2 = [],[],[]
	for key in overlap:
		loci_nums = [len(indexer1[key]),len(indexer2[key])]
		if loci_nums==2:
			loci1+=indexer1[key]
			loci2+=indexer2[key]
			taxids.append(key)
		else:
			# distance matrix between locus numbers
			a = [ float(filter(lambda x: x,[ filter(str.isdigit,frag) for frag in locus.split('_') ])[-1]) for locus in indexer1[key] ]
			b = [ float(filter(lambda x: x,[ filter(str.isdigit,frag) for frag in locus.split('_') ])[-1]) for locus in indexer2[key] ]
			distmat = -abs(np.subtract.outer(a,b))
			for count in range(min(loci_nums)): # limit to 1:1 matches
				i,j = np.unravel_index(distmat.argmax(),distmat.shape)
				# remove the selected row/column
				distmat[i,:] = np.NINF
				distmat[:,j] = np.NINF
				loci1.append(indexer1[key][i])
				loci2.append(indexer2[key][j])
				taxids.append(key)
	return loci1, loci2, taxids

def spec_sim_mat(msa_mat):
	'''
	Computes a species similarity matrix for the provided multiple sequence alignment. For simplicity and to preserve linearity with respect to positions, gaps are considered to be the 21st amino acid.

		 **Arguments:**
		   -  `msa_mat`  =  Numerical multiple sequence alignment, code doesn't matter

		 **Example:**
		   >>> sim_mat  =  spec_sim_mat(msa)
	'''
	if len(msa_mat.shape) == 1: # binarize_msa requires msa_mat to be 2-dimensional
		msa_mat = np.array([msa_mat]).T
	msa_bin = binarize_msa(msa_mat,Naa=21)
	Npos = msa_mat.shape[1]
	return msa_bin.dot(msa_bin.T)/float(Npos)

def specw_calc(sim_mat, delta=0.1): # no clue what delta should be
	'''
	Compute species weights based on a cutoff detla. Species with similarity greater than (1-delta) are combined into a single effective species.

		 **Arguments:**
		   -  `sim_mat`  =  Species similarity matrix

		 **Keyword Arguments:**
		   -  `delta`  =  Cutoff for sequence similarity

		 **Example:**
		   >>> specw = specw_calc(sim_mat, delta=0.1)
	'''
	return 1./(sim_mat>(1.-delta)).sum(axis=0)

##	Mirror tree

def pos_mirror_tree(ps, phylo, delta=0.1, thread_lim = 16):
	'''
	Compute positional mirror-tree for all proteins and protein pairs. Uses parallel computing via pool. Returns results as one list for Prots and one for ProtPairs.

		 **Arguments:**
		   -  `ps`  =  ProtSet object
		   -  `phylo`  =  PhyloSet object

		 **Keyword Arguments:**
		   -  `thread_lim`  =  Maximum number of simultaneous threads (I think its 48 on a 256GB WebGUI)

		 **Example:**
		   >>> protmats, pairmats = pos_mirror_tree(ps, phylo1, thread_lim=48)
	'''
	Nprot = len(ps.prots)
	prot_inputs = []
	full_sigma = lambda msa,phylo_vec,pwsw: wdev(phylo_correct(triu_flatten(spec_sim_mat(msa)), phylo_vec, pwsw),pwsw)
	for pt in ps.prots: # should probably move these inside their own functions
		msa,tids = pt.filter_aln(phylo.index)
		phylo_vec,specw,pwsw,Meff = phylo.info(tids,delta)
		#pt.specw = specw
		prot_inputs.append([pt.name,msa,phylo_vec,pwsw,Meff,full_sigma(msa,phylo_vec,pwsw)**2])
	pair_inputs = []
	for pair in ps.pairs:
		msa1,msa2,tids = pair.filter_alns(phylo.index)
		phylo_vec, specw, pwsw, Meff = phylo.info(tids,delta)
		#pair.specw = specw
		norm_factor = full_sigma(msa1,phylo_vec,pwsw)*full_sigma(msa2,phylo_vec,pwsw)
		pair_inputs.append([pair.names,msa1,msa2,phylo_vec,pwsw,Meff,norm_factor])
	pool = Pool(thread_lim)
	prot_results = pool.map_async(pt_pos_mt,prot_inputs)
	# can I free up memory after initiating these?
	pair_results = pool.map_async(pair_pos_mt,pair_inputs)
	protmats = prot_results.get()
	pairmats = pair_results.get()
	pool.terminate()
	return protmats,pairmats

def pt_pos_mt(args):
	'''
	Calculates positional mirror-tree within a protein. Arguments are provided as a list for the sake of parallelizing using pool.

		 **Arguments:**
		   -  `args`  =   (list) msa, phylo_vec, pwsw, Meff, norm_factor
		   -  `msa`  =   Multiple sequence alignment of the protein in numeric form
		   -  `phylo_vec`  =  A corresponding phylogeny vector with norm = 1 and the same species dimensions
		   -  `pwsw`  =  Corresponding species pair weights
		   -  `Meff`  =  Effective number of species pairs (pwsw.sum())
		   -  `norm_factor`  =  Product of the two relevant full length standard deviations

		 **Example:**
		   >>> pos_mt1 = pt_pos_mt(args)
	'''
        name, msa, phylo_vec, pwsw, Meff, norm_factor = args
        Npos = msa.shape[1]
        sim_list = [ triu_flatten(spec_sim_mat(msa[:,k])) for k in range(Npos) ]
        pos_mt_mat = np.zeros((Npos,Npos))
        for (i,sim1),(j,sim2) in it.combinations(enumerate(sim_list), 2):
                pos_mt_mat[i,j] = pos_mt_calc([sim1,sim2],phylo_vec,pwsw,Meff,norm_factor)
        print('Completed calculations for prot: '+ name)
        return symmetrize(pos_mt_mat)

def pair_pos_mt(args):
	'''
	Calculates positional mirror-tree between a pair of proteins. Arguments are provided as a list for the sake of parallelizing using pool.

		 **Arguments:**
		   -  `args`  =   (list) msa1, msa2, phylo_vec, pwsw, Meff, norm_factor
		   -  `msa1`  =   Multiple sequence alignment of protein 1 in numeric form
		   -  `msa2`  =   Multiple sequence alignment of protein 2 in numeric form
		   -  `phylo_vec`  =  A corresponding phylogeny vector with norm = 1 and the same species dimensions
		   -  `pwsw`  =  Corresponding species pair weights
		   -  `Meff`  =  Effective number of species pairs (pwsw.sum())
		   -  `norm_factor`  =  Product of the two relevant full length standard deviations

		 **Example:**
		   >>> pos_mt1 = pair_pos_mt(args)
	'''
        names, msa1, msa2, phylo_vec, pwsw, Meff, norm_factor = args
        Npos1,Npos2 = msa1.shape[1], msa2.shape[1]
        sim_list1 = [ triu_flatten(spec_sim_mat(msa1[:,k])) for k in range(Npos1) ]
        sim_list2 = [ triu_flatten(spec_sim_mat(msa2[:,k])) for k in range(Npos2) ]
        pos_mt_mat = np.zeros((Npos1,Npos2))
        for (i,sim1),(j,sim2) in it.product(*map(enumerate,[sim_list1,sim_list2])):
                pos_mt_mat[i,j] = pos_mt_calc([sim1,sim2],phylo_vec,pwsw,Meff,norm_factor)
        print('Completed calculations for prots: '+names[0]+', '+names[1])
        return pos_mt_mat

def pos_mt_calc(sims, phylo_vec, pwsw, Meff, norm_factor):
	'''
	Compute the positional mirror tree score. Only differs from the full length score in that it divides through by the full length standard deviations rather than positional dev.

		 **Arguments:**
		   -  `sims`  =  A list of two vectorized similarity matrices
		   -  `phylo_vec`  =  A corresponding phylogeny vector with norm = 1 and the same species dimensions
		   -  `pwsw`  =  Corresponding species pair weights
		   -  `Meff`  =  Effective number of species pairs (pwsw.sum())
		   -  `norm_factor`  =  Product of the two relevant full length standard deviations

		 **Example:**
		   >>> pos_score = pos_mt_calc(sims, phylo_vec, pwsw, Meff, norm_factor)
	'''
	eps = [ phylo_correct(sim, phylo_vec, pwsw) for sim in sims ]
	eps_centered = [ wcenter(ep,pwsw) for ep in eps ]
	return wdot(*eps_centered, w=pwsw)/(Meff*norm_factor)

def mirror_tree(ps, phylo, delta=0.1):
	'''
	Compute the full length mirror tree matrix over a set of proteins given a phylogenetic model.

		 **Arguments:**
		   -  `ps`  =  A ProtSet object
		   -  `phylo`  =  A PhyloSet object
                   -  `delta` = distance cutoff for sequence weights

		 **Example:**
		   >>> mt_mat = mirror_tree(ps, phylo)
	'''
	Nprot = len(ps.prots)
	mt_mat = np.zeros((Nprot,Nprot))
	for pair in ps.pairs:
		msa1,msa2,tids = pair.filter_alns(phylo.index)
		sims = [ triu_flatten(spec_sim_mat(msa)) for msa in [msa1,msa2] ]
		phylo_vec, specw, pwsw, Meff = phylo.info(tids,delta)
		#pair.specw = specw # should I do this?
		mt_mat[pair.ij] = full_mt_calc(sims, phylo_vec, pwsw, Meff)
	return symmetrize(mt_mat)


def scamat_calc(msa_num, specw, posw, lbda=0.03):
        '''
        A function for computing the SCA matrix.

        **Arguments:**
        -  `msa_num` = An msa, converted to numerical encoding
        -  `specw` = species weights
        -  `posw` = position weights
        -  `lbda` = lambda for setting pseudocounting frequency

        **Example:**
        >>> scamat = scamt_calc(msa,specw,posw)

        '''
        freq1,freq2 = freq(msa_num, specw, lbda=lbda)[:2]
        cmat = freq2 - np.outer(freq1,freq1)
        wmat = np.outer(posw,posw)
        tildeC = np.multiply(cmat,wmat)
        Naa = 20; Npos = msa_num.shape[1]
        Cfrob = np.zeros((Npos,Npos))
        for i,j in it.combinations_with_replacement(range(Npos),2):
            Cfrob[i,j] = np.sqrt(np.square(tildeC[Naa*i:Naa*(i+1), Naa*j:Naa*(j+1)]).sum())
        Cfrob += np.triu(Cfrob,1).T
        return Cfrob

def sca(ps, phylo, delta=0.1):
        '''
        Compute sca matrices over a set of proteins.
        Even though we aren't really using the phylogenetic model here, we include one so that the
        sequence filtering and species weights are identical to mirror tree.

        Return a list of sca matrices for single proteins (protSCAmats), and protein pairs (pairSCAmats).
        Also computes an "interaction score" which is the average over all positional correlations
        (more specifically SCA couplings) between proteins

        **Arguments:**
        -  `ps` = A ProtSet object
        -  `phylo`  =  A PhyloSet object
        -  `delta` = distance cutoff for sequence weights

        **Example:**
        >>> protSCAmats, pairSCAmats, scaScoreMat = sca(ps, phylo1, delta=0.1)

        '''
        Nprot = len(ps.prots)
        sca_intrxn_score = np.zeros((Nprot,Nprot))
        sca_intrxn_score2 = np.zeros((Nprot,Nprot))
        protSCAmats = []
        pairSCAmats = []
        pairSCArandmats = []

        for i,pt in enumerate(ps.prots):
                msa,tids = pt.filter_aln(phylo.index)
                phylo_vec,specw,pwsw,Meff = phylo.info(tids,delta)
                posw = posw_calc(msa,specw)
                pt.posw = posw

                #recalculate the species weights for just this alignment to be
                #more consistent with the original SCA. specw above are computed
                #via the phylogenetic model.
                X2d = binarize_msa(msa,Naa=21)
                simMat = (X2d.dot(X2d.T))/msa.shape[1]
                seqw = np.array(1./(simMat>(1-delta)).sum(axis=0))

                protSCAmats.append(scamat_calc(msa,seqw,posw))
                sca_intrxn_score[i,i] = np.mean(protSCAmats[-1])
                #print('Completed calcs for: '+pt.name)

        for pair in ps.pairs:
                msa1, msa2, tids = pair.filter_alns(phylo.index)
                phylo_vec,specw,pwsw,Meff = phylo.info(tids,delta)
                cat = np.concatenate([msa1,msa2],axis=1)
                posw = np.hstack([pt.posw for pt in pair.prots])

                #recalculate the species weights for just this alignment to be
                #more consistent with the original SCA. specw above are computed
                #via the phylogenetic model.
                X2d = binarize_msa(cat,Naa=21)
                simMat = (X2d.dot(X2d.T))/cat.shape[1]
                seqw = np.array(1./(simMat>(1-delta)).sum(axis=0))

                sca_cat = scamat_calc(cat, seqw, posw)
                sca_cat = sca_cat[0:msa1.shape[1],msa1.shape[1]:]
                pairSCAmats.append(sca_cat)
                #print('Completed calcs for: '+pair.names[0]+', '+pair.names[1])

                #also compute a randomized matrix
                randSample = random.sample(range(msa2.shape[0]),msa2.shape[0])
                m2shuffle= np.zeros(msa2.shape).astype(int)
                for i,randN in enumerate(randSample):
                        m2shuffle[i,:] = msa2[randN,:]

                mcombRand = np.concatenate((msa1,m2shuffle),axis=1)
                X2d = binarize_msa(mcombRand,Naa=21)
                simMat = (X2d.dot(X2d.T))/mcombRand.shape[1]
                seqw = np.array(1./(simMat>(1-delta)).sum(axis=0))

                sca_cat_rand = scamat_calc(mcombRand, seqw, posw)
                sca_cat_rand = sca_cat_rand[0:msa1.shape[1],msa1.shape[1]:]
                pairSCArandmats.append(sca_cat-sca_cat_rand)

                #compute interaction
                sca_intrxn_score[pair.ij[0],pair.ij[1]] = np.mean(sca_cat)
                sca_intrxn_score2[pair.ij[0],pair.ij[1]] = np.mean(sca_cat-sca_cat_rand)

        return protSCAmats,pairSCAmats,pairSCArandmats,symmetrize(sca_intrxn_score),symmetrize(sca_intrxn_score2)



def full_mt_calc(sims, phylo_vec, pwsw, Meff):
	'''
	Computes the mirror tree score (partial correlation) given a pair of similarity vectors and corresponding phylogenetic information (which can be easily obtained through the PhyloSet bound method .info())

		 **Arguments:**
		   -  `sims`  =  A list of two vectorized similarity matrices
		   -  `phylo_vec`  =  A corresponding phylogeny vector with norm = 1 and the same dimensions
		   -  `pwsw`  =  Corresponding species pair weights
		   -  `Meff`  =  Effective number of species pairs (pwsw.sum())

		 **Example:**
		   >>> r = full_mt_calc(sims, phylo_vec, pwsw, Meff)
	'''
	eps = [ phylo_correct(sim, phylo_vec, pwsw) for sim in sims ]
	scored = [ wscore(ep,pwsw) for ep in eps ]
	return wdot(*scored, w=pwsw)/Meff

def phylo_correct(sim, phylo_vec, pwsw):
	'''
	Computes the phylogenetically orthogonal component of a similarity matrix using the dot pro

		 **Arguments:**
		   -  `sim`  =  Similarity vector where the elements are pairs of species
		   -  `phylo_vec`  =  A corresponding phylogeny vector with norm = 1 and the same dimensions
		   -  `pwsw`  =  Corresponding species pair weights

		 **Example:**
		   >>> eps = phylo_correct(sim, phylo_vec, pwsw)
	'''
	return sim - phylo_vec*wdot(sim, phylo_vec, pwsw)

def wdot(x, y, w=None):
	'''
	Weighted dot product

		 **Arguments:**
		   -  `x`  =  Numpy array
		   -  `y`  =  Numpy array

		 **Keyword Arguments**
		   -  `w`  =  Weights in a numpy array

		 **Example:**
		   >>> proj = wdot(sim, phylo_vec, w=pwsw)
	'''
	if w is None:
		w = np.ones(x.shape)
	return float((w*x).dot(y.T))

def wcenter(x, w=None):
	'''
	Weighted centring a vector around its mean. Affine translation of the values such that the mean is zero. Useful for calculating a weighted standard score.

		 **Arguments:**
		   -  `x` = numpy array

		 **Keyword Arguments**
		   -  `w` = weights in a numpy array

		 **Example:**
		   >>> eps_centered = wcenter(eps, w=pwsw)
	'''
	return x - np.average(x,weights=w)

def wdev(x, w=None):
	'''
	Calculates a weighted standard deviation

		 **Arguments:**
		   -  `msas`  =  A 1D numpy array

		 **Keyword Arguments:**
		   -  `w`  =  Weights for the dimension of x

		 **Example:**
		   >>> eps_dev = wscore(eps, w=specw)
	'''
	return np.sqrt((w*(wcenter(x)**2)).sum()/(w.sum()-1))

def wscore(x, w=None):
	'''
	Calculates a weighted standard score (sets the mean to zero and divides by the standard deviation).

		 **Arguments:**
		   -  `msas`  =  A 1D numpy array

		 **Keyword Arguments:**
		   -  `w`  =  Weights for the dimension of x

		 **Example:**
		   >>> eps_scored = wscore(eps, w=specw)
	'''
	ctr = wcenter(x)
	return ctr/np.sqrt((w*(ctr**2)).sum()/(w.sum()-1))

## Reading in alignments

def prots_from_scratch(names, path2alg='eggNOG/'):
	'''
	Create a list of Prot objects based on the .fasta files specified by names and path2alg. Filters out overly gapped positions (>20%) and overly gapped sequences (>30%).

		 **Arguments:**
		   -  `names`  =  List of names. These strings must be at the beginning of the filenames of the corresponding alignments, separated by an underscore.

		 **Keyword Arguments:**
		   -  `path2alg`  =  Path to the directory containing the .fasta files

		 **Example:**
		   >>> prots = prots_from_scratch(['trpA','trpB'], path2alg='alignments/eggNOG/')
	'''
	msas, indexers, speciesmaps = read_algs(names, path=path2alg)
	msas, indexers, speciesmaps = filter_seqs(msas, indexers, speciesmaps)
	prots = []
	for name,msa,indexer,speciesmap in it.izip(names,msas,indexers,speciesmaps): prots.append(Prot(name,msa,indexer,speciesmap))
	return prots

def unpack_file(path):
	'''
	Unpack a pickle created by create_db, returns the contents in a list of prot objects.

		 **Arguments:**
		   -  `path`  =  Path to the pickle

		 **Example:**
		   >>> prots, names, extras = unpack_file('example_prots.db')
	'''
	with open(path, "rb") as p:
		db = pickle.load(p)
		names, msas, indexers, speciesmaps = db['names'], db['msas'], db['indexers'], db['speciesmaps']
		if 'extras' in db.keys():extras = db['extras']
		else: extras = {}
		prots = []
		for name,msa,indexer,speciesmap in it.izip(names,msas,indexers,speciesmaps): prots.append(Prot(name,msa,indexer,speciesmap))
	return prots, names, extras

def create_db(names, filename, path2alg='eggNOG/'): # somewhat outdated since I mostly use numpy to save
	'''
	Create a .db pickle that contains the alignments and indexing information for the list of proteins specified by `names`

		 **Arguments:**
		   -  `names`  =  List of names. These strings must be at the beginning of the filenames of the corresponding alignments, separated by an underscore.
		   -  `filename`  =  Filename without a file extension

		 **Keyword Arguments**
		   -  `path2alg`  =  Path to the directory containing the .fasta files

		 **Example:**
		   >>> create_db(names, 'example_prots.db', path2alg='alignments/eggNOG/')
	'''
	# i should make a pipeline directly from alignment reading to a .npy file
	msas, indexers, speciesmaps = read_algs(names, path=path2alg)
	msas, indexers, speciesmaps = filter_seqs(msas, indexers, speciesmaps)
	ProtsInfo = {}
	ProtsInfo['names'] = names
	ProtsInfo['msas'] = msas
	ProtsInfo['indexers'] = indexers
	ProtsInfo['speciesmaps'] = speciesmaps
	pickle.dump(ProtsInfo, open(filename+'.db','wb'))

def read_algs(names, path='eggNOG/'):
	'''
	Read eggNOG alignments, return their sequences, a dictionary mapping taxonomic ids to loci called `indexer`, and a pandas DataFrame mapping loci to taxids called speciesmap

		 **Arguments:**
		   -  `names`  =  List of names. These strings must be at the beginning of the filenames of the corresponding alignments, separated by an underscore.

		 **Keyword Arguments**
		   -  `path`  =  Path to the directory containing the .fasta files

		 **Example:**
		   >>> read_algs(names, path='alignments/eggNOG/')
	'''
	msas, indexers, speciesmaps = [],[],[]
	for name in names:
		taxid2loci = dict()
		loci = []
		taxids = []
		seqs = []
		filepath = glob.glob(path+name+'*')[0]
		# sort sequences, taxids, and locus ids
		with open(filepath,'rb') as f:
			for line in f:
				if line[0] is '>':
					temp = line.split('.')
					taxid, locus =  temp[0][1:], temp[1].rstrip('\n')
					taxids.append(taxid)
					loci.append(locus)
					if taxid in taxid2loci.keys(): taxid2loci[taxid].append(locus)
					else: taxid2loci[taxid] = list([locus])
				else: seqs.append(line.rstrip('\n').upper())
		# index alignments
		temp, pos = filter_pos(seqs)
		msas.append(pd.DataFrame(lett2num(temp),index=loci)) # I want to keep this line the same
		indexers.append(taxid2loci)
		speciesmaps.append(pd.DataFrame(taxids,index=loci))
	return msas, indexers, speciesmaps

def filter_seqs(msas, indexers, speciesmaps, thresh=0.3):
	'''
	Filter sequences from a multiple sequence alignment that are composed of more than a threshold percentage of gaps.

		 **Arguments:**
		   -  `msas`  =  (list) Multiple sequence alignments as a list of strings (single letter amino acid code)
		   -  `indexers`  =  (list) Mappings from taxonomic IDs to gene IDs in dictionary form
		   -  `speciesmaps`  =  (list) Mappings from gene IDs to taxonomic IDs in DataFrame form, 1:1 with no duplicate indices

		 **Keyword Arguments:**
		   -  `thresh`  =  Maximum tolerated fraction of gaps

		 **Example:**
		   >>> filtered_algs = filter_seqs(msas, indexers, speciesmaps)
	'''
	thresh = 0.3
	for i,(msa,indexer,speciesmap) in enumerate(it.izip(msas,indexers,speciesmaps)):
		gap_bool = msa.values==0
		seq_bool = gap_bool.mean(axis=1) > thresh
		temp = speciesmap[seq_bool]
		for j in range(temp.shape[0]):
			locus, tid = temp.index[j], temp.values[j][0]
			indexer[tid].remove(locus)
			if not indexer[tid]: indexer.pop(tid)
		msas[i] = msa[~seq_bool]
		indexers[i] = indexer
		speciesmaps[i] = speciesmap[~seq_bool]
	return msas, indexers, speciesmaps

def filter_pos(alg, seqw=[1], max_fracgaps=.2):
	'''
	Code adapted from scaTools: Filter positions from a multiple sequence alignment that are composed of more than a threshold percentage of gaps.

		 **Arguments:**
		   -  `alg`  =  Multiple sequence alignment as a list of strings (single letter amino acid code)

		 **Keyword Arguments:**
		   -  `seqw`  =  Sequence weighting
		   -  `max_fracgaps`  =  Maximum tolerated fraction of gaps

		 **Example:**
		   >>> alg_filtered = filter_pos(alg)
	'''
        Nseq, Npos = len(alg), len(alg[0])
        if len(seqw) == 1: seqw = np.tile(1.0, (1, Nseq))
        # Fraction of gaps, taking into account sequence weights:
        gapsMat = np.array([[int(alg[s][i]=='-') for i in range(Npos)] for s in range(Nseq)])
        seqwn = seqw/seqw.sum()
        gapsperpos = seqwn.dot(gapsMat)[0]
        # Selected positions:
        selpos = [i for i in range(Npos) if gapsperpos[i] < max_fracgaps]
        # Truncation:
        alg_tr = [''.join([alg[s][i] for i in selpos]) for s in range(Nseq)]
        return alg_tr, selpos

def freq(msa_num, specw, Naa=20, lbda=0, freq0=np.ones(20)/21):
	'''
	Compute amino acid frequencies for a given alignment.

	**Arguments:**
		-  `msa_num` = a MxL sequence alignment (converted using lett2num)

	**Keyword Arguments:**
		- `seqw` = a vector of sequence weights (1xM)
		- `Naa` = the number of amino acids
		- `lbda` = lambda parameter for setting the frequency of pseudo-counts (0 for no pseudo counts)
		- `freq0` = expected average frequency of amino acids at all positions

	**Returns:**
		-  `freq1` = the frequencies of amino acids at each position taken independently (Naa*L)
		-  `freq2` = the joint frequencies of amino acids at pairs of positions (freq2, Naa*L * Naa*L)
		-  `freq0` = the average frequency of amino acids at all positions (Naa)

	:Example:
		>>> freq1, freq2, freq0 = freq(msa_num, seqw, lbda=lbda)
	'''
	Nseq, Npos = msa_num.shape
	msa_bin = binarize_msa(msa_num, Naa)
	Meff = specw.sum()
	freq1 = specw.dot(msa_bin)/Meff
	#freq2 = msa_bin.T.dot(np.diag(seqwn)).dot(msa_bin)
	# Background:
	block = np.outer(freq0,freq0)
	freq2_bkg = np.zeros((Npos*Naa, Npos*Naa))
	for i in range(Npos): freq2_bkg[Naa*i:Naa*(i+1),Naa*i:Naa*(i+1)] = block
	# Regularizations:
	freq1_reg = (1-lbda)*freq1 + lbda*np.tile(freq0,Npos)
	#freq2_reg = (1-lbda)*freq2 + lbda*freq2_bkg # doing the calculation all at once actually saves memory
	freq2_reg = (1-lbda)*msa_bin.T.dot(np.diag(specw)).dot(msa_bin)/Meff + lbda*freq2_bkg
	freq0_reg = freq1_reg.reshape(Npos, Naa).mean(axis=0) # Is this actually what I want for doing DCA?
	return freq1_reg, freq2_reg, freq0_reg

def posw_calc(msa_num, specw, lbda=0.03, freq0 = np.array([.073, .025, .050, .061, .042, .072,\
    .023, .053, .064, .089,.023, .043, .052, .040, .052, .073, .056, .063, .013, .033])):
	'''
	Code adapted from scaTools: Calculates the position and amino acid specific conservation weights (KL entropy) from SCA. In contrast to the scaTools version, this only computes Wia and not the other terms.

		 **Arguments:**
		   -  `msa_num`  =  Multiple sequence alignment with a numeric code in a numpy array
		   -  `specw`  =  Array of species weights

		 **Keyword Arguments**
		   -  `lbda`  =  Lambda parameter for setting the frequency of pseudo-counts (0 for no pseudo counts)
		   -  `freq0`  =  Expected frequency of each amino acid

		 **Example:**
		   >>> posw = posw_calc(msa_num, specw)
	'''

	Nseq, Npos = msa_num.shape; Naa = 20
	msa_bin = binarize_msa(msa_num, Naa)
	specwn = specw/specw.sum()
	freq0v = np.tile(freq0,Npos)
	freq1 = specwn.dot(msa_bin)*(1-lbda) + lbda*freq0v
	iok = [i for i in range(Npos*Naa) if (freq1[i]>0 and freq1[i]<1)]
	Wia = np.zeros(Npos*Naa)
	Wia[iok] = abs(np.log((freq1[iok]*(1-freq0v[iok]))/((1-freq1[iok])*freq0v[iok])))
	return Wia

def lett2num(msa_lett, code='ACDEFGHIKLMNPQRSTVWY'):
	'''
	Convert a multiple sequence alignment from a list of amino acid letters to a numeric code.

		 **Arguments:**
		   -  `msa_lett`  =  Multiple sequence alignment as a list of string sequences in single letter amino acid code format

		 **Keyword Arguments**
		   -  `code`  =  The position of each letter (+1) in the string represents the mapping from the letter code to numeric code, gaps map to zero.

		 **Example:**
		   >>> msa_num = lett2num(seqs)
	'''
	[Nseq, Npos] = [len(msa_lett), len(msa_lett[0])]
	msa_num = np.zeros((Nseq, Npos)).astype(int)
	for s, seq in enumerate(msa_lett):
		for i, lett in enumerate(seq):
			if lett in code:
				msa_num[s, i] = code.index(lett)+1
	return msa_num

def num2lett(msa_num, code='-ACDEFGHIKLMNPQRSTVWY'):
	'''
	Convert a numeric multiple sequence alignment to a list of string sequences with a single letter amino acid code.

		 **Arguments:**
		   -  `msa_num`  =  Multiple sequence alignment with a numeric code in a numpy array

		 **Keyword Arguments**
		   -  `code`  =  The position of each letter in the string represents the mapping from numeric code to single letter amino acid code

		 **Example:**
		   >>> seqs = num2lett(msa_num)
	'''
        Nseq = msa_num.shape[0]
        msa_lett = []
        for k, row in enumerate(msa_num):
                seq = ''
                for num in row: seq = seq+code[num]
                msa_lett.append(seq)
        return msa_lett

def binarize_msa(msa_num, Naa=21): # Reshaping 12 million elements is really inefficient
	'''
	Convert a multiple sequence alignment from a numeric format to a binary one in which the positional axis is length Naa x Npos.

		 **Arguments:**
		   -  `msa_num`  =  Multiple sequence alignment with a numeric code (zero is assumed to be a gap)

		 **Keyword Arguments**
		   -  `Naa`  =  Length of the numeric code

		 **Example:**
		   >>> msa_bin = binarize_msa(msa_num)
	'''
	Nseq,Npos = msa_num.shape
	msa_bin = np.zeros((Nseq,Naa*Npos))
	spacer = np.arange(Npos)*Naa
	for k,row in enumerate(msa_num):
		if Naa is 20:
			nonzero = np.where(row>0)[0]
			idx = row[nonzero] + spacer[nonzero] - 1
		else: idx = row + spacer
		msa_bin[k,idx] = 1
	return msa_bin

def symmetrize(mat, eye=False):
	'''
	Create a symmetric matrix from an upper triangular one.

		 **Arguments:**
		   -  `mat`  =  Upper triangle matrix

		 **Keyword Arguments**
		   -  `eye`  =  (boolean) Fill the diagonal with ones

		 **Example:**
		   >>> mt_mat = symmetrize(triu_rmat, eye=True)
	'''
	smat = mat + mat.T
	if eye:
		smat += np.eye(mat.shape[0])
	return smat

def triu_flatten(mat): ## ambiguous if its columns or rows first
	'''
	Returns the a vectorized form of the upper triangle for a provided matrix.

		 **Arguments:**
		   -  matrix

		 **Example:**
		   >>> vec = triu_flatten(mat)
	'''
	idx = np.triu_indices(mat.shape[0],1)
	return np.array(np.squeeze(mat[idx])).flatten()

def subsample_mat(mat, index):
	'''
	Obtain a square subsample of a matrix with new dimensions described by index. Can be used to re-order rows and columns.

		 **Arguments:**
		   -  `mat`  =  Matrix to be subsampled
		   -  `index`  =  A list of indices to subsample

		 **Example:**
		   >>> phylo_sub = subsample_mat(phylo_mat,indices)
	'''
	arr = np.array([index])
	return mat[arr.T,arr]

def mats_avg(mats):
	'''
	Create a matrix by averaging across a list of matrices

		 **Arguments:**
		   -  `mats`  =  A list of matrices with the same dimensions

		 **Example:**
		   >>> avg_mat = mats_avg(list_of_mats)
	'''
	return reduce(lambda x,y: x+y, mats)/float(len(mats))

def makeATS(sequences, refpos, refseq, iref=0, truncate=False):
	'''
	A function directly adapted from scaTools: If specified, truncate the alignment to the structure (assumes MSAsearch_ has already been run
	to identify the reference sequence (iref)) and produce a mapping (ats) between alignment positions and the positions in the reference sequence (refpos).

	**Arguments:**
	   -  `sequences`  =  Sequences
	   -  `refpos`  =  Reference positions
	   -  `refseq`  =  Reference sequence

	:Keyword Arguments:
	   -  `iref`  =  The index of the sequence in the alignment with the highest identity to the reference
	   -  `truncate`  =  (boolean) Truncate the alignment to the structure

	:Example:
	   >>> sequences_trun, ats_new = ce.makeATS(sequences_full, ats_pdb, seq_pdb, i_ref)

	'''
	from Bio import pairwise2
	if truncate == True:
		print("truncating to reference sequence...")
		# Removing gaps:
		pos_ref = [i for i,a in enumerate(refseq) if a != '-']
		seq_ref = ''.join([refseq[i] for i in pos_ref])
		ats_ref = [refpos[i] for i in pos_ref]
		pos_alg = [i for i,a in enumerate(sequences[iref]) if a != '-']
		seq_tr = [''.join([sq[i] for i in pos_alg]) for sq in sequences]
		# Positions to keep in the alignment and pbd sequences
		# (no gap in any of them after co-alignment):
		seqal_ref, seqal_alg, _, _, _ = pairwise2.align.globalms(seq_ref, seq_tr[iref],\
										2, -1, -.5, -.1)[0]
		keep_ref, keep_alg = list(), list()
		j_ref, j_alg = 0, 0
		for i in range(len(seqal_ref)):
			if seqal_ref[i] != '-' and seqal_alg[i] != '-':
				keep_ref.append(j_ref)
				keep_alg.append(j_alg)
			if seqal_ref[i] != '-': j_ref += 1
			if seqal_alg[i] != '-': j_alg += 1
		sequences_out = [''.join([sq[i] for i in keep_alg]) for sq in seq_tr]
		ats_out = [ats_ref[i] for i in keep_ref]
	else:
		tmp = sequences[iref].replace('-','.')
		refseq = refseq.replace('-','');
		seqal_ref, seqal_alg, _, _, _ = pairwise2.align.globalms(refseq, tmp,\
										2, -1, -.5, -.1)[0]
		print ('Len refseq %i, len refpos %i, Len alg seq %i, len pairalg %i, len gloalg %i' %\
                       (len(refseq),len(refpos), len(tmp),len(seqal_alg),len(sequences[0])))
		#print seqal_ref
		#print seqal_alg
		ats_out = list()
		j_ref = 0
		j_pdb = 0
		for i in range(len(seqal_alg)):
			if seqal_alg[i] == '.' and seqal_ref[i] == '-':
				ats_out.insert(j_ref,'-')
				j_ref += 1
			elif seqal_alg[i] != '.' and seqal_alg[i] != '-':
				if seqal_ref[i] != '-':
					ats_out.insert(j_ref,refpos[j_pdb])
					j_ref += 1
					j_pdb += 1
				else:
					ats_out.insert(j_ref,'-')
					j_ref += 1
			elif seqal_alg[i] == '.' and seqal_ref[i] != '-':
				ats_out.insert(j_ref, refpos[j_pdb])
				j_ref += 1
				j_pdb += 1
			elif seqal_alg[i] == '-':
				j_pdb += 1
		sequences_out = sequences
	return sequences_out, ats_out

def _ppair_list(prots):
	'''
	Create a list of all possible ProtPairs from an input list of prots

		 **Arguments:**
		   -  `prots`  =  A list of Prot objects

		 **Example:**
		   >>> pairs, ijs = _ppair_list(prots)
	'''
	pairs = list()
	ijs = [[],[]]
	for (i,p1),(j,p2) in it.combinations(enumerate(prots), 2):
		pairs.append(ProtPair(p1, p2, (i,j)))
		ijs[0].append(i)
		ijs[1].append(j)
	return pairs, ijs
