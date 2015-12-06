from eden.util import configure_logging
import logging
configure_logging(logging.getLogger(),verbosity=1)


'''
GET RNA DATA
'''
from eden.converter.fasta import fasta_to_sequence
import itertools
from eden.util import random_bipartition_iter
import random

def rfam_uri(family_id):
    return 'http://rfam.xfam.org/family/%s/alignment?acc=%s&format=fastau&download=0'%(family_id,family_id)
def rfam_uri(family_id):
    return '%s.fa'%(family_id)
 
    
    
RFAM="RF01725"
#cutoff 162 (44.0)
#cutoff 1725 (38.0)
#cutoff rest (29)
sizes=[50,100,200,400]
repeats=3




def get_sequences(size=9999,rand=False):
    sequences = get_sequences_with_names(size=size,rand=rand)
    return [ b for (a,b) in sequences ]

def get_sequences_with_names(size=9999, rand=0):
    if rand>0:
        sequences , boring = random_bipartition_iter(fasta_to_sequence("../toolsdata/%s.fa" % RFAM),.9,random_state=random.random()*rand)
        sequences = itertools.islice( sequences , size)
    else:
        sequences = itertools.islice( fasta_to_sequence("../toolsdata/%s.fa" % RFAM), size)
    return sequences





import random
import graphlearn.abstract_graphs.RNA as rna
from  graphlearn.feasibility import FeasibilityChecker as Checker
from graphlearn.estimator import Wrapper as estimatorwrapper
import graphlearn.utils.draw as draw
from graphlearn.graphlearn import Sampler as GLS
import itertools



def fit_sample(graphs, random_state=random.random()):
    '''
    graphs -> more graphs
    '''
    graphs = list(graphs)
    estimator=estimatorwrapper( nu=.33, cv=2, n_jobs=-1)
    sampler=rna.AbstractSampler(radius_list=[0,1],
                                thickness_list=[2], 
                                min_cip_count=1, 
                                min_interface_count=2, 
                                preprocessor=rna.PreProcessor(base_thickness_list=[1],ignore_inserts=True), 
                                postprocessor=rna.PostProcessor(),
                                estimator=estimator
                                #feasibility_checker=feasibility
                               )
    sampler.fit(graphs,grammar_n_jobs=4,grammar_batch_size=1)
    


    
    
    
    #logger.info('graph grammar stats:')
    dataset_size, interface_counts, core_counts, cip_counts = sampler.grammar().size()
    #logger.info('#instances:%d   #interfaces: %d   #cores: %d   #core-interface-pairs: %d' % (dataset_size, interface_counts, core_counts, cip_counts))
    
    
    graphs = [ b for a ,b in graphs  ]
    
    graphs = sampler.sample(graphs,
                            n_samples=3,
                            batch_size=1,
                            n_steps=50,
                            n_jobs=4,
                            quick_skip_orig_cip=True,
                            probabilistic_core_choice=True,
                            burnin=10,
                            improving_threshold=0.9,
                            improving_linear_start=0.3,
                            max_size_diff=20,
                            accept_min_similarity=0.65,
                            select_cip_max_tries=30,
                            keep_duplicates=False,
                            include_seed=False,
                            backtrack=10,
                            monitor=False)
    result=[]
    for graphlist in graphs:
        result+=graphlist
    # note that this is a list [('',sequ),..]
    return result
    
    
    
    
    
    
    

def eval(repeats,size):
    result=[]
    for i in range(repeats):
        graphs=get_sequences_with_names(size=size, rand=(i+3)*10)
        zz=fit_sample(graphs)
        z=[b for a ,b in zz]
        cmpath='../%s.cm' % RFAM
        result+=rna.infernal_checker(z,cmfile=cmpath, cmsearchbinarypath='../toolsdata/cmsearch')
        
    a = numpy.array(res)
    mean = numpy.mean(a, axis=0)
    std = numpy.std(a, axis=0)
    
    print 'size:%d mean:%f std:%f' % (size,mean,std)
    return [mean,std]
    


import numpy as np
import matplotlib.pyplot as plt
def make_bar_plot(labels=('G1', 'G2', 'G3', 'G4', 'G5'),means=(20, 35, 30, 35, 27),stds=(2, 3, 4, 1, 2)):
    N = len(labels)
    ind = np.arange(N) 
    width = .5 #0.35
    p1 = plt.bar(ind, means, width, color='r', yerr=stds)
    plt.ylabel('Scores')
    plt.title('Scores by training size')
    plt.xticks(ind + width/2, labels )
    plt.yticks(np.arange(0, 100, 10))
    plt.show()
    



means=[]
stds=[]
for size in sizes:
    m,s=eval(repeats,size)
    means.append(m)
    stds.append(s)
    
    
print 'size: '+str(sizes)
print 'means: '+str(means)
print 'stds: '+ str(stds)

#make_bar_plot(sizes,means,stds)

