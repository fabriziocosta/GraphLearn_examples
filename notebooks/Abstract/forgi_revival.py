#from eden.converter.fasta import fasta_to_sequence
from eden_rna.io.fasta import load 
import itertools


def rfam_uri(family_id):
    return 'http://rfam.xfam.org/family/%s/alignment?acc=%s&format=fastau&download=0'%(family_id,family_id)
def rfam_uri(family_id):
    return '%s.fa'%(family_id)

def get_sequences(size=9999,withoutnames=False):
    sequences = itertools.islice( load("../../toolsdata/RF00005.fa"), size)
    if withoutnames:
        return [ b for (a,b) in sequences ]
    return sequences


import graphlearn.minor.rna.infernal as infernal
from  graphlearn.feasibility import FeasibilityChecker as Checker
import graphlearn.estimate as estimate
from graphlearn.minor.rna import forgitransform as forgitransform
#from graphlearn.minor.decompose import MinorDecomposer
from graphlearn.minor.rna.rnadecomposer import RnaDecomposer
# not really needed since after refolding we get an RNA
#feasibility=Checker()
#feasibility.checklist.append(rna.is_rna)
graphs = get_sequences(size=100)

estimator=estimate.OneClassEstimator( nu=.33, cv=2, n_jobs=-1)

sampler=infernal.AbstractSampler(
                            #radius_list=[0,1],
                            #thickness_list=[2],
                            #min_cip_count=1,
                            #min_interface_count=2,
                            graphtransformer=forgitransform.GraphTransformerForgi(),
                            decomposer=RnaDecomposer(output_sequence=True,
                                                     pre_vectorizer_rm_f=True)
                            #estimator=estimator
                            #feasibility_checker=feasibility
                           )

sampler.fit(graphs)
graphs = get_sequences(size=5,withoutnames=True)
r= list( sampler.transform(graphs))


