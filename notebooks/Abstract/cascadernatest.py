from graphlearn.minor.rna.rnadecomposer import RnaDecomposer
import graphlearn.minor.rna.infernal as infernal
from graphlearn.minor.rna import forgitransform as forgitransform
from graphlearn.learnedlayer import cascade


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

graphs = get_sequences(size=100)


sampler=infernal.AbstractSampler(
                                #radius_list=[0,1],
                                #thickness_list=[2],
                                #min_cip_count=1,
                                #min_interface_count=2,
                                #graphtransformer= learntrans.GraphMinorTransformer(),
                                graphtransformer= cascade.RNACascade(num_classes=1,debug=False, debug_rna=True),
                                decomposer=RnaDecomposer(output_sequence=True,pre_vectorizer_rm_f=True,calc_contracted_edge_nodes=False),
                                #estimator=estimator
                                #feasibility_checker=feasibility
                                include_seed=False
                               )

sampler.fit(graphs)
graphs = get_sequences(size=5,withoutnames=True)
r= list( sampler.transform(graphs))
