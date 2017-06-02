



# get some graphs


from eden.io.gspan import  gspan_to_eden
from itertools import  islice

from eden.graph import _edge_to_vertex_transform as etvt
def get_graphs(dataset_fname= '../toolsdata/bursi.pos.gspan', size=100):
        return  map (etvt, list(islice(gspan_to_eden(dataset_fname),size)))



# fastest way to get a grammar is by using the sampler class.. 
# but here is how to do it manualy


from graphlearn.graphlearn import decompose
from graphlearn.graphlearn import LocalSubstitutableGraphGrammar as LSGG
from graphlearn.utils.ascii import nx_to_ascii 
from graphlearn.utils.neighbors import getallneighbors
# make graphs decomposeable
deco = decompose.Decomposer()
decomps = [deco.make_new_decomposer(g) for g in get_graphs()]

# train a grammar
grammar=LSGG()
grammar.fit(decomps)



# here happens the magic

graphs=get_graphs(size=10)
graph=graphs[1]
graphs= getallneighbors(deco.make_new_decomposer(graph), grammar)


# print system:


for e in graphs:
    print nx_to_ascii(e, xmax=40,ymax=20)







