from eden.io.gspan import gspan_to_eden
from itertools import islice
def get_graphs(dataset_fname='../../toolsdata/bursi.pos.gspan', size=100):
    return  list(islice(gspan_to_eden(dataset_fname),size))
import warnings
warnings.filterwarnings('ignore')



from graphlearn.utils import draw
import graphlearn.minor.molecule.transform_cycle as mole
import graphlearn.minor.decompose as decompose
from graphlearn.graphlearn import Sampler as GLS
from eden.graph import Vectorizer


from graphlearn.graphlearn import Sampler as graphlearn_sampler
from graphlearn.learnedlayer import transform
graphs = get_graphs(size=200)
import graphlearn.minor.decompose as decompose

sampler=graphlearn_sampler(
            decomposer=decompose.MinorDecomposer(),
            graphtransformer=transform.GraphMinorTransformer(group_score_threshold=0.4,num_classes=1,debug=False),
            n_samples=5,
     
            batch_size=1,
            n_steps=50,
            n_jobs=1,
            quick_skip_orig_cip=False,
            core_choice_byfrequency=True,
            burnin=0,
            improving_threshold_fraction=0.5,
            select_cip_max_tries=100,
            keep_duplicates=True,
            monitor=True,
            include_seed=True)

sampler.fit(graphs)
