"""
This whole file here is really not used for simprepper at all.
I fixed the problem before I noticed that, though...
"""
import logging
import numpy as np

# cython-optimized pairwise distance function
# profiled to run in ~50% of the time of pdist
# much less memory hungry. No storage of all distances.
try:
    # NOTE PAQ: This is borrowed from 
    # https://github.com/JoaoRodrigues/openmm_scripts/blob/master/openmm/amberff/setPeriodicBox.py
    """
    # add to PYTHONPATH current workdir and script dir
    #NOTE: PAQ: I don't know why we would do that...?
    # sys.path.insert(0, os.curdir)
    # sys.path.insert(0, os.path.dirname(__file__))
    from _pwdistance import pw_dist
    """
    # PAQ: get from under the hood of openmm
    from simprepper.lib._pwdistance import pw_dist
    # logging.info("\nActually managed to import the cython based pwdist routine.\n")
    print("\nActually managed to import the cython based pwdist routine.\n")
except ImportError:
    logging.warning("\nUsing numpy/scipy (slower) pwdist routine for simulation setup.\n")
    from scipy.spatial.distance import pdist

    def pw_dist(xyz_array):
        return np.amax(pdist(xyz_array, "euclidean"))
