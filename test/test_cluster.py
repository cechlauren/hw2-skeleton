from hw2skeleton import cluster
from hw2skeleton import io
import os

def test_similarity():
    filename_a = os.path.join("data", "276.pdb")
    filename_b = os.path.join("data", "4629.pdb")
    filename_c = os.path.join("data", "82238.pdb")
    filename_d = os.path.join("data", "7674.pdb")

    activesite_c = io.read_active_site(filename_c)
    activesite_d = io.read_active_site(filename_d)
    activesite_a = io.read_active_site(filename_a)
    activesite_b = io.read_active_site(filename_b)

    # update this assertion
    assert cluster.compute_similarity(activesite_c, activesite_c) == 1.0
    #assert cluster.compute_similarity(activesite_a, activesite_a) == 1.0

#def test_distance():
   # filename_a = os.path.join("data", "276.pdb")
   # filename_b = os.path.join("data", "4629.pdb")

   # activesite_a = io.read_active_site(filename_a)
   # activesite_b = io.read_active_site(filename_b)

    # update this assertion
   # simulation1 = test distance(activesite_a, activesite_a)
   # simulation2 = test distance(activesite_a, activesite_a)
   # simulation3 = test distance(activesite_a, activesite_b)
   # simulation4 = test distance(activesite_b, activesite_a)
    
   # assert simulation1 == simulation2
   # assert simulation3 == simulation4
    
    
    
def test_partition_clustering():
    # tractable subset
    pdb_ids = [276, 4629, 10701]

    active_sites = []
    for id in pdb_ids:
        filepath = os.path.join("data", "%i.pdb"%id)
        active_sites.append(io.read_active_site(filepath))

    # update this assertion
    assert cluster.cluster_by_partitioning(active_sites) == []

def test_hierarchical_clustering():
    # tractable subset
    pdb_ids = [276, 4629, 10701]

    active_sites = []
    for id in pdb_ids:
        filepath = os.path.join("data", "%i.pdb"%id)
        active_sites.append(io.read_active_site(filepath))

    # update this assertion
    assert cluster.cluster_hierarchically(active_sites) == []
