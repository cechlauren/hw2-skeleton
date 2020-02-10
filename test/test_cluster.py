from hw2skeleton import cluster
from hw2skeleton import io
import os

def test_similarity():
    filename_a = os.path.join("data", "276.pdb")
    filename_b = os.path.join("data", "4629.pdb")
    filename_c = os.path.join("data", "82238.pdb")
    filename_d = os.path.join("data", "7674.pdb")

    activesite_c = 82238
    activesite_d = 7674
    activesite_a = io.read_active_site(filename_a)
    activesite_b = io.read_active_site(filename_b)
    
    convert = {'ARG': 'B', 'LYS': 'B', 'ASP': 'A', 'GLN': 'P', 'ASN': 'P', 'GLU': 'A', 'HIS': 'B', 'SER': 'P', 'THR': 'P', 'PRO': 'N', 'TYR': 'P', 'CYS': 'P', 'GLY': 'N', 'ALA': 'F', 'MET': 'N', 'TRP': 'L', 'LEU': 'F', 'VAL': 'A', 'PHE': 'L', 'ILE': 'F'}
    
    My_active_Site = read_active_sites("../data/")
    temp = []
    for site in My_active_Site:
        temp_physio = []
        amino_acid_list = site.residues
        for amino_acid in amino_acid_list:
            temp_physio.append(amino_acid.type)
        temp.append(temp_physio)
    new_temp = []
    for item in temp:
        new_temp_physio = []
        for aa in item:
            new_temp_physio.append(convert[aa])
        new_temp.append(new_temp_physio)
    new_super_temp = []
    for almost in new_temp:
        new_super_temp.append(''.join(almost))
    my_res_list = [ [x] for x in new_super_temp]
    
    counter_list=list(map(lambda x: list(map(lambda y: Counter(y),x)),my_res_list))
    list_of_dfs=list(map(lambda z: pd.DataFrame.from_dict(z),counter_list))
    dfObj = pd.DataFrame()
    a_df = dfObj.append(list_of_dfs).fillna(0)
    
    active_site_index = [46495,23812,41729,91911,82212,15813,85232,20856,3458,13052,82993,32088,64392,29047,42633,53272,97218,69893,96099,82238,50018,68578,55996,93456,33838,71389,43878,57602,58445,39939,81697,17622,81859,6040,57370,22711,52235,46042,34563,78796,9776,83741,7674,91796,63064,24307,63703,8208,62186,28919,56394,81816,34088,37224,40084,82886,23319,42296,93168,42269,17526,38472,83227,29209,83394,10701,94372,93192,98170,85492,73462,38846,91426,27031,28672,46975,25551,91194,18773,37237,25196,61242,56029,42074,70005,35014,19267,1806,57644,26095,64258,84035,45127,34047,14181,57481,29773,54203,36257,4629,26246,37438,98797,81563,65815,38181,63634,23760,49624,39117,24634,94652,7780,73624,3733,73183,42202,32054,50362,276,70919,94719,10814,25878,39299,27312,88042,38031,52954,20326,8304,72058,34958,41719,97612,47023]
    a_df['']= active_site_index
    b_df = a_df.set_index('')
    
    
    b_df["total"]=b_df.sum(axis=1) #cant use a_df sum bc will include index values...
    
    similarity_df = b_df.T.corr()
    

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
