import pandas as pd
from .utils import Atom, Residue, ActiveSite
from .io import read_active_sites, read_active_site, write_clustering, write_mult_clusterings

# I consider each residue in my active sites to embody 4 different physiochemical properties--
#this means each of the 20 amino acids is now condensed to an alphabet of just 6 letters.
#these properties are: L = large, A= acidic, B= basic, N= neutral, F= aliphatic, P= polar
#By reducing the features of my active sites from 20 to 6, I hope to find trends in composition across active sites.
#Similar composition and similar size would suggest similar ligands (the things that bind to active sites), and thus 
#similar biological function.

#Here's a dictionary of those residues with their corresponding physiochemical definitions:
convert = {'ARG': 'B', 'LYS': 'B', 'ASP': 'A', 'GLN': 'P', 'ASN': 'P', 'GLU': 'A',
'HIS': 'B', 'SER': 'P', 'THR': 'P', 'PRO': 'N', 'TYR': 'P', 'CYS': 'P', 'GLY': 'N', 'ALA': 'F',
'MET': 'N', 'TRP': 'L', 'LEU': 'F', 'VAL': 'A', 'PHE': 'L', 'ILE': 'F'}

#Now I need to extract all of those active sites from the data folder, using those utility functions.
#First call the function that reads all active sites in a folder.

#So first I extract the residues from each active site, make them a list, then add them to some temporary list of all active sites.
#Next I need to change each of those residue names into my own kind.
#The residues are still separated, and I need them in one string, so join them together.
#Finally, I get them into the format I used to start making my similarity matrix. 
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
my_lame_list = [ [x] for x in new_super_temp]
# print(my_lame_list) #if you want


def compute_similarity(site_a, site_b):
    """
    Compute the similarity between two given ActiveSite instances.

    Input: two ActiveSite instances
    Output: the similarity between them (a floating point number)
    def levenshteinDistance(str1, str2):
    m = len(str1)
    n = len(str2)
    lensum = float(m + n)
    d = []           
    for i in range(m+1):
        d.append([i])        
    del d[0][0]    
    for j in range(n+1):
        d[0].append(j)       
    for j in range(1,n+1):
        for i in range(1,m+1):
            if str1[i-1] == str2[j-1]:
                d[i].insert(j,d[i-1][j-1])           
            else:
                minimum = min(d[i-1][j]+1, d[i][j-1]+1, d[i-1][j-1]+2)         
                d[i].insert(j, minimum)
    ldist = d[-1][-1]
    ratio = (lensum - ldist)/lensum
    return {'distance':ldist, 'ratio':ratio}
 
print(levenshteinDistance("kitten","sitting"))   
print(levenshteinDistance("rosettacode","raisethysword"))
    """

    similarity = 0.0

    # Fill in your code here!

    return similarity


def cluster_by_partitioning(active_sites):
    """
    Cluster a given set of ActiveSite instances using a partitioning method.

    Input: a list of ActiveSite instances
    Output: a clustering of ActiveSite instances
            (this is really a list of clusters, each of which is list of
            ActiveSite instances)
    """
    # Fill in your code here!

    return []


def cluster_hierarchically(active_sites):
    """
    Cluster the given set of ActiveSite instances using a hierarchical algorithm.                                                                  #

    Input: a list of ActiveSite instances
    Output: a list of clusterings
            (each clustering is a list of lists of Sequence objects)
    """

    # Fill in your code here!

    return []
