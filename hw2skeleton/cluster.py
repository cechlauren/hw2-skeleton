import pandas as pd
from .utils import Atom, Residue, ActiveSite
from .io import read_active_sites, read_active_site, write_clustering, write_mult_clusterings

###########################################################################################################
#PROVIDE CONTEXT TO MY METHODS:
# I consider each residue in my active sites to embody 4 different physiochemical properties--
#this means each of the 20 amino acids is now condensed to an alphabet of just 6 letters.
#these properties are: L = large, A= acidic, B= basic, N= neutral, F= aliphatic, P= polar
#By reducing the features of my active sites from 20 to 6, I hope to find trends in composition across active sites.
#Similar composition and similar size would suggest similar ligands (the things that bind to active sites), and thus 
#similar biological function.
###########################################################################################################
#MAKE A DICTIONARY:
#Here's a dictionary of those residues with their corresponding physiochemical definitions:
convert = {'ARG': 'B', 'LYS': 'B', 'ASP': 'A', 'GLN': 'P', 'ASN': 'P', 'GLU': 'A',
'HIS': 'B', 'SER': 'P', 'THR': 'P', 'PRO': 'N', 'TYR': 'P', 'CYS': 'P', 'GLY': 'N', 'ALA': 'F',
'MET': 'N', 'TRP': 'L', 'LEU': 'F', 'VAL': 'A', 'PHE': 'L', 'ILE': 'F'}
###########################################################################################################
#PREP MY ACTIVE SITE DATA:
#Now I need to extract all of those active sites from the data folder, using those utility functions.
#First call the function that reads all active sites in a folder.

#I extract the residues from each active site, make them a list, then add them to some temporary list of all active sites.
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
my_res_list = [ [x] for x in new_super_temp]
# print(my_lame_list) #if you want
#now I should have a list of all my active sites' residues, composed of only those 6 specified features.
###########################################################################################################

#MAKE A DATAFRAME OF MY ACTIVESITES X FEATURES:

#now I need to make my dataframe with all my final values of the counts of those residues in my active sites.

#(Miriam helped with this):
#Using map allows us to apply a function to directly to a list.
counter_list=list(map(lambda x: list(map(lambda y: Counter(y),x)),my_res_list)) #This is a counter for each of the features
list_of_dfs=list(map(lambda z: pd.DataFrame.from_dict(z),counter_list)) #I think this makes a data frame for each of the counts
dfObj = pd.DataFrame() #make an empty data frame to add our count data frame to
a_df = dfObj.append(list_of_dfs).fillna(0) #join those two data frames together and fill in all NAs with zero
#This outputs a 136xfeature dataframe with each feature having some number of counts based on the activesite. 

#My activesites are not labeled, so need to indexed them
a_df.index=range(len(my_res_list)) #index them based on length of the list of active sites I put in

#I also included a total count list that can inform me about the size of the active site itself
a_df["total"]=a_df.sum(axis=1)

#Determine the similarity between all given ActiveSite instances, where the input is all active site's features
similarity_df = a_df.T.corr()
#(residues and total) counts, and the output will be the "distance" or dissimilarity between them (since I do 1-correlation).

distance_matrix = 1-a_df.T.corr()

###########################################################################################################


#I'll consider my similarity metric to be the "correlation" between active site instances that I produced above.
#I'll copy and paste all of those matrix construction machinations here for the similarity function.
def compute_similarity(site_a, site_b):
    """
    Compute the similarity between two given ActiveSite instances.

    Input: two ActiveSite instances
    Output: the similarity between them (a floating point number)
  
    """

    
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
    
    a_df.index=range(len(my_res_list))
    a_df["total"]=a_df.sum(axis=1)
    similarity_df = a_df.T.corr()
    
    
    #now I have a matrix containing all of the correlations between activesites. 
    #The higher the score (0:1), the higher the similarity between them. 
    #To get the similarity, find the instance in the similarity matrix when the two activesites intersect.
    similarity = similarity_df[site_a][site_b]
    return similarity

def Euclidean_distance(x,y):
    """
    Given two vectors, calculate the euclidean distance between them.
    Input: Two vectors
    Output: Distance
    
    """
    D = math.sqrt(sum([(a-b)**2 for a,b in zip(x,y)]))
    return D

def cluster_by_partitioning(active_sites):
    """
    Cluster a given set of ActiveSite instances using a partitioning method.

    Input: a list of ActiveSite instances
    Output: a clustering of ActiveSite instances
            (this is really a list of clusters, each of which is list of
            ActiveSite instances)
    """
    # Fill in your code here!
    
    """
    K means summary:
    A good general purpose way to think about discovering groups in data. However, there are several cons. First, 
    one has to specify the number of clusters in advance, or do some selection of clusters post hoc.
    Second, the rules behind belonging to a group are sometimes too simplistic: one point belongs to cluster k if it 
    is closer to the kth center than any other center.
    Third, the solution kmeans finds is dependent on the initialization, which often has some randomized component.
    What happens if our data is more complicated than some partitioning method can describe?
    ##########################################################################################################
    #WHAT A SINGLE ITERATION MIGHT LOOK LIKE:
    Choose number of iterations that might result in convergence
    m=X.shape[0] #number of training examples to use, given the shape of our data frame
    n=X.shape[1] #number of features
    n_iter=100
    
    Choose the number of clusters:
    K=5
    
    #start implementing clustering steps
    #Initialize centroids randomly from the data
    Centroids=np.array([]).reshape(n,0) 
    #will be an n*K dimensional matrix where each column is a centroid for one cluster(5)
    
    for i in range(K):
        rand=rd.randint(0,m-1) #generate some random points to use for each cluster 
        Centroids=np.c_[Centroids,X[rand]] #now we have randomly placed centroids
        
    #compute euc distance from centroid and assign cluster based on minimum distance:
    #initialize our 'output' dictionary: 
    Output={}
    #find euc from each point to the centroids and store it in that matrix
    EuclidianDistance=np.array([]).reshape(m,0) #every row in matrix will have dist of a point to all centroids given
    
    for k in range(K):
       tempDist=np.sum((X-Centroids[:,k])**2,axis=1) #basic definition of euclidean distance, where ,k is ref each centroid location
        EuclidianDistance=np.c_[EuclidianDistance,tempDist] 
        #this will concatenate those temporary distances to our euc distances, which helps us as we (re)assign data to clusters later
    C=np.argmin(EuclidianDistance,axis=1)+1 #find minimum dist and store index of column in C vector
    
  
    #Regroup data points based on cluster index C and store in Output. Compute mean of sep. clusters and assign 
    #new centroids. Y is now temp dictionary that stores solution for some iteration.
    Y={}
    for k in range(K):
        Y[k+1]=np.array([]).reshape(2,0) #this the array where all the solutions will be stored
    for i in range(m):
        Y[C[i]]=np.c_[Y[C[i]],X[i]] #now the minimum euclidean distance for each activesite to the centroid is chosen
    for k in range(K):
        Y[k+1]=Y[k+1].T #transpose bc will need to recalculate euc distance to new centroids
    for k in range(K):
        Centroids[:,k]=np.mean(Y[k+1],axis=0) #gives us the new centroid location based off mean dist. values
    ############################################################################################################
    #WHAT A TRUE KMEANS ALG MIGHT LOOK LIKE:
    #repeat as described prev. until the convergence is achieved, basically looping over n_int and repeating 
    
    for i in range(n_iter):
     
      EuclidianDistance=np.array([]).reshape(m,0)
        for k in range(K):
            tempDist=np.sum((X-Centroids[:,k])**2,axis=1)
            EuclidianDistance=np.c_[EuclidianDistance,tempDist]
        C=np.argmin(EuclidianDistance,axis=1)+1
     
     
      Y={}
        for k in range(K):
            Y[k+1]=np.array([]).reshape(2,0)
        for i in range(m):
            Y[C[i]]=np.c_[Y[C[i]],X[i]]
     
        for k in range(K):
            Y[k+1]=Y[k+1].T
    
        for k in range(K):
            Centroids[:,k]=np.mean(Y[k+1],axis=0)
      Output=Y
      
      ####################################################
      #VISUALIZATION
      #1. scatter original unclustered data
      plt.scatter(X[:,0],X[:,1],c='black',label='unclustered data')
      plt.xlabel('')
      plt.ylabel('')
      plt.legend()
      plt.title('Plot of data points')
      plt.show()
      
      
      #2. now plot clustered data:
      color=['red','blue','green','cyan','magenta']
      labels=['cluster1','cluster2','cluster3','cluster4','cluster5']
      for k in range(K):
        plt.scatter(Output[k+1][:,0],Output[k+1][:,1],c=color[k],label=labels[k])
      plt.scatter(Centroids[0,:],Centroids[1,:],s=300,c='yellow',label='Centroids')
      plt.xlabel('')
      plt.ylabel('')
      plt.legend()
      plt.show()

    """

    return []

##################################################################################################
#SOME BACKGROUND:
#So far, we have taken each active site and reduced it down to a set of four features.
#Each feature has a count for each active site, and each active site has a total sum of all feature counts. 
#Our dataframe then is composed of floating values from [0:1] that inform us about how correlated our active site features are to each other. 
#Then if we do 1-corr, we can get the 'distance' each of our active sites are from each other.

#The end product of hierarchical clustering usually is a binary tree that covers the data. The root of this tree
# is a single cluster that contains all the data, while its leaves would represent individual data points.
#To make the tree, one needs to make clusters of clusters.

#Agglomerative approach:
#Agglomerative clustering starts with each data as its own cluster, then merges groups together. Essentially, 
# an algorithm describing this would take and maintain an active set of clusters that would be decided upon to
# merge or not at each stage of distance calculations. When two clusters are merged, they are each removed from 
# the active set and their union is added to the active set. 
#One iterates this until there is only one cluster in the active set.
#To get the hierarchical tree, one should keep track of which clusters were merged. 
#The main decision one needs to make when using agglomerative hierarchical clustering is the distance criteria
# between groups. While Kmeans does distances between data items, agglomerative will do distances between 
# groups of data items. Common ones are single linkage (minimum between all inter-group pairs), 
# complete linkage (maximum), average, and centroid. 

################################################################################################
#EXAMPLE OF SINGLE LINKAGE:

#def distance_single_linkage({Xn},{Ym}) 
    #min||Xn-Ym||
################################################################################################

def cluster_hierarchically(active_sites):
    """
    Cluster the given set of ActiveSite instances using a hierarchical algorithm.                                                                  #

    Input: a list of ActiveSite instances
    Output: a list of clusterings
            (each clustering is a list of lists of Sequence objects)
    """

    # Fill in your code here!
    """
    #Will use agglomerative clustering
    #pseudocode:
    #Input: (1)some active sites, (2) recalculated similarity matrix between those active sites
    #output:  some clusterings of  activesites
    #while TRUE do
        #a single cluster <-- each original active site
        #find a pair of clusters (a) and (b) such that their similarity s[(a),(b)] = max(s[(m),(n)]);
        #if (s[(a),(b)]) < or = to some threshold, then
            #terminate the while loop;
        #else
            #merge (a) and (b) into a new cluster(a+b);
            #update similarity matrix by deleting both the row and the column corresponding to (a) and (b);
        #end
    #end
    #output current clusters containing similar activesites
    """
    

    return []

################################################################################################
#CLUSTER QUALITY:

#Silhoutte score to check quality of clustering assignment.
#The range of a SS is from [-1,1].
#Closer to +1 means that the sample is far away from its neighboring cluster
#0 means that the sample is on/close to the decision boundary separating two neighboring clusters
#-1 menas the samples were assigned to the wrong clusters


"""
    #Since:
    silhouette_score = (p-q)/max(p,q)
    #where p = mean distance to the points in the nearest cluster
    #where q = mean intra-cluster distance to all the points
    #This indicates that we need some functions that can extract the nearest cluster and intra-
    #cluster distances.
        
"""
################################################################################################
#CLUSTER QUALITY PREP: DEFINE OUR DISTANCES (INTER/INTRA-CLUSTER):

def _intra_cluster_dist_singular(subsetX, metric):
    """
    Must find the pairwise distances for all of the points in a cluster, and then find the average distance for that point
    to the others.
    """
    distances = pairwise_distances(subsetX, metric = metric)
    return distances.sum(axis=1)/ (distances.shape[0]-1)

#Now we can find the mean intra-cluster distance for all of our points in the cluster:
def _intra_cluster_dist(X, labels, metric, n_jobs = 1):
    """
    Calculate mean intra-cluster distance for some sample i where:
    X is an array containing some number of samples with some number of features
    labels is an array that contains the labels for each sample
    metric can be euclidean or your favorite distance metric, or if X is the dist. array itself, 
       just assume 'precomputed.'
    
    The output should be an array containing the intra-cluster distance
    #So this will compute the distance of points within a cluster
    """
    intra_dist = np.zeros(labels.size, dtype=float)
    values = Parallel(n_jobs = n_jobs)(
            delayed(_intra_cluster_dist_singular) #use the singular av distance from above
                (X[np.where(labels == label)[0]], metric)
                for label in np.unique(labels))
    for label, values_ in zip(np.unique(labels), values):
        intra_dist[np.where(labels == label)[0]] = values_
    return intra_dist

#That completes intra-, or within-, cluster distances. Now we need to compare distances between clusters:

def _inter_cluster_dist_single(subsetX_a, subsetX_b, metric):
    """
    Now we find the mean distance for a pair of clusters within all of our clustering assignments.
    Should probably add a case for when the cluster doesn't exist?
    """
    dist = pairwise_distances(subX_a, subX_b, metric=metric)
    dist_a = dist.mean(axis=1)
    dist_b = dist.mean(axis=0)
    return dist_a, dist_b
    
#Now we can use our single intercluster distances function and apply that to all of our clusters.    
def _inter_cluster_dist(X, labels, metric, n_jobs = 1):
    """
    Calc the mean nearest cluster distance for sample i compared to all other clusters where:
    X is an array containing some number of samples with some number of features
    labels is an array that contains the labels for each sample
    metric can be euclidean or your favorite distance metric, or if X is the dist. array itself, 
       just assume 'precomputed.'
       
    The output will be a float, a 'mean-nearest' cluster distance for our sample i
    
    """
    inter_dist = np.empty(labels.size, dtype = float)
    inter_dist.fill(np.inf)
    #will compute the cluster distances between a pair of clusters
    
    some_lab = np.unique(labels)
    
    values = Parallel(n_jobs=n_jobs)(
            delayed(_inter_cluster_dist_single)(
                X[np.where(labels == label_a)[0]],
                X[np.where(labels == label_b)[0]],
                metric)
            for label_a, label_b in combinations(some_lab, 2))
    for (label_a, label_b), (valeus_a, values_b) in \
            zip(combinations(some_lab, 2), values):
            
            indices_a = np.where(labels == label_a)[0]
            inter_dist[indices_a] = np.minimum(values_a, inter_dist[indices_a])
            del indices_a
            indices_a = np.where(labels == label_b)[0]
            inter_dist[indices_b] = np.minimum(values_b, inter_dist[indices_b])
            del indices_b

    return inter_dist

################################################################################################
#USE SILHOUETTE SCORE:

def cluster_quality(X, labels, metric='precomputed', n_jobs=1):
    """
    compute silhouette coefficient for each sample. 
    Models with high Silhouette coefficient are dense because samples in the same cluster are
    similar to each other.
    Returns: silhouette coefficient for each sample, where best is 1 and worst is -1, 0 means
    overlapping clusters
    
    X is a feature array
    labels could be the label values for each sample
    metric is the string that we describe for measuring distances between instances in a feature array
    """
    B = _inter_cluster_dist(X, labels, metric, n_jobs=n_jobs)
    A = _intra_cluster_dist(X, labels, metric, n_jobs =n_jobs)
    
    sil = (B-A)/np.maximum(A, B)
    #some clusters might have a size of 1 when they should be zero so:
    return np.nan_to_num(sil)


################################################################################################
#COMPARE YOUR TWO CLUSTERINGS:

#rand index to check clustering assignment, is agnostic to size/dimensionality
def cluster_compare(cluster_kmeans, cluster_agglomerative):
    """
    for i in cluster:
        for j in cluster:
            #for some cluster with i, is it the same as the cluster of j? and do this for each point
            #if they are the same:
            cluster(i) == cluster(j)?
            if True then TRUE #our clustering was good
            cluster2(i) == cluster2(j)?
            if False then FALSE #our clustering was bad
    """
    return[]
        
