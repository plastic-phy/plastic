import sys
# sys.path.insert(0, 'D:/University/RA/Celluloid/celluloid-master/')

import numpy as np
from kmodes.kmodes import KModes
from kmodes.util.dissim import matching_dissim


from sklearn.model_selection import train_test_split
from sklearn.linear_model import Lasso, LogisticRegression
from sklearn.feature_selection import SelectFromModel
from sklearn.preprocessing import StandardScaler

import pandas as pd       

# the conflict dissimilarity measure
def conflict_dissim(a, b, **_) :
    v = np.vectorize(lambda ai, bi : ai != 2 and bi != 2 and ai != bi)
    return np.sum(v(a,b), axis = 1)

import argparse, os, sys, errno              

# parser = argparse.ArgumentParser(description='K-means clustering')
# parser.add_argument('-f', '--file', type=str, required=True,
#                     help='SCS matrix')
# parser.add_argument('-k', type=int, required=True,
#                     help='K value of celluloid')
# parser.add_argument('-n', type=int, default=10,
#                     help='n_init')
# parser.add_argument('-c', '--cluster', required=True,
#                     choices=['cells', 'mutations', 'both'],
#                     help='Cluster either cells or mutations')
# parser.add_argument('-o', '--outdir', type=str, required=True,
#                     help='output path')

# args = parser.parse_args()



def cluster_and_output(k, matrix, clust_type, inputpath, outdir):
    km = KModes(
    n_clusters = k,
    cat_dissim = conflict_dissim,
    init = 'huang',
    n_init = 10,
    verbose = 1)
    km.fit(matrix)

    from collections import defaultdict
    cluster_groups = defaultdict(list)

    for j in range(matrix.shape[0]):
        cluster_groups[km.labels_[j]].append(j)
    
    tot_rows = 0
    for cluster in cluster_groups:
        tot_rows += len(cluster_groups[cluster])

    filename = os.path.splitext(os.path.basename(inputpath))[0]
    outfile = os.path.join(outdir, filename)

    centroids  = km.cluster_centroids_
    out_matrix = list()
    for ix_c, c in enumerate(centroids):
        if ix_c in cluster_groups:
            x = list(map(int, list(map(round, c))))
            out_matrix.append(x)

    out_matrix = np.transpose(np.array(out_matrix))

    print(out_matrix.shape)
    print(len(cluster_groups))
    np.savetxt('{}_ridge-reducedCells-celluloid.matrix'.format(outfile), out_matrix, fmt='%d', delimiter=' ')

    with open('{}_ridge-reducedCells-celluloid_clusters.txt'.format(outfile), 'w+') as file_out:
        for cluster in sorted(cluster_groups):
            file_out.write('{0}\t"{1}"\n'.format(
                cluster, ','.join([ str(x+1) for x in cluster_groups[cluster]])
            ))

    with open('{}_ridge-reducedCells-celluloid.mutations'.format(outfile), 'w+') as file_out:
        for cluster in sorted(cluster_groups):
            file_out.write('{0}\n'.format(
                ','.join([ str(x+1) for x in cluster_groups[cluster]])
            ))

    print('Done.')

################################

exp_name = "exp3"

for i in range(1,51):
    print("i = ",i)
    file_num = str(i)
    # datafile = "data_for_lsh.csv"
    datafile = "sim_" + file_num + "_scs.txt"
    
    # output_directory = "D:/University/RA/Celluloid/celluloid-master/data/" + exp_name + "/reduced/JLL/kmodes/"
    # input_file_name = "D:/University/RA/Celluloid/celluloid-master/data/" + exp_name + "/JLL/" + datafile
    input_file_name = "D:/University/RA/Celluloid/celluloid-master/data/" + exp_name + "/ground/" + datafile
    output_directory = "D:/University/RA/Celluloid/celluloid-master/data/" + exp_name + "/ridge-reducedCells-celluloid/"
    
    clust_type = "mutations"
    num_clust = 100
    num_init = 10
        
    scs_matrix_input = np.loadtxt(input_file_name, dtype='d', delimiter=' ')
    check_data = np.transpose(scs_matrix_input)

    file_num = str(i)
    datafile = "sim_" + file_num + "_log_cluster_mutations.txt"
    label_file_name = "D:/University/RA/Celluloid/celluloid-master/data/" + exp_name + "/ground/" + datafile
    
    lines = []
    with open(label_file_name) as f:
        lines = f.readlines()
    
    df_0 = [ [ 0 for i in range(1000) ] for j in range(1) ]
    aa = df_0[0]
    
    for clust_num in range(20):
        clust_lab = clust_num
        str_temp = lines[clust_lab]
        new_string = str_temp.replace(str(clust_lab) + "\t", "")
        new_string = new_string.replace("\n", "")
        new_string = new_string.replace("\"", "")
        new_string_2 = new_string.split(",")
    
        for clust_ind in range(len(new_string_2)):
            val = int(new_string_2[clust_ind])
            aa[(val - 1)] = clust_lab
    
    
    input_labels = aa    
    check_data2 = pd.DataFrame(check_data)
    scaler = StandardScaler()
    scaler.fit(check_data2.fillna(0))

    # L1 = Lasso, L2 = Ridge
    sel_ = SelectFromModel(LogisticRegression(C=1, penalty='l2', solver='liblinear'))
    sel_.fit(scaler.transform(check_data2.fillna(0)), input_labels)
    
    sel_.get_support()
    
    selected_feat = check_data2.columns[(sel_.get_support())]
    print('total features: {}'.format((check_data2.shape[1])))
    print('selected features: {}'.format(len(selected_feat)))
    print('features with coefficients shrank to zero: {}'.format(
          np.sum(sel_.estimator_.coef_ == 0)))
    
    X_train_selected = sel_.transform(check_data2.fillna(0))
    
    check_data_final = np.transpose(X_train_selected)
    scs_matrix_input = check_data_final
    if clust_type == 'cells':
        # print('Clustering cells')
        # cluster_and_output(num_clust, scs_matrix_input, clust_type, input_file_name, output_directory)
        print('not fully implemented yet')
    elif clust_type == 'mutations':
        scs_matrix_input = np.transpose(scs_matrix_input)
        print('Clustering mutations')
        cluster_and_output(num_clust, scs_matrix_input, clust_type, input_file_name, output_directory)
    elif clust_type == 'both':
        print('not implemented yet')
    else:
        sys.exit('Something very wrong happened.')
