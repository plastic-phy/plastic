import sys
# sys.path.insert(0, 'D:/University/RA/Celluloid/celluloid-master/')

import numpy as np
from kmodes.kmodes import KModes
from kmodes.util.dissim import matching_dissim

import argparse, os, sys, errno

from sklearn.model_selection import train_test_split
from sklearn.linear_model import Lasso, LogisticRegression
from sklearn.feature_selection import SelectFromModel
from sklearn.preprocessing import StandardScaler

import pandas as pd  

from hyperopt import Trials, STATUS_OK, tpe
from keras.datasets import mnist
from keras.layers.core import Dense, Dropout, Activation
from keras.models import Sequential
from keras.utils import np_utils

from tensorflow.keras.models import Model, Sequential
from tensorflow.keras.layers import Conv2D, MaxPooling2D, UpSampling2D, Dropout, Flatten, Dense, Reshape

from keras.layers import Input,Dense
                   

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



def cluster_and_output(k, matrix, clust_type, inputpath, outdir,num_init2):
    km = KModes(
    n_clusters = k,
    cat_dissim = matching_dissim,
    init = 'huang',
    n_init = num_init2,
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
    np.savetxt('{}_autoencoder-reducedCells-kmodes.matrix'.format(outfile), out_matrix, fmt='%d', delimiter=' ')

    with open('{}_autoencoder-reducedCells-kmodes_clusters.txt'.format(outfile), 'w+') as file_out:
        for cluster in sorted(cluster_groups):
            file_out.write('{0}\t"{1}"\n'.format(
                cluster, ','.join([ str(x+1) for x in cluster_groups[cluster]])
            ))

    with open('{}_autoencoder-reducedCells-kmodes.mutations'.format(outfile), 'w+') as file_out:
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
    output_directory = "D:/University/RA/Celluloid/celluloid-master/data/" + exp_name + "/autoencoder-reducedCells-kmodes/"
    
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
    
    # reduce to 20 features
    encoding_dim = 50
    
    orig_shape = len(check_data2.columns)
    input_df = Input(shape=(orig_shape,))
    encoded = Dense(encoding_dim, activation='relu')(input_df)
    decoded = Dense(orig_shape, activation='sigmoid')(encoded)
    
    # encoder
    autoencoder = Model(input_df, decoded)
    
    # intermediate result
    encoder = Model(input_df, encoded)
    
    autoencoder.compile(optimizer='adadelta', loss='mean_squared_error')
    
    X_train, X_test, y_train, y_test = train_test_split(
    check_data2,input_labels,
    test_size=0.3,
    random_state=0)
    
    autoencoder.fit(X_train, X_train,
                    epochs=20,
                    batch_size=256,
                    shuffle=True,
                    validation_data=(X_test, X_test))
    
    X_train_encoder = encoder.predict(check_data2)
    check_data_final = np.transpose(X_train_encoder)
    scs_matrix_input = check_data_final
    if clust_type == 'cells':
        # print('Clustering cells')
        # cluster_and_output(num_clust, scs_matrix_input, clust_type, input_file_name, output_directory)
        print('not fully implemented yet')
    elif clust_type == 'mutations':
        scs_matrix_input = np.transpose(scs_matrix_input)
        print('Clustering mutations')
        cluster_and_output(num_clust, scs_matrix_input, clust_type, input_file_name, output_directory,num_init)
    elif clust_type == 'both':
        print('not implemented yet')
    else:
        sys.exit('Something very wrong happened.')
