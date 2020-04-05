#!/usr/bin/env python
# coding: utf-8

# In[1]:


from node2vec import Node2Vec
import networkx as nx
#import igraph
import matplotlib as plt
import csv
import pandas as pd
from sklearn.neighbors import kneighbors_graph
import numpy as np
import scanpy as sc
from annoy import AnnoyIndex
from sklearn.manifold import TSNE
import seaborn as sns 
import umap.umap_ as umap
from sklearn import preprocessing
from sklearn.decomposition import PCA
import subprocess
import os
import shutil
import sys

# # DEFAULT PARAMETER SETTING


sc.settings.verbosity = 2
np.set_printoptions(precision=2)
reducer = umap.UMAP(random_state=42)
G_normalized = nx.Graph()
adata = None;

# # READ 10X DATA


def read_data(inp_path,out_path):
    global adata;
    adata = sc.read_10x_mtx(
    inp_path,  # the directory with the `.mtx` file
    var_names='gene_symbols',                  # use gene symbols for the variable names (variables-axis index)
    cache=True)


# # FILTER CELL AND GENES
def compute(pca_n,cN,gN,out_path):
    t=sc.pp.filter_cells(adata,min_counts=3)
    sc.pp.filter_genes(adata,min_cells=3)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata,n_top_genes=1000)

    variable_genes=adata.var['highly_variable']
    variable_gene_list = []
    for ind in variable_genes.index:
        if variable_genes[ind]==False:
            variable_gene_list.append(ind)


# # NORMALIZATION
    sc.pp.normalize_total(adata, target_sum=1, inplace=True)
    data_df = adata.to_df()
    read_data_df = adata.to_df()
    column_names = data_df.columns


# # LOUVAIN CLUSTERING

    sc.pp.neighbors(adata,n_neighbors=25)
    sc.tl.louvain(adata)
    p=np.array(adata.uns['neighbors']['distances'].todense())
    neighbor_df = pd.DataFrame(p,columns=list(read_data_df.index),index=list(read_data_df.index))
    louvain_cluster_df=sc.get.obs_df(adata,keys=["louvain"])
    n_cluster=louvain_cluster_df['louvain'].nunique()
    size_per_cluster = louvain_cluster_df['louvain'].value_counts()


# # SUBSAMPLING

    pl=0.1
    pu=0.9
    k=500
    prob_per_cluster=[]
    for i in range(0,n_cluster):
        prob_per_cluster.append(pl-np.exp(-(size_per_cluster[i]/k))*(pl-pu))
    sample_per_cluster=[]
    for i in range(0,n_cluster):
        sample_per_cluster.append(int(prob_per_cluster[i]*size_per_cluster[i]))
    sampled_df = pd.DataFrame(columns=data_df.columns)
    #print(sample_per_cluster)
    #print(size_per_cluster)
    left_out_cells=[]
    for i in range(0,n_cluster):
        cells=list(louvain_cluster_df.loc[louvain_cluster_df['louvain'] == str(i)].index)
        sample = np.random.choice(cells,sample_per_cluster[i],replace=False)
        remainder = list(set(cells)-set(sample))
        left_out_cells=left_out_cells+remainder
        for smp in sample:
            sampled_df=sampled_df.append(data_df.loc[smp])
    data_df=sampled_df
    column_names = data_df.columns
    #print(len(data_df))


# # TOP PCA GENES

    pca = PCA(n_components=pca_n+200, svd_solver='full')
    pca.fit(data_df)
    rotation_mat = np.transpose(abs(pca.components_))
    rotation_matrix_df = pd.DataFrame(rotation_mat,index=column_names)
    gene_max = rotation_matrix_df.max(axis=1)
    gene_max=gene_max.sort_values(ascending=False)
    pca_genes = list(gene_max[:pca_n].index)
    pca_genes_df = data_df[pca_genes]
    cell_names = list(pca_genes_df.index.values)
    column = list(pca_genes_df.columns)
    pca_genes_df.to_csv(out_path+'/pca_genes_df.csv')


# # VISUALISE SELECTED GENES

# In[23]:


#df=reducer.fit_transform(pca_genes_df)
#df=pd.DataFrame(df)
#ax = sns.scatterplot(x=0, y=1, data=df)


# # TOP VARIABLE GENES INSTEAD OF PCA GENES

    #pca_genes_df = data_df
    #pca_genes_df=pca_genes_df.drop(variable_gene_list,axis=1)
    #cell_names = list(pca_genes_df.index.values)
    #column = list(pca_genes_df.columns)


# # CELL CELL GRAPH BY ANNOY


    f = len(pca_genes_df.columns)
    t = AnnoyIndex(f, 'angular')  # Length of item vector that will be indexed
    for i in range(len(pca_genes_df.index)):
        v  = pca_genes_df.iloc[i]
        t.add_item(i, v)

    t.build(30) # 10 trees
#t.save('test.ann')



    cell_cell_df = pd.DataFrame(index=[0,1,2,3,4,5])
    for i in range(len(pca_genes_df.index)):
               knn_values = list(t.get_nns_by_item(i,7,include_distances=True))
               #cells = knn_values[0]
               cell_cell_df["item."+str(i)]=knn_values[0][1:]
               cell_cell_df["distance."+str(i)] = knn_values[1][1:]
    cell_graph_column = list(cell_cell_df.columns)



# # CELL-GENE CONNECTIONS

    cell_gene_edge_list = []
    for col in column:
    	n_cells = pca_genes_df.nlargest(gN,col)[col].index.values
    
    	for cell in n_cells:
        	temp=[]
        	temp.append(col)
        	#temp.append(cell_names.index(cell)+1)
        	temp.append(cell)
        	temp.append(float(1))
        	cell_gene_edge_list.append(tuple(temp))
        	i=5
        	while(i>=1):
            		distance="distance."+str(cell_names.index(cell))
            		if cell_cell_df.iloc[i][distance]!=999:
                		cell_cell_df.at[i,distance]=999
                		break
            		else:
                		i=i-1
    G_normalized.add_weighted_edges_from(cell_gene_edge_list)

# # CELL -CELL WITHOUT DROPPING CELL CONNECTIONS

    cell_cell_edge_list=[]
    c=0
    g=0
    for col in cell_graph_column:
        if "distance" in col or "Unnamed" in col:
            continue
        else:
            if col=="item":
                distance="distance"
            else:
                distance="distance."+str(col[5:])
            for i in range(1,6):
                temp_edge=[]
                #if cell_cell_df.iloc[i][distance]!=999:
                temp = list(cell_cell_df["distance."+str(int(cell_cell_df.iloc[i][col]))])
                    #if col[5:] not in temp or (col[5:] in temp and col[5:]!=999):
                #temp_edge.append(str(int(col[5:])+1))
                #temp_edge.append(str(int(cell_cell_df.iloc[i][col])))
                temp_edge.append(pca_genes_df.index[int(col[5:])])
                temp_edge.append(pca_genes_df.index[int(cell_cell_df.iloc[i][col])])
                temp_edge.append(float(1))
                    #if len(temp_edge) !=0:
                cell_cell_edge_list.append(tuple(temp_edge))
                c=c+1
                        #print(temp_edge)
                #else:
                #    g=g+1
    G_normalized.add_weighted_edges_from(cell_cell_edge_list)
    print(c,g)


# # CONNECTION CHECK

    cc=0
    cg=0
    gg=0
    for (node1,node2) in G_normalized.edges:
        node1= str(node1)
        node2=str(node2)
        if node1.isnumeric() and node2.isnumeric():
            cc=cc+1
        if (node1.isnumeric() and not node2.isnumeric()) or (not node1.isnumeric() and node2.isnumeric()):
            cg = cg+1
        if not node1.isnumeric() and not node2.isnumeric():
            gg=gg+1
    print(cc,cg,gg)


# # WRITE GRAPH

    nx.write_gexf(G_normalized, out_path+"/graph_normalized_own_cc.gexf")
    nx.write_gexf(G_normalized, "graph_normalized_own_cc.gexf")
    neighbor_df=neighbor_df.drop(list(sampled_df.index),axis=0)
    neighbor_df=neighbor_df.drop(left_out_cells,axis=1)
    neighbor_df.to_csv(out_path+'/neighbor.csv')

    subprocess.call(['java', '-jar', 'openOrd.jar',out_path+'graph_normalized_own_cc.gexf',out_path+'output_normalized_own_cc.csv'])


input_path=""
output_path=""
file=""
arguments=sys.argv[1:]
if len(arguments)>=4 and arguments[0]=="-i" and arguments[2]=="-o" and len(arguments)<=10:
    input_path=arguments[1]
    output_path=arguments[3]
    temp = arguments[4:]
    pca_n = 300
    cN = 7
    gN = 20
    #file=open(output_file,"w")
    if len(temp)>0:
        if temp[0] == 'pca_n':
            pca_n = arguments[5]
        if temp[0] == 'cN':
            cN = temp[1]
        elif len(temp)>2 and temp[2] == 'cN':
            cN = temp[3]
        if temp[0] == 'gN':
            gN = temp[1]
        elif len(temp)>2 and temp[2] == 'gN':
            gN = temp[3]
        elif len(temp)>4 and temp[4] == 'gN':
            gN = temp[5]
    read_data(input_path,output_path);
    compute(pca_n,cN,gN,output_path);
else:
        print ("\nPlease follow the command: <file>.py -i inputpath -o outputpath -pca_n integer -cN integer -gN integer\n")
        print('inputpath : path for read 10x data\n')
        print('outputpath : path where the output is to be stored\n')
        print('pca_n : number of pca components\n');
        print('cN : number of connections per cell\n');
        print('gN : number of connections per gene\n');
        sys.exit()



