#!/usr/bin/env python
# coding: utf-8

# In[17]:


import json
import matplotlib as plt
import csv
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from decimal import Decimal
import seaborn as sns 
import pandas as pd
import networkx as nx
from sklearn.cluster import DBSCAN
from sklearn.cluster import KMeans
import operator
import numpy as np
import random
import sys

# In[18]:


import seaborn
seaborn_colors=[]
for k in seaborn.xkcd_rgb.keys():
    seaborn_colors.append(k)


# In[19]:

def compute(inp_dataset,input_path,output_path):

    import json
    import matplotlib as plt
    import csv 
    from sklearn.manifold import TSNE
    import matplotlib.pyplot as plt
    from sklearn.decomposition import PCA
    from decimal import Decimal
    import seaborn as sns
    import pandas as pd
    import networkx as nx
    from sklearn.cluster import DBSCAN
    from sklearn.cluster import KMeans
    import operator
    import numpy as np
    import random
    import sys

    

    #csvData=[['data','x','y','type']]
    print("Processing the input data into datafames....")
    csvData=[]
    count=0
    #filename = "G:/Thesis/Dropclust/plots/output_normalized_own_cc.csv"
    #filename = "G:/Thesis/Dropclust/plots/PCA_GENES/output_normalized_own_cc.csv"
    #filename = "G:/Thesis/Dropclust/output_normalized_zscore_cc1.csv"
    #filename = "C:/Users/Swagatam/IdeaProjects/openOrd/output_normalized_own_cc.csv"
    filename = input_path+"/output_normalized_own_cc.csv"
    coord_data = pd.read_csv(filename,names=['data','x','y'])
    coord_data.set_index('data',inplace=True)
    data=[]
    data_outlier=[]
    with open(filename, 'r') as csvfile:
         csvreader = csv.reader(csvfile)
         for row in csvreader:
            #f=0
             #row=[float(i) for i in row]
            data.append(row)
            temp_outlier=[]
            temp_outlier.append(row[1])
            temp_outlier.append(row[2])
            data_outlier.append(temp_outlier)
            temp=row
            #if row[0].isnumeric():
            #    temp.append('cell')
            if len(row[0]) ==16:
                temp.append('cell')
            else:
                temp.append('gene')
                count=count+1
            csvData.append(temp)


    # # DB SCAN

    # In[20]:


    noise =[]
    print("Performing clustering....")
    db = DBSCAN(eps=35,min_samples=15).fit_predict(data_outlier)
    final_data=[]
    csvData=[['data','x','y','type']]
    for i in range(0,len(list(db))):
        if db[i]!=-1:
            final_data.append(data[i])
            csvData.append(data[i])
        if db[i]==-1:
            noise.append(data[i][0])
    data=final_data
    print(len(csvData))
    n_clusters = len(set(db)) - (1 if -1 in list(db) else 0)
    print("Clustering done. the number of obtained clusters: ",n_clusters)


    # In[13]:



    # # OUTLIER VISUALIZATION

    # In[21]:

    print("Starting outlier detection....")
    data_type = []
    c=0
    g=0
    for i in range(0,len(coord_data)):
        if db[i]!=-1:
            data_type.append("data")
        else:
            if len(coord_data.index[i])==16:
                data_type.append("cell_outliers")
            else:
                data_type.append("gene_outliers")
    coord_data["data_type"] = data_type
    data_colors = ["lightblue"]
    noise_colors = ['blue','red']
    coord_data["alpha"] = np.where(coord_data['data_type'] == 'data', 0.5,1.0)
    plt.figure(figsize=(15,10))
    ax = sns.scatterplot(x="x", y="y", data=coord_data[coord_data['alpha']==0.5],hue="data_type",palette=sns.xkcd_palette(data_colors),sizes=(50,100),size="data_type",alpha=0.3)
    sns.scatterplot(x="x", y="y", data=coord_data[coord_data['alpha']==1.0],hue="data_type",palette=sns.xkcd_palette(noise_colors),sizes=(50,100),size="data_type",marker="^",alpha=1.0,ax=ax)
    plt.savefig(output_path+'outliers_visualization.pdf',dpi=300)
    print("Outliers removed from the dataset....")

    # # POST-HOC CLUSTER ASSIGNMENT

    # In[23]:


    print("Starting post hoc clustering....")
    neighbor_df = pd.read_csv(input_path+'/neighbor.csv')
    neighbor_df.set_index('Unnamed: 0',inplace=True)
    col = list(neighbor_df.columns)
    index = list(neighbor_df.index)
    cell_dict = dict()
    column_dict = dict()
    for i in range(len(col)):
        column_dict[i]=col[i]
    for i in range(len(list(neighbor_df.index))):
        row = neighbor_df.iloc[i]
        col_ind = list(row.to_numpy().nonzero())[0]
        for ind in col_ind:
            if index[i] in cell_dict.keys():
                cell_dict[index[i]].append(column_dict[ind])
            else:
                temp=[]
                temp.append(column_dict[ind])
                cell_dict[index[i]]=temp
    cluster_assign = []
    for key_cell in cell_dict.keys():
        clust = dict()
        cells = cell_dict[key_cell]
        for cell in cells:
            cluster = db[list(coord_data.index).index(cell)]
            if cluster not in clust.keys():
                clust[cluster]=1
            else:
                clust[cluster]=clust[cluster]+1
        max_cluster=max(clust.items(), key=operator.itemgetter(1))[0]
        if max_cluster == -1:
            continue
        cluster_assign.append(max_cluster)
        x_total = 0
        y_total = 0
        count = 0
        for cell in cells:
            if db[list(coord_data.index).index(cell)] == max_cluster:
                count=count+1
                x_total=x_total+coord_data.loc[cell]['x']
                y_total=y_total+coord_data.loc[cell]['y']
        temp=[]
        temp.append(key_cell)
        temp.append(x_total/count)
        temp.append(y_total/count)
        temp.append('cell')
        csvData.append(temp)
    print("Post hoc clustering done....")


    # In[24]:



    with open(output_path+'data.csv', 'w') as csvFile:
        writer = csv.writer(csvFile)
        writer.writerows(csvData)
    csvFile.close()
    data_df = pd.read_csv(output_path+"data.csv",delimiter=",",index_col=False)
    clusters_info = [x for x in db if x!=-1]
    clusters_info = clusters_info + cluster_assign
    data_df['cluster'] = clusters_info
    data_df.to_csv(output_path+'data.csv')


    # In[26]:


    #colors = ["red","blue","green","yellow","pink","black","orange"]
    #colors = ["red","blue","green","yellow","pink","black","orange","brown"]
    #colors = ["red","blue","green","yellow","pink","black","orange","brown","purple"]
    #colors = ["red","blue","green","yellow","pink","black","orange","brown","purple","magenta"]
    #colors = ["red","blue","green","yellow","pink","black","orange","brown","purple","magenta","grey","silver"]
    colors = random.sample(seaborn_colors,n_clusters)
    plt.figure(figsize=(15,10))
    #cmap = sns.cubehelix_palette(dark=.3, light=.8, as_cmap=True)
    ax = sns.scatterplot(x="x", y="y", data=data_df,hue="cluster",palette=sns.xkcd_palette(colors))
    plt.savefig(output_path+"cluster_visualization.pdf",dpi=300);

    # # CONVEX HULL CLUSTER REPRESENTATION

    # In[27]:

    print("Convex hull representation of the clustering started....")
    from scipy.spatial import ConvexHull, convex_hull_plot_2d
    colors = random.sample(seaborn_colors,n_clusters)
    plt.figure(figsize=(15,10))
    ax = sns.scatterplot(x="x", y="y", data=data_df,hue="cluster",palette=sns.xkcd_palette(colors))
    for c in range(n_clusters):
        p=data_df[data_df["cluster"]==c]
        p=p[['x','y']]
        points = p.values
        hull = ConvexHull(points)
        for simplex in hull.simplices:
            sns.lineplot(points[simplex, 0], points[simplex, 1])
    plt.savefig(output_path+'convex_hull_representation.pdf',dpi=300)
    print("Convex hull representation of the clustering done....")


    # # CELL GENE MARKER

    # In[28]:

    import matplotlib.pyplot as plt, mpld3
    from scipy.spatial import ConvexHull, convex_hull_plot_2d
    data_df["alpha"] = np.where(data_df['type'] == 'gene', 1.0, 0.5)
    color_gene = ["lightblue"]
    color_cell = ["red"]
    #fig,ax1 = plt.subplots()
    plt.figure(figsize=(20,15))

    ax = sns.scatterplot(x="x", y="y", data=data_df[data_df['alpha']==0.5],hue="type",palette=sns.xkcd_palette(color_gene),sizes=(70,50),size="type",alpha=0.3)
    sns.scatterplot(x="x", y="y", data=data_df[data_df['alpha']==1.0],hue="type",palette=sns.xkcd_palette(color_cell),sizes=(70,50),size="type",marker="^",alpha=1.0,ax=ax)
    for c in range(n_clusters):
        p=data_df[data_df["cluster"]==c]
        p=p[['x','y']]
        points = p.values
        hull = ConvexHull(points)
        for simplex in hull.simplices:
            sns.lineplot(points[simplex, 0], points[simplex, 1])
    plt.savefig(output_path+'cell_gene_marker.pdf',dpi=300)





    # # GRAPH VISUALIZE

    # In[30]:

    print("graph network started loading....")
    position_dict={}
    for d in data:
        temp=[]
        temp.append(float(d[1]))
        temp.append(float(d[2]))
        position_dict[str(d[0])]=np.asarray(temp)


    # In[31]:


    G_normalized = nx.read_gexf(input_path+"/graph_normalized_own_cc.gexf")
    for node in noise:
        G_normalized.remove_node(node)


    # In[32]:


    nodes = list(G_normalized.nodes())
    cell_nodes=[]
    gene_nodes=[]
    for i in range(0,len(nodes)):
        if nodes[i].isnumeric():
            cell_nodes.append(i)
        else:
            gene_nodes.append(i)


    # In[33]:


    color_map=[]
    size_map=[]
    shape = []
    nodes = list(G_normalized.nodes())
    for node in nodes:
        if len(node)==16:
            color_map.append('blue')
            size_map.append(15)
            shape.append("o")
        else:
            color_map.append('red')
            size_map.append(100)
            shape.append("^")


    # In[34]:


    plt.figure(figsize=(20,15))
    nx.draw(G_normalized,pos=position_dict,node_color=color_map,node_size=size_map,edge_color='lightgrey',edge_width=0.015,node_shape="^")
    plt.savefig(output_path+'network_visualization.pdf',dpi=300)

    # # cells attached clusterwise for input gene

    # In[49]:


    #earch_gene = "GAPDH"
    #search_gene = "TXNIP"
    #search_gene = "LGALS2"
    #earch_gene = "CD79A"
    #search_gene = "CD79B"
    ##############MARKER GENES ########
    search_gene = "GNLY"
    #search_gene = "S100A8"
    nodes = G_normalized.nodes()
    cells = []
    for e in G_normalized.edges():
        if e[0]==search_gene and len(e[1]) ==16:
            cells.append(e[1])
        elif e[1]==search_gene and(len(e[0]))==16:
            cells.append(e[0])

    import matplotlib.pyplot as plt, mpld3
    from scipy.spatial import ConvexHull, convex_hull_plot_2d
    data_df = pd.read_csv(output_path+"data.csv",delimiter=",",index_col=False)
    cell_gene_df = data_df
    data = cell_gene_df["data"]
    alpha=[]
    for d in data:
        if d in cells:
            alpha.append(1)
        else:
            alpha.append(0.5)
    cell_gene_df["alpha"] = alpha
    color_cells = ["lightblue"]
    color_rest = ["red"]
    #fig,ax1 = plt.subplots()
    plt.figure(figsize=(20,15))

    ax = sns.scatterplot(x="x", y="y", data=cell_gene_df[cell_gene_df['alpha']==0.5],hue="alpha",palette=sns.xkcd_palette(color_cells),sizes=(70,50),size="alpha",alpha=0.3)
    sns.scatterplot(x="x", y="y", data=cell_gene_df[cell_gene_df['alpha']==1.0],hue="alpha",palette=sns.xkcd_palette(color_rest),sizes=(70,50),size="alpha",alpha=1.0,ax=ax)
    sns.scatterplot(x="x",y="y",data = cell_gene_df[cell_gene_df['data']==search_gene],s=500,color="blue",marker="^")
    plt.savefig(output_path+'gene_specific_cells.pdf',dpi=300)
    #for c in range(8):
    #    p=data_df[data_df["cluster"]==c]
    #    p=p[['x','y']]
    #    points = p.values
    #    hull = ConvexHull(points)
    #    for simplex in hull.simplices:
    #        sns.lineplot(points[simplex, 0], points[simplex, 1])


    # In[38]:


    import scanpy as sc
    sc.settings.verbosity = 2
    np.set_printoptions(precision=2)
    adata = sc.read_10x_mtx(
    inp_dataset,  # the directory with the `.mtx` file
    var_names='gene_symbols',                  # use gene symbols for the variable names (variables-axis index)
    cache=True)
    sc.pp.log1p(adata)
    orig_data_df = adata.to_df()
    read_data_df = adata.to_df()
    column_names = data_df.columns


    # In[50]:


    adata = sc.read_10x_mtx(
    inp_dataset,  # the directory with the `.mtx` file
    var_names='gene_symbols',                  # use gene symbols for the variable names (variables-axis index)
    cache=True)
    #earch_gene = "GAPDH"
    #search_gene = "TXNIP"
    #search_gene = "LGALS2"
    #search_gene = "CD79A"
    #search_gene = "CD79B"
    ####marker genes###########
    #search_gene = "PF4"
    search_gene = "GNLY"
    #search_gene = "S100A8"
    nodes = G_normalized.nodes()
    cells = []
    for e in G_normalized.edges():
        if e[0]==search_gene and len(e[1]) ==16:
            cells.append(e[1])
        elif e[1]==search_gene and(len(e[0]))==16:
            cells.append(e[0])

    import matplotlib.pyplot as plt, mpld3
    from scipy.spatial import ConvexHull, convex_hull_plot_2d
    data_df = pd.read_csv(output_path+"data.csv",delimiter=",",index_col=False)
    cell_gene_df = data_df
    cell_data = list(cell_gene_df[cell_gene_df["type"]=='cell']['data'])
    value_1 = []
    value_2 = []
    x_val_1 =[]
    x_val_2 = []
    y_val_1 = []
    y_val_2 = []
    cell_2 = []
    for d in cell_data:

        x_val_1.append(cell_gene_df[cell_gene_df['data']==d]['x'])
        y_val_1.append(cell_gene_df[cell_gene_df['data']==d]['y'])
        value_1.append(orig_data_df.loc[d][search_gene])
        #if d in cells:
        #    x_val_1.append(cell_gene_df[cell_gene_df['data']==d]['x'])
        #    y_val_1.append(cell_gene_df[cell_gene_df['data']==d]['y'])
        #    value_1.append(orig_data_df.loc[d][search_gene])
        #    #print(orig_data_df.loc[d][search_gene])
        #else:
        #    cell_2.append(d)
        #    x_val_2.append(cell_gene_df[cell_gene_df['data']==d]['x'])
        #    y_val_2.append(cell_gene_df[cell_gene_df['data']==d]['y'])
        #    value_2.append(0)

    plt.figure(figsize=(20,15))
    final_data_1 = pd.DataFrame();
    final_data_2 = pd.DataFrame();
    final_data_1['x'] = x_val_1
    final_data_1['y'] = y_val_1
    final_data_1['value'] = value_1
    final_data_1['data'] = cell_data
    final_data_2['x'] = x_val_2
    final_data_2['y'] = y_val_2
    final_data_2['value'] = value_2
    final_data_2['data'] = cell_2

    #ax = sns.scatterplot(x="x", y="y", data=final_data_2,hue="value",alpha=0.5)
    ax = sns.scatterplot(x="x", y="y", data=final_data_1,hue="value",alpha=1.0,size="value")
    plt.savefig(output_path+'gene_specific_expression_heatmap.pdf',dpi=300)
    #sns.scatterplot(x="x", y="y", data=final_data[final_data['value']!=0],hue="value",alpha=1.0,ax=ax)
    #sns.scatterplot(x="x", y="y", data=cell_gene_df[cell_gene_df['alpha']==1.0],hue="alpha",palette=sns.xkcd_palette(color_rest),sizes=(70,50),size="alpha",ax=ax)
    #sns.scatterplot(x="x",y="y",data = cell_gene_df[cell_gene_df['data']==search_gene],s=500,color="blue",marker="^")
    #for c in range(8):
    #    p=data_df[data_df["cluster"]==c]
    #    p=p[['x','y']]
    #    points = p.values
    #    hull = ConvexHull(points)
    #    for simplex in hull.simplices:
    #        sns.lineplot(points[simplex, 0], points[simplex, 1])


    # In[53]:


    #search_gene = "S100A8"
    #search_gene = "PF4"
    search_gene = "GNLY"
    import matplotlib.pyplot as plt, mpld3
    data_df = pd.read_csv(output_path+"data.csv",delimiter=",",index_col=False)
    info_data_df = data_df
    info_data_df.set_index('data',inplace=True)
    pca_df = pd.read_csv(input_path+"pca_genes_df.csv",delimiter=",",index_col=False)
    pca_df.set_index('Unnamed: 0',inplace=True)
    cells_data = pca_df.loc[:,search_gene]
    format_df = cells_data.to_frame()
    cells = format_df.index
    cluster = []
    for cell in cells:
        if cell in info_data_df.index :
            cluster.append(info_data_df.loc[cell]['cluster'])
        else:
            format_df=format_df.drop(cell)
    format_df['cluster'] = cluster
    sx = sns.violinplot(x="cluster",y=search_gene,data=format_df)
    plt.savefig(output_path+'gene_specific_violinplot.pdf',dpi=300)

    # In[43]:


    heat_map_df = pd.read_csv(input_path+"pca_genes_df.csv",delimiter=",",index_col=False)
    heat_map_df.fillna(0)
    data_df = pd.read_csv(output_path+"data.csv",delimiter=",",index_col=False)
    data_df.set_index('data',inplace=True)
    data_df.drop('Unnamed: 0',axis=1,inplace=True)
    cells = list(data_df.index)
    heat_map_df = pd.read_csv(input_path+"pca_genes_df.csv",delimiter=",",index_col=False)
    heat_map_df.set_index('Unnamed: 0',inplace=True)
    cluster=[]
    for ind in range(len(list(heat_map_df.index))):
        if data_df.iloc[ind]['cluster']==-1:
            heat_map_df.drop(ind,inplace=True)
        else:
            cluster.append(data_df.iloc[ind]['cluster'])
    heat_map_df['cluster'] = cluster
    colors = {0:"red",1:"blue",2:"green",3:"yellow",4:"pink",5:"black",6:"orange",7:"brown"}
    map_colors = heat_map_df.cluster.map(colors)
    heat_map_df=heat_map_df.transpose()
    #heat_map_df.columns = heat_map_df.iloc[0]
    #heat_map_df =heat_map_df[1:]
    heat_map_df = heat_map_df[heat_map_df.columns].astype(float)
    #heat_map_df = heat_map_df[heat_map_df.columns].astype(float)
    #map_colors = heat_map_df.cluster.map(colors)
    plt.figure(figsize=(20,15))
    sns.clustermap(heat_map_df,metric="correlation",col_colors=map_colors)
    plt.savefig(output_path+'clustermap.pdf',dpi=300)


# In[44]:


input_path=""
output_path=""
file=""
arguments=sys.argv[1:]
if len(arguments)==6 and arguments[0]=="-inp_data" and arguments[2]=="-i" and arguments[4]=="-o":
    input_dataset=arguments[1]
    input_path=arguments[3]
    output_path=arguments[5]
    #file=open(output_file,"w")
    compute(input_dataset,input_path,output_path);
else:
        print ("\nPlease follow the command: <file>.py -inp_data inp_dataset -i inputpath -o outputpath \n")
        print('inputpath : path for read 10x data\n')
        print('outputpath : path where the output is to be stored\n')
        sys.exit()

