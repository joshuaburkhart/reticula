#!/usr/bin/env python
# coding: utf-8

# In[2]:


get_ipython().system('pip install -q torch-scatter -f https://pytorch-geometric.com/whl/torch-1.8.0+cu101.html')
get_ipython().system('pip install -q torch-sparse -f https://pytorch-geometric.com/whl/torch-1.8.0+cu101.html')
get_ipython().system('pip install -q torch-geometric')
get_ipython().system('pip install -q sklearn')
get_ipython().system('pip install -q captum')


# In[4]:


import time
import torch
import numpy
import random
import sklearn
import numpy as np
import networkx as nx
import torch.nn.functional as F
from torch.nn import Linear
from sklearn import preprocessing
from collections import defaultdict
from IPython.display import Javascript
from torch_geometric.utils import to_networkx
from torch_geometric.datasets import TUDataset
from torch_geometric.data import Data, DataLoader
from captum.attr import Saliency, IntegratedGradients
from torch_geometric.nn import GraphConv, global_mean_pool


# In[5]:


random.seed = 88888888
device = cuda0 = torch.device('cuda:0')
cpu = torch.device('cpu')
INPUT_CHANNELS = 1
OUTPUT_CHANNELS = 51
HIDDEN_CHANNELS = 64

edges_fn = '/content/gdrive/My Drive/Academia/OHSU/Proposal/Aim 2/Edge Sets/pathway_hierarchy_edges.txt'
node_features_fn = '/content/gdrive/My Drive/Academia/OHSU/Proposal/Aim 2/Node Features/pathway_hierarchy_node_features.txt'
graph_targets_fn = '/content/gdrive/My Drive/Academia/OHSU/Proposal/Aim 2/Target Sets/pathway_hierarchy_graph_targets.txt'

transformed_targets_path = '/content/gdrive/My Drive/Academia/OHSU/Proposal/Aim 2/pathway_hierarchy_transformed_targets.txt'
inverted_targets_path = '/content/gdrive/My Drive/Academia/OHSU/Proposal/Aim 2/pathway_hierarchy_inverted_targets.txt'

model_save_name = 'pathway_hierarchy_trained_pytorch_model_fold_full_dataset.pt'
model_save_path = F"/content/gdrive/My Drive/Academia/OHSU/Proposal/Aim 2/Pytorch Models/PathwayHierarchy_10FoldCV_Full_RandFeat_RandTarget/{model_save_name}"

drive.mount('/content/gdrive')

get_ipython().system('ls "{edges_fn}"')
get_ipython().system('ls "{node_features_fn}"')
get_ipython().system('ls "{graph_targets_fn}"')


# In[15]:


class GNN(torch.nn.Module):
    def __init__(self, hidden_channels):
        super(GNN, self).__init__()

        self.conv1 = GraphConv(INPUT_CHANNELS, hidden_channels)
        self.conv2 = GraphConv(hidden_channels,hidden_channels)
        self.conv3 = GraphConv(hidden_channels,hidden_channels)
        self.lin = Linear(hidden_channels, OUTPUT_CHANNELS)

    def forward(self, x, edge_index, batch, edge_weight=None):
        # 1. Obtain node embeddings 
        x = self.conv1(x, edge_index, edge_weight)
        x = x.relu()
        x = self.conv2(x, edge_index, edge_weight)
        x = x.relu()
        x = self.conv3(x, edge_index, edge_weight)

        # 2. Readout layer
        x = global_mean_pool(x, batch)  # [batch_size, hidden_channels]

        # 3. Apply a final classifier
        x = F.dropout(x, training=self.training)
        x = self.lin(x)
        
        return x

def read_reactome_graph(edges_fn, node_features_fn):
    edge_v1 = []
    edge_v2 = []

    for line in open(edges_fn, 'r'):
        data = line.split()
        node1 = int(data[0]) - 1 #subtracting to convert R idx to python idx
        node2 = int(data[1]) - 1 # " "
        edge_v1.append( node1 )
        edge_v2.append( node2 )

    return edge_v1, edge_v2

def build_reactome_graph_datalist(edge_v1, edge_v2, node_features_fn, graph_targets_fn):
    edge_index = torch.tensor([edge_v1, edge_v2], dtype = torch.long)
    feature_v = numpy.loadtxt(node_features_fn)
    target_v = numpy.loadtxt(graph_targets_fn,dtype=str,delimiter="\n")
    
    target_encoder = sklearn.preprocessing.LabelEncoder()
    target_v = target_encoder.fit_transform(target_v)

    data_list = []
    for row_idx in range(len(feature_v)):
      x = torch.tensor(feature_v[row_idx,:],dtype=torch.float)
      x = x.unsqueeze(1)
      y = torch.tensor([target_v[row_idx]])
      data_list.append(Data(x = x, y = y, edge_index = edge_index))

    return data_list

def explain(method, data, target=0):
    input_mask = torch.ones(data.edge_index.shape[1]).requires_grad_(True).to(device)
    if method == 'ig':
        ig = IntegratedGradients(model_forward)
        mask = ig.attribute(input_mask, target=target,
                            additional_forward_args=(data,),
                            internal_batch_size=data.edge_index.shape[1])
    elif method == 'saliency':
        saliency = Saliency(model_forward)
        mask = saliency.attribute(input_mask, target=target,
                                  additional_forward_args=(data,))
    else:
        raise Exception('Unknown explanation method')

    edge_mask = np.abs(mask.cpu().detach().numpy())
    if edge_mask.max() > 0:  # avoid division by zero
        edge_mask = edge_mask / edge_mask.max()
    return edge_mask

def aggregate_edge_directions(edge_mask, data):
  edge_mask_dict = defaultdict(float)
  for val, u, v in list(zip(edge_mask, *data.edge_index)):
    u, v = u.item(), v.item()
    if u > v:
      u, v = v, u
    edge_mask_dict[(u, v)] += val
  return edge_mask_dict

def model_forward(edge_mask, data):
  batch = torch.zeros(data.x.shape[0], dtype=int).to(device)
  out = model(data.x,
              data.edge_index, 
              batch,
              edge_mask)
  return out


# In[16]:


(edge_v1, edge_v2) = read_reactome_graph(edges_fn, node_features_fn)

data_list = build_reactome_graph_datalist(edge_v1, edge_v2, node_features_fn, graph_targets_fn)
data_loader = DataLoader(data_list)


# In[17]:


# rebuild label encoder to invert numerical transformation
target_v = numpy.loadtxt(graph_targets_fn,dtype=str,delimiter="\n")
target_encoder = sklearn.preprocessing.LabelEncoder()


# In[18]:


target_v = target_encoder.fit_transform(target_v)
path = transformed_targets_path
numpy.savetxt(path, target_v, delimiter=",", fmt="%.0f")
print(F"target_v saved as {path}")


# In[19]:


target_l = target_encoder.inverse_transform(target_v)
path = inverted_targets_path 
numpy.savetxt(path, target_l, delimiter=",", fmt="%s")
print(F"target_l saved as {path}")


# In[20]:


model = GNN(hidden_channels=HIDDEN_CHANNELS)
model = model.to(device)
path = model_save_path
model.load_state_dict(torch.load(path))
model.eval()


# In[21]:


d = data_loader.dataset[0]
d.edge_index.shape[1]


# In[22]:


data = data_loader.dataset[0]


# In[23]:


for target_tissue in range(51):
  for title, method in [('Integrated Gradients', 'ig'), ('Saliency', 'saliency')]:
    data.to(device)
    print(F"processing tissue {target_tissue} with {title}, a.k.a. {method}")
    edge_mask = explain(method, data, target=target_tissue)
    #edge_mask_dict = aggregate_edge_directions(edge_mask, data)
    path = F"/content/gdrive/My Drive/Academia/OHSU/Proposal/{method}_{target_tissue}.txt"
    numpy.savetxt(path, edge_mask, delimiter=",")
    print(F"{method} {target_tissue} edges saved as {path}")

