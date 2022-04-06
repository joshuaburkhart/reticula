#!/usr/bin/env python
# coding: utf-8


import os.path as op
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

random.seed = 88888888

node_features_fn = '/home/jgburk/PycharmProjects/reticula/data/SRP049593/input/node_features.txt'
graph_targets_fn = '/home/jgburk/PycharmProjects/reticula/data/SRP049593/input/graph_targets.txt'
edges_fn = '/home/jgburk/PycharmProjects/reticula/data/gtex/Copy of edges.txt'
model_fn = '/home/jgburk/PycharmProjects/reticula/data/SRP049593/output/500_epoch_tuned_pytorch_srp049593_model.pt'
transformed_targets_fn = '/home/jgburk/PycharmProjects/reticula/data/SRP049593/output/transformed_targets.txt'
inverted_targets_fn = '/home/jgburk/PycharmProjects/reticula/data/SRP049593/output/inverted_targets.txt'

# test graph_targets.txt, node_features.txt and edges.txt
features_exist = op.exists(node_features_fn)
targets_exist = op.exists(graph_targets_fn)
edges_exist = op.exists(edges_fn)
model_exists = op.exists(model_fn)

print(f'features exist: {features_exist},'
      f' targets exist: {targets_exist},'
      f' edges exist: {edges_exist}',
      f' model exists: {model_exists}')
assert features_exist
assert targets_exist
assert edges_exist
assert model_exists

# magic numbers
INPUT_CHANNELS = 1
OUTPUT_CHANNELS = 2
HIDDEN_CHANNELS = 64

class GNN(torch.nn.Module):
    def __init__(self, hidden_channels):
        super(GNN, self).__init__()

        self.conv1 = GraphConv(INPUT_CHANNELS, hidden_channels)
        self.conv2 = GraphConv(hidden_channels, hidden_channels)
        self.conv3 = GraphConv(hidden_channels, hidden_channels)
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


def read_reactome_graph(e_fn):
    e_v1 = []
    e_v2 = []

    for line in open(e_fn, 'r'):
        dt = line.split()
        node1 = int(dt[0]) - 1  # subtracting to convert R idx to python idx
        node2 = int(dt[1]) - 1  # " "
        e_v1.append(node1)
        e_v2.append(node2)

    return e_v1, e_v2


def build_reactome_graph_datalist(e_v1, e_v2, n_features_fn, g_targets_fn):
    edge_index = torch.tensor([e_v1, e_v2], dtype=torch.long)
    feature_v = numpy.loadtxt(n_features_fn)
    t_v = numpy.loadtxt(g_targets_fn, dtype=str, delimiter="\n")

    t_encoder = sklearn.preprocessing.LabelEncoder()
    t_v = t_encoder.fit_transform(t_v)

    d_list = []
    for row_idx in range(len(feature_v)):
        x = torch.tensor(feature_v[row_idx, :], dtype=torch.float)
        x = x.unsqueeze(1)
        y = torch.tensor([t_v[row_idx]])
        d_list.append(Data(x=x, y=y, edge_index=edge_index))

    return d_list


def explain(m, dt, target):
    input_mask = torch.ones(dt.edge_index.shape[1]).requires_grad_(True).to(device)
    if m == 'ig':
        ig = IntegratedGradients(model_forward)
        mask = ig.attribute(input_mask, target=target,
                            additional_forward_args=(dt,),
                            internal_batch_size=dt.edge_index.shape[1])
    else:
        raise Exception('Unknown explanation method')

    e_mask = np.abs(mask.cpu().detach().numpy())
    if e_mask.max() > 0:  # avoid division by zero
        e_mask = e_mask / e_mask.max()
    return e_mask


def aggregate_edge_directions(e_mask, dt):
    edge_mask_dict = defaultdict(float)
    for val, u, v in list(zip(e_mask, *dt.edge_index)):
        u, v = u.item(), v.item()
        if u > v:
            u, v = v, u
        edge_mask_dict[(u, v)] += val
    return edge_mask_dict


def model_forward(e_mask, dt):
    batch = torch.zeros(dt.x.shape[0], dtype=int).to(device)
    out = model(dt.x,
                dt.edge_index,
                batch,
                e_mask)
    return out


(edge_v1, edge_v2) = read_reactome_graph(edges_fn)

data_list = build_reactome_graph_datalist(edge_v1, edge_v2, node_features_fn, graph_targets_fn)
data_loader = DataLoader(data_list)

# rebuild label encoder to invert numerical transformation
target_v = numpy.loadtxt(graph_targets_fn, dtype=str, delimiter="\n")
target_encoder = sklearn.preprocessing.LabelEncoder()

target_v = target_encoder.fit_transform(target_v)
path = transformed_targets_fn
numpy.savetxt(path, target_v, delimiter=",", fmt="%.0f")
print(F"target_v saved as {path}")

target_l = target_encoder.inverse_transform(target_v)
path = inverted_targets_fn
numpy.savetxt(path, target_l, delimiter=",", fmt="%s")
print(F"target_l saved as {path}")

model = GNN(hidden_channels=HIDDEN_CHANNELS)
device = cpu = torch.device('cpu')
model = model.to(device)
path = model_fn
model.load_state_dict(torch.load(path, map_location=device))
model.eval()

d = data_loader.dataset[0]
d.edge_index.shape[1]

data = data_loader.dataset[0]

for target_tissue in range(2):
    title = 'Integrated Gradients'
    method = 'ig'
    data.to(device)
    print(F"processing tissue {target_tissue} with {title}, a.k.a. {method}")
    edge_mask = explain(method, data, target=target_tissue)
    # edge_mask_dict = aggregate_edge_directions(edge_mask, data)
    path = F"/home/jgburk/PycharmProjects/reticula/data/SRP049593/output/{method}_{target_tissue}.txt"
    numpy.savetxt(path, edge_mask, delimiter=",")
    print(F"{method} {target_tissue} edges saved as {path}")
