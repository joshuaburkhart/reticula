#!/usr/bin/env python
# coding: utf-8

import os.path as op
import random
import time

import matplotlib.pyplot as plt
import numpy
import sklearn
import torch
import torch.nn.functional as nn_func
from sklearn import preprocessing
from sklearn.metrics import adjusted_rand_score
from torch.nn import Linear
from torch_geometric.data import Data
from torch_geometric.loader import DataListLoader
from torch_geometric.nn import GraphConv, global_mean_pool, DataParallel

random.seed = 88888888

node_features_fn = '/home/jgburk/zsl_validation/test/zsl_tcga_node_features.txt'
graph_targets_fn = '/home/jgburk/zsl_validation/test/zsl_tcga_graph_targets.txt'
edges_fn = '/home/jgburk/zsl_validation/test/edges.txt'
model_fn = '/home/jgburk/zsl_validation/test/zsl_fully_trained_pytorch_gtex_gnn_model.pt'
output_fn = '/home/jgburk/zsl_validation/test/zsl_gnn_predictions.tsv'

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
OUTPUT_CHANNELS = 51
HIDDEN_CHANNELS = 64
BATCH_SIZE = 64
BENCHMARKING = False
EPOCHS = 500


# from https://colab.research.google.com/drive/1I8a0DfQ3fI7Njc62__mVXUlcAleUclnb?usp=sharing#scrollTo=CN3sRVuaQ88l
class GNN(torch.nn.Module):
    def __init__(self, hidden_channels):
        super(GNN, self).__init__()

        self.conv1 = GraphConv(INPUT_CHANNELS, hidden_channels)
        self.conv2 = GraphConv(hidden_channels, hidden_channels)
        self.conv3 = GraphConv(hidden_channels, hidden_channels)
        self.lin = Linear(hidden_channels, OUTPUT_CHANNELS)

    def forward(self, data, edge_weight=None):
        x = data.x
        edge_index = data.edge_index
        batch = data.batch
        # 1. Obtain node embeddings
        x = self.conv1(x, edge_index, edge_weight)
        x = x.relu()
        x = self.conv2(x, edge_index, edge_weight)
        x = x.relu()
        x = self.conv3(x, edge_index, edge_weight)

        # 2. Readout layer
        x = global_mean_pool(x, batch)  # [batch_size, hidden_channels]

        # 3. Apply a final classifier
        x = nn_func.dropout(x, training=self.training)
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


def build_reactome_graph_datalist(e_v1, e_v2, n_fn, g_fn):
    edge_index = torch.tensor([e_v1, e_v2], dtype=torch.long)
    feature_v = numpy.loadtxt(n_fn)
    target_v = numpy.loadtxt(g_fn, dtype=str, delimiter="\n")

    target_encoder = sklearn.preprocessing.LabelEncoder()
    target_v = target_encoder.fit_transform(target_v)

    d_list = []
    for row_idx in range(len(feature_v)):
        features = feature_v[row_idx, :]
        x = torch.tensor(features, dtype=torch.float)
        x = x.unsqueeze(1)
        y = torch.tensor([target_v[row_idx]])
        d = Data(x=x,y=y,edge_index=edge_index)
        d_list.append(d)

    return d_list


def build_reactome_graph_loader(d_list, batch_size):
    loader = DataListLoader(d_list, batch_size=batch_size, shuffle=False)#True)

    return loader

def test(loader, dv):
    model.eval()

    targets = []
    predictions = []
    for data in loader:  # Iterate in batches over the test dataset.
        #data.to(dv)
        y = torch.cat([d.y for d in data]).to(dv)
        targets += torch.Tensor.tolist(y)
        out = model(data)  # Perform a single forward pass.
        pred = out.argmax(dim=1)  # Use the class with highest probability.
        predictions += torch.Tensor.tolist(pred)
    print(targets)
    print(predictions)
    numpy.savetxt(output_fn, numpy.transpose([targets, predictions]),
                  fmt='%d', delimiter='\t', header='target\tprediction')
    ari = adjusted_rand_score(targets, predictions)
    print(f'ari: {ari}')
    return ari


# from https://stackoverflow.com/questions/12150872/change-key-in-ordereddict-without-losing-order
def change_key(self, old, new):
    for _ in range(len(self)):
        k, v = self.popitem(False)
        self[new if old == k else k] = v


(edge_v1, edge_v2) = read_reactome_graph(edges_fn)
model = GNN(hidden_channels=HIDDEN_CHANNELS)
model = DataParallel(model)
device = cuda0 = torch.device('cuda:0')
model.to(device)
sd = torch.load(model_fn, map_location=device)
#change_key(sd, 'conv1.lin_l.weight', 'conv1.lin_rel.weight')
#change_key(sd, 'conv1.lin_l.bias', 'conv1.lin_rel.bias')
#change_key(sd, 'conv1.lin_r.weight', 'conv1.lin_root.weight')
#change_key(sd, 'conv2.lin_l.weight', 'conv2.lin_rel.weight')
#change_key(sd, 'conv2.lin_l.bias', 'conv2.lin_rel.bias')
#change_key(sd, 'conv2.lin_r.weight', 'conv2.lin_root.weight')
#change_key(sd, 'conv3.lin_l.weight', 'conv3.lin_rel.weight')
#change_key(sd, 'conv3.lin_l.bias', 'conv3.lin_rel.bias')
#change_key(sd, 'conv3.lin_r.weight', 'conv3.lin_root.weight')
#change_key(sd, 'lin.weight', 'lin.weight')
#change_key(sd, 'lin.bias', 'lin.bias')
model.load_state_dict(sd)
model.eval()

data_list = build_reactome_graph_datalist(edge_v1, edge_v2, node_features_fn, graph_targets_fn)
print(len(data_list))

test_data_list = data_list
print(len(test_data_list))
print(f'Number of test graphs: {len(test_data_list)}')

test_data_loader = build_reactome_graph_loader(test_data_list, BATCH_SIZE)
test_ari = test(test_data_loader, device)
print(f'test_ari: {test_ari}')

# real network gets to 0.8417
