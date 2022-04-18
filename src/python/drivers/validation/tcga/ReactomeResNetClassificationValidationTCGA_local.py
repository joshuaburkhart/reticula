#!/usr/bin/env python
# coding: utf-8

import os.path as op
import random
import time

import matplotlib.pyplot as plt
import numpy
import sklearn
import torch.nn.functional as nn_func
from sklearn import preprocessing
from sklearn.metrics import adjusted_rand_score
from torch.nn import Linear
from torch_geometric.data import Data, DataLoader
import torch
from torch import nn
import torchvision.models as models

random.seed = 88888888

node_features_fn = '/home/jgburk/PycharmProjects/reticula/data/tcga/input/node_features.txt'
graph_targets_fn = '/home/jgburk/PycharmProjects/reticula/data/tcga/input/graph_targets.txt'
model_fn = '/home/jgburk/PycharmProjects/reticula/data/gtex/resnet_trained_pytorch_model_fold_full_dataset.pt'
output_fn = '/home/jgburk/PycharmProjects/reticula/data/tcga/output/resnet_predictions.tsv'

# test graph_targets.txt, node_features.txt and edges.txt
features_exist = op.exists(node_features_fn)
targets_exist = op.exists(graph_targets_fn)
model_exists = op.exists(model_fn)

print(f'features exist: {features_exist},'
      f' targets exist: {targets_exist},'
      f' model exists: {model_exists}')
assert features_exist
assert targets_exist
assert model_exists

# magic numbers
INPUT_CHANNELS = 1
OUTPUT_CHANNELS = 51
NEW_CHANNELS = 20
HIDDEN_CHANNELS = 64
BATCH_SIZE = 64
BENCHMARKING = False
EPOCHS = 500


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


def build_resnet_datalist(n_features_fn, g_targets_fn):
    feature_v = numpy.loadtxt(n_features_fn)
    target_v = numpy.loadtxt(g_targets_fn, dtype=str, delimiter="\n")

    target_encoder = sklearn.preprocessing.LabelEncoder()
    target_v = target_encoder.fit_transform(target_v)

    d_list = []
    for row_idx in range(len(feature_v)):
        x = torch.tensor(feature_v[row_idx, :], dtype=torch.float)
        x = x.reshape(2, 31, 173)
        y = torch.tensor([target_v[row_idx]])
        d_list.append({'x': x, 'y': y})

    return d_list


def build_reactome_graph_loader(d_list, batch_size):
    loader = DataLoader(d_list, batch_size=batch_size, shuffle=False)  # True)

    return loader


def train(loader, dv):
    model.train()

    correct = 0
    for batch in loader:  # Iterate in batches over the training dataset.
        x = batch['x'].to(dv)
        y = batch['y'].to(dv)
        out = model(x)  # Perform a single forward pass.
        y = torch.squeeze(y)
        loss = criterion(out, y)  # Compute the loss.
        loss.backward()  # Derive gradients.
        optimizer.step()  # Update parameters based on gradients.
        optimizer.zero_grad()  # Clear gradients.
        pred = out.argmax(dim=1)  # Use the class with highest probability.
        correct += int((pred == y).sum())  # Check against ground-truth labels.
    return correct / len(loader.dataset)  # Derive ratio of correct predictions.


def test(loader, dv):
    model.eval()

    targets = []
    predictions = []
    for batch in loader:  # Iterate in batches over the test dataset.
        x = batch['x'].to(dv)
        y = batch['y'].to(dv)
        targets += torch.Tensor.tolist(torch.squeeze(y))
        out = model(x)  # Perform a single forward pass.
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


model = models.resnet18(num_classes=51)
conv1 = model.conv1
model.conv1 = nn.Conv2d(2,
                        conv1.out_channels,
                        conv1.kernel_size,
                        conv1.stride,
                        conv1.padding,
                        conv1.dilation,
                        conv1.groups,
                        conv1.bias)
device = cpu = torch.device('cpu')

sd = torch.load(model_fn, map_location=device)
# change_key(sd, 'conv1.lin_l.weight', 'conv1.lin_rel.weight')
# change_key(sd, 'conv1.lin_l.bias', 'conv1.lin_rel.bias')
# change_key(sd, 'conv1.lin_r.weight', 'conv1.lin_root.weight')
# change_key(sd, 'conv2.lin_l.weight', 'conv2.lin_rel.weight')
# change_key(sd, 'conv2.lin_l.bias', 'conv2.lin_rel.bias')
# change_key(sd, 'conv2.lin_r.weight', 'conv2.lin_root.weight')
# change_key(sd, 'conv3.lin_l.weight', 'conv3.lin_rel.weight')
# change_key(sd, 'conv3.lin_l.bias', 'conv3.lin_rel.bias')
# change_key(sd, 'conv3.lin_r.weight', 'conv3.lin_root.weight')
# change_key(sd, 'lin.weight', 'lin.weight')
# change_key(sd, 'lin.bias', 'lin.bias')

model.load_state_dict(sd)
model.eval()

# replace final layer with new shape matching new dataset
#print(f'in_features: {model.fc.in_features}, out_features:{model.fc.out_features}, bias: {model.fc.bias}')
model.fc = Linear(in_features=model.fc.in_features,
                  out_features=NEW_CHANNELS,
                  bias=True)

for param in model.parameters():
    param.requires_grad = False

model.fc.weight.requires_grad = True
model.fc.bias.requires_grad = True

optimizer = torch.optim.AdamW(model.parameters())
criterion = torch.nn.CrossEntropyLoss()

data_list = build_resnet_datalist(node_features_fn, graph_targets_fn)
print(len(data_list))
# retrain model for fine tuning transfer learning
train_data_list = data_list[:370]
print(len(train_data_list))
print(f'Number of training graphs: {len(train_data_list)}')
train_data_loader = build_reactome_graph_loader(train_data_list, BATCH_SIZE)
for epoch in range(EPOCHS):
    train(train_data_loader, device)
    train_acc = train(train_data_loader, device)
    print(f'Epoch: {epoch}, Train Acc: {train_acc}')

test_data_list = data_list[370:]
print(len(test_data_list))
print(f'Number of test graphs: {len(test_data_list)}')

test_data_loader = build_reactome_graph_loader(test_data_list, BATCH_SIZE)
test_ari = test(test_data_loader, device)
print(f'test_ari: {test_ari}')

model_save_name = f'tuned_pytorch_tcga_resnet_model.pt'
path = f'/home/jgburk/PycharmProjects/reticula/data/tcga/output/{model_save_name}'
torch.save(model.state_dict(), path)
print(f'model saved as {path}')

# real network gets to 0.8417
