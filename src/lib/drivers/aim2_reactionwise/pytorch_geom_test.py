# testing install from https://pytorch-geometric.readthedocs.io/en/latest/notes/installation.html
# https://pytorch.org/get-started/locally/#macos-version
# https://www.jetbrains.com/help/pycharm/conda-support-creating-conda-virtual-environment.html

from __future__ import print_function
from torch_geometric.datasets import Planetoid
import torch
import torch.nn.functional as F
from torch_geometric.data import Data
from torch_geometric.nn import GCNConv

x = torch.rand(5, 3)
print(x)

edge_index = torch.tensor([[0, 1],
                           [1, 0],
                           [1, 2],
                           [2, 1]],
                          dtype=torch.long)
x_list = [torch.tensor([[-1.1, 2.3], [0.5, 21.1], [1.1, 3.9]], dtype=torch.float),
          torch.tensor([[-2.1, 4.3], [0.2, 21.3], [2.1, 3.2]], dtype=torch.float),
          torch.tensor([[2.1, 3.3], [3.5, 61.1], [5.1, 3.5]], dtype=torch.float),
          torch.tensor([[2.1, 3.3], [3.5, 61.1], [5.1, 3.5]], dtype=torch.float),
          torch.tensor([[2.1, 3.3], [1.5, 91.1], [1.1, 9.5]], dtype=torch.float),
          torch.tensor([[2.1, 3.8], [2.5, 81.1], [2.1, 8.5]], dtype=torch.float),
          torch.tensor([[2.1, 3.7], [3.5, 71.1], [3.1, 7.5]], dtype=torch.float),
          torch.tensor([[2.1, 3.6], [4.5, 61.1], [4.1, 6.5]], dtype=torch.float),
          torch.tensor([[2.1, 3.5], [5.5, 51.1], [5.1, 5.5]], dtype=torch.float),
          torch.tensor([[2.1, 3.4], [6.5, 41.1], [5.6, 4.5]], dtype=torch.float),
          torch.tensor([[2.1, 3.3], [3.5, 11.1], [5.7, 3.5]], dtype=torch.float),
          torch.tensor([[2.1, 3.3], [3.5, 21.1], [5.8, 2.5]], dtype=torch.float),
          torch.tensor([[2.1, 2.1], [1.5, 7.1], [7.1, 3.2]], dtype=torch.float)
          ]

y_list = [torch.tensor([1]),
          torch.tensor([3]),
          torch.tensor([2]),
          torch.tensor([2]),
          torch.tensor([3]),
          torch.tensor([4]),
          torch.tensor([1]),
          torch.tensor([2]),
          torch.tensor([3]),
          torch.tensor([2]),
          torch.tensor([3]),
          torch.tensor([4]),
          torch.tensor([4])]

data_list = [Data(x=x_list[0],
                  edge_index=edge_index.t().contiguous(),
                  y=y_list[0]),
             Data(x=x_list[1],
                  edge_index=edge_index.t().contiguous(),
                  y=y_list[1]),
             Data(x=x_list[2],
                  edge_index=edge_index.t().contiguous(),
                  y=y_list[2]),
             Data(x=x_list[3],
                  edge_index=edge_index.t().contiguous(),
                  y=y_list[3]),
             Data(x=x_list[4],
                  edge_index=edge_index.t().contiguous(),
                  y=y_list[4]),
             Data(x=x_list[5],
                  edge_index=edge_index.t().contiguous(),
                  y=y_list[5]),
             Data(x=x_list[6],
                  edge_index=edge_index.t().contiguous(),
                  y=y_list[6]),
             Data(x=x_list[7],
                  edge_index=edge_index.t().contiguous(),
                  y=y_list[7]),
             Data(x=x_list[8],
                  edge_index=edge_index.t().contiguous(),
                  y=y_list[8]),
             Data(x=x_list[9],
                  edge_index=edge_index.t().contiguous(),
                  y=y_list[9]),
             Data(x=x_list[10],
                  edge_index=edge_index.t().contiguous(),
                  y=y_list[10]),
             Data(x=x_list[11],
                  edge_index=edge_index.t().contiguous(),
                  y=y_list[11]),
             Data(x=x_list[12],
                  edge_index=edge_index.t().contiguous(),
                  y=y_list[12])]


class Net(torch.nn.Module):
    def __init__(self):
        super(Net, self).__init__()
        self.conv1 = GCNConv(dataset.num_node_features, 16)  # 2 features per node, 13 examples of the graph
        self.conv2 = GCNConv(16, dataset.num_classes)  # 13 examples of the graph, 4 classes of graph

    def forward(self, data):
        x, edge_index = data.x, data.edge_index

        x = self.conv1(x, edge_index)
        x = F.relu(x)
        x = F.dropout(x, training=self.training)
        x = self.conv2(x, edge_index)

        return F.log_softmax(x, dim=1)


dataset = Planetoid(root='/tmp/Cora', name='Cora')

device = torch.device('cpu')
model = Net().to(device)
data = dataset[0].to(device)
optimizer = torch.optim.Adam(model.parameters(), lr=0.01, weight_decay=5e-4)

model.train()
for epoch in range(200):
    optimizer.zero_grad()
    out = model(data)
    loss = F.nll_loss(out[data.train_mask], data.y[data.train_mask])
    loss.backward()
    optimizer.step()

model.eval()
_, pred = model(data).max(dim=1)
correct = float(pred[data.test_mask].eq(data.y[data.test_mask]).sum().item())
acc = correct / data.test_mask.sum().item()
print('Accuracy: {:.4f}'.format(acc))
