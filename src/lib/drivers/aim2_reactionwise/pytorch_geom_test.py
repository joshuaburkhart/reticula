# testing install from https://pytorch-geometric.readthedocs.io/en/latest/notes/installation.html
# https://pytorch.org/get-started/locally/#macos-version
# https://www.jetbrains.com/help/pycharm/conda-support-creating-conda-virtual-environment.html

from __future__ import print_function

import torch
import torch.nn.functional as F
from torch_geometric.data import Data, DataLoader
from torch_geometric.nn import GCNConv

edge_index = torch.tensor([[0, 1],
                           [1, 0],
                           [1, 2],
                           [2, 1]],
                          dtype=torch.long)

x_list = [torch.tensor([[2.3], [21.1], [3.9]], dtype=torch.float),
          torch.tensor([[4.3], [21.3], [3.2]], dtype=torch.float),
          torch.tensor([[3.3], [61.1], [3.5]], dtype=torch.float),
          torch.tensor([[3.3], [61.1], [3.5]], dtype=torch.float),
          torch.tensor([[3.3], [91.1], [9.5]], dtype=torch.float),
          torch.tensor([[3.8], [81.1], [8.5]], dtype=torch.float),
          torch.tensor([[3.7], [71.1], [7.5]], dtype=torch.float),
          torch.tensor([[3.6], [61.1], [6.5]], dtype=torch.float),
          torch.tensor([[3.5], [51.1], [5.5]], dtype=torch.float),
          torch.tensor([[3.4], [41.1], [4.5]], dtype=torch.float),
          torch.tensor([[3.3], [11.1], [3.5]], dtype=torch.float),
          torch.tensor([[3.3], [21.1], [2.5]], dtype=torch.float),
          torch.tensor([[2.1], [17.1], [3.2]], dtype=torch.float)
          ]

y_list = [1,
          2,
          2,
          2,
          3,
          4,
          1,
          2,
          3,
          2,
          3,
          4,
          4]

train_datalist = [Data(x=x_list[0],
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
                       y=y_list[4])]
val_datalist = [Data(x=x_list[5],
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
                     y=y_list[8])]
test_datalist = [Data(x=x_list[9],
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

batch_size = 3
train_datalaoader = DataLoader(train_datalist, batch_size=batch_size,shuffle=False)
val_dataloader = DataLoader(val_datalist, batch_size=batch_size,shuffle=False)
test_dataloader = DataLoader(test_datalist, batch_size=batch_size,shuffle=False)


class JNet(torch.nn.Module):
    def __init__(self):
        super(JNet, self).__init__()
        self.conv1 = GCNConv(3, 3)  # dataset.num_node_features, 16)  # 3 nodes with 2 features per node?
        self.conv2 = GCNConv(3, 4)  # dataset.num_classes)  # 4 classes of graph

    def forward(self, data):
        x, edge_index = data.x, data.edge_index
        x = self.conv1(x, edge_index)
        x = F.relu(x)
        x = F.dropout(x, training=self.training)
        x = self.conv2(x, edge_index)

        return F.log_softmax(x, dim=1)


def train(model, optimizer, loss_fn, train_loader, val_loader, epochs=20, device="cpu"):
    model = model.to(device)
    for epoch in range(epochs):
        training_loss = 0.0
        valid_loss = 0.0
        model.train()

        for batch in train_loader:
            optimizer.zero_grad()
            print("batch.num_graphs:",batch.num_graphs)
            print("batch.batch:",batch.batch)
            inputs,targets = batch
            inputs = inputs.to(device)
            targets = targets.to(device)
            output = model(inputs)
            loss = loss_fn(output, targets)
            loss.backward()
            optimizer.step()
            training_loss += loss.data.item()
        training_loss /= len(train_loader)

        model.eval()
        num_correct = 0
        num_examples = 0
        for batch in val_loader:
            inputs,targets = batch
            inputs = inputs.to(device)
            output = model(inputs)
            targets = targets.to(device)
            loss = loss_fn(output, targets)
            valid_loss += loss.data.item()
            correct = torch.eq(torch.max(F.softmax(output), dim=1)[1], targets).view(-1)
            num_correct += torch.sum(correct).item()
            num_examples += correct.shape[0]
        valid_loss /= len(val_loader)

        print('Epoch: {},'
              'Training Loss: {:.2f},'
              'Validation Loss: {:.2f},'
              'accuracy = {:.2f'.format(epoch,
                                        training_loss,
                                        valid_loss,
                                        num_correct / num_examples))


model = JNet()
train(model=model,
      optimizer=torch.optim.Adam(model.parameters(), lr=0.01, weight_decay=5e-4),
      loss_fn=torch.nn.NLLLoss,
      train_loader=train_datalaoader,
      val_loader=val_dataloader)
