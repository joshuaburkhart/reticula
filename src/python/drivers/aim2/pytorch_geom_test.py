# testing install from https://pytorch-geometric.readthedocs.io/en/latest/notes/installation.html
# https://pytorch.org/get-started/locally/#macos-version
# https://www.jetbrains.com/help/pycharm/conda-support-creating-conda-virtual-environment.html

from __future__ import print_function

import torch
import torch.nn.functional as F
from torch.nn import Sequential, Linear, ReLU
from torch_geometric.data import Data, DataLoader
from torch_geometric.nn import GCNConv
from torch_geometric.nn import GINConv

edge_index = torch.tensor([[0, 1],
                           [1, 0],
                           [1, 2],
                           [2, 1]],
                          dtype=torch.long)

x_list = [torch.tensor([[1.3,2], [1.1,2], [2.9,2]], dtype=torch.float), # 1
          torch.tensor([[2.3,1], [2.3,2], [2.2,2]], dtype=torch.float), # 0
          torch.tensor([[2.3,2], [2.1,1], [2.5,1]], dtype=torch.float), # 0
          torch.tensor([[2.3,4], [2.1,2], [2.5,2]], dtype=torch.float), # 0
          torch.tensor([[1.3,2], [1.1,2], [2.5,2]], dtype=torch.float), # 1
          torch.tensor([[2.8,2], [2.1,1], [2.5,1]], dtype=torch.float), # 0
          torch.tensor([[1.7,6], [1.1,2], [2.5,2]], dtype=torch.float), # 1
          torch.tensor([[2.6,2], [2.1,2], [2.5,2]], dtype=torch.float), # 0
          torch.tensor([[1.5,2], [1.1,2], [2.5,1]], dtype=torch.float), # 1
          torch.tensor([[2.4,1], [2.1,2], [2.5,2]], dtype=torch.float), # 0
          torch.tensor([[1.3,2], [1.1,2], [2.5,2]], dtype=torch.float), # 1
          torch.tensor([[1.3,2], [1.1,2], [2.5,1]], dtype=torch.float), # 1
          torch.tensor([[1.1,8], [1.1,2], [2.2,2]], dtype=torch.float), # 1
          torch.tensor([[2.1,2], [2.1,4], [2.2,2]], dtype=torch.float)  # 0
          ]

y_list = [torch.tensor([1,1,1], dtype=torch.long), # 1
          torch.tensor([0,0,0], dtype=torch.long), # 0
          torch.tensor([0,0,0], dtype=torch.long), # 0
          torch.tensor([0,0,0], dtype=torch.long), # 0
          torch.tensor([1,1,1], dtype=torch.long), # 1
          torch.tensor([0,0,0], dtype=torch.long), # 0
          torch.tensor([1,1,1], dtype=torch.long), # 1
          torch.tensor([0,0,0], dtype=torch.long), # 0
          torch.tensor([1,1,1], dtype=torch.long), # 1
          torch.tensor([0,0,0], dtype=torch.long), # 0
          torch.tensor([1,1,1], dtype=torch.long), # 1
          torch.tensor([1,1,1], dtype=torch.long), # 1
          torch.tensor([1,1,1], dtype=torch.long), # 1
          torch.tensor([0,0,0], dtype=torch.long)] # 0

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
                      y=y_list[12]),
                 Data(x=x_list[13],
                      edge_index=edge_index.t().contiguous(),
                      y=y_list[13])]

batch_size = 1
train_datalaoader = DataLoader(train_datalist, batch_size=batch_size, shuffle=False)
val_dataloader = DataLoader(val_datalist, batch_size=batch_size, shuffle=False)
test_dataloader = DataLoader(test_datalist, batch_size=batch_size, shuffle=False)

dataset_num_node_features = 2
dataset_num_classes = 2

class JNet(torch.nn.Module):
    def __init__(self):
        super(JNet, self).__init__()
        nn1 = Sequential(Linear(dataset_num_node_features, 2), ReLU(), Linear(2, 16))
        nn2 = Sequential(Linear(16, 2), ReLU(), Linear(2, dataset_num_classes))
        self.conv1 = GINConv(nn1)
        self.conv2 = GINConv(nn2)

    def forward(self, data):
        x, edge_index = data.x, data.edge_index
        x = self.conv1(x, edge_index)
        x = F.relu(x)
        x = F.dropout(x, training=self.training)
        x = self.conv2(x, edge_index)
        ret = F.log_softmax(x, dim=1)
        return ret


def train(model, optimizer, loss_fn, train_loader, val_loader, epochs=100, device="cpu"):
    model = model.to(device)
    for epoch in range(epochs):
        training_loss = 0.0
        valid_loss = 0.0
        model.train()
        for batch in train_loader:
            optimizer.zero_grad()
            output = model(batch)
            loss = loss_fn(output,batch.y)
            loss.backward()
            optimizer.step()
            training_loss += loss.data.item()
        training_loss /= len(train_loader.dataset)

        model.eval()
        num_correct = 0
        num_examples = 0
        for batch in val_loader:
            output = model(batch)
            loss = loss_fn(output, batch.y)
            valid_loss += loss.data.item() * batch.num_graphs
            t_max = torch.argmax(output.exp(),dim=1)
            t_eq = torch.eq(t_max, batch.y)
            t_view = t_eq.view(-1)
            correct = t_view
            num_correct += torch.sum(correct).item()
            num_examples += correct.shape[0]
        valid_loss /= len(val_loader.dataset)

        print('Epoch: {}, '
              'n correct: {},'
              'n examples: {},'
              'Training Loss: {:.2f}, '
              'Validation Loss: {:.2f}, '
              'accuracy = {:.2f}'.format(epoch,
                                         num_correct,
                                         num_examples,
                                        training_loss,
                                        valid_loss,
                                        num_correct / num_examples))


model = JNet()
train(model=model,
      optimizer=torch.optim.Adam(model.parameters(),lr=0.03),
      loss_fn=torch.nn.NLLLoss(),
      train_loader=train_datalaoader,
      val_loader=val_dataloader)

for batch in test_dataloader:
    output = model(batch)
    #print("output.exp():",output.exp())
    print("predct:",torch.argmax(output.exp(),dim=1))
    print("labels:",batch.y)
