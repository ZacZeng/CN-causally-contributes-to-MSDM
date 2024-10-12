import torch.nn.init as init
import torch
import torch.nn as nn
import os
class MainNet(nn.Module):
    def __init__(self, main_list, dev, gap_connections=[]):
        super(MainNet, self).__init__()
        self.layers = nn.ModuleList([nn.ModuleList(subnets).to(dev) for subnets in main_list])
        self.activation = nn.ReLU()
        self.gap_connections = gap_connections  # list of tuples
        self.input_dims = [[subnet.input_dim for subnet in subnets] for subnets in main_list]
        self.output_dims = [[subnet.output_dim for subnet in subnets] for subnets in main_list]
        self.layer_activations = {}
        self.layer_hiddens = {}
        self.device = dev
    def forward(self, x):
        gap_outputs = [None] * len(self.gap_connections)
            # to save outputs for gap connections
        # gidx = 0
        for list_index, subnets in enumerate(self.layers):
            adjusted_input_dims = self.input_dims[list_index].copy()
            input_left = None
            total_input_dim = sum(adjusted_input_dims)
            # adjust input_dims based on gap_connections
            for source_list_index, source_net_index, target_list_index, target_net_index in self.gap_connections:
                if target_list_index == list_index:
                    adjusted_input_dims[target_net_index] -= self.output_dims[source_list_index][source_net_index]
                    total_input_dim = sum(adjusted_input_dims)
            if total_input_dim < x.shape[-1]:
                input_left = x[:, :, total_input_dim:]
                x = x[:, :, :total_input_dim]
            split_x = list(torch.split(x, adjusted_input_dims, dim=-1))
            for gap_index,(source_list_index, source_net_index, target_list_index, target_net_index) in enumerate(self.gap_connections):
                if target_list_index == list_index:
                    split_x[target_net_index] = torch.cat([split_x[target_net_index], gap_outputs[gap_index]], dim=-1)
            # if list_index>0:
            #     split_x = [self.activation(sx) for sx in split_x]  # activation input
            output = [subnet(sx) for subnet, sx in zip(subnets, split_x)]#compute output
            for subnet_index, subnet in enumerate(subnets):
                if isinstance(subnet, CustomRNN):
                    self.layer_hiddens[(list_index, subnet_index)] = subnet.hiddens  # Save the hidden states of the RNNs
            if input_left is not None:
                output.append(input_left)
            self.layer_activations[list_index] = output  # Save the activations of the current layer
            # save output for gap connections
            for gap_index,(source_list_index, source_net_index, target_list_index, target_net_index) in enumerate(self.gap_connections):
                if source_list_index == list_index:
                    gap_outputs[gap_index] = output[source_net_index]
                    # gidx += 1
            x = torch.cat(output, dim=-1)
        return x
    def patiral_forward(self, x, li,gap_outputs):
        selected_layers = self.layers[li[0]:li[1]]
        # gidx = len(gap_outputs)
        for list_index, subnets in enumerate(selected_layers,start=li[0]):
            adjusted_input_dims = self.input_dims[list_index].copy()
            input_left = None
            total_input_dim = sum(adjusted_input_dims)
            # adjust input_dims based on gap_connections
            for source_list_index, source_net_index, target_list_index, target_net_index in self.gap_connections:
                if target_list_index == list_index:
                    adjusted_input_dims[target_net_index] -= self.output_dims[source_list_index][source_net_index]
                    total_input_dim = sum(adjusted_input_dims)
            if total_input_dim < x.shape[-1]:
                input_left = x[:, :, total_input_dim:]
                x = x[:, :, :total_input_dim]
            split_x = list(torch.split(x, adjusted_input_dims, dim=-1))
            for gap_index,(source_list_index, source_net_index, target_list_index, target_net_index) in enumerate(self.gap_connections):
                if target_list_index == list_index:
                    split_x[target_net_index] = torch.cat([split_x[target_net_index], gap_outputs[gap_index]], dim=-1)
            output = [subnet(sx) for subnet, sx in zip(subnets, split_x)]
            if input_left is not None:
                output.append(input_left)
            for gap_index,(source_list_index, source_net_index, target_list_index, target_net_index) in enumerate(self.gap_connections):
                if source_list_index == list_index:
                    gap_outputs[gap_index] = output[source_net_index]
                    # gidx += 1
            x = torch.cat(output, dim=-1)
        return x, gap_outputs
    def set_all_hidden_reset_flag(self, reset_flag):
        for subnets in self.layers:
            for subnet in subnets:
                if isinstance(subnet, CustomRNN):
                    subnet.reset_hiddens = reset_flag

    def get_layer_activations(self, layer_index=None):
        if layer_index is None:
            return self.layer_activations
        else:
            return self.layer_activations.get(layer_index, None)
def load_model_weights(model, model_path):
    if os.path.exists(model_path):
        try:
            model.load_state_dict(torch.load(model_path), strict=True)
            print("Loaded all weights successfully.")
        except RuntimeError as e:
            print("Error in loading all weights: ", str(e))
            print("Loading matching weights...")
            model.load_state_dict(torch.load(model_path), strict=False)
            print("Loaded matching weights successfully.")
    else:
        print("No weights found at specified path.")
class MLP(nn.Module):
    def __init__(self, input_size, hidden_sizes, output_size, activation=nn.ReLU(),
                 dropout=0, use_bn=False, device='cuda'):
        super(MLP, self).__init__()
        layers = []
        layer_sizes = [input_size] + hidden_sizes + [output_size]
        self.input_dim = input_size
        self.output_dim = output_size
        self.activation = activation
        for i in range(len(layer_sizes) - 1):
            layers.append(nn.Linear(layer_sizes[i], layer_sizes[i + 1]))
            if use_bn:
                layers.append(nn.BatchNorm1d(layer_sizes[i + 1]))
            if i < len(layer_sizes) - 2:
                layers.append(self.activation)
                if dropout > 0:
                    layers.append(nn.Dropout(p=dropout))
        self.net = nn.Sequential(*layers)
        self.activations = []  # List to save the activations of each layer
        self.device = device
    def _init_weights(self):
        for m in self.modules():
            if isinstance(m, nn.Linear):
                # 使用Kaiming/He初始化进行权重初始化
                init.kaiming_normal_(m.weight, mode='fan_out', nonlinearity='relu')
                # 或者使用Xavier/Glorot初始化
                # init.xavier_normal_(m.weight)
                if m.bias is not None:
                    # 将偏置初始化为0
                    init.constant_(m.bias, 0)
    def set_weights(self, layer_index, weights, biases=None):
        with torch.no_grad():
            self.net[layer_index].weight.copy_(torch.from_numpy(weights))
            if biases is not None:
                self.net[layer_index].bias.copy_(torch.from_numpy(biases))
    def forward(self, x):
        self.activations = []  # Reset the activations list
        for layer in self.net:
            x = layer(x)
            if isinstance(layer, nn.Linear):
                self.activations.append(x.detach())  # Save the output of the linear layer
        return x
    def get_hidden_activations(self, layer_index=None):
        """Return the activations of a specific hidden layer or all layers if layer_index is None."""
        if layer_index is not None:
            if layer_index < len(self.activations):
                return self.activations[layer_index]
            else:
                return None
        else:
            return self.activations
class LSTModel(nn.Module):
    def __init__(self, input_size, hidden_size, output_size):
        super(LSTModel, self).__init__()
        self.input_dim = input_size
        self.output_dim = output_size
        self.hidden_size = hidden_size
        self.num_layers = 1
        self.lstm = nn.LSTM(input_size, hidden_size, 1, batch_first=True)
        self.fc = nn.Linear(hidden_size, output_size)
    def forward(self, x):
        if x.dim() == 3:  # Batch mode
            batch_size = x.size(0)
            h0 = torch.zeros(self.num_layers, batch_size, self.hidden_size).to(x.device)
            c0 = torch.zeros(self.num_layers, batch_size, self.hidden_size).to(x.device)
            time_dim = 1
        else:  # Non-batch mode
            h0 = torch.zeros(self.num_layers, self.hidden_size).to(x.device)
            c0 = torch.zeros(self.num_layers, self.hidden_size).to(x.device)
            time_dim = 0
        out, _ = self.lstm(x, (h0.detach(), c0.detach()))
        out_sequence = torch.stack([self.fc(out_timestep) for out_timestep in out.unbind(time_dim)], dim=time_dim)
        return out_sequence
class CustomRNN(nn.Module):
    def __init__(self, input_size, hidden_size, noise_stddev = 0, alpha=0.1, reset_hiddens=True, input_units = None, output_units = None):
        super(CustomRNN, self).__init__()
        self.input_dim = input_size
        self.hidden_size = hidden_size
        self.output_dim = hidden_size
        self.noise_stddev = noise_stddev
        # Define the RNN parameters (weights and biases)
        self.Wxh = nn.Linear(input_size, hidden_size, bias=False)
        init.kaiming_normal_(self.Wxh.weight, mode='fan_out', nonlinearity='tanh')
        self.Whh = nn.Linear(hidden_size, hidden_size)
        init.kaiming_normal_(self.Whh.weight, mode='fan_out', nonlinearity='tanh')
        self.alpha = alpha
        self.hiddens = None  # Store all hidden states
        self.h0 = None
        self.reset_hiddens = reset_hiddens
        if input_units is None:
            input_units = hidden_size
        self.input_units = input_units
        self.output_units = hidden_size
        self.Wxh = nn.Linear(input_size, self.input_units, bias=False)
    def forward(self, x, h0 = None):
        #take input with shape (batch_size, seq_len, input_size)
        if h0 is None:
            if self.reset_hiddens:
                h0 = torch.zeros(x.size(0), self.hidden_size, device=self.Wxh.weight.device)
            else:
                h0 = self.h0
        h = h0
        hall = []
        # Forward pass through the RNN layers, return all hidden state
        for i in range(x.size(1)):
            input_step = x[:, i, :]
            noise = torch.randn_like(h) * self.noise_stddev
            h = torch.tanh(self.Wxh(input_step) + self.Whh(h) + noise)
            hall.append(h)
        self.h0 = h
        hiddens = torch.stack(hall,axis = 1)
        self.hiddens = hiddens  # Store all hidden states
        return hiddens
    def stable_loss(self):
        alpha = self.alpha
        stability_loss = torch.sum(torch.pow(self.hiddens[:, 1:, :] - self.hiddens[:, :-1, :], 2))
        return alpha * stability_loss
    def activity_loss(self):
        return torch.sum(torch.pow(self.hiddens, 2))
class IdenticalMapping(nn.Module):
    def __init__(self, input_dim):
        super(IdenticalMapping, self).__init__()
        self.input_dim = input_dim
        self.output_dim = input_dim
        self.hidden_size = input_dim

    def forward(self, x):
        return x

def loss_calculation(data, model, loss_functions, output_slice):
    stability_losses = 0
    activity_losses = 0
    inputs = data[0]
    targets_weights = data[1:]
    outputs = model(inputs)
    total_loss = 0
    for i in range(0, len(targets_weights)-1, 2):  # Step size is 2 because each target has a corresponding weight
        target = targets_weights[i]
        weight = targets_weights[i+1]
        output = outputs[...,output_slice[i // 2]]
        weight = weight.reshape(-1, )
        output = output.reshape(-1, output.shape[-1])
        target = target.reshape(-1, target.shape[-1])
        if len(weight.unique()) <= 2:
            weight = weight.bool()
            loss = torch.sum(loss_functions[i//2](output[weight,], target[weight,]))
        else:
            weight = weight.unsqueeze(1)
            loss = torch.sum(loss_functions[i//2](output, target) * weight)
        total_loss += loss
    items = torch.sum(weight)
    if hasattr(model, 'layers'):
        for subnets in model.layers:
            for subnet in subnets:
                if isinstance(subnet, CustomRNN):
                    stability_losses += subnet.stable_loss()
                    activity_losses += subnet.activity_loss()
    return total_loss, stability_losses, activity_losses, items
def loss_calculation_back(data, model, loss_functions, device='cpu'):
    stability_losses = 0
    if len(data) == 3:  # Check if sample weights are included
        inputs, targets, weights = data
    else:
        inputs, targets = data
        weights = torch.ones(*targets.shape[:-1], 1).to(device)
    outputs = model(inputs)
    total_loss = 0
    for loss_function, slice_obj in loss_functions:
        output_slice = outputs[..., slice_obj]
        target_slice = targets[..., slice_obj]
        output_slice = output_slice.reshape(-1, outputs.shape[-1])
        target_slice = target_slice.reshape(-1, targets.shape[-1])
        weights = weights.reshape(-1, )
        if len(weights.unique()) == 2:
            weights = weights.bool()
            loss = torch.sum(loss_function(output_slice[weights,], target_slice[weights,]))
        else:
            weights = weights.unsqueeze(1)
            loss = torch.sum(loss_function(output_slice, target_slice) * weights)
        total_loss += loss
    items = torch.sum(weights)
    for subnets in model.layers:
        for subnet in subnets:
            if isinstance(subnet, CustomRNN):
                stability_loss = subnet.stable_loss()  # Calculate the stability loss
                stability_losses += stability_loss  # Add the stability loss to the total loss
    return total_loss, stability_losses, items
def train_model(model, train_loader, test_loader, HP, log_name):
    # log_dir = 'logs/' + log_name
    # writer = SummaryWriter(log_dir)
    num_epochs = HP['num_epochs']
    loss_function = HP['loss_function']
    optimizer = HP['optimizer']
    L2_lambda = HP['L2_lambda']
    L1_lambda = HP['L1_lambda']
    la = HP['la']
    device = HP['device']
    output_slice = HP['output_slice']
    fixnetindexs = HP.get('fixnetindexs', None)
    zeroindices = HP.get('zeroindices', None)
    model.to(device)
    for epoch in range(num_epochs):
        epoch_loss = 0.0
        test_loss = 0.0
        item_num = 0
        for data in train_loader:
            if fixnetindexs is not None:
                for index, indices in zip(fixnetindexs, zeroindices):
                    set_weights_to_zero(model, index, indices)
            loss,stability_losses,activity_losses,items = loss_calculation(data,model,loss_function,output_slice)
            epoch_loss += loss.item()
            item_num += items
            L1_reg = 0
            L2_reg = 0
            for subnets in model.layers:
                for subnet in subnets:
                    if isinstance(subnet, CustomRNN):
                        L1_reg += sum(torch.norm(param, 1) for param in subnet.parameters())
                        L2_reg += sum(torch.norm(param, 2) for param in subnet.parameters())
            loss += L2_lambda * L2_reg * len(train_loader.dataset)+stability_losses + L1_reg*L1_lambda + activity_losses*la
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
        avg_epoch_loss = epoch_loss / item_num
        # writer.add_scalar('Loss/Train', avg_epoch_loss, epoch + 1)
        print('Epoch %d/%d, Loss: %f' % (epoch + 1, num_epochs, avg_epoch_loss))
        item_num = 0
        for data in test_loader:
            loss,stability_losses,activity_losses,items = loss_calculation(data,model,loss_function,output_slice)
            test_loss += loss.item()
            item_num += items
        avg_test_loss = test_loss / item_num
        print('Epoch %d/%d, Test Loss: %f' % (epoch + 1, num_epochs, avg_test_loss))
    #     writer.add_scalar('Loss/Test', avg_test_loss, epoch + 1)
    #     for name, param in model.named_parameters():
    #         writer.add_histogram(name, param.clone().cpu().data.numpy(), epoch)
    # writer.close()
    # model_name = 'model_' + log_name + '.pth'
    #torch.save(model.state_dict(), model_name)

    return avg_test_loss, model

def set_weights_to_zero(model, index, indices):
    # 获取指定子网络的权重
    weight = model.layers[index[0]][index[1]].weight.data

    # 将指定的权重设置为零
    weight[indices] = 0.0
def load_model_weights(model, model_path):
    if os.path.exists(model_path):
        try:
            model.load_state_dict(torch.load(model_path), strict=True)
            print("Loaded all weights successfully.")
        except RuntimeError as e:
            print("Error in loading all weights: ", str(e))
            print("Loading matching weights...")
            # Create a new state dict where layers dimensions match with the model
            state_dict = torch.load(model_path)
            new_state_dict = {k: v for k, v in state_dict.items() if
                              k in model.state_dict() and v.shape == model.state_dict()[k].shape}

            # Update the current model state dict with the new one
            model_dict = model.state_dict()
            model_dict.update(new_state_dict)

            # Load the new state dict
            model.load_state_dict(model_dict)
            print("Loaded matching weights successfully.")
    else:
        print("No weights found at specified path.")
    return model


def freeze_model_weights(model):
    # Freeze weights in net1
    for net in model.layers[0]:
        for param in net.parameters():
            param.requires_grad = False
    # Freeze weights in net2
    for net in model.layers[1]:
        for param in net.parameters():
            param.requires_grad = False
