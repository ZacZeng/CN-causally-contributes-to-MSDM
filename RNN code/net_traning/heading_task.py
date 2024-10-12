import os
import numpy as np
from scipy.io import savemat
from tool.visulization_function import plot_model_predictions
from net_traning.generate_dataset import create_dataset_heading2
import torch
import torch.nn as nn
from simulationNN import train_model, MainNet, CustomRNN, MLP
from tool.data_processing import generate_dataset_multi_target
import torch.optim as optim
from Input_fucntion import a_profile


def su_simulation_train(model, HP, data, Hp_log, file_name,model_path):
    RNNweights = []
    Accs = []
    sigs = []
    # data prepare
    targets = data['data_target'].copy()
    input_data = data['input_data']
    headings = data['headings']
    modalities = data['modalities']
    train_loader, valid_loader, test_loader = generate_dataset_multi_target(input_data, targets, HP)
    log_name = file_name + '_' + Hp_log
    # model training
    for sub_epoch in range(300):
        avg_test_loss, model = train_model(model, train_loader, test_loader, HP, log_name)
        model_path_full = model_path[:-4] + '_full.pth'
        torch.save(model.state_dict(), model_path)
        torch.save(model, model_path_full)
        data_show = next(iter(test_loader))
        input_data = data_show[0]
        end_time = SeqSize-1
        data_indice = data_show[-1].cpu().numpy().astype(int)
        heading = headings[data_indice,0,0] # target heading
        indices = modalities[data_indice,0].astype(int)
        label_names = ['vestibular', 'visual', 'combined']
        acc,props = plot_model_predictions(model, input_data, end_time, heading, labels=indices, label_names=label_names,
                               bandwidth=None, choice_type='1d')
        RNNweights.append(model.layers[0][0].Whh.cpu().weight.detach().numpy())
        Accs.append(acc)
        sigs.append(props[:,1])
    return model,RNNweights,Accs,sigs
# 超参数定义
device = torch.device('cuda' if torch.cuda.is_available else 'cpu')
#input def
StimTime = 1.5
ChoiceTime = 0.5
dt = 0.1
dt = 0.05
SeqSize = int(np.ceil((StimTime + ChoiceTime)/dt))
StimSeqSize = int(np.ceil(StimTime/dt))
ChoiceSeqSize = int(np.ceil(ChoiceTime/dt))
SampSize = 6000
#herperparameters
time_offset = 0 # offset between input and target,>0 means prediction
lr = 1e-5
# Define the loss function and optimizer
HP = {
    'time_offset': time_offset,
    'num_epochs': 1,
    'batch_size': 1024,
    'L1_lambda': 0,
    'L2_lambda': 0,
    'device': device,
    'use_hist_activity': False,
    'use_pop_activity': False,
    'time_range': (2, 0),
    'use_rnn': True,
    'la':0,
    'time_set':[0],
    'output_format':'1dchoice',
}
#traning_dataset def
filename = '../decoding/m18c38r2_Spike.mat'
a = a_profile('heading', 1, StimSeqSize, ChoiceSeqSize)
noisesigset = np.array([[0.4,0.4]])
vis_delay = 5
data, input_size, output_size, output_slices = create_dataset_heading2(a, noisesigset, StimSeqSize, ChoiceSeqSize, dt, SampSize,HP,'psyco',vis_delay=vis_delay)
targets = data['data_target'].copy()
input_data = data['input_data']
headings = data['headings']
modality = data['modalities']
#net structure define
neu_num = 256
hidden_sizes = [32]
net1 = CustomRNN(4, neu_num,0.1,0)
net2 = MLP(neu_num, [], output_size)
model = MainNet(([net1],[net2]), device)
model_name = 'RNN_NMI_A4S_Vdelay_N0.40.4_H_TPNM_N256n0.1dt0.05_C1d'
load_model_path = '../PreTrainPolarN256n0.01.pth'
model_path = model_name + '.pth'
model_path_full = model_name+ '_full.pth'
HP['loss_function'] = [nn.BCEWithLogitsLoss(reduction='none')]
HP['optimizer'] = optim.Adam(model.parameters(), lr=lr)
HP['output_slice'] = output_slices
model,RNNweights,Accs,sigs = su_simulation_train(model, HP, data, 'test', 'file_name',model_path)
data_dict = {'Accs':Accs}
file_name = model_name + 'Training.npz'
file_name_mat = model_name + 'Training.mat'
save_dir = './Training_results'
sigs = np.array(sigs)
colors = ['b','r','g']
if not os.path.exists(save_dir):
    os.makedirs(save_dir)
np.savez(save_dir + '/' + file_name, **data_dict)
savemat(os.path.join(save_dir, file_name_mat), data_dict)
file_name = model_name + 'Trainingsigs.npz'
file_name_mat = model_name + 'Trainingsigs.mat'
data_dict = {'sigs':sigs}
np.savez(save_dir + '/' + file_name, **data_dict)
savemat(os.path.join(save_dir, file_name_mat), data_dict)
