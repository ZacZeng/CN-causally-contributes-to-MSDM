import torch
import os
import numpy as np
from tool.net_validation import predict_choice, predict_1dchoice
from visulization.plot_function import plot_psychophysical_function
os.environ['CUDA_LAUNCH_BLOCKING'] = "1"

def plot_model_predictions(NNmodel, input, time, variable, labels=None, label_names=None, label_color=None, bandwidth=None,refvar=0,choice_type = 'norm'):
    # 使用神经网络模型进行预测
    if choice_type == '1d':
        predicted_choice = predict_1dchoice(NNmodel, input, time)
    else:
        predicted_choice = predict_choice(NNmodel, input, time)
    output = NNmodel(input).cpu().detach().numpy()
    if output.shape[2] > 2:
        variable_encoded_t = np.zeros((variable.shape[0],))
        variable_encoded = output[:,:,:-2]
        if isinstance(time, np.ndarray):
            for i in range(variable_encoded.shape[0]):
                variable_encoded_t[i] = variable_encoded[i, time[i], -1]
        else:
            variable_encoded_t = variable_encoded[:, time, -1]
    # 将预测结果转换为numpy数组，以便在绘图函数中使用
    # 调用plot_psychophysical_function函数进行绘图
    variable = variable.cpu().numpy() if isinstance(variable, torch.Tensor) else variable
    _,props,acc = plot_psychophysical_function(variable, predicted_choice, labels = labels, label_names=label_names, label_color=label_color, bandwidth=bandwidth,refvar=refvar)
    return acc,props
