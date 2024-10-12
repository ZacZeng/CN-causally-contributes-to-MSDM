import numpy as np
import torch


def predict_choice(NNmodel, input, time=None):
    output = NNmodel(input)
    if output.size(2) >= 2:
        labelout = output[:, :, -2:]
        probabilities = torch.nn.functional.softmax(labelout, dim=2)
        predicted_choice = torch.argmax(probabilities, dim=2)
    else:
        predicted_choice = output>0

    # Check if time is an array
    if isinstance(time, torch.Tensor) or isinstance(time, np.ndarray):
        if input.size(0) != len(time):
            raise ValueError("Length of input and time must be equal when time is an array.")
        output_choice = np.zeros(len(time))
        for i in range(output.size(0)):
            output_choice[i] = predicted_choice[i, time[i]].cpu()
    elif time is None:
        output_choice = predicted_choice.detach().cpu().numpy()
    else:  # If time is a single value, use the original code
        output_choice = predicted_choice[:, time].detach().cpu().numpy()
    output_choice = output_choice.astype(int).squeeze()
    return output_choice
def predict_1dchoice(NNmodel, input, time=None):
    output = NNmodel(input)
    predicted_choice = output[:, :, -1:]>0
    # Check if time is an array
    if isinstance(time, torch.Tensor) or isinstance(time, np.ndarray):
        if input.size(0) != len(time):
            raise ValueError("Length of input and time must be equal when time is an array.")
        output_choice = np.zeros(len(time))
        for i in range(output.size(0)):
            output_choice[i] = predicted_choice[i, time[i]].cpu()
    elif time is None:
        output_choice = predicted_choice.detach().cpu().numpy()
    else:  # If time is a single value, use the original code
        output_choice = predicted_choice[:, time].detach().cpu().numpy()
    output_choice = output_choice.astype(int).squeeze()
    return output_choice

