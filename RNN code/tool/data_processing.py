import numpy as np
import torch
from torch.utils.data import TensorDataset, DataLoader
from sklearn.model_selection import train_test_split


def one_hot_coding(labels, num_classes=None):
    """Convert label array with integer values to one-hot array.
    Args:
        labels: Array of integer type label values
        num_classes: Number of classes
    Returns:
        one_hot: Numpy array with one-hot encoding
    """
    labels = np.array(labels, dtype=np.int32)
    if num_classes is None:
        num_classes = np.max(labels) + 1
    else:
        labels[labels>=num_classes] = num_classes-1
    one_hot = np.zeros((labels.shape[0], num_classes))
    for idx, label in enumerate(labels):
        one_hot[idx, label] = 1
    return one_hot
def expand_label(output, label, trial_end=None, method='trial_end'):
    # 获取 trial, timesteps, labeldim
    trial, timesteps, labeldim = output.shape[0], output.shape[1], label.shape[1]
    # 初始化 target
    target = np.zeros((trial, timesteps, labeldim))
    if trial_end is None:
        trial_end = timesteps
    if method == 'trial_end':
        # 把trial_end标记的时间点的值设定的和原标签相同
        if isinstance(trial_end, int):
            target[:, trial_end:,:] = label[:, np.newaxis, :]
        else:
            for i in range(trial):
                target[i, trial_end[i],:] = label[i]
    elif method == 'all_same':
        # 把所有时间点的值设定的和原标签相同
        target = np.repeat(label[:, np.newaxis, :], timesteps, axis=1)
    elif method == 'linear_change':
        # 从在0-trial_end时间内从值从1/labeldim 向原标签均匀变化
        for i in range(trial):
            target[i, :trial_end[i]] = np.linspace(1/labeldim, label[i], trial_end[i]).T
    elif method == 'stimlen':
        if isinstance(trial_end, int):
            target[:, :trial_end, :] = label[:, np.newaxis, :]
        else:
            for i in range(trial):
                target[i, :trial_end[i],:] = label[i]
    else:
        print("Invalid method. Please choose from 'trial_end', 'all_same', or 'linear_change'.")
    return target
def pad_sequence(sequences, max_length=None):
    """
    Pad a list of variable length Tensors with zero
    sequences:  (list of 2D array) (time * feature_dim)
    Padded Sequence:  (3D PyTorch tensor) (sample_size, max_length, feature_dim)
    Mask:  (3D PyTorch tensor) (sample_size, max_length, 1) 1 represent valid elements, and 0 represent padding elements
    """
    sequences = [torch.tensor(data, dtype=torch.float32) for data in sequences]
    lengths = torch.tensor([seq.shape[0] for seq in sequences])
    if max_length is None:
        max_len = lengths.max()
    else:
        max_len = max_length
    Padded_Sequence = torch.zeros(lengths.size(0), max_len, sequences[0].shape[1])
    mask = torch.ones(lengths.size(0), max_len, 1)
    for i, seq in enumerate(sequences):
        end = lengths[i]
        Padded_Sequence[i, :end, :] = seq
        mask[i, end:, 0] = 0
    return Padded_Sequence, mask, max_len



def seg_data(input_data, target_data, additional_input, time_range):
    """
    Segment the input data and target data for each time slice, and concatenate additional info. For CNN-like net

    Args:
        input_data: (2d array) Input data of shape (time_steps, input_features).
        target_data: (2d array) Target data of shape (time_steps, target_features).
        additional_input: (1d array) Additional input data to be concatenated.
        time_range: (tuple) Time range for segmenting the data, specified as (start_time, end_time). input from -start_time to end_time for each time slice.`

    Returns:
        segmented_input: (list) List of segmented input data arrays(trials) of shape (input_features * time_range + additional_features,).
        segmented_target: (list) List of segmented target data arrays(trials) of shape (target_features,).
    """
    segmented_input = []
    segmented_target = []
    start_time, end_time = time_range

    for time in range(start_time, input_data.shape[0] - end_time):
        input_segment = input_data[time - start_time: time + end_time, ...]

        if input_segment.ndim == 2:
            input_segment = input_segment.reshape(-1)

        target_segment = target_data[time, ...]

        # Check for NaN values
        if np.isnan(input_segment).any() or np.isnan(target_segment).any():
            continue

        if additional_input is not None:
            input_segment = np.concatenate((input_segment, additional_input), axis=0)
        segmented_input.append(input_segment)
        segmented_target.append(target_segment)

    return segmented_input, segmented_target


def prepare_sequential_data(input_data, target_data, additional_input, time_offset, stim_start=0, from_stim_start=True):
    """
    Prepare sequential data for RNN-like net.
    Args:
        input_data: (2d array) Input data of shape (time_steps, input_features).
        target_data: (2d array)Target data of shape (time_steps, target_features).
        additional_input: (1d array) Additional input data of shape (additional_features,).
        time_offset: (int) Time offset for target data. time_shift between input and target.
        stim_start: (int) Start time for stimulus (default: 0).
        from_stim_start: (bool) Flag to indicate whether to start from stim_start (default: True).
    Returns:
        prepared_input: (2d array) Prepared input data of shape (selected_steps, input_features + additional_features).
        prepared_target: (2d array) Prepared target data of shape (selected_steps, target_features). Selected_steps is variable here!!!
    """
    if additional_input is not None:
        additional_input_resized = np.tile(additional_input[np.newaxis, :], (input_data.shape[0], 1))
        # Concatenate input and resized additional_input
        input_data = np.concatenate((input_data, additional_input_resized), 1)

    # Determine valid time indices
    time_range = ~np.isnan(input_data).any(axis=1)
    if from_stim_start:
        time_range[:stim_start] = False
    selected_time = np.where(time_range)[0]

    # Calculate target indices and filter out-of-bounds indices
    target_selected_time = selected_time + time_offset
    max_time = target_data.shape[0]
    valid_indices = (target_selected_time >= 0) & (target_selected_time < max_time)
    target_selected_time = target_selected_time[valid_indices]
    selected_time = selected_time[valid_indices]

    # Select input and target data based on valid indices
    prepared_input = input_data[selected_time, :]
    prepared_target = target_data[target_selected_time, :]
    valid_indices = ~np.isnan(prepared_target).any(axis=1)
    prepared_input = prepared_input[valid_indices, :]
    prepared_target = prepared_target[valid_indices, :]
    return prepared_input, prepared_target


def task_data_align(data, HP):
    """
    for my task, get full_seq_input for prepare_data
    """
    pop_activity = data['population_activity'].copy()[0]
    neural_activity = data['neural_activity'].copy()[0]
    variable_data_3d = data['variables'].copy()
    condition = data['condition'].copy()
    use_pop_activity = HP['use_pop_activity']
    use_hist_activity = HP['use_hist_activity']
    full_seq_input = []
    for trial in range(neural_activity.shape[0]):
        conditionT = condition[trial, 0] - 1
        if use_pop_activity:
            datain = np.concatenate((pop_activity[trial, ...], variable_data_3d[conditionT, ...]), axis=-1)
        elif use_hist_activity:
            datain = np.concatenate((neural_activity[trial, ...], variable_data_3d[conditionT, ...]), axis=-1)
        else:
            datain = variable_data_3d[:, conditionT, :]
        full_seq_input.append(datain)
    full_seq_input = np.stack(full_seq_input, axis=0)
    return full_seq_input


def prepare_data(data, HP, condition_sep=False):
    """
    selection and alignment of data.
    Args:
        data: (dict) Dictionary of data.
        HP: (dict) Dictionary of hyperparameters.
    Returns:
        input_data: input data array      (trials) * (selected_steps, input_features + additional_features)
                                                         /(trials*time_steps) * (input_features + additional_features)
        target_data:target data array/dict of target and mask array    (trials) * (selected_steps, target_features)
                                                           / (trials*time_steps) * (output_features)
    """
    input_data = []
    target_data = []
    mask = None
    if data['additional_input'] is not None:
        additional_inputs = data['additional_input'].copy()
    else:
        additional_inputs = None
        additional_input = None
    stim_start = data['stim_start'].copy()
    seq_input = data['seq_input'].copy()
    data_target = data['data_target'].copy()
    use_rnn = HP['use_rnn']
    if isinstance(data_target, dict):
        # Prepare the combined target data
        combined_target_data = np.concatenate([value[0] for value in data_target.values()], axis=-1)
        # Prepare the length of each target
        target_lengths = [value[0].shape[-1] for value in data_target.values()]
    else:
        combined_target_data = data_target.copy()
    for trial in range(combined_target_data.shape[0]):
        if additional_inputs is not None:
            additional_input = additional_inputs[trial, :]
        data_in = seq_input[trial, ...].copy()
        data_tar = combined_target_data[trial, ...].copy()
        if use_rnn:
            trial_in, trial_tar = prepare_sequential_data(data_in, data_tar, additional_input,
                                                          HP['time_offset'], stim_start, True)
            input_data.append(trial_in)
            target_data.append(trial_tar)
        else:
            trial_in, trial_tar = seg_data(data_in, data_tar, additional_input,
                                           HP['time_range'])
            if condition_sep:
                if len(trial_in) > 0:
                    input_data.append(np.stack(trial_in, axis=0))
                    target_data.append(np.stack(trial_tar, axis=0))
            else:
                input_data.extend(trial_in)
                target_data.extend(trial_tar)
    if use_rnn:
        input_data, mask, max_len = pad_sequence(input_data)
        target_data = pad_sequence(target_data,max_len)[0]
    else:
        input_data = np.stack(input_data, axis=0)
        target_data = np.stack(target_data, axis=0)

    if isinstance(data_target, dict):
        # Slice the combined target data back into individual targets and weights
        target_data_dict = data_target.copy()
        start = 0
        for key,length in zip(data_target.keys(), target_lengths):
            end = start + length
            target_data_dict[key] = (target_data[...,start:end],)+target_data_dict[key][1:]
            start = end
        target_data = target_data_dict
    return input_data, target_data, mask


def prepare_conditional_averaged_data(raw_data, condition_labels_list, trial_dim):
    """
    Prepare averaged data for training and testing.
    Args:
        raw_data (nd np array): Raw data with arbitrary dimensions. The dimension of trials is specified by trial_dim.
        condition_labels_list (list of 1d np arrays): List of arrays containing condition labels for each trial.
        trial_dim (integer): The dimension of trials in raw_data.
    Returns:
        averaged_condition_data (nd np array): Averaged data with the same shape as raw_data, but the trial dimension replaced by conditions.
        condition_indices (2d np array): Array containing condition indices for each condition.
    """
    trial_dim = trial_dim - 1  # Convert to 0-based indexing
    num_conditions = [np.max(labels) + 1 for labels in condition_labels_list]

    total_conditions = np.prod(num_conditions)

    avg_dims = list(raw_data.shape)
    avg_dims[trial_dim] = total_conditions
    averaged_condition_data = np.full(avg_dims, np.nan)

    condition_indices = np.full((total_conditions, len(num_conditions)), np.nan)
    curr_condition = np.zeros(len(num_conditions), dtype=int)
    condition_mask = np.zeros((len(num_conditions), raw_data.shape[trial_dim]), dtype=bool)

    # Average data by condition
    for i in range(total_conditions):
        resi = i
        for j in range(len(num_conditions)):
            curr_condition[j] = resi % num_conditions[j]
            condition_mask[j, :] = curr_condition[j] == condition_labels_list[j]
            resi = resi // num_conditions[j]

        condition_masks_combined = np.all(condition_mask, axis=0)
        #select the trials that satisfy all conditions
        avg_slice = [slice(None)] * raw_data.ndim
        avg_slice[trial_dim] = condition_masks_combined
        avg_slice = tuple(avg_slice)
        #select the final index
        avg_cond_slice = [slice(None)] * raw_data.ndim
        avg_cond_slice[trial_dim] = i
        avg_cond_slice = tuple(avg_cond_slice)

        averaged_condition_data[avg_cond_slice] = np.mean(raw_data[avg_slice], axis=trial_dim)

        condition_indices[i, :] = curr_condition

    return averaged_condition_data, condition_indices

def handle_weights(weight, data, device):
    if weight is not None:
        if len(weight.shape) == 1:
            new_shape = (*data.shape[:-1], 1)
            weight = weight*np.ones(new_shape)
        return torch.tensor(weight, dtype=torch.float32).to(device)
    else:
        new_shape = (*data.shape[:-1], 1)
        return torch.ones(new_shape, dtype=torch.float32, device=device)

def generate_dataset_multi_target(input_data, targets, HP):
    device = HP['device']
    batch_size = HP['batch_size']
    test_size = 0.4
    random_state = 42

    # Keep input data as numpy array
    all_arrays = [input_data]
    for target_name, (target, weight) in targets.items():
        weight = handle_weights(weight, target, device)
        all_arrays.extend([target, weight])  # Extend the list with target and weight

    # Create index array
    indices = np.arange(len(input_data))

    # Function to split data into train, valid and test sets
    def split_data(indices, test_size=0.1,valid_size=0.1, random_state=random_state):
        train_valid_indices, test_indices = train_test_split(indices, test_size=test_size, random_state=random_state)
        valid_size_adjusted = valid_size / (1 - test_size)
        train_indices, valid_indices = train_test_split(train_valid_indices, test_size=valid_size_adjusted, random_state=random_state)
        return train_indices, valid_indices, test_indices

    # Split the indices
    train_indices, valid_indices, test_indices = split_data(indices, test_size=test_size, random_state=random_state)

    # Function to create TensorDatasets and DataLoaders
    def generate_dataloader(indices, arrays):
        data = [array[indices] for array in arrays]  # Get data by indices
        data.append(indices)  # Append indices to the data
        tensors = [torch.tensor(d, dtype=torch.float32).to(device) for d in data]
        dataset = TensorDataset(*tensors)
        return DataLoader(dataset, batch_size=batch_size, shuffle=True)

    # Create DataLoaders
    train_loader = generate_dataloader(train_indices, all_arrays)
    valid_loader = generate_dataloader(valid_indices, all_arrays)
    test_loader = generate_dataloader(test_indices, all_arrays)

    return train_loader, valid_loader, test_loader

def generate_dataset(input_data, target_data, HP, weight=None):
    """
    Generate training, validation and test datasets.
    Args:
        input_data: (numpy.ndarray) Input data of shape (sample_size, input_features).
        target_data: (numpy.ndarray) Target data of shape (sample_size, target_features).
        weight: (numpy.ndarray) Mask of shape (sample_size, 1) to indicate valid elements.
    Returns:
        train_loader: (torch.utils.data.DataLoader) Training dataset loader.
        valid_loader: (torch.utils.data.DataLoader) Validation dataset loader.
        test_loader: (torch.utils.data.DataLoader) Test dataset loader.
    """
    device = HP['device']
    batch_size = HP['batch_size']
    random_state = 42
    input_data = torch.tensor(input_data, dtype=torch.float32).clone().detach().to(device)
    target_data = torch.tensor(target_data, dtype=torch.float32).clone().detach().to(device)
    if weight is not None:
        weight = torch.tensor(weight, dtype=torch.float32).clone().detach().to(device)
    else:
        input_data_shape = input_data.shape
        new_shape = (*input_data_shape[:-1], 1)
        weight = torch.ones(new_shape, dtype=input_data.dtype, device=device)
    x_train, x_test, y_train, y_test, mask_train, mask_test = train_test_split(
        input_data, target_data, weight, test_size=0.1, random_state=random_state
    )
    x_train, x_valid, y_train, y_valid, mask_train, mask_valid = train_test_split(
        x_train, y_train, mask_train, test_size=0.1, random_state=random_state
    )
    train_dataset = TensorDataset(x_train, y_train, mask_train)
    valid_dataset = TensorDataset(x_valid, y_valid, mask_valid)
    test_dataset = TensorDataset(x_test, y_test, mask_test)
    train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
    valid_loader = DataLoader(valid_dataset, batch_size=batch_size, shuffle=False)
    test_loader = DataLoader(test_dataset, batch_size=batch_size, shuffle=False)
    return train_loader, valid_loader, test_loader

def loader2array(loader):
    data_array = []
    labels_array = []
    mask_array = []
    if len(loader.dataset.tensors) == 3:  # Check if mask are included
        for batch_data, batch_labels, batch_mask in loader:
            data_array.append(batch_data.cpu().numpy())
            labels_array.append(batch_labels.cpu().numpy())
            mask_array.append(batch_mask.cpu().numpy())
        mask_array = np.concatenate(mask_array, axis=0)
    else:
        for batch_data, batch_labels in loader:
            data_array.append(batch_data.cpu().numpy())
            labels_array.append(batch_labels.cpu().numpy())
    # Concatenate the lists to get arrays
    data_array = np.concatenate(data_array, axis=0)
    labels_array = np.concatenate(labels_array, axis=0)
    return data_array,labels_array,mask_array

