import numpy as np
from tool.data_processing import prepare_data, one_hot_coding, expand_label
from Input_fucntion import GenerateTaskData

def process_noise(noisesig, trial_num):
    if isinstance(noisesig, np.ndarray):
        noise_sigi = np.repeat(noisesig[:,:,np.newaxis],trial_num,2) #sig_num*2*trial_num
        noise_sigi = np.swapaxes(noise_sigi,0,1)#2*sig_num*trial_num
        noise_sigi = np.reshape(noise_sigi,(2,-1))#2*(trial_num*sig_num)
        return noise_sigi
    else:
        return noisesig
def generate_mask(RT_seqs,max_time):
    mask = np.zeros((RT_seqs.shape[0],max_time))
    for i in range(RT_seqs.shape[0]):
        RT = RT_seqs[i,0]
        mask[i,:RT] = 1
    return mask

def create_dataset_heading2(a, noisesig, StimSeqSize, ChoiceSeqSize, dt, SampSize,HP, mode = 'rand',hsf=0,norm=False,vis_delay=0):
    stim_start = np.ceil(0).astype(int)
    x_datas = []
    y_datas = []
    labels = []
    modalities = []
    noisesigs = []
    headings = []
    StimSeqSizes = []
    H_all = []
    RT_seqs = []
    contexts = []
    output_format = HP['output_format']
    #Generate heading
    if mode == 'rand':
        heading = np.random.rand(a.shape[0], 1) * 2 * np.pi
        #Extend to time dimension
        heading = np.repeat(heading, a.shape[1], axis=1)
        a_all = a
    elif (mode == 'psyco') | (mode == 'timeshift'):
        Hrange = np.arange(-8, 10, 2)
        H_num = len(Hrange)
        a_all = np.repeat(a, H_num, 0)
        for H in Hrange:
            heading = np.ones_like(a) * H / 180 * np.pi
            H_all.append(heading)
        heading = np.concatenate(H_all, 0)
        if mode == 'psyco':
            modality_num = 3
        else:
            modality_num = 2
    #Generate data
    for M in range(modality_num):
        x_data, y_data, label, refvar, SS = GenerateTaskData(a_all, heading, noisesig, StimSeqSize, ChoiceSeqSize, M, dt,
                                                            SampSize, stim_start, 0, 5, heading_shift=hsf,vis_delay=vis_delay)
        if M == 0:
            RT_seq = np.ones((x_data.shape[0], 1))*16
            context = np.concatenate((np.ones((x_data.shape[0],x_data.shape[1], 1)),np.zeros((x_data.shape[0],x_data.shape[1], 1))),axis=2)
        elif M==1:
            RT_seq = np.ones((x_data.shape[0], 1))*26
            context = np.concatenate((np.zeros((x_data.shape[0],x_data.shape[1], 1)),np.ones((x_data.shape[0],x_data.shape[1], 1))),axis=2)
        elif M==2:
            RT_seq = np.ones((x_data.shape[0], 1))*18
            context = np.concatenate((np.ones((x_data.shape[0],x_data.shape[1], 1)),np.ones((x_data.shape[0],x_data.shape[1], 1))),axis=2)
        noise_rep = a_all.shape[0]*SampSize
        noisesigs.append(process_noise(noisesig, noise_rep))
        x_datas.append(x_data)
        y_datas.append(y_data)
        labels.append(label)
        RT_seqs.append(RT_seq)
        contexts.append(context)
        headrep = noisesig.shape[0]*SampSize
        headings.append(np.tile(heading, (headrep, 1)))
        modality = np.repeat(np.ones_like(a_all) * M, headrep, 0)
        modalities.append(modality)
        StimSeqSizes.append(SS)
    if mode == 'timeshift':
        time_set = HP['time_set']
        for i,Ts in enumerate(time_set):
            x_data, y_data, label, refvar, SS = GenerateTaskData(a_all, heading, noisesig, StimSeqSize, ChoiceSeqSize, 3, dt, SampSize,
                                                                stim_start, 0, 5, heading_shift=hsf)#regenerate data from different noise
            if Ts > 0:#delay acceleration
                time_shift = Ts
                shifted = np.roll(x_data[:, :, 0:2], time_shift, axis=1)
                shifted[:, :time_shift, :] = 0
                add_seq = x_data[:]#avoid change original data
                add_seq[:, :, 0:2] = shifted
            else:#delay velocity
                time_shift = -Ts
                shifted = np.roll(x_data[:, :, 2:4], time_shift, axis=1)
                shifted[:, :time_shift, :] = 0
                add_seq = x_data[:]
                add_seq[:, :, 2:4] = shifted
            noisesigs.append(process_noise(noisesig, noise_rep))
            x_datas.append(add_seq)
            y_datas.append(y_data)
            labels.append(label)
            headrep = noisesig.shape[0]*SampSize
            headings.append(np.tile(heading, (headrep, 1)))
            modality = np.repeat(np.ones_like(a_all) * (2 + i), headrep, 0)
            modalities.append(modality)
            StimSeqSizes.append(SS)
    x_datas = np.concatenate(x_datas, axis=0)
    y_datas = np.concatenate(y_datas, axis=0)
    labels = np.concatenate(labels, axis=0)
    modalities = np.concatenate(modalities, axis=0)
    contexts = np.concatenate(contexts, axis=0)
    contextves = contexts[:, :, 0:1]
    contextvis = contexts[:, :, 1:2]
    StimSeqSizes = np.concatenate(StimSeqSizes, axis=0)
    RT_seqs = np.concatenate(RT_seqs, axis=0).astype(int)
    if noisesigs and (len(noisesigs)>2):
        noisesigs = np.concatenate(noisesigs, axis=1)
    headings = np.concatenate(headings, axis=0)
    headings = np.expand_dims(headings, 2)
    y_datas = y_datas[:, :, :6]
    y_datas = np.concatenate((y_datas,headings),axis=2)
    one_hot_modality = one_hot_coding(modalities,3)
    additional_input = None
    input_seq = x_datas
    label_one_hot = one_hot_coding(labels, 2)
    label_one_hot = expand_label(y_datas, label_one_hot, StimSeqSize, 'all_same')
    label_expanded = expand_label(y_datas, labels, StimSeqSize, 'all_same')
    if output_format == 'Choice':
        targets = {
            'target1': (label_one_hot, None),
        }
    elif output_format == 'All':
        targets = {
            'target1': (y_datas, None),
            'target2': (label_one_hot, None),
        }
    elif output_format =='1dchoice':
        targets = {
            'target1': (label_expanded, None),
            # 'target1': (label_expanded_RT, None),
        }
    elif output_format == 'ChoiceContext':
        targets = {
            'target1': (contexts, None),
            # 'target1': (contextves, None),
            # 'target2': (contextvis, None),
            'target3': (label_expanded, None),
        }
    data = {
        'stim_start': stim_start,
        'seq_input': input_seq.copy(),
        'data_target': targets.copy(),
        'additional_input': additional_input,
        'headings': headings,
        'modalities': modalities,
        'noisesigs': noisesigs,
        'refvar': refvar,
        'StimSeqSizes': StimSeqSizes,
    }
    input_data, target_data, mask = prepare_data(data, HP, condition_sep=False)
    mask_RT=generate_mask(RT_seqs,mask.shape[1])
    data['input_data'] = input_data
    data['data_target'] = target_data
    input_size = input_data.shape[2]  # RNN
    output_size = 0
    output_slices = []
    for key in targets:
        # 获取当前 target 的最后一个维度
        current_size = targets[key][0].shape[-1]
        # 计算输出的切片
        current_slice = slice(output_size, output_size + current_size)
        # 将当前切片添加到切片列表中
        output_slices.append(current_slice)
        # 更新 output_size
        output_size += current_size
    return data, input_size, output_size, output_slices
def create_dataset_all(a, noisesig, StimSeqSize, ChoiceSeqSize, dt, SampSize,HP, mode = 'rand',hsf=0,norm=False,vis_delay=0):
    stim_start = np.ceil(0).astype(int)
    x_datas = []
    y_datas = []
    labels = []
    modalities = []
    noisesigs = []
    headings = []
    StimSeqSizes = []
    H_all = []
    output_format = HP['output_format']
    ref_var = HP['ref_var']
    ref_dim = HP['ref_dim']
    modality_num = 3
    #Generate heading
    if mode == 'rand':
        heading = np.random.rand(a.shape[0], 1) * 2 * np.pi
        #Extend to time dimension
        heading = np.repeat(heading, a.shape[1], axis=1)
        a_all = a
    elif (mode == 'psyco') | (mode == 'timeshift'):
        Hrange = np.arange(-8, 10, 2)
        H_num = len(Hrange)
        a_all = np.repeat(a, H_num, 0)
        for H in Hrange:
            heading = np.ones_like(a) * H / 180 * np.pi
            H_all.append(heading)
        heading = np.concatenate(H_all, 0)
        if mode == 'timeshift':
            modality_num = 2
    elif mode == 'distance':
        a_all = a
        heading = np.random.rand(a.shape[0], 1) * 0
        heading = np.repeat(heading, a.shape[1], axis=1)
    #Generate data
    for M in range(modality_num):
        x_data, y_data, label, refvar, SS = GenerateTaskData(a_all, heading, noisesig, StimSeqSize, ChoiceSeqSize, M, dt,
                                                            SampSize, stim_start, ref_var, ref_dim, heading_shift=hsf,vis_delay=vis_delay)
        noise_rep = a_all.shape[0]*SampSize
        noisesigs.append(process_noise(noisesig, noise_rep))
        x_datas.append(x_data)
        y_datas.append(y_data)
        labels.append(label)
        headrep = noisesig.shape[0]*SampSize
        headings.append(np.tile(heading, (headrep, 1)))
        modality = np.repeat(np.ones_like(a_all) * M, headrep, 0)
        modalities.append(modality)
        StimSeqSizes.append(SS)
    if mode == 'timeshift':#shift inout in cominecondition
        time_set = HP['time_set']
        for i,Ts in enumerate(time_set):
            x_data, y_data, label, refvar, SS = GenerateTaskData(a_all, heading, noisesig, StimSeqSize, ChoiceSeqSize, 3, dt, SampSize,
                                                                stim_start, ref_var, ref_dim, heading_shift=hsf)#regenerate data from different noise
            if Ts > 0:#delay acceleration
                time_shift = Ts
                shifted = np.roll(x_data[:, :, 0:2], time_shift, axis=1)
                shifted[:, :time_shift, :] = 0
                add_seq = x_data[:]#avoid change original data
                add_seq[:, :, 0:2] = shifted
            else:#delay velocity
                time_shift = -Ts
                shifted = np.roll(x_data[:, :, 2:4], time_shift, axis=1)
                shifted[:, :time_shift, :] = 0
                add_seq = x_data[:]
                add_seq[:, :, 2:4] = shifted
            noisesigs.append(process_noise(noisesig, noise_rep))
            x_datas.append(add_seq[:])
            y_datas.append(y_data[:])
            labels.append(label[:])
            headrep = noisesig.shape[0]*SampSize
            headings.append(np.tile(heading, (headrep, 1)))
            modality = np.repeat(np.ones_like(a_all) * (2 + i), headrep, 0)
            modalities.append(modality)
            StimSeqSizes.append(SS)
    x_datas = np.concatenate(x_datas, axis=0)
    y_datas = np.concatenate(y_datas, axis=0)
    labels = np.concatenate(labels, axis=0)
    modalities = np.concatenate(modalities, axis=0)
    StimSeqSizes = np.concatenate(StimSeqSizes, axis=0)
    if noisesigs and (len(noisesigs)>2):
        noisesigs = np.concatenate(noisesigs, axis=1)
    headings = np.concatenate(headings, axis=0)
    headings = np.expand_dims(headings, 2)
    y_datas = y_datas[:, :, :6]
    if ref_var == 0:#heading task
        y_datas = np.concatenate((y_datas,headings),axis=2)
    # one_hot_modality = one_hot_coding(modalities,3)#for modality input
    # additional_input = one_hot_modality.copy()
    additional_input = None
    input_seq = x_datas
    label_one_hot = one_hot_coding(labels, 2)
    label_one_hot = expand_label(y_datas, label_one_hot, StimSeqSize, 'all_same')
    if output_format == 'Choice':
        targets = {
            'target2': (label_one_hot, None),
        }
    elif output_format == 'All':
        targets = {
            'target1': (y_datas, None),
            'target2': (label_one_hot, None),
        }
    data = {
        'stim_start': stim_start,
        'seq_input': input_seq.copy(),
        'data_target': targets.copy(),
        'additional_input': additional_input,
        'headings': headings,
        'modalities': modalities,
        'noisesigs': noisesigs,
        'refvar': refvar,
        'StimSeqSizes': StimSeqSizes,
    }
    input_data, target_data, mask = prepare_data(data, HP, condition_sep=False)
    data['input_data'] = input_data
    data['data_target'] = target_data
    input_size = input_data.shape[2]  # RNN
    output_size = 0
    output_slices = []
    for key in targets:
        # 获取当前 target 的最后一个维度
        current_size = targets[key][0].shape[-1]
        # 计算输出的切片
        current_slice = slice(output_size, output_size + current_size)
        # 将当前切片添加到切片列表中
        output_slices.append(current_slice)
        # 更新 output_size
        output_size += current_size
    return data, input_size, output_size, output_slices
