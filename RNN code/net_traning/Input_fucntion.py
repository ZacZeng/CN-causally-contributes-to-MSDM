import numpy as np
import math
from scipy.special import erf
def bool_to_01(labels1):
    l = np.zeros(labels1.shape)
    l[labels1] = 1
    return l
def AN(inp, InputChanNum, noise_sig=1e-2, trial_rep=1):
    repeat_times = (trial_rep,) + (1,) * (inp.ndim - 1)
    # trial repetition
    inp = np.tile(inp, repeat_times)
    # multi-channel input
    inp_ext = np.repeat(inp[..., None], InputChanNum, axis=-1)
    noise = np.random.normal(0, noise_sig, inp_ext.shape)
    #    noise = np.random.normal(0, noise_sig * (1+np.sqrt(np.abs(inp_ext))), inp_ext.shape)
    return inp_ext + noise
def calculate_pos(t, ampl, duration, num_sigs):
    pos = ampl * 0.5 * (erf((t - duration/2) / ((duration/2/num_sigs) * math.sqrt(2))) + 1)
    return pos
def generate_acc_matrix(trial_num, time_steps):
    matrix = []
    for _ in range(trial_num):  # 对于每个trial
        while True:
            # 分别随机确定加速、匀速和减速三个阶段的时间步数
            accelerate_time = np.random.randint(1, time_steps - 3)
            decelerate_time = np.random.randint(1, time_steps - accelerate_time)
            steady_time = time_steps - accelerate_time - decelerate_time

            # 生成加速阶段的加速度，这里采用均匀分布，这里的1和100都可以替换成你需要的固定值或者范围
            accelerate_value = np.random.uniform(3, 4, size=accelerate_time)

            # 因为匀速阶段的加速度为0，所以可以直接生成一个全0的数组
            steady_value = np.zeros(steady_time)

            # 随机生成减速阶段的加速度，这里也采用均匀分布
            decelerate_value = np.random.uniform(3, 4, size=decelerate_time)

            # 计算gain值，确保所有阶段的加速度总和为零
            gain = -np.sum(accelerate_value) / np.sum(decelerate_value)

            # 应用gain值
            decelerate_value = -np.sqrt(np.abs(gain)) * decelerate_value
            accelerate_value = accelerate_value /np.sqrt(np.abs(gain))

            if np.all(decelerate_value < 0):  # 确保减速阶段的加速度为负
                break  # 如果满足条件，则跳出循环

        # 将三个阶段的加速度值连起来并添加到matrix中
        matrix.append(np.concatenate((accelerate_value, steady_value, decelerate_value)))

    return np.array(matrix)  # 将list转化为np.array

def a_profile(mode, SampSize, StimSeqSize, ChoiceSeqSize, *args):
    af = np.zeros((SampSize, StimSeqSize+ChoiceSeqSize))
    hfs = int(StimSeqSize/2)
    if mode == 'random':
        a = generate_acc_matrix(SampSize, StimSeqSize)
        # a = np.abs(np.random.rand(SampSize, hfs))
        # a = np.concatenate((a, np.flip(-a, axis=1)), axis=1)
        #SMOOTH
        # for i in range(SampSize):
        #     af[i, :StimSeqSize] = np.convolve(a[i, :], np.ones(5) / 5, 'same')
        #     af[i,StimSeqSize] = -np.sum(af[i,:StimSeqSize])
        af[:, :StimSeqSize] = a
    elif mode == 'trapezium':
        a0 = np.random.rand(SampSize, 1) * 2 + 8
        t1 = np.random.randint(1, int(StimSeqSize / 2), (SampSize,))
        t2 = StimSeqSize - 2 * t1
        for i in range(SampSize):
            af[i, :StimSeqSize] = np.concatenate(
                (np.ones(t1[i]) * a0[i], np.zeros(t2[i]), -np.ones(t1[i]) * a0[i]), axis=0)
    elif mode == 'astrapezium':
        a0 = np.random.rand(SampSize, 1) * 2 + 8
        t10 = np.random.randint(1, int(StimSeqSize / 4), (SampSize,))
        t11 = t10 * 3
        t2 = StimSeqSize - 4 * t10
        for i in range(0, SampSize, 2):
            af[i, :StimSeqSize] = np.concatenate(
                (a0[i] * 2 * np.ones(t10[i]), np.zeros(t2[i]), -a0[i] * 2 / 3 * np.ones(t11[i])), axis=0)
            af[i + 1, :StimSeqSize] = np.concatenate(
                (a0[i] * 2 / 3 * np.ones(t11[i]), np.zeros(t2[i]), -a0[i] * 2 * np.ones(t10[i])), axis=0)
    elif mode == 'input':
        a0 = args[0]
        t1 = int(args[1])
        t2 = int(args[2]) + 1
        ss = int(SampSize / 2)
        StimSeqSize = int(np.ceil(t1 / 2)) + int(np.ceil(t1 * 1.5)) + t2
        t10 = int(np.ceil(t1 / 2))
        t11 = int(np.ceil(t1 * 1.5))
        a10 = 0.64 / t10
        a11 = 0.64 / t11
        if StimSeqSize<22:
            t2 = 20 - t10 - t11
            StimSeqSize = 20
        else:
            t2 = 25 - t10 - t11
            StimSeqSize = 25
        af[0:ss, 0:StimSeqSize] = np.concatenate(
            (a10 * np.ones((ss,t10)), np.zeros((ss, t2)), -a11 * np.ones((ss, t11))), axis=1)
        af[ss:, 0:StimSeqSize] = np.concatenate(
            (a11 * np.ones((ss, t11)), np.zeros((ss, t2)), -a10 * np.ones((ss, t10))), axis=1)
    elif mode == 'heading':
        p = calculate_pos(np.arange(0, StimSeqSize+ChoiceSeqSize), 1, StimSeqSize, 3.5)
        p[StimSeqSize:] = p[StimSeqSize-1]
        v = np.concatenate(([0], np.diff(p)))
        a = np.concatenate(([0], np.diff(v)))
        a = np.repeat(a[None, :], SampSize, axis=0)
        af = a/np.max(a)

    return af

def clean_stim_data(i, stim_end, dt, FixPoint, ChoiceTarget, ax, ay, vx, vy, sx, sy, ax_shift, ay_shift, vx_shift, vy_shift):
    FixPoint[i, 0:stim_end - 1] = 1
    ChoiceTarget[i, stim_end:] = 1

    # ax[i, stim_end-1] += -vx[i, stim_end-1] / dt
    # ay[i, stim_end-1] += -vy[i, stim_end-1] / dt
    #
    # ax_shift[i, stim_end-1] += -vx_shift[i, stim_end-1] / dt
    # ay_shift[i, stim_end-1] += -vy_shift[i, stim_end-1] / dt

    sx[i, stim_end-1:] = sx[i, stim_end-1].reshape(-1, 1) - vx[i, stim_end-1].reshape(-1, 1)  * dt / 2
    sy[i, stim_end-1:] = sy[i, stim_end-1].reshape(-1, 1) - vy[i, stim_end-1].reshape(-1, 1)  * dt / 2

    vx[i, stim_end-1:] = 0
    vy[i, stim_end-1:] = 0

    vx_shift[i, stim_end-1:] = 0
    vy_shift[i, stim_end-1:] = 0

    return FixPoint, ChoiceTarget, ax, ay, vx, vy, sx, sy, ax_shift, ay_shift, vx_shift, vy_shift
def GenerateTaskData(a, heading, noise_sig, StimSeqSize, ChoiceSeqSize, modality, dt, trial_rep = 1, stim_start=0, ref_var=None, ref_dim = 5, heading_shift=0,vis_delay = 0):
    InputChanNum = 1
    a = a / np.max(a)
    v = np.nancumsum(a, axis=1) * dt
    s = np.nancumsum(v, axis=1) * dt
    v = v / np.max(v)/4
    # a = v.copy()    #same profile
    s = s / np.max(s)
    astrength = np.mean(np.abs(a))
    vstrength = np.mean(np.abs(v))
    shifted = np.roll(v, vis_delay, axis=1)
    shifted[:, :vis_delay] = 0
    v = shifted
    ax = a*np.cos(heading)
    ay = a*np.sin(heading)
    vx = v*np.cos(heading)
    vy = v*np.sin(heading)
    sx = s*np.cos(heading)
    sy = s*np.sin(heading)
    heading_shift = np.ones(heading.shape)* heading_shift/2
    ax_shift = a*np.cos(heading-heading_shift)
    ay_shift = a*np.sin(heading-heading_shift)
    vx_shift = v*np.cos(heading+heading_shift)
    vy_shift = v*np.sin(heading+heading_shift)
    time = np.arange(0, a.shape[1]) * dt*1000
    if isinstance(StimSeqSize, np.ndarray):
        max_stim_size = int(max(StimSeqSize) + ChoiceSeqSize + stim_start)
    else:
        max_stim_size = int(StimSeqSize + ChoiceSeqSize + stim_start)
    FixPoint = np.zeros((ax.shape[0], max_stim_size, 1))
    ChoiceTarget = np.zeros((ax.shape[0], max_stim_size, 1))
    if isinstance(StimSeqSize, np.ndarray):
        #variable stimulus size
        for i, stim_size in enumerate(StimSeqSize):
            stim_end = stim_size + stim_start
            # FixPoint, ChoiceTarget, ax, ay, vx, vy, sx, sy, ax_shift, ay_shift, vx_shift, vy_shift = clean_stim_data(i,stim_end,dt,FixPoint,ChoiceTarget,
            #                                                                                                      ax,ay,vx,vy,sx,sy,ax_shift,ay_shift,vx_shift,vy_shift)
    else:
        #fixed stimulus size
        stim_end = StimSeqSize + stim_start
        i = np.arange(FixPoint.shape[0])
        # FixPoint, ChoiceTarget, ax, ay, vx, vy, sx, sy, ax_shift, ay_shift, vx_shift, vy_shift = clean_stim_data(i,stim_end,dt,FixPoint,ChoiceTarget,
        #                                                                                                          ax,ay,vx,vy,sx,sy,ax_shift,ay_shift,vx_shift,vy_shift)
        StimSeqSize = np.repeat(StimSeqSize, ax.shape[0])
    # plt.plot(ax[0])
    # plt.plot(vx[0])
    # plt.show()
    y_data = np.stack((ax, ay, vx, vy, sx, sy,a,v,s), axis=2)
    if modality == 0:  # ves
        vx.fill(0)
        vy.fill(0)
        vx_shift.fill(0)
        vy_shift.fill(0)
        vstrength = 0
    elif modality == 1:  # vis
        ax.fill(0)
        ay.fill(0)
        ax_shift.fill(0)
        ay_shift.fill(0)
        astrength = 0
    noise_siga = noise_sig[:, 0]*astrength
    noise_sigv = noise_sig[:, 1]*vstrength
    x_datas = []
    y_datas = []
    for i in range(noise_sig.shape[0]):
        x_data_i = np.concatenate((AN(ax_shift, InputChanNum, noise_siga[i], trial_rep=trial_rep),
                                   AN(ay_shift, InputChanNum, noise_siga[i], trial_rep=trial_rep),
                                   AN(vx_shift, InputChanNum, noise_sigv[i], trial_rep=trial_rep),
                                   AN(vy_shift, InputChanNum, noise_sigv[i], trial_rep=trial_rep)),
                                  axis=2)  # AN(-ax_shift, inputdim, noise_siga[i],trial_rep=trial_rep),
        repeat_times = (trial_rep,) + (1,) * (y_data.ndim - 1)
        # trial repetition
        y_data_i = np.tile(y_data, repeat_times)
        x_datas.append(x_data_i)
        y_datas.append(y_data_i)
    x_data = np.concatenate(x_datas, axis=0)
    y_data = np.concatenate(y_datas, axis=0)
    StimSeqSize = np.repeat(StimSeqSize, trial_rep * noise_sig.shape[0])
    if ref_var>0: #distance task
        var_ref = y_data[:,:,ref_dim]
        labels = var_ref > ref_var
    else: #heading task
        var_ref = y_data[:, -1, ref_dim]
        if ref_var == None:
            ref_var = var_ref.mean()
            labels = np.expand_dims(bool_to_01(var_ref > var_ref.mean()), axis=1)
        else:
            var_ref[var_ref == 0] = np.random.choice([-1, 1], size=sum(var_ref == 0))
            labels = np.expand_dims(bool_to_01(var_ref > ref_var), axis=1)

    return x_data, y_data, labels, ref_var, StimSeqSize


