import numpy as np



def tuning_function_bined(variable, fr, variable_range=None):
    # Compute variable range
    if variable_range is None:
        variable_range = np.linspace(np.min(variable), np.max(variable), 11)
    uv = np.unique(variable)
    if len(uv)<20:
        variable_range = uv
    # put variable into bins
    bins = np.digitize(variable, variable_range)
    bins[bins==len(variable_range)] = len(variable_range)-1
    bins[bins==len(variable_range)] = len(variable_range)-1
    # Compute mean firing rate and std and mean variable for each bin
    tuning_function = np.zeros((len(variable_range)-1))
    std = np.zeros((len(variable_range)-1))
    m_variable = np.zeros((len(variable_range)-1))
    for i in range(len(variable_range)-1):
        tuning_function[i] = np.mean(fr[bins==i+1], axis=0)
        std[i] = np.std(fr[bins==i+1], axis=0)/np.sqrt(np.sum(bins==i+1))
        m_variable[i] = np.mean(variable[bins==i+1])
    return tuning_function, m_variable, std

