import os
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import pearsonr, norm
from tool.tuning_function import tuning_function_bined


def plot_psychophysical_function(variable, choice, labels=None, label_names=None, bandwidth=None, label_color=None, marker_style='o', line_style='-', ax=None,refvar=0):
    if bandwidth is None:
        bandwidth = 0.001 * (np.max(variable) - np.min(variable))
    if labels is None:
        labels = np.zeros_like(variable)
    Vr = np.max(np.abs(variable-refvar))
    variable_range = np.linspace(refvar-Vr, refvar+Vr, 11)
    variable_range_func = np.linspace(refvar-Vr, refvar+Vr, 20)
    unique_labels = np.unique(labels)
    unique_variable = np.unique(variable)
    if len(unique_variable)<20:
        variable_range = unique_variable
        bandwidth = 1e-10
    if refvar>0:
        distance_flag = True
    else:
        distance_flag = False

    # Default colors
    default_colors = ['blue', 'red', 'green']

    # If label_color is None, assign default colors
    if label_color is None:
        label_color = default_colors[:len(unique_labels)]

    if ax is None:
        fig, ax = plt.subplots()
    props = []
    accuracy = np.zeros(len(unique_labels))

    #fit direction
    corr, p_value = pearsonr(variable, choice)
    if corr<0:
        choice = 1-choice
        ylab = 'P(Left)'
        if distance_flag:
            ylab = 'P(Short)'
    else:
        ylab = 'P(Right)'
        if distance_flag:
            ylab = 'P(Long)'
    for idx, label in enumerate(unique_labels):
        mask = labels == label
        variable_label = variable[mask]
        choice_label = choice[mask]
        correct_choice = variable_label>refvar
        accuracy[idx] = np.mean(correct_choice==choice_label)
        cl = choice_label[variable_label>1]
        # Compute weighted averages
        weighted_average, variable_range, error = tuning_function_bined(variable_label, choice_label, variable_range)
        # Fit a cumulative Gaussian model
        initial_guess = [np.mean(variable_label), np.std(variable_label)]
        popt, _ = curve_fit(cum_gauss, variable_label, choice_label, p0=initial_guess, maxfev=10000, bounds=([np.min(variable_label), 0], [np.max(variable_label), np.std(variable_label)*20]))
        # Generate fitted curve
        fit = cum_gauss(variable_range_func, *popt)
        # Plot the weighted averages and the fitted line
        fit_results = f"mu: {popt[0]:.4f}, std: {popt[1]:.4f}"  # Assuming result.x contains mu and std
        ax.errorbar(variable_range, weighted_average, yerr=error, fmt=marker_style,
                    color=label_color[idx])  #label=f'Data {label_names[label] if label_names else label}',
        ax.plot(variable_range_func, fit, linestyle=line_style, label=f'Fit {label_names[label] if label_names else label} ({fit_results})',
                color=label_color[idx])
        props.append(popt)
    props = np.array(props)
    ax.set_xlabel('Variable Value')
    ax.set_ylabel(ylab)
    ax.legend(loc='upper left', fontsize='small')
    return ax,props,accuracy  # Optionally return ax, so it can be used further


def save_plot(plot_name, directory='test', plot_flag= False):
    # Modify the directory to be under 'img/'
    directory = os.path.join('img', directory)
    # Check if the directory exists, if not, create it
    if not os.path.exists(directory):
        os.makedirs(directory)
    # Save the plot
    filename = f'{plot_name}.png'
    plt.savefig(os.path.join(directory, filename))
    if plot_flag:
        plt.show()
    plt.close()  # Close the plot to free up memory


def cum_gauss(x, mu, sigma):
    return norm.cdf(x, mu, sigma)
