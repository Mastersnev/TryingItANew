import matplotlib.pyplot as plt
import numpy as np


def Plotlog(x_data, y_data, show=True, **kwargs):
    """
    :param x_data:
    :param y_data:
    :param show: If true, graph is plotted
    :param kwargs: additional arguments [title, x_limit, y_limit, x_label, y_label, labels]
    :return:
    """
    if 'title' in kwargs:
        title = kwargs['title']
        plt.title(title)

    if 'x_limit' in kwargs:
        x_limit = kwargs['x_limit']
        plt.xlim(x_limit)

    if 'y_limit' in kwargs:
        y_limit = kwargs['y_limit']
        plt.ylim(y_limit)

    if 'x_label' in kwargs:
        x_label = kwargs['x_label']
        plt.xlabel(x_label)

    if 'y_label' in kwargs:
        y_label = kwargs['y_label']
        plt.ylabel(y_label)

    if 'labels' in kwargs:
        labels = np.array(kwargs['labels'])
        legend = True
    else:
        legend = False

    if 'linetype' in kwargs:
        lt = np.array(kwargs['linetype'])

    plt.yscale('log')
    plt.xscale('log')

    if type(x_data) == tuple:
        num_plots = len(x_data)
    else:
        num_plots = 1
    x_data = np.array(x_data)
    y_data = np.array(y_data)


    if num_plots > 1:
        if legend:
            for i in range(num_plots):
                plt.plot(x_data[i], y_data[i], label=labels[i])

            plt.legend()
        else:
            for i in range(num_plots):
                plt.plot(x_data[i], y_data[i])

    else:
        if legend:
            plt.plot(x_data, y_data, label=labels)

            plt.legend()
        else:
            plt.plot(x_data, y_data)

    if show:
        plt.show()
