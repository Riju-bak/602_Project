import matplotlib.pyplot as plt


def plot_2D(y, x, **kwargs):
    figure_name = kwargs.get("figure_name", "")
    title = kwargs.get("title", "")
    xlabel = kwargs.get("xlabel", "")
    ylabel = kwargs.get("ylabel", "")
    plt.figure(figure_name)
    plt.title(title)
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)
    plt.plot(x, y)
    plt.draw()
