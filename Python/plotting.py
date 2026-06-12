# Functions to be used with Matplotlib.pyplot

def plot(
    data = [], # Nx2xM array [[[x1],[y1]],...], where xi and yi are lists of size M
    title = 'Title',
    ylabel = 'ylabel',
    xlabel = 'xlabel',
    scatter = False, # that the plot will be a scatterplot
    line = True, # that the plot will be linear interpolant (default = True)
):

    # plot each data set
    