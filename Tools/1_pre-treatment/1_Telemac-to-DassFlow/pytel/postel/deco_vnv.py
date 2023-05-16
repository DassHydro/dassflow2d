from cycler import cycler

decoVNV = {
    'font.size': 14.0,
    'font.serif': 'Helvetica',
    'text.color': 'black',
    'text.usetex': False,
    'axes.grid': True,
    'axes.facecolor': 'white',
    'axes.edgecolor': 'black',
    'axes.linewidth': 0.3,
    'axes.xmargin': .0,
    'axes.ymargin': .0,
    'axes.prop_cycle': cycler(
        'color',
        ['#002d74', '#e85113', '#1fa12e', '#c9d200', '#f49e00', '#006ab3',
         '#002d74', '#e85113', '#1fa12e', '#c9d200', '#f49e00', '#006ab3',
         '#002d74', '#e85113', '#1fa12e', '#c9d200', '#f49e00', '#006ab3'])
        + cycler('marker',
                 ['', '', '', '', '', '',
                  '+', '1', '2', '3', '4', 'x',
                  '+', '1', '2', '3', '4', 'x'])
        + cycler('linestyle',
                 ['-', '-', '-', '-', '-', '-',
                  '-', '-', '-', '-', '-', '-',
                  ':', ':', ':', ':', ':', ':']),
    'axes.labelsize': 14.,
    'axes.labelpad': 4.0,
    'axes.labelcolor': 'black',
    'axes.titlesize': 14,
    'lines.markeredgewidth': 1.0,
    'lines.markersize': 4,
    'grid.color': 'b0b0b0',
    'grid.linestyle': ':',
    'grid.linewidth': 0.7,
    'grid.alpha': 1.0,
    'xtick.minor.visible': False,
    'xtick.major.size': 3.5,
    'xtick.minor.size': 1.5,
    'ytick.minor.visible': False,
    'ytick.major.size': 3.5,
    'ytick.minor.size': 1.5,
    'legend.loc': 'best',
    'legend.fontsize': 12.,
    'legend.frameon': True,
    'legend.fancybox': True,
    'figure.dpi': 100,
    'figure.autolayout': True,
    'figure.max_open_warning': 200,
    'image.aspect': 'equal',
    }

decoVNV_bar = {
    'axes.xmargin': .05,
    'axes.ymargin': .05,
    }

decoVNV_1d = {
    'xtick.minor.visible': True,
    'ytick.minor.visible': True,
    }

decoVNV_markers = {
    'axes.prop_cycle': cycler(
        'color',
        ['#002d74', '#e85113', '#1fa12e', '#c9d200', '#f49e00', '#006ab3',
         '#002d74', '#e85113', '#1fa12e', '#c9d200', '#f49e00', '#006ab3',
         '#002d74', '#e85113', '#1fa12e', '#c9d200', '#f49e00', '#006ab3'])
        + cycler('marker',
                 ['+', '1', '2', '3', '4', 'x',
                  'o', 'v', '^', '<', '>', 's',
                  '+', '1', '2', '3', '4', 'x'])
        + cycler('linestyle',
                 ['-', '-', '-', '-', '-', '-',
                  '--', '--', '--', '--', '--', '--',
                  ':', ':', ':', ':', ':', ':']),
    }
