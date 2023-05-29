import numpy as np
import matplotlib.pylab as pl

pl.close('all')

def line_hover(event):
    ax = pl.gca()
    for line in ax.get_lines():
        if line.contains(event)[0]:
            print(line.get_label())

labels = ['line 1','line 2','line 3']

fig = pl.figure()
for i in range(len(labels)):
    pl.plot(np.arange(10), np.random.random(10), label=labels[i])
pl.legend(frameon=False)

fig.canvas.mpl_connect('motion_notify_event', line_hover)
pl.show()