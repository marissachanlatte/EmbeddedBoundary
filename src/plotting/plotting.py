import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.style as style
from matplotlib import collections  as mc
import seaborn as sns
import pandas as pd
style.use('ggplot')

filename = sys.argv[1]
array = np.loadtxt(filename)
data = {'x': array[:, 0], 'y': array[:, 1], 'area': array[:, 2]}
df = pd.DataFrame(data=data)
pivot = df.pivot(index='y', columns='x', values='area')
ax = sns.heatmap(pivot)
ax.invert_yaxis()
plt.savefig("boundary_plot")
