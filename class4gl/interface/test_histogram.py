

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

data = [220,14.2,150,400,420]
error = [10, 1, 20, 60, 10]
x = [i + .5 for i in range(5)]

fig, ax = plt.subplots()
bar = ax.bar(x, data,0.1 +np.arange(len(data))*0.9/len(data), align="center", yerr=error)
plot = ax.plot(x, data)
ax.set_xticks(x)
ax.set_xticklabels(('wt', 'N23PP', 'N23PP/PEKN', 'PEKN', 'N23PP/PEKN/L28F'))
ax.set_title(r"Everything in the document can use m$\alpha$th language",
             y=1.05)
ax.set_ylabel(r"Rate (s$^{-1}$)", labelpad=10)
ax.set_xlabel("Mutant",labelpad=10)
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
plt.savefig('test.png')
plt.show()


