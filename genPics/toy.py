import kCyl as kc
import pyStorm as pS
import os
import numpy as np
import datetime
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec
import matplotlib.dates as mdates
import lfmViz as lfmv

lfmv.ppInit()

x = np.linspace(0,100)
y = np.sin(2*np.pi*x)

plt.plot(x,y,'r')
plt.ylabel(r'\textcolor{blue}{Value}')