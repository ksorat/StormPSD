#!/usr/bin/env python
from pylab import *
import datetime

LW = 0.75

if sys.argv[1]:
    filename=sys.argv[1]
else:
    'Enter file name!'

if sys.argv[2]:
    start_time_str = sys.argv[2]
else:
    'Enter start time (YYYY-MM-DD-HH-MM-SS)!'

if sys.argv[3]:
    if sys.argv[3] == '1': 
        PUBLICATION = True
    else:
        PUBLICATION = False
else:
    PUBLICATION = True

if not PUBLICATION:
    rcdefaults()
else:
##### DEFINE FIGURE SIZE ETC FOR PUBLICATION QUALITY ######
    golden_mean = (sqrt(5)-1.0)/2.0
    fig_width = 20./6. ## 20 pica (pc) = 20/6 inches
    fig_height = 5.#golden_mean*fig_width
    fig_size = [fig_width,fig_height]
    
    params = {'axes.labelsize': 7,
              'text.fontsize': 7,
              'legend.fontsize': 7,
              'xtick.labelsize': 6,
              'ytick.labelsize': 6,
              'figure.figsize': fig_size}
    rcParams.update(params)

rc('mathtext',fontset='stixsans',default='regular')
############################################################


(year,month,day,hours,minutes,seconds) = [int(s) for s in start_time_str.split('-')]
start_time = datetime.datetime(year,month,day,hours,minutes,seconds)
print('Start time: ',start_time)
print('PUBLICATION: ',PUBLICATION)

#start_time=datetime.datetime(2012,7,14,15)
#start_time=datetime.datetime(2012,1,24,13)

nstart = where(['DATA:' in l for l in open(filename)])[0]
swdata = loadtxt(filename,skiprows=nstart+1)

time=[]
[time.append(start_time+datetime.timedelta(minutes=dt)) for dt in swdata[:,0]]

fig,ax=subplots(3,1,sharex=True)

ax[0].plot(time,swdata[:,1],'k',linewidth=LW)
ax[0].set_ylabel('n, cm$^\mathrm{-3}$')

ax[1].plot(time,swdata[:,2],label='V$_x$',linewidth=2*LW)
ax[1].plot(time,swdata[:,3],label='V$_y$',linewidth=LW)
ax[1].plot(time,swdata[:,4],label='V$_z$',linewidth=LW)
ax[1].legend(loc='lower left')
ax[1].set_ylabel("V, km/s")

#ax[2].plot(time,swdata[:,5],'k')
#ax[2].set_ylabel('C$_sound$, km/s')

ax[2].plot(time,swdata[:,6],label='B$_x$',linewidth=LW)
ax[2].plot(time,swdata[:,7],label='B$_y$',linewidth=LW)
ax[2].plot(time,swdata[:,8],label='B$_z$',linewidth=2*LW)
#ax[2].plot(time,sqrt(swdata[:,6]**2+swdata[:,7]**2+swdata[:,8]**2),'k',label='B$_t$')
ax[2].legend(loc='upper left')
ax[2].axhline(y=0,linestyle='--',color='k')
ax[2].set_ylabel("Magnetic field, nT")
ax[2].set_xlabel('Time, h')

ax[2].xaxis.set_major_formatter(DateFormatter('%H:%M'))
fig.tight_layout()

for a in ax: a.grid(True)


if not PUBLICATION:
    show()
else:
    savefig('sw_imf.pdf')


# ax1=subplot(4,1,1)
# plot(time,swdata[:,1],'k')
# setp( ax1.get_xticklabels(), visible=False)
# ylabel('n, cm$^\mathrm{-3}$')

# ax2=subplot(4,1,2)
# plot(time,swdata[:,2],label='V$_x$')
# plot(time,swdata[:,3],label='V$_y$')
# plot(time,swdata[:,4],label='V$_z$')
# legend(loc='upper left')
# setp( ax2.get_xticklabels(), visible=False)
# ax2.yaxis.tick_right()
# ax2.yaxis.set_label_position("right")
# ylabel("V, km/s")

# ax3=subplot(4,1,3)
# plot(time,swdata[:,5],'k')
# setp( ax3.get_xticklabels(), visible=False)
# ylabel('Cs, km/s')

# ax4=subplot(4,1,4)
# plot(time,swdata[:,6],label='B$_x$')
# plot(time,swdata[:,7],label='B$_y$')
# plot(time,swdata[:,8],label='B$_z$')
# plot(time,sqrt(swdata[:,6]**2+swdata[:,7]**2+swdata[:,8]**2),'k',label='B$_t$')
# legend(loc='upper left')
# axhline(y=0,linestyle='--',color='k')
# #setp( ax4.get_xticklabels(), visible=False)
# ax4.yaxis.tick_right()
# ax4.yaxis.set_label_position("right")
# ylabel("B, nT")

# ax5=subplot(5,1,5)
# plot(time,sqrt(swdata[:,6]**2+swdata[:,7]**2+swdata[:,8]**2),'k')
# ylabel("B$_{total}$, nT")
# xlabel('Time, h')

# subplots_adjust(hspace=0)
# show()
