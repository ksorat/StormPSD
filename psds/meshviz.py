#Toy to play with various meshes
from pylab import *


Np = 48
Nr = 30

Rin = 2.05
Rout = 18

pI = np.linspace(0,2*np.pi,Np+1)

dP = pI[1]-pI[0]

rILog = np.logspace(np.log10(Rin),np.log10(Rout),Nr+1)
rILin = np.linspace(Rin,Rout,Nr+1)

xxLin = np.zeros((Nr+1,Np+1))
yyLin = np.zeros((Nr+1,Np+1))

xxLog = np.zeros((Nr+1,Np+1))
yyLog = np.zeros((Nr+1,Np+1))

for i in range(Nr+1):
	for j in range(Np+1):
		xxLin[i,j] = rILin[i]*np.cos(pI[j])
		yyLin[i,j] = rILin[i]*np.sin(pI[j])
		xxLog[i,j] = rILog[i]*np.cos(pI[j])
		yyLog[i,j] = rILog[i]*np.sin(pI[j])

#plt.plot(xxLog,yyLog,'r-',xxLog.T,yyLog.T,'r-')

#plt.plot(xxLin,yyLin,'b-',xxLin.T,yyLin.T,'b-')
#plt.show()

dLoLin = np.zeros(Nr)
dLoLog = np.zeros(Nr)
dLog = np.zeros(Nr)
for i in range(Nr):
	dLog[i] = rILog[i+1]-rILog[i]
	L = 0.5*(rILog[i+1]+rILog[i])
	dLoLog[i] = dLog[i]/L

	dL = rILin[i+1]-rILin[i]
	L = 0.5*(rILin[i+1]+rILin[i])
	dLoLin[i] = dL/L

rCLog = 0.5*(rILog[1:]+rILog[0:-1])
rCLin = 0.5*(rILin[1:]+rILin[0:-1])

#plt.plot(rCLin,dLoLin,'r',rCLog,dLoLog,'b')
#plt.show()

aRat = dP/dLoLog[0]
print("Aspect ratio (dp/(dL/L)) = %f"%(aRat))