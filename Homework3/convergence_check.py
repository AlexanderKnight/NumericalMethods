import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def power_function(y, a, b, c):
    return a*(y**b)+c


method = "RungeKutta3-SineWave"
cf = "0.500000"
exts = np.array((100,200,400,800,1600))
x = np.linspace(100,1600,1000)
maxdiff = np.zeros(len(exts))
extents = np.array(("100","200","400","800","1600"))
for i in range(len(exts)):
    dat = np.genfromtxt("data_files/"+method+"-cf-"+str(cf)+"-"+str(exts[i])+".dat")
    sin = np.genfromtxt("data_files/Sinewave-cf-"+str(cf)+"-"+str(exts[i])+".dat")
    diff = np.abs(dat[:,2:-2]-sin);
    maxdiff[i]=np.amax(diff)

opt,cov = curve_fit(power_function,exts,maxdiff,[1,-3,0])
plt.plot(exts,maxdiff)
plt.plot(x,power_function(x,opt[0],opt[1],opt[2]),label="ax^b+c:\na={0}\nb={1}\nc={2}".format(opt[0],opt[1],opt[2]))
plt.legend(bbox_to_anchor=(1,1))
#plt.show()
plt.savefig(method+"_ConvergenceTest.png",dpi=200)
