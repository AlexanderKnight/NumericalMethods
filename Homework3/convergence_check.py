import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def power_function(y, a, b, c):
    return a*(y**b)+c

def linear_func(y,a,b):
    return a*y+b


method = "RungeKutta3-SineWave"
cf = "0.100000"
exts = np.array((100,200,400,800))
dx = np.array((2*np.pi)/exts)
x = np.linspace(100,800,1000)
maxdiff = np.zeros(len(exts))
extents = np.array(("100","200","400","800"))
for i in range(len(exts)):
    dat = np.genfromtxt("data_files/"+method+"-cf-"+str(cf)+"-"+str(exts[i])+".dat")
    sin = np.genfromtxt("data_files/Sinewave-cf-"+str(cf)+"-"+str(exts[i])+".dat")
    diff = np.abs(dat[-1,2:-2]-sin[-1,:]);
    maxdiff[i]=np.sum(diff)*dx[i]
    

opt,cov = curve_fit(power_function,exts,maxdiff,[1.,-3.,0.])
print("The convergence order is {0:.2f}. See accompanied plot.".format(np.abs(opt[1])))
plt.title("Convergence Check")
plt.plot(exts,maxdiff, label="Error vs Extents")
plt.plot(x,power_function(x,opt[0],opt[1],opt[2]),label="a*(x**b)+c:\na={0}\nb={1}\nc={2}".format(opt[0],opt[1],opt[2]))
plt.legend(bbox_to_anchor=(1,1))
#plt.show()
plt.savefig(method+"_ConvergenceTest.png",dpi=200)
