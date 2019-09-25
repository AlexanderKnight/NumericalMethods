import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def power_function(x, a, b):
    return a*(x**b)

Dt = ["0.100000", "0.200000","0.300000","0.400000","0.500000"]
Diff = ["centered", "downstream", "upstream"]

x = np.linspace(0,0.0175,500)
prev_best = np.zeros(len(Dt)-1)
next_best = np.zeros(len(Dt)-1)

for diff in Diff:
    plt.figure(figsize=(5,5))
    plt.title("Forward Euler, "+diff)
    plt.xlim(0,0.0175)
    for i,dt in enumerate(Dt):
        sine_data = np.genfromtxt("data_files/Sinewave-dt-"+dt+".dat")
        data = np.genfromtxt("data_files/ForwardEuler-"+diff+"-SineWave-dt-"+dt+".dat")
        data = data[:,1:-1]
        diff_data = np.abs(data-sine_data)
        if(i==0):
            best = np.max(diff_data)
            prev_best[i]=best
        else:
            next_best[i-1]=np.max(diff_data)
            #dt = float(dt)*(2*np.pi)/(300)
            plt.scatter(best,np.max(diff_data))
            if(i!=len(Dt)-1):
                prev_best[i]=np.max(diff_data)
            best=np.max(diff_data)
    opt,cov = curve_fit(power_function,prev_best,next_best)
    print("Convergence order of {0} is {1}".format(diff,opt[1]))
    plt.savefig("figures/ForwardEuler/convergence_plots/"+diff+".png",dpi=200)
    plt.close()
