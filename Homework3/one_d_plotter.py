import numpy as np
import matplotlib.pyplot as plt



sine_data = np.genfromtxt("data_files/Sinewave-cf-0.500000.dat")
data = np.genfromtxt("data_files/RungeKutta3-SineWave-cf-0.500000.dat")
for i in range(len(data))[0::50]:
    print("            "+str(i))
    plt.figure(figsize=(3,3))
    plt.title("RK3,3Upstream,0.5")
    plt.plot(data[i,1:-1])
    plt.plot(sine_data[i,:])
    plt.savefig("figures/RungeKutta3_upstream_0.5_"+str(i)+".png", dpi=200)
    plt.close()

