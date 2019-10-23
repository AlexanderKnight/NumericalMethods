import numpy as np
import matplotlib.pyplot as plt



extents = ["200"]
for j in range(len(extents)):
    #sine_data = np.genfromtxt("data_files/Sinewave-cf-0.100000-{0}.dat".format(extents[j]))
    data = np.genfromtxt("data_files/LLF-SineWave-cf-0.100000-{0}.dat".format(extents[j]))
    for i in range(len(data)-1)[0::50]:
        print("            "+str(i))
        plt.figure(figsize=(3,3))
        plt.title("LLF,3Upstream,0.5\n{0},{1}".format(extents[j],i))
        plt.plot(data[i,2:-2])
        #plt.plot(sine_data[i+1,:])
        plt.savefig("figures/RungeKutta3_upstream_0.5_{0}_{1}.png".format(extents[j],i), dpi=200)
        plt.close()

