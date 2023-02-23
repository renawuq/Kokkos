import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
#rc('text', usetex=True)

size=np.array([50,550,1050,1550,2050])

v0=np.array([0.0112829,1.17746,6.00018,17.7317,26.3108])

v1=np.array([2.19471e-314,6.9516e-310,2.15692e-314,2.19475e-314,6.9516e-310])

v2=np.array([6.9516e-310,2.15693e-314,6.9516e-310,6.9516e-310,nan])

v3=np.array([nan,4.94016e-320,6.9516e-310,0,0])

v4=np.array([0,6.92455e-310,1.18576e-322,1.4822e-322,0])

v5=np.array([0,0,2.15826e-314,6.9516e-310,2.15689e-314])

plt.plot(size,v0, label='serial')
plt.plot(size,v1, label='kokkod with struct')
plt.plot(size,v2, label='kokkod half flat half struct')
plt.plot(size,v3, label='kokkod using struct with 2d range')
plt.plot(size,v4, label='kokkod flat without 2d range')
plt.plot(size,v5, label='kokkod flat with 2d range')
plt.grid(True)
plt.title('Size vs Time')
plt.xlabel('Domain size')
plt.ylabel(r'Time')
plt.legend()
plt.show()
