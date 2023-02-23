import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
#rc('text', usetex=True)

size=np.array([0,1,2])

v0=np.array([2.35704,5.55185,7.88122])

v1=np.array([4.06752,3.61922,10.7565])

v2=np.array([3.89305,6.05328,17.9082])

v3=np.array([7.54227,7.06954,7.02248])

plt.plot(size,v0, label='serial')
plt.plot(size,v1, label='kokkod 0')
plt.plot(size,v2, label='kokkod 1')
plt.plot(size,v3, label='kokkod 2')
plt.grid(True)
plt.title('draft')
plt.xlabel('Time Stamp')
plt.ylabel(r'Time')
plt.legend()
plt.show()
