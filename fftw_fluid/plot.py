import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("dye_step1.txt")
plt.imshow(data, cmap='viridis')
plt.colorbar()
plt.show()
