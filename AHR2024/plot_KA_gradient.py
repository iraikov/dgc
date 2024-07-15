import numpy as np
import matplotlib.pyplot as plt
x = np.arange(100, 250, 0.1)
prox_gmax_KA = 0.0017
gmax_KA_scale = -0.03
gmax_KA_shift = 70
gmax_KA = prox_gmax_KA / (np.exp(gmax_KA_scale*(x+gmax_KA_shift)) + 1.0)
plt.plot(x, gmax_KA)

plt.show()
gmax_KA = prox_gmax_KA / (1.0 + np.exp(gmax_KA_scale*(x - gmax_KA_shift)))
plt.plot(x, gmax_KA)

plt.show()
x = np.arange(60, 250, 0.1)
gmax_KA = prox_gmax_KA / (1.0 + np.exp(gmax_KA_scale*(x - gmax_KA_shift)))
plt.plot(x, gmax_KA)

plt.show()
x = np.arange(60, 350, 0.1)
gmax_KA = prox_gmax_KA / (1.0 + np.exp(gmax_KA_scale*(x - gmax_KA_shift)))
plt.plot(x, gmax_KA)

plt.show()
gmax_KA_scale = -0.02
gmax_KA = prox_gmax_KA / (1.0 + np.exp(gmax_KA_scale*(x - gmax_KA_shift)))
plt.plot(x, gmax_KA)
