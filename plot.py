import matplotlib.pyplot as plt
import numpy as np

# Read file02
data = np.loadtxt("run/file02", skiprows=1)

x       = data[:, 0]
u_init  = data[:, 1]
u_final = data[:, 2]
u_anal  = data[:, 3]

# Plot
plt.figure(figsize=(9, 5))
plt.plot(x, u_init,  'b--', label='Initial distribution')
plt.plot(x, u_final, 'r-',  label='Final numerical')
plt.plot(x, u_anal,  'go',  markersize=6, label='Analytical steady state')

plt.xlabel('x (m)')
plt.ylabel('Temperature u')
plt.title('1D Convection-Diffusion — Backward Euler')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("run/result.png", dpi=150)
plt.show()