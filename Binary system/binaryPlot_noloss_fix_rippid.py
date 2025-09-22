import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import LineString
import matplotlib.animation as animation
import math

# Constants
G = (6.6743 * (10**-11))*(1.9891*(10**30))*(1/((1.495978707*10**11)**3))*((365.25*24*3600)**2) # (AU^3)*(year^-2)*(solarmass^-1)
M = 36 # solar mass
m = 29 # solar mass
l = math.sqrt(0.5)  # solarmass*AU^2/year
e = 0 # Eccentricity
eps = 1e-10

M_total = M+m
M_025 = M_total**0.25
tE = 60/365.25
tB = (tE/M_025)*6
a = (M_total*(tB**2))**(1/3)

# Set to find distance between star (r)
theta_values = 0
r_values = a / (1 + e * np.cos(theta_values))
x1 = []
y1 = []
x2 = []
y2 = []

# Calculate positions
for t in np.arange(0, tB+tB, (1/365.25)):
    omaga = (2*math.pi)/tB #day^-1
    theta_values = omaga*t
    r_values = a / (1 + e * np.cos(theta_values))
    r1 = M/(M+m)*r_values
    r2 = (-1*(m/(M+m)))*r_values
    x1.append(r1 * np.cos(theta_values))
    y1.append(r1 * np.sin(theta_values))
    x2.append(r2 * np.cos(theta_values))
    y2.append(r2 * np.sin(theta_values))

# Animation
fig, ax = plt.subplots(figsize=(6, 6))
scat = ax.scatter([], [], c="indigo", s=5, label='Object A')
SCAT = ax.scatter([], [], c="r", s=5, label='Object B')
#ax.set(xlim=[-3, 3], ylim=[-3, 3])
ax.legend()

def update(frame):
    # Update data for scatter plots
    data_min = np.array([x1[frame], y1[frame]]).reshape(1, -1)
    data_max = np.array([x2[frame], y2[frame]]).reshape(1, -1)
    
    # Update scatter plot positions without clearing
    scat.set_offsets(data_min)
    SCAT.set_offsets(data_max)
    return scat, SCAT

# Create animation
ani = animation.FuncAnimation(fig=fig, func=update, frames=len(x1), interval=500, blit=True)

writer = animation.PillowWriter(fps=15,
                                metadata=dict(artist='Me'),
                                bitrate=1800)

plt.scatter(x1, y1, s=2, c='m', marker = 'o', label='Object A')
plt.scatter(x2, y2, s=2, c='pink', marker = 'o', label='Object B')
plt.scatter(0,0, s=10, c='r', marker = '+')
ani.save('C:\\Users\\LeZen_e595\\Desktop\\blackhole_program\\BinaryOrbit\\Figure\\binary_orbit.gif', writer=writer)
plt.show()