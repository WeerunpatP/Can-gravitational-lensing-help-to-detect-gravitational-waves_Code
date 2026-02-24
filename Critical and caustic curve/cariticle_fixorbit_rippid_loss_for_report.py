import numpy as np
import matplotlib.pyplot as plt
import math
from sympy import symbols, Eq, solve, sqrt
import asyncio
import pandas as pd

# Basic unit AU, year, solar mass, times of tE
# Constants
G = (6.6743 * (10**-11)) * (1.9891 * (10**30)) * (1 / ((1.495978707 * 10**11)**3)) * ((24 * 3600 * 365.25)**2)  # (AU^3)*(year^-2)*(solar mass^-1)
e = 0  # Eccentricity
c = (299792458) * (1 / (1.495978707 * 10**11)) * (24 * 3600 * 365.25)  # AU/year
ly = 9.4605284 * (10**15)  # m
AU = 149597871*(10**3) #m
H0 = 69.8 * (10**3) / ((10**6) * 3.26 * ly)  # m/s/m
Ds = (2.537*(10**6)*ly)/AU #AU, andromeda galaxy
Dl = (5*ly)/AU
M_solar = 1.9891*(10**30) #kg
M = 36 # solar mass
m = 29 # solar mass
q = m/M
print("q: ",q, " solar mass")
M_A = M+m
E1 = M/M_A
E2 = m/M_A

rs = 2 * G * (M + m) / (c**2)  # Schwarzschild radius unit AU
Re = math.sqrt((4*G*M_A*Dl*(Ds-Dl))/((c**2)*Ds))
tE = 60/365.25 #year
Vt = Re/tE

print("tE: ", tE, " day")
print("Vt: ", Vt, " AU/day")

Dls = Ds-Dl
print("Dl: ",Dl, " AU")
print("Ds: ", Ds, " AU")

M_025 = M_A**0.25
tB = (tE/M_025)*6
a = (M_A*(tB**2))**(1/3)
eps = 1e-10

Ai_data = []
A_data = []

print("a: ", a, " AU")

# Parameters for lens
scale = 40
chrip_time = a**4/((32/5)*((m*M)/(m+M))*(rs**3)*c)
start_time = chrip_time/1000 - (scale/2)
print(f"start time : {start_time}")
end_time = start_time+scale
time = [-4*tE, -2*tE, eps,2*tE,4*tE,6*tE,10*tE]
print(f"time: {time}")

r_list = []
omega_list = []
omega = math.sqrt((G*(M+m))/(a**3))
    
for t in time:
    a_t = (a**4 - ((32/5)*((m*M)/(m+M))*(rs**3)*c*(t+start_time)))**(1/4)
    if not(isinstance(a_t, complex)) and a_t > 0:
        #++++++++++++++++++++++++++++++++++++++++++++++++++++++
        r_valuse = a / (1 + e * np.cos(omega*(t+start_time))) # This updates r_valuse based on that omega and time
        r_list.append(r_valuse)
        omega = math.sqrt((G*(M+m))/(r_valuse**3))
        omega_list.append(omega)
    else:
         break

r_list = np.array(r_list)
omega_list = np.array(omega_list)
time = np.array(time)

print(len(r_list))
print(len(omega_list))

r1 = M / (M + m) * r_list
r2 = -1 * (m / (M + m)) * r_list
theta_values = omega_list*(time+start_time)
i_unitv = np.cos(theta_values)
j_unitv = np.sin(theta_values)

# Using the calculated values for X1, Y1, X2, Y2
X1=(r1*i_unitv) #AU
Y1=(r1*j_unitv) #AU
X2=(r2*i_unitv) #AU
Y2=(r2*j_unitv) #AU

star_1 = X1 + Y1*1j
star_2 = X2 + Y2*1j

source_i_position = -1*(Vt * ((scale/2)*tE))
x_max = 5
x_min = -5
y_max = 5
y_min = -5
x = np.linspace(x_min, x_max, 2000)
y = np.linspace(y_min, y_max, 2000)
X, Y = np.meshgrid(x, y)

n = 0

# Calculate positions
for i in range(len(star_1)):
    Z = X + (1j * Y)
    z1 = (1/Re)*star_1[i] #no unit !!
    z2 = (1/Re)*star_2[i] #no unit
    dwz = (E1/((Z-z1)**2))+(E2/((Z-z2)**2))

    # Take the magnitude or real part for plotting
    dwz_magnitude = np.abs(dwz)  # Use magnitude for critical curve

    # Plot where the magnitude hits a level — e.g., level=1
    fig = plt.figure()
    contour = plt.contour(X, Y, dwz_magnitude, levels=[1], colors='blue')

    light = []
    for collection in contour.collections:
        for path in collection.get_paths():
            vertices = path.vertices
            z_image = vertices[:, 0] + 1j * vertices[:, 1]
            lens_term = (E1 / (np.conj(z_image) - np.conj(z1))) + (E2 / (np.conj(z_image) - np.conj(z2)))
            z_source = z_image - lens_term
            light.append(z_source)#find light intersection (image is len plane for criticle, source is source plane for causetic)
            plt.scatter(np.real(z_source), np.imag(z_source), color='red', s=1)
    
    #path of source
    m_source = math.tan(math.radians(0))

    b_source_a = 0 #pass though both criticle and caustic curve
    b_source_b = 1.5 #pass though only criticle curve
    b_source_c = 4 #not pass though any criticle or caustic curve

    x_source = source_i_position + Vt * (time)

    y_source_a = m_source*x_source + b_source_a
    y_source_b = m_source*x_source + b_source_b
    y_source_c = m_source*x_source + b_source_c

    w_array_a = (x_source + (y_source_a*1j))*(1/Re)
    w_array_b = (x_source + (y_source_b*1j))*(1/Re)
    w_array_c = (x_source + (y_source_c*1j))*(1/Re)
    #plt.scatter(x_source[i],y_source[i], color='cyan', s=7)
    arrow_length = x_source[len(x_source)-1]-x_source[0]
    plt.arrow(-4, b_source_a, 8,0,  width=0.05, head_width=0.2, head_length=0.3, fc='#006182', ec='#006182')
    plt.arrow(-4, b_source_b, 8,0,  width=0.05, head_width=0.2, head_length=0.3, fc='#fd6b1b', ec='#fd6b1b')
    plt.arrow(-4, b_source_c, 8,0,  width=0.05, head_width=0.2, head_length=0.3, fc='#006e2b', ec='#006e2b')
    
    # plt.plot([-4, 4], [b_source_a, b_source_a], linestyle='--', color='#ab298c')
    # plt.plot([-4, 4], [b_source_b, b_source_b], linestyle='--', color='#13a837')
    # plt.plot([-4, 4], [b_source_c, b_source_c], linestyle='--', color='#00304e')

    #Plot critical + caustic
    plt.contour(X, Y, dwz_magnitude, levels=[1], colors='blue')
    #plt.title(f"tE = {(time[i]/tE+start_time):.2f}")
    plt.xlabel("x (AU)")
    plt.ylabel("y (AU)")
    # plt.figtext(0.5, 0.001, f"tE = {(time[i]/tE)-(end_time/2):.2f}", ha="center", fontsize=10)
    plt.gca().set_aspect("equal")
    plt.grid(True)
    #-(end_time/2)
    plt.savefig(f"D:\\blackhole_program\\graphmag\\CC_GW_fix_loss\\{(time[i]/tE):.4f}_{M}{m}.png")
    print(f"save {r_list[i]}")
    # plt.show()
    n += 1
    plt.close()