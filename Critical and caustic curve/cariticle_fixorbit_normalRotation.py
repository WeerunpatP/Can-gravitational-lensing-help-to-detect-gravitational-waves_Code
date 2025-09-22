import numpy as np
import matplotlib.pyplot as plt
import math
from sympy import symbols, Eq, solve, sqrt

# move only object
# Constants
G = (6.6743 * (10**-11)) * (1.9891 * (10**30)) * (1 / ((1.495978707 * 10**11)**3)) * ((24 * 3600)**2)  # (AU^3)*(day^-2)*(solar mass^-1)
l = math.sqrt(0.5)  # solar mass*AU^2/day
e = 0  # Eccentricity
c = (299792458) * (1 / (1.495978707 * 10**11)) * (24 * 3600)  # AU/day
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

Re = math.sqrt((4*G*M_A*Dl*(Ds-Dl))/((c**2)*Ds))
#tE = Re/Vt #year
tE = 60/365.25 #year
Vt = Re/tE

# ds = symbols('x')
# Ds_equation = Eq(sqrt((4*G*M_A*Dl*(ds-Dl))/((c**2)*ds)), Re)
# Ds = solve(Ds_equation, ds)
# Ds = Ds[0]
print("tE: ", tE, " year")
print("Vt: ", Vt, " AU/year")

Dls = Ds-Dl
print("Dl: ",Dl, " AU")
print("Ds: ", Ds, " AU")

a = l**2 / (G * M * m**2)
c = (3*(10**8))*(1/(1.495978707*10**11))*(24*3600)  # AU/day
g = 6.6743*(10**-11) #m3 kg-1 s-2
eps = 1e-10
scale = 1e5

Ai_data = []
A_data = []

# Derived quantities
mu = (M * m) / (M + m)  # solar mass
k = G * M * m  # N(AU)^2
r0 = l**2 / (G * M * m**2)  # Initial radius unit AU
print("a: ", a, " AU")
rc = 2 * G * (M + m) / (c**2)  # Schwarzschild radius unit AU
E = ((e**2) - 1) * mu * ((-k)**2) * (1 / (2 * (l**2)))
T = math.sqrt((a**3))  # year # change to day

# Parameters for lens
time = np.arange(0,4*tE, tE*0.1)
u0 = 206265/Dls 
omega = (2*math.pi)/T
print("T/tE: ", T/tE, " year")

r_list = []
omega_list = []

#conserved energy
for t in time:
    r_valuse = a / (1 + e * np.cos(omega*t)) # This updates r_valuse based on that omega and time
    r_list.append(r_valuse)
    omega_list.append(omega)

r_list = np.array(r_list)
omega_list = np.array(omega_list)

r1 = M / (M + m) * r_list
r2 = -1 * (m / (M + m)) * r_list
theta_values = omega_list*time
i_unitv = np.cos(theta_values)
j_unitv = np.sin(theta_values)

# Using the calculated values for X1, Y1, X2, Y2
X1=(r1*i_unitv) #AU
Y1=(r1*j_unitv) #AU
X2=(r2*i_unitv) #AU
Y2=(r2*j_unitv) #AU

star_1 = X1 + Y1*1j
star_2 = X2 + Y2*1j

#Relative line
Impact_parameter = u0  #acrscond
c_plot = (Impact_parameter/206265)*Dl

# Set to find distance between star (r)
theta_values = 0
n=0
R1 = []
R2 = []
X_cautic = []
Y_cautic = []
time_re = []
x1 = 0
y1 = 0
x2 = 0
y2 = 0
source_i_position = -2
x_max = 5
x_min = -5
y_max = 5
y_min = -5
x = np.linspace(x_min, x_max, 2000)
y = np.linspace(y_min, y_max, 2000)
X, Y = np.meshgrid(x, y)

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
    
    m_source = math.tan(math.radians(0.247))
    b_source = 0
    x_source = source_i_position + Vt * (time)
    y_source = m_source*x_source + b_source
    w_array = (x_source + (y_source*1j))*(1/Re)
    plt.scatter(x_source[i],y_source[i], color='cyan', s=7)
    #plt.plot(x_source,y_source,color='g')

    #Plot critical + caustic
    plt.contour(X, Y, dwz_magnitude, levels=[1], colors='blue')
    plt.title(f"Critical and Cuastic Curve of q: {q} system at tE = {time[i]/tE}")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.gca().set_aspect("equal")
    plt.grid(True)
    plt.savefig(f"C:\\Users\\LeZen_e595\\Desktop\\blackhole_program\\graphmag\\CC_GW_fix\\CC_noloss_at_{time[i]:.4f}_3629.png")
    print(f"save {time[i]}")
    #plt.show()
    n += 1
    plt.close()