import numpy as np
import matplotlib.pyplot as plt
import math
from sympy import symbols, Eq, solve, sqrt
import asyncio
import pandas as pd

# move only object
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

Re = math.sqrt((4*G*M_A*Dl*(Ds-Dl))/((c**2)*Ds))
#tE = Re/Vt #year
tE = 60/365.25 #year
Vt = Re/tE

print("tE: ", tE, " year")
print("Vt: ", Vt, " AU/year")

Dls = Ds-Dl
print("Dl: ",Dl, " AU")
print("Ds: ", Ds, " AU")

M_025 = M_A**0.25
tB = (tE/M_025)*6
a = (M_A*(tB**2))**(1/3)
eps = 1e-10

Ai_data = []
A_data = []

# Derived quantities
print("a: ", a, " AU")
rs = 2 * G * (M + m) / (c**2)  # Schwarzschild radius unit AU
T = math.sqrt((a**3))  # year # change to day

# Parameters for lens
scale = 40
chrip_time = a**4/((32/5)*((m*M)/(m+M))*(rs**3)*c)
far_from_chirp = 10**8
start_time = chrip_time - scale
print(f"start time : {start_time}")
end_time = chrip_time
time = np.arange(eps,(end_time-start_time), (10**6)*tE)
print(f"time: {time}")

time_pd = pd.Series(time)  
num_parts = len(time)
step = int(num_parts/8)
print(step)
parts_t = [
    time_pd[i:i+step]
    for i in range(0, num_parts, step)
]
print(len(parts_t))

#loss energy
async def binary_sim(time_part):
    r_list = []
    omega_list = []
    omega = math.sqrt((G*(M+m))/(a**3))
    
    for i, value in time_part.items():
        a_t = (a**4 - ((32/5)*((m*M)/(m+M))*(rs**3)*c*(value+start_time)))**(1/4)
        if not(isinstance(a_t, complex)) and a_t > 0:
            r_valuse = a_t / (1 + e * np.cos(omega*(value+start_time))) # This updates r_valuse based on that omega and time
        else:
            break
        r_list.append(r_valuse)
        if (G*(M+m))/(r_valuse**3) > 0:
            omega = math.sqrt((G*(M+m))/(r_valuse**3))
        else:
            break
        omega_list.append(omega)

    r_list = np.array(r_list)
    omega_list = np.array(omega_list)
    return [r_list, omega_list]

async def main_bs(parts):
    tasks = [asyncio.create_task(binary_sim(part))for part in parts]
    results = await asyncio.gather(*tasks)
    return results

result_bs = asyncio.run(main_bs(parts_t))

r_list = [res[0] for res in result_bs]
omega_list = [res[1] for res in result_bs]

r_list = np.concatenate(r_list)
omega_list = np.concatenate(omega_list)

print(len(r_list))
print(len(omega_list))

r1 = M / (M + m) * r_list
r2 = -1 * (m / (M + m)) * r_list
theta_values = omega_list*time[:len(omega_list)]
i_unitv = np.cos(theta_values)
j_unitv = np.sin(theta_values)

# Using the calculated values for X1, Y1, X2, Y2
X1=(r1*i_unitv) #AU
Y1=(r1*j_unitv) #AU
X2=(r2*i_unitv) #AU
Y2=(r2*j_unitv) #AU

star_1 = X1 + Y1*1j
star_2 = X2 + Y2*1j

cta = 2 * omega_list*(time[:len(omega_list)]+start_time)
eta = (M*m)/((M+m)**2)
h_strain = ((eta*((G*(M+m))**(5/3))*(omega_list**(2/3))) / (2*(c**4)*r_list)) *  np.cos(cta)
# M_chirp = ((M*m)**(3/5))/((M+m)**(1/5))
# f_gw = (((8*math.pi)**(8/3))/5)*(((G*M_chirp)/(c**3))**5/3)*time[:len(omega_list)]
# omega_strain = 2*math.pi*f_gw
# h0 = (eta*((G*(M+m))*(5/3))*(omega_list*(2/3)))/(2*(c**4)*r_list)
# h_strain = h0*np.cos(2*omega_list*time[:len(omega_list)])

# async def strain_signal(time_part):
#     h_strain = []
#     for i, t in enumerate(time_part):
#         cta = 2 * omega_list[i]*(t+start_time)
#         nta = 0.247
#         try:
#             ht = ((nta*((G*(M+m))**(5/3))*(omega_list[i]**(2/3))) / (2*(c**4)*r_list[i])) *  math.cos(cta)
#         except ValueError:
#             ht = None
#         # Mc = (m*M)**(3/5) / (m*M)**(1/5)
#         # f = omega_list[i]/(math.pi*2) # year-1
#         # ft = (96/5) * (math.pi**(8/3)) * (((G*Mc)/(c**3))**(5/3)) * f
#         # cta = omega_list[i]*(t+start_time)
#         # try:
#         #     ht = ((4*((G*Mc)**(5/3)))/((c**4)*r_list[i])) * ((math.pi*ft)**(2/3)) * math.cos(cta)
#         # except ValueError:
#         #     ht = None
#         #     print("Can't find")
#         # h_strain.append(ht)
#     #print(f"{(i/len(time))*100:.2f}%")
#     return h_strain

# async def main(parts):
#     tasks = [asyncio.create_task(strain_signal(part))for part in parts]
#     results = await asyncio.gather(*tasks)
#     return results

# h_strain= asyncio.run(main(parts_t))
# h_strain= [item for sublist in h_strain for item in sublist]
# h_strain = np.array(h_strain)
print(len(h_strain))

#Plot critical + caustic
fig, ax = plt.subplots(figsize=(10, 5))
ax.set_xscale('log')
plt.title(f"Strain signal from mass {M} and {m} M\u2609") #fast rotation binary system with
plt.xlabel("time (tE)")
plt.ylabel("strain signal")
ax.plot(time[:len(h_strain)], h_strain)
plt.show()
plt.close()