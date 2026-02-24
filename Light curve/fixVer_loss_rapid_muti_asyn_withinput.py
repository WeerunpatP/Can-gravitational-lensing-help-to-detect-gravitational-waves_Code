import numpy as np
import matplotlib.pyplot as plt
import math
import asyncio
import pandas as pd
from scipy.optimize import fsolve
from scipy import stats
import time as time_libery

time1 = time_libery.perf_counter()

def binary_lens_equation_system(vars, w, z1_c, z2_c):
    # Lens masses normalized
    """
    The system of real equations derived from the complex binary lens equation.
    f(z) = z - m1/(conj(z) - conj(z1)) - m2/(conj(z) - conj(z2)) - w = 0
    """
    x, y = vars
    z = complex(x, y)

    # Avoid division by zero by adding a small epsilon to the denominator if it's too close to zero
    epsilon = 1e-10

    denom1 = (z - z1_c)
    if denom1.real < epsilon and denom1.img < epsilon:
        term1 = 0 # Or handle as a singularity (for numerical stability, treat as zero contribution)
    else:
        term1 = (E1 / (z - z1_c))

    denom2 = (z - z2_c)
    if denom2.real < epsilon and denom2.img < epsilon:
        term2 = 0
    else:
        term2 = (E2 / (z - z2_c))
    
    #equation_complex = w - z + (E1 / (z - z1_c)) + (E2 / (z - z2_c))
    equation_complex = w - z + term1 + term2

    return [equation_complex.real, equation_complex.imag]

#input value
M = float(input("Enter M (big star mass) in solar mass: "))
m = float(input("Enter m (small star mass) in solar mass: "))
b_source_a = float(input("Enter source a impact parameter: "))
b_source_b = float(input("Enter source b impact parameter: "))
b_source_c = float(input("Enter source c impact parameter: "))
cheak_loss = str(input("Loss energy system or not (y/n): "))
# from colorama import Fore, Back, Style
# while not(cheak_loss == 'y' or cheak_loss == 'n'):
#     print("Please enter ", end='')
#     print(Fore.RED + "y", end='')
#     print(Style.RESET_ALL + " for yes it is binary system with loss energy or ", end='')
#     print(Fore.RED + "n", end='')
#     print(Style.RESET_ALL + " for no it the conserved energy system.")
#     cheak_loss = str(input("Loss energy system or not (y/n): "))

# move only object
# Constants
G = (6.6743 * (10**-11)) * (1.9891 * (10**30)) * (1 / ((1.495978707 * 10**11)**3)) * ((365.35 * 24 * 3600)**2)  # (AU^3)*(year^-2)*(solar mass^-1)
e = 0  # Eccentricity
c = (299792458) * (1 / (1.495978707 * 10**11)) * (365.25 * 24 * 3600)  # AU/year
ly = 9.4605284 * (10**15)  # m
AU = 149597871*(10**3) #m
Ds = (2.537*(10**6)*ly)/AU #AU, andromeda galaxy
Dl = (5*ly)/AU
M_solar = 1.9891*(10**30) #kg
q = m/M
print("q: ",q, " solar mass")
M_A = M+m
E1 = M/M_A
E2 = m/M_A

Re = math.sqrt((4*G*M_A*Dl*(Ds-Dl))/((c**2)*Ds))
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
print("a: ", a, "AU")
eps = 1e-10

# Derived quantities
rs = 2 * G * (M + m) / (c**2)  # Schwarzschild radius unit AU

# Parameters for lens
scale = 40
chrip_time = a**4/((32/5)*((m*M)/(m+M))*(rs**3)*c)
print(f"chrip at {chrip_time:,}")
far_from_chirp = float(input("Enter time before chirp in power of ten (year): "))
far_from_chirp = 10**(far_from_chirp)
start_time = chrip_time/far_from_chirp - (scale/2)
print(f"strat time from chrip : {chrip_time - start_time} year")
end_time = start_time+scale
time = np.arange(eps,scale*tE, 0.01*tE)

print(f"input in --> M = {M}, m = {m}.")
print(f"setting source impact parameter as --> {b_source_a}, {b_source_b}, {b_source_c}")

time_pd = pd.Series(time)  
num_parts = len(time)
step = int(num_parts/8)
parts_t = [
    time_pd[i:i+step]
    for i in range(0, num_parts, step)
]

#loss energy
async def binary_sim(time_part):
    r_list = []
    omega_list = []
    omega = math.sqrt((G*(M+m))/(a**3))
    
    for i, value in time_part.items():
        if cheak_loss == 'y':
            a_t = (a**4 - ((32/5)*((m*M)/(m+M))*(rs**3)*c*(value+start_time)))**(1/4)
        elif cheak_loss == 'n':
            a_t = a
        #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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

z1 = (X1 + Y1*1j)*(1/Re)
z2 = (X2 + Y2*1j)*(1/Re)

# ---------------------------------------------Sorce moving in linear motion--------------------------------------------------------------
# Set the source start position (AU)
source_i_position = -1*(Vt * ((scale/2)*tE))
print(f"source start at : {source_i_position}")

m_source = math.tan(math.radians(0))

x_source = source_i_position + (Vt * (time))

y_source_a = m_source*x_source + b_source_a
y_source_b = m_source*x_source + b_source_b
y_source_c = m_source*x_source + b_source_c

w_array_a = (x_source + (y_source_a*1j))*(1/Re)
w_array_b = (x_source + (y_source_b*1j))*(1/Re)
w_array_c = (x_source + (y_source_c*1j))*(1/Re)

# Define the search space for initial guesses
x_min, x_max = -7.0, 7.0
y_min, y_max = -7.0, 7.0
num_points = 30 # Increase density

x_guesses = np.linspace(x_min, x_max, num_points)
y_guesses = np.linspace(y_min, y_max, num_points)

#-----------------------------------megnification function------------------------------------------------------------
async def megnification_find(w_array):
    A_data = []
    for i, value in w_array.items():
        all_solutions = []
        for x_guess in x_guesses:
            for y_guess in y_guesses:
                initial_guess = [x_guess, y_guess]
                try:
                    # fsolve can return a warning for non-convergence, suppress for cleaner output
                    solution = fsolve(binary_lens_equation_system, initial_guess, args=(value, z1[i], z2[i]), xtol=1e-8, maxfev=2000)
                    
                    # Check for convergence and if the solution is valid (not NaN or inf)
                    if not np.any(np.isnan(solution)) and not np.any(np.isinf(solution)):
                        # Convert to complex number
                        sol_complex = complex(solution[0], solution[1])
                        all_solutions.append(sol_complex)

                except Exception as e:
                    # Catch other potential errors, but fsolve often handles non-convergence gracefully
                    pass

        # Filter out duplicate solutions
        unique_solutions = []
        tolerance = 1e-4 # Define a tolerance for comparing complex numbers
        for sol in all_solutions:
            is_unique = True
            for unique_sol in unique_solutions:
                if abs(sol - unique_sol) < tolerance:
                    is_unique = False
                    break
            if is_unique:
                unique_solutions.append(sol)

        A_list = []
        #calculated magnification
        #Path a
        if len(unique_solutions) > 0:
            for i, sol in enumerate(unique_solutions):
                    #print(f"Solution {i+1}: z = {sol:.6f}")
                    dwz = (E1/((sol-z1[i])**2))+(E2/((sol-z2[i])**2))
                    J = 1 - abs(dwz**2)
                    A = abs(1/J)
                    A_list.append(A)
        else:
            dwz = 0
            J = 1 - abs(dwz**2)
            A = abs(1/J)
            A_list.append(A)

        A_total = sum(A_list)
        A_total = abs(A_total)
        A_data.append(A_total)
        
    return pd.Series(A_data)

#--------------------asynchronous data input----------------------
w_array_a_pd = pd.Series(w_array_a)
w_array_b_pd = pd.Series(w_array_b)
w_array_c_pd = pd.Series(w_array_c)

num_parts = len(time)
step = int(num_parts/8)
parts_a = [w_array_a_pd[i:i+step] for i in range(0, num_parts, step)]
parts_b = [w_array_b_pd[i:i+step] for i in range(0, num_parts, step)]
parts_c = [w_array_c_pd[i:i+step] for i in range(0, num_parts, step)]

async def main(parts):
    tasks = [asyncio.create_task(megnification_find(part))for part in parts]
    results = await asyncio.gather(*tasks)
    return results

time_tranform = time/tE-(scale/2) # time that use to plot x-axis can be negative

#------------------------Run simulation-------------------------------------------
A_data_a = asyncio.run(main(parts_a))
time2 = time_libery.perf_counter()
print(f"Source A {b_source_a} Finished in {(time2 - time1)/60:.2f} minminutes")

A_data_b = asyncio.run(main(parts_b))
time3 = time_libery.perf_counter()
print(f"Source B {b_source_b} Finished in {(time3 - time2)/60:.2f} minutes")

A_data_c = asyncio.run(main(parts_c))
time4 = time_libery.perf_counter()
print(f"Source C {b_source_c} Finished in {(time4 - time3)/60:.2f} minutes")

time5 = time_libery.perf_counter()
print(f"Finished in {(time5 - time1)/60:.2f} minutes")

#--------------------------Trainform simulation data----------------------------------------
A_data_full_a = [item for sublist in A_data_a for item in sublist]
A_data_full_a = np.array(A_data_full_a)
A_data_full_b = [item for sublist in A_data_b for item in sublist]
A_data_full_b = np.array(A_data_full_b)
A_data_full_c = [item for sublist in A_data_c for item in sublist]
A_data_full_c = np.array(A_data_full_c)

#style define
if cheak_loss == 'y':
    lable = 'loss'
    color_lable = 'red'
elif cheak_loss == 'n':
    lable = "conserved"
    color_lable = 'green'
    
#graph plot
fig, ax = plt.subplots(figsize=(10, 6))
#plt.xlim(-20, -10)
ax.set_xlabel("Time (tE)", fontsize=16)
ax.set_ylabel("Magnification", fontsize=16)
ax.set_title(f"Light Curve of binary system of mass {M} and {m} M\u2609 with {lable} energy and fast rotation at -10^{int(math.log10(far_from_chirp))}", color=color_lable, fontsize=14)
ax.plot(time_tranform[:len(A_data_full_a)], A_data_full_a, label=f"impact parmeter: {b_source_a}", color="red", lw = '3', ls = '-') #Path a
ax.plot(time_tranform[:len(A_data_full_b)], A_data_full_b, label=f"impact parmeter: {b_source_b}", color='green', lw = '2', ls = '--') #Path b
ax.plot(time_tranform[:len(A_data_full_c)], A_data_full_c, label=f"impact parmeter: {b_source_c}", color="blue", lw = '1', ls = ':') #Path c
# ax.plot(time_tranform[:len(r_list)], r_list)
plt.legend() 
plt.grid()
plt.show()

#writing csv
data = {
    'time(tE)': list(time_tranform[:len(A_data_full_a)]),
    'Magnification_A': list(A_data_full_a),
    'Magnification_B': list(A_data_full_b),
    'Magnification_C': list(A_data_full_c),
}
df = pd.DataFrame(data)

df.to_csv('D:\\blackhole_program\\graphmag\\csv\\simu_data\\' + f'{M}{m}_{lable}_{far_from_chirp}_magdata.csv')

# from playsound import playsound
# from gtts import gTTS
# time2 = time_libery.perf_counter()
# voice_text = f"Program complete time usad {time2 - time1:.2f}"
# language = 'en'
# file_path = "D:\\blackhole_program\\graphmag\\Other\\complete_voice.mp3"
# print(voice_text)
# voice_obj = gTTS(text=voice_text, lang=language, slow=False)
# voice_obj.save(file_path)
# playsound(file_path)