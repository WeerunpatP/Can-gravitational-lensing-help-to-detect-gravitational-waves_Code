import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.optimize import fsolve

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
M = 30 # solar mass
m = M*0.25 # solar mass
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
time = np.arange(0,10*tE, (tE*0.1)+eps)
u0 = 206265/Dls 
omega = (2*math.pi)/tB

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

z1 = (X1 + Y1*1j)*(1/Re)
z2 = (X2 + Y2*1j)*(1/Re)

#Relative line
Impact_parameter = u0  #acrscond
c_plot = (Impact_parameter/206265)*Dl

# Set to find distance between star (r)
source_i_position = -5
time= np.arange(0,4*tE, (tE*0.001)+eps)

#plot light curve----------------------------------------------------------------

# Slightly perturb z1 and z2 to avoid perfect symmetry issues for some solvers
m_source = math.tan(math.radians(0.247))
b_source = 100
x_source = source_i_position + Vt * (time)
y_source = m_source*x_source + b_source
w_array = x_source + (y_source*1j)

# Define the search space for initial guesses
x_min, x_max = -5.0, 5.0
y_min, y_max = -5.0, 5.0
num_points = 30 # Increase density

x_guesses = np.linspace(x_min, x_max, num_points)
y_guesses = np.linspace(y_min, y_max, num_points)

for i in range(len(w_array)):
    all_solutions = []
    for x_guess in x_guesses:
        for y_guess in y_guesses:
            initial_guess = [x_guess, y_guess]
            try:
                # fsolve can return a warning for non-convergence, suppress for cleaner output
                solution = fsolve(binary_lens_equation_system, initial_guess, args=(w_array[i], z1[i], z2[i]), xtol=1e-8, maxfev=2000)
                
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
    tolerance = 1e-2 # Define a tolerance for comparing complex numbers
    for sol in all_solutions:
        is_unique = True
        for unique_sol in unique_solutions:
            if abs(sol - unique_sol) < tolerance:
                is_unique = False
                break
        if is_unique:
            unique_solutions.append(sol)

    A_list = []
    #print("\nAll possible complex roots z:")
    if len(unique_solutions) > 0:
        for i, sol in enumerate(unique_solutions):
            print(f"Solution {i+1}: z = {sol:.6f}")
            dwz = (E1/((sol-z1[i])**2))+(E2/((sol-z2[i])**2))
            abs_dwz = math.sqrt((dwz.real**2)+(dwz.imag**2))
            J = 1 - (abs_dwz**2)
            A = 1/J
            A_list.append(A)
    else:
        print(f"Solution {i+1}: z = {sol:.6f}")
        dwz = (E1/((sol-z1[i])**2))+(E2/((sol-z2[i])**2))
        abs_dwz = math.sqrt((dwz.real**2)+(dwz.imag**2))
        J = 1 - (abs_dwz**2)
        A = 1/J
        A_list.append(A)

    A_total = sum(A_list)
    A_total = abs(A_total)
    A_data.append(A_total)
    #print(f"\nFound {len(unique_solutions)} unique solutions. Whan source at {w_array[i]}")
    print(A_total)

    # # Verify solutions by plugging back into the equation (optional)
    # print("\nVerification (residuals should be close to zero):")
    # for i, sol in enumerate(unique_solutions):
    #     x_sol, y_sol = sol.real, sol.imag
    #     residual_real, residual_imag = binary_lens_equation_system((x_sol, y_sol), w_c, z1_c, z2_c)
    #     #print(f"Solution {i+1} residual: Real={residual_real:.2e}, Imag={residual_imag:.2e}")

# Plotting
plt.figure(figsize=(10, 6))
#min_length = min(len(time), len(A_data[i]))  # Find the smallest length
plt.plot(time/tE-5, A_data, label=f"Magnification_{M}")
plt.xlabel("Time (tE)")
plt.ylabel("Magnification")
plt.title(f"Light Curve (No Loss Energy)")
plt.legend()
plt.grid()
plt.show()
plt.figure(figsize=(10, 6))