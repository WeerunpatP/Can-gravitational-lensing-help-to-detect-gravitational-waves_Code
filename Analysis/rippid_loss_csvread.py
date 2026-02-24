import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks

G = (6.6743 * (10**-11)) * (1.9891 * (10**30)) * (1 / ((1.495978707 * 10**11)**3)) * ((365.35 * 24 * 3600)**2)
c = (299792458) * (1 / (1.495978707 * 10**11)) * (365.25 * 24 * 3600)  # AU/year

M = float(input("Enter M (big star mass) in solar mass: "))
m = float(input("Enter m (small star mass) in solar mass: "))
b_source_a = float(input("Enter source a impact parameter: "))
b_source_b = float(input("Enter source b impact parameter: "))
b_source_c = float(input("Enter source c impact parameter: "))
cheak_loss = str(input("Loss energy system or not (y/n): "))

M_A = M+m
tE = 60/365.25 #year
rs = 2 * G * (M + m) / (c**2)  # Schwarzschild radius unit AU
M_025 = M_A**0.25
tB = (tE/M_025)*6
a = (M_A*(tB**2))**(1/3)

chrip_time = a**4/((32/5)*((m*M)/(m+M))*(rs**3)*c)
print(f"chrip at {chrip_time:,}")
far_from_chirp = float(input("Enter time before chirp in power of ten (year): "))
far_from_chirp = 10**(far_from_chirp)
end_lim = int(input("End limit of graph : "))

if cheak_loss == 'y':
    lable = 'loss'
    color_lable = 'red'
elif cheak_loss == 'n':
    lable = "conserved"
    color_lable = 'green'

#improt csv
file_path = "D:\\blackhole_program\\graphmag\\csv\\simu_data\\" + f"{M}{m}_{lable}_{far_from_chirp}_magdata.csv"
print(file_path)
data = pd.read_csv(file_path)

Mag_A = data["Magnification_A"]
Mag_B = data["Magnification_B"]
Mag_C = data["Magnification_C"]

peaks_a, _ = find_peaks(Mag_A, height=1.05)
peak_values_a = Mag_A[peaks_a]
peaks_b, _ = find_peaks(Mag_B, height=1.05)
peak_values_b = Mag_B[peaks_b]
peaks_c, _ = find_peaks(Mag_C, height=1.05)
peak_values_c = Mag_C[peaks_c]

print(f"number of peak in line A : {len(peaks_a)}")
print(f"number of peak in line B: {len(peaks_b)}")
print(f"number of peak in line C : {len(peaks_c)}")
    
#graph plot
fig, ax = plt.subplots(figsize=(10, 6))
plt.xlim(-20, end_lim)
ax.set_yscale('log')
ax.set_xlabel("Time (tE)", fontsize=16)
ax.set_ylabel("Magnification", fontsize=16)
ax.set_title(f"Light Curve of binary system of mass {M} and {m} M\u2609 with {lable} energy and fast rotation", color=color_lable, fontsize=14)
# ax.scatter(data["time(tE)"][peaks_a], peak_values_a, color="orange")
# ax.scatter(data["time(tE)"][peaks_b], peak_values_b, color="orange")
# ax.scatter(data["time(tE)"][peaks_c], peak_values_c, color="orange")
ax.plot(data["time(tE)"], Mag_A, label=f"impact parmeter: {b_source_a}", color="red", lw = '3', ls = '-') #Path a
ax.plot(data["time(tE)"], Mag_B, label=f"impact parmeter: {b_source_b}", color='green', lw = '2', ls = '--') #Path b
ax.plot(data["time(tE)"], Mag_C, label=f"impact parmeter: {b_source_c}", color="blue", lw = '1', ls = ':') #Path c
# ax.plot(time_tranform[:len(r_list)], r_list)
plt.legend() 
plt.grid()
plt.show()