# -*- coding: utf-8 -*-
""" 
FINAL PROJECT
Jack Olson
""" 
import numpy as np 
import scipy as sp  
import os
import pandas as pd
import matplotlib.pyplot as plt
#Givens
L =(0.030) #(m) length 
A = (25e-6) #(m^2) cross-sectional area
P = 0.022 #(m) perimeter 
T_Amb = 298.15 #(K) Ambient Temperature  
T_b = 793.15 #(K) base temperature
#inportig and reading excel
def Import_Excel (): 
   file_path = 'Cutting_Tool_data.xlsx'
   df1 =pd.read_excel(file_path, sheet_name='k_vs_phi')
   df2 =pd.read_excel(file_path,sheet_name='h_vs_V')
   phi = df1['Carbon_fraction_wt_percent'].to_numpy()
   k = df1['Thermal_conductivity_W_mK'].to_numpy()
   Air_Speed = df2['Air_speed_m_s'].to_numpy() 
   h = df2['h_W_m2K'].to_numpy() 
   return phi,k,Air_Speed, h
phi,k,Air_Speed,h = Import_Excel ()

def g_func (count,xf): 
    if count == 0 :
        g = 1.0 
    elif count == 1 :
        g = xf 
    elif count == 2: 
        g = xf**2 
    #safe call
    else: 
        g = 0.0
    return g
#-------------------Part 1b--------------------------------------------------
def mat_prop_fit(phi, k): 
    m =2 
    n = len(k)-1 
    A = np.zeros((m+1,m+1))
    B = np.zeros(m+1) 
    #A
    for K in range(m+1): 
        for j in range(m+1): 
            for i in range(n+1): 
                A[K,j] += g_func(K,phi[i])*g_func(j,phi[i])
                
    for K_B in range(m+1): 
        for i_B in range (n+1): 
            B[K_B] += g_func(K_B,phi[i_B])*(k[i_B]) 
    a_bar = sp.linalg.solve(A,B) 
    return a_bar
a_bar = mat_prop_fit(phi, k)  
c,b,a = a_bar 


#evaluate to fit 
def k_of_phi(phi: float, a: np.ndarray) -> float:
    # k = a0 + a1*phi + a2*phi^2
    return a[0] + a[1]*phi + a[2]*phi**2

# user input for (phi) 
phi_user = float(input("Enter phi: "))
k_chosen = k_of_phi(phi_user, a_bar)
print(f"fitted quad:(a,b,c)(phi)={c:.3f}+{b:.3f}*phi+{a:.3f}*phi^2")


#---------------------------#part 1b-----------------------------------------
#convert to log
logV = np.log(Air_Speed)
logh = np.log(h)
#was told through email polyfit is allowed on final 
m, log_alpha = np.polyfit(logV, logh, 1)   # y=logh, x=logV 
#convert back from log
alpha = np.exp(log_alpha)

#user input
V_chosen = float(input("enter a speed : "))

#puts user entered value into equation
h_chosen = alpha * V_chosen**m
print("m=",m)
print("alpha =",alpha)
print("h =", h_chosen)

#-------------------Part 2---------------------------------------------------
# .005 step size

def solve_part2_temperature(k, h, L, Area, Perim, Tb, Tinf, dx=None, N=None):
    # grid setup and ensure constraints
    if dx is None and N is None:
        N = 101  # default choice
    if dx is not None and N is not None:
        raise ValueError("Specify either dx or N, not both.")
    if dx is not None:
        N_float = L/dx + 1.0
        if abs(N_float - round(N_float)) > 1e-12:
            raise ValueError("dx must divide L exactly so that L/dx is an integer.")
        N = int(round(N_float))
    else:
        dx = L/(N-1)

    x = np.linspace(0.0, L, N)

    #build linear system M*T = b
    M = np.zeros((N, N), dtype=float)
    b = np.zeros(N, dtype=float)

    # Dirichlet: T1 = Tb
    M[0, 0] = 1.0
    b[0] = Tb

    c = k * Area / dx**2

    # Interior nodes i = 2..N-1  -> indices 1..N-2
    for i in range(1, N-1):
        M[i, i-1] = c
        M[i, i]   = -2.0*c - h*Perim
        M[i, i+1] = c
        b[i]      = -h*Perim*Tinf

    # (8) Tipat x=L 
    M[N-1, N-2] = k*Area/dx
    M[N-1, N-1] = -(k*Area/dx) - h*Perim
    b[N-1]      = -h*Perim*Tinf

    T = sp.linalg.solve(M, b)

    # --- required outputs: max T and max slope (central differences on data) ---
    Tmax = float(np.max(T))

    # central difference slopes on interior nodes only 
    dTdx_interior = (T[2:] - T[:-2])/(2.0*dx)   # corresponds to nodes 2-N-1
    max_slope = float(np.max(np.abs(dTdx_interior)))

    # --- diagnostics---
    lin_res = float(np.linalg.norm(M @ T - b, ord=np.inf))
    tip_res = float((k*Area)*(T[-2]-T[-1])/dx - (h*Perim)*(T[-1]-Tinf))  # should be ~0 if Eq. (8) satisfied

    return x, T, dx, Tmax, max_slope, lin_res, tip_res


# ---------------- Part 2 driver ---------------------------------------------

x, T, dx, Tmax, max_slope, lin_res, tip_res = solve_part2_temperature(
    k=k_chosen, h=h_chosen,
    L=L, Area=A, Perim=P,
    Tb=T_b, Tinf=T_Amb,
    dx=0.005  #chosen 
)

print(f"dx = {dx} m")
print(f"Max T = {Tmax} K")
print(f"Max |dT/dx| (central, interior) = {max_slope} K/m")
plt.figure()
plt.plot(x, T)
plt.xlabel("x (m)")
plt.ylabel("T (K)")
plt.title("Temperature distribution T(x)")
plt.grid(True)
plt.show()
#-------------------------------part 3--------------------------------------- 
out_xlsx = "part3_results.xlsx"  #excel file for data out

rows = []

for phi_val in phi:                 # phi array from excel()
    k_val = k_of_phi(phi_val, a_bar)

    for V_val in Air_Speed:         # Air Speed array from excel
        h_val = alpha * (V_val ** m)

        x, T, dx, Tmax, max_slope, lin_res, tip_res = solve_part2_temperature(
            k=k_val, h=h_val,
            L=L, Area=A, Perim=P,
            Tb=T_b, Tinf=T_Amb,
            dx=0.005
        )

        rows.append([phi_val, V_val, k_val, h_val, Tmax, max_slope])

df = pd.DataFrame(rows, columns=["phi", "V", "k", "h", "Tmax", "max_abs_dTdx"])
df.to_excel(out_xlsx, index=False)

print("Saved:", os.path.abspath(out_xlsx))

print(os.getcwd())

