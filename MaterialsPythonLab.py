# -*- coding: utf-8 -*-
"""
Welcome to the Shear, Moment, and Deflection Calculator for a Simply Supported Beam!
Code by Joe Lubin
CEE 3010 - Materials and Structures Lab
"""

#%%imports
import numpy as np
import matplotlib.pyplot as plt

#%%inputs

# Loads in pounds
P1 = float(input("Enter load P1 (lb): "))
P2 = float(input("Enter load P2 (lb): "))
P3 = float(input("Enter load P3 (lb): "))

# Distances (in inches)
a = float(input("Enter distance a (in) for load P1: "))
b = float(input("Enter distance b (in) for load P2: "))
c = float(input("Enter distance c (in) for load P3: "))

# Beam span in inches
L = float(input("Enter span L (in): "))

# Point along the beam at which deflection and rotation are calculated (in inches)
X = float(input("Enter point X (in) where deflection and rotation are to be calculated: "))

# Material properties
E = float(input("Enter modulus of elasticity E (psi): "))
I = float(input("Enter moment of inertia I (in^4): "))

#%%rxn forces
# Reaction Forces Calculation
# 
Rb = (P1 * a + P2 * b + P3 * c) / L  # Reaction at right support
Ra = (P1 + P2 + P3) - Rb              # Reaction at left support

#%%shear&moment calcs
dx = L / 10000.0                 #   Dividing beams into 1000 sections so that we can get smooth curve.
x = np.arange(0, L + dx, dx)    #   np.arrange (start, stop, step) - creates 1D array with evenly spaced values
M = np.zeros_like(x)            #   prepare array to store bending moment
V = np.zeros_like(x)            #   prepare array to store shear force

for i, xi in enumerate(x):
    if xi <= a:
        M[i] = Ra * xi
        V[i] = Ra
    elif xi <= b:
        M[i] = Ra * xi - P1 * (xi - a)
        V[i] = Ra - P1
    elif xi <= c:
        M[i] = Ra * xi - P1 * (xi - a) - P2 * (xi - b)
        V[i] = Ra - P1 - P2
    else:
        M[i] = Ra * xi - P1 * (xi - a) - P2 * (xi - b) - P3 * (xi - c)
        V[i] = Ra - P1 - P2 - P3

#%%shear&moment plots
plt.figure(figsize=(10, 8))         #create new figure space of wifth x height

# Bending Moment Diagram
plt.subplot(2, 1, 1)                #divides the figure into a grid of 2 rows and 1 column
plt.plot(x, M, linewidth=3)
plt.grid(True)
plt.title("Bending Moment Diagram")
plt.xlabel("Distance from Left Support [in]")
plt.ylabel("Bending Moment [lb-in]")

# Shear Force Diagram
plt.subplot(2, 1, 2)
plt.plot(x, V, 'r', linewidth=3)
plt.grid(True)
plt.title("Shear Force Diagram")
plt.xlabel("Distance from Left Support [in]")
plt.ylabel("Shear Force [lb]")

plt.tight_layout()                  # Adjusts the spacing between subplots so they dont overlap
plt.show()

#%%virtualwork deflection&rotation

# Compute the virtual moment diagram m_t and its derivative
m_t = np.zeros_like(x)
m_t_deriv = np.zeros_like(x)
for i, xi in enumerate(x):
    if xi <= X:
        m_t[i] = xi * (L - X) / L
        m_t_deriv[i] = (L - X) / L  # derivative when xi <= X
    else:
        m_t[i] = X * (L - xi) / L
        m_t_deriv[i] = -X / L       # derivative when xi > X

# Use the virtual work theorem to compute deflection and rotation at X
delta = np.trapz(M * m_t, x) / (E * I)              # x = np.arange(0, L + dx, dx) knows to take range 0 to L
theta = np.trapz(M * m_t_deriv, x) / (E * I)

# 
# Display Results
# 
print("\n--- Results at x = {:.2f} in ---".format(X))
print("Deflection, Δ = {:.6e} in".format(delta))
print("Rotation, θ  = {:.6e} rad".format(theta))