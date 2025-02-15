import numpy as np
from scipy.optimize import fsolve
import CoolProp.CoolProp as CP

# Constants
epsilon = 1e-6  # Convergence criterion
R = 8.314  # Universal gas constant in J/(mol·K)
P_total = 400 * 101325 / 14.696  # Total pressure in Pa (400 psia converted to Pa)

# Antoine constants for components (A, B, C)
antoine_constants = {
    'C1': {'A': 4.22057, 'B': 516.689, 'C': -14.170},  # Methane
    'C2': {'A': 4.50798, 'B': 944.615, 'C': -16.719},  # Ethane
    'C3': {'A': 4.53678, 'B': 1149.36, 'C': -24.906},  # Propane
    'nC4': {'A': 4.35597, 'B': 1305.198, 'C': -32.445},  # n-Butane
    'nC5': {'A': 4.29152, 'B': 1436.014, 'C': -41.485},  # n-Pentane
}

# Feed and absorbent oil data (lbmol/h)
feed = {
    'C1': 160.0,
    'C2': 370.0,
    'C3': 240.0,
    'nC4': 25.0,
    'nC5': 5.0
}

# Number of stages and components
n_stages = 13
components = ['C1', 'C2', 'C3', 'nC4', 'nC5']

# Initial guesses for stage temperatures, vapor, and liquid flow rates
initial_guesses = {
    'stage1': {'T': 150, 'V': 530, 'L': 700},
    'stage13': {'T': 350, 'V': 600, 'L': 770}
}

# Function to calculate vapor pressure using Antoine equation
def vapor_pressure(A, B, C, T):
    return 10 ** (A - B / (T + C)) * 133.322  # Convert mmHg to Pa

# Function to calculate enthalpy using CoolProp
def calculate_enthalpy(component, T, P, phase):
    fluid = {
        'C1': 'Methane',
        'C2': 'Ethane',
        'C3': 'Propane',
        'nC4': 'Butane',
        'nC5': 'Pentane'
    }[component]
    if phase == 'liquid':
        return CP.PropsSI('H', 'T', T + 273.15, 'P', P, fluid)  # Liquid enthalpy
    elif phase == 'vapor':
        return CP.PropsSI('H', 'T', T + 273.15, 'P', P, fluid)  # Vapor enthalpy

# Material balance equation
def material_balance(V, L, F, gamma_V, gamma_L):
    return (1 + gamma_V) * V + (1 + gamma_L) * L - F

# Equilibrium balance equation
def equilibrium_balance(K, L, V):
    return K * (L / np.sum(L)) - (V / np.sum(V))

# Enthalpy balance equation (Fixed F handling)
def enthalpy_balance(h_L, h_V, h_F, L, V, F, Q, gamma_L, gamma_V):
    return h_L * (1 + gamma_L) * np.sum(L) + h_V * (1 + gamma_V) * np.sum(V) - h_V * np.sum(V) - h_L * np.sum(L) - h_F * sum(F.values()) + Q

# Main function to solve the SC method
def solve_sc_method():
    # Initialize variables
    T = np.zeros(n_stages)  # Stage temperatures
    V = np.zeros((n_stages, len(components)))  # Vapor flow rates
    L = np.zeros((n_stages, len(components)))  # Liquid flow rates
    gamma_V = np.zeros(n_stages)  # Vapor draw ratio
    gamma_L = np.zeros(n_stages)  # Liquid draw ratio

    # Set initial guesses
    T[0] = initial_guesses['stage1']['T']
    T[-1] = initial_guesses['stage13']['T']
    V[0, :] = initial_guesses['stage1']['V']
    V[-1, :] = initial_guesses['stage13']['V']
    L[0, :] = initial_guesses['stage1']['L']
    L[-1, :] = initial_guesses['stage13']['L']

    # Iterative solution
    while True:
        # Calculate vapor pressures
        for i, comp in enumerate(components):
            A, B, C = antoine_constants[comp].values()
            for j in range(n_stages):
                V[j, i] = vapor_pressure(A, B, C, T[j])

        # Solve material balance
        for j in range(n_stages):
            for i, comp in enumerate(components):
                F_ij = feed[comp] if j == 0 else 0
                M_ij = material_balance(V[j, i], L[j, i], F_ij, gamma_V[j], gamma_L[j])
                L[j, i] += 0.01 * M_ij  # Adjust L based on M_ij
                V[j, i] -= 0.01 * M_ij  # Adjust V based on M_ij

        # Solve equilibrium balance
        for j in range(n_stages):
            for i, comp in enumerate(components):
                if L[j, i] > 0:
                    K_ij = V[j, i] / L[j, i]
                    E_ij = equilibrium_balance(K_ij, L[j, i], V[j, i])
                    L[j, i] += 0.01 * E_ij  # Adjust L based on E_ij
                    V[j, i] -= 0.01 * E_ij  # Adjust V based on E_ij

        # Solve enthalpy balance
        for j in range(n_stages):
            h_L = sum(calculate_enthalpy(comp, T[j], P_total, 'liquid') for comp in components)
            h_V = sum(calculate_enthalpy(comp, T[j], P_total, 'vapor') for comp in components)
            h_F = sum(calculate_enthalpy(comp, T[j], P_total, 'liquid') for comp in components)  # Feed is liquid
            Q_j = 0  # Assuming no heat duty for now
            H_j = enthalpy_balance(h_L, h_V, h_F, L[j], V[j], feed, Q_j, gamma_L[j], gamma_V[j])
            T[j] += 0.01 * H_j  # Adjust T based on H_j

        # Check for convergence
        tau = np.sum(H_j*2 + np.sum(M_ij2 + E_ij*2))
        beta = sum(feed.values()) * n_stages * (2 * len(components) + 1)
        if tau / beta < epsilon:
            break

    return T, V, L

# Run the solver
T, V, L = solve_sc_method()

# Output results
print("Stage Temperatures (K):", T)
print("Vapor Flow Rates (lbmol/h):", V)
print("Liquid Flow Rates (lbmol/h):", L)