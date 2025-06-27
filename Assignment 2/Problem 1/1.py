import numpy as np
from scipy.optimize import fsolve

# Given data
components = ['IC4', 'nC4', 'IC5', 'nC5', 'C6', 'C7', 'C8', 'C9']
feed_flow = np.array([12, 448, 36, 15, 23, 39.1, 272.2, 31.0])  # kgmol/h

# Updated K-values at 80 psia (ensure these are accurate for the given conditions)
K_values = np.array([1.3, 1.05, 0.69, 0.6, 0.3, 0.2, 0.1, 0.05])

# Part a: Estimate N_min using the Fenske Equation
def fenske_equation(alpha_LK_HK, x_LK_D, x_HK_D, x_LK_B, x_HK_B):
    return np.log((x_LK_D / x_HK_D) / (x_LK_B / x_HK_B)) / np.log(alpha_LK_HK)

# Relative volatility for LK (nC4) and HK (IC5)
alpha_LK_HK = K_values[1] / K_values[2]

# Mole fractions in distillate and bottoms (example values, replace with actual calculations)
x_LK_D = 442 / 468
x_HK_D = 13 / 468
x_LK_B = 6 / 408.3
x_HK_B = 23 / 408.3

N_min = fenske_equation(alpha_LK_HK, x_LK_D, x_HK_D, x_LK_B, x_HK_B)
print(f"Minimum number of theoretical stages (N_min): {N_min:.2f}")

# Part b: Estimate product distributions using the Fenske Equation
def component_distribution(alpha_i, N_min, x_LK_D, x_HK_D, x_LK_B, x_HK_B):
    return alpha_i ** N_min * (x_LK_D / x_HK_D) / (x_LK_B / x_HK_B)

D_i = []
B_i = []
for i, alpha_i in enumerate(K_values):
    ratio = component_distribution(alpha_i, N_min, x_LK_D, x_HK_D, x_LK_B, x_HK_B)
    D_i.append(feed_flow[i] * ratio / (1 + ratio))
    B_i.append(feed_flow[i] / (1 + ratio))

print("\nDistillate and Bottoms Flow Rates:")
print("{:<10} {:<15} {:<15}".format("Component", "Distillate (D_i)", "Bottoms (B_i)"))
for i, comp in enumerate(components):
    print("{:<10} {:<15.2f} {:<15.2f}".format(comp, D_i[i], B_i[i]))

# Part c: Calculate minimum reflux using the Underwood Equation
def underwood_equation(theta, alpha, feed_flow, q=1):
    return sum((alpha * feed_flow) / (alpha - theta)) - (1 - q)

theta_initial_guess = 1.1
theta = fsolve(underwood_equation, theta_initial_guess, args=(K_values, feed_flow))[0]

R_min = sum((K_values * D_i) / (K_values - theta)) - 1
print(f"\nMinimum reflux (R_min): {R_min:.2f}")

# Part d: Estimate number of equilibrium stages using the Gilliland Chart
R = 1.3 * R_min
X = (R - R_min) / (R + 1)

# Example Gilliland correlation (simplified)
Y = 1 - np.exp((1 + 54.4 * X) / (11 + 117.2 * X) * (X - 1) / np.sqrt(X))
N = N_min + Y * (1 + N_min) / (1 - Y)

print(f"Number of equilibrium stages (N): {N:.2f}")

# Recommendations
print("\nRecommendations:")
print("1. Ensure that the K-values used are accurate for the given operating conditions (80 psia).")
print("2. Verify the implementation of the Underwood equation to ensure correct calculation of theta and R_min.")
print("3. Validate the feed composition and separation requirements to ensure realistic results.")