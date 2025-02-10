import numpy as np

def rachford_rice(vapor_fraction, equilibrium_constants, mole_fractions):
    """Compute Rachford-Rice function value."""
    return np.sum((equilibrium_constants - 1) * mole_fractions / (vapor_fraction * (equilibrium_constants - 1) + 1))

def rachford_rice_derivative(vapor_fraction, equilibrium_constants, mole_fractions):
    """Compute derivative of the Rachford-Rice function."""
    denominator = (vapor_fraction * (equilibrium_constants - 1) + 1) ** 2
    return -np.sum((equilibrium_constants - 1) ** 2 * mole_fractions / denominator)

def solve_vapor_fraction(equilibrium_constants, mole_fractions, tol=1e-6, max_iter=100):
    """Solve for vapor fraction using Newton-Raphson method."""
    vapor_fraction = 0.5  # Initial guess
    for _ in range(max_iter):
        f_val = rachford_rice(vapor_fraction, equilibrium_constants, mole_fractions)
        df_val = rachford_rice_derivative(vapor_fraction, equilibrium_constants, mole_fractions)
        
        if abs(df_val) < tol:
            return None  # Derivative too small, avoid division by zero
        
        new_vapor_fraction = vapor_fraction - f_val / df_val
        if abs(new_vapor_fraction - vapor_fraction) < tol:
            return new_vapor_fraction
        
        vapor_fraction = new_vapor_fraction
    
    return None  # No convergence

def compute_phase_composition(vapor_fraction, equilibrium_constants, mole_fractions):
    """Compute liquid and vapor phase compositions."""
    liquid_composition = mole_fractions / (vapor_fraction * (equilibrium_constants - 1) + 1)
    vapor_composition = equilibrium_constants * liquid_composition
    return liquid_composition, vapor_composition

def compute_equilibrium_constants(temp_K, pressure_bar):
    """Compute equilibrium constants using Antoine equation."""
    A = np.array([4.01814, 4.14157, 4.12928])  # Benzene, Toluene, o-Xylene
    B = np.array([1203.835, 1377.578, 1478.244])
    C = np.array([-53.226, -50.507, -59.076])
    
    log_p_sat = A - (B / (C + temp_K))
    p_sat = 10 ** log_p_sat
    return p_sat / pressure_bar

def optimize_temperature(feed_flow=100, pressure_bar=1.01325, temp_range=(100, 105), step=0.0001):
    """Find the optimal temperature for maximum benzene recovery."""
    mole_fractions = np.array([0.6, 0.25, 0.15])  # Benzene, Toluene, o-Xylene
    best_temp = None
    max_recovery = -1
    
    for temp_C in np.arange(temp_range[0], temp_range[1], step):
        temp_K = temp_C + 273.18  # Convert to Kelvin
        equilibrium_constants = compute_equilibrium_constants(temp_K, pressure_bar)
        vapor_fraction = solve_vapor_fraction(equilibrium_constants, mole_fractions)
        
        if vapor_fraction is None:
            continue
        
        vapor_flow = vapor_fraction * feed_flow
        liquid_flow = feed_flow - vapor_flow
        _, vapor_composition = compute_phase_composition(vapor_fraction, equilibrium_constants, mole_fractions)
        
        recovery_benzene = (vapor_flow * vapor_composition[0]) / (feed_flow * mole_fractions[0])
        
        if recovery_benzene > max_recovery:
            max_recovery = recovery_benzene
            best_temp = temp_C

    return best_temp, max_recovery * 100

if __name__ == "__main__":
    optimal_temp, max_recovery = optimize_temperature()
    print(f"Optimal temperature for maximum benzene recovery: {optimal_temp:.4f}°C")
    print(f"Maximum benzene recovery: {max_recovery:.2f}%")

