import math

def equilibrium_function(vapor_fraction, equilibrium_constant, feed_composition):
    return (equilibrium_constant - 1) * feed_composition / (vapor_fraction * (equilibrium_constant - 1) + 1)

def equilibrium_function_derivative(vapor_fraction, equilibrium_constant, feed_composition):
    return -(equilibrium_constant - 1) ** 2 * feed_composition / ((vapor_fraction * (equilibrium_constant - 1) + 1) ** 2)

def liquid_phase_composition(vapor_fraction, equilibrium_constant, feed_composition):
    return feed_composition / (vapor_fraction * (equilibrium_constant - 1) + 1)

def vapor_phase_composition(vapor_fraction, equilibrium_constant, feed_composition):
    return equilibrium_constant * feed_composition / (vapor_fraction * (equilibrium_constant - 1) + 1)

def calculate_equilibrium():
    feed_flow_rate = 100
    feed_compositions = [0.6, 0.25, 0.15]
    total_pressure = 1.01325  # bar
    
    temperature_celsius = float(input("Input temperature in Celsius: "))
    temperature_kelvin = temperature_celsius + 273.18
    
    # Antoine equation coefficients for Benzene, Toluene, and o-Xylene
    antoine_A = [4.01814, 4.14157, 4.12928]
    antoine_B = [1203.835, 1377.578, 1478.244]
    antoine_C = [-53.226, -50.507, -59.076]
    
    equilibrium_constants = []
    for A, B, C in zip(antoine_A, antoine_B, antoine_C):
        vapor_pressure = 10 ** (A - (B / (C + temperature_kelvin)))
        equilibrium_constants.append(vapor_pressure / total_pressure)
    
    print("Equilibrium constants:")
    for K in equilibrium_constants:
        print(K)
    
    vapor_fraction = float(input("Choose initial guess for vapor fraction: "))
    prev_value = 0
    tolerance = 1e-6
    max_iterations = 100
    
    # Checking valid range
    fa, fb = 0, 0
    for K, z in zip(equilibrium_constants, feed_compositions):
        fa += equilibrium_function(0.0, K, z)
        fb += equilibrium_function(1.0, K, z)
    
    if fa * fb > 0:
        print("No valid root found in [0,1]. Check equilibrium constants and try a different temperature.")
        return
    
    for iteration in range(max_iterations):
        sum_function = sum(equilibrium_function(vapor_fraction, K, z) for K, z in zip(equilibrium_constants, feed_compositions))
        sum_derivative = sum(equilibrium_function_derivative(vapor_fraction, K, z) for K, z in zip(equilibrium_constants, feed_compositions))
        
        if sum_derivative == 0:
            print("Derivative is zero. Choose another initial guess.")
            return
        
        if abs(vapor_fraction - prev_value) < tolerance:
            print(f"Vapor fraction found in {iteration} iterations: {vapor_fraction:.6f}")
            break
        
        prev_value = vapor_fraction
        vapor_fraction -= sum_function / sum_derivative
    else:
        print("Max iterations reached. Solution may not be accurate.")
    
    vapor_flow_rate = vapor_fraction * feed_flow_rate
    liquid_flow_rate = feed_flow_rate - vapor_flow_rate
    
    liquid_compositions = [liquid_phase_composition(vapor_fraction, K, z) for K, z in zip(equilibrium_constants, feed_compositions)]
    vapor_compositions = [vapor_phase_composition(vapor_fraction, K, z) for K, z in zip(equilibrium_constants, feed_compositions)]
    
    print(f"Vapor flow rate: {vapor_flow_rate:.2f}")
    print(f"Liquid flow rate: {liquid_flow_rate:.2f}")
    
    components = ["Benzene", "Toluene", "o-Xylene"]
    for i, component in enumerate(components):
        print(f"{component} - Liquid Phase: {liquid_compositions[i]:.6f}, Vapor Phase: {vapor_compositions[i]:.6f}")

if __name__ == "__main__":
    calculate_equilibrium()
