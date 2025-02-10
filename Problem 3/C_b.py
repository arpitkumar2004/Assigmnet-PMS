# Q3: BUBBLE POINT METHOD – MULTICOMPONENT MULTISTAGE DISTILLATION
""" 
1000 kmol/h of a saturated-liquid mixture of 60% methanol (normal boiling point 65 C), 20 mol% ethanol 
(normal boiling point 98 C), and 20 mol% n-propanol (normal boiling point 97 C) is fed to the middle 
stage of a distillation column having three equilibrium stages, a total condenser, a partial reboiler, all 
operated at 1 atm. The distillate rate is 600 kmol/h, and the external reflux rate is 2,000 kmol/h of 
saturated liquid.
 """
# Write your own code to obtain converged solutions of the temperature profiles.

import numpy as np
from typing import Tuple

class Distillation:
    def __init__(self, feed_rate: int = 1000, distillate_rate: int = 600, reflux_rate: int = 2000, stages: int = 3, pressure: int = 1, components: list = ['methanol', 'ethanol', 'n-propanol'], feed_composition: list = [0.6, 0.2, 0.2], boiling_points: list = [65, 98, 97], convergence_tolerance: float = 1e-6, max_iterations: int = 1000):
        """ 
        Initialize the distillation column parameters.
        
        Parameters:
        feed_rate : int
            Feed rate of the mixture in kmol/h.
        distillate_rate : int
            Distillate rate in kmol/h.
        reflux_rate : int
            Reflux rate in kmol/h.
        stages : int
            Number of equilibrium stages.
        pressure : int
            Operating pressure in atm.
        components : list
            List of component names.
        feed_composition : list 
            List of feed composition in mole fractions.
        boiling_points : list
            List of boiling points for each component.
        temperature_profile : np.array
            Temperature profile for the distillation column.
        convergence_tolerance : float
            Convergence tolerance for the solution.
        max_iterations : int
            Maximum number of iterations for the solution.

        Returns:
        None
        """

        self.feed_rate = feed_rate
        self.distillate_rate = distillate_rate
        self.bottom_pdt_rate = feed_rate - distillate_rate
        self.reflux_rate = reflux_rate
        self.stages = stages
        self.pressure = pressure
        self.components = components
        self.num_components = len(components)
        self.feed_composition = feed_composition
        self.boiling_points = boiling_points
        self.convergence_tolerance = convergence_tolerance
        self.max_iterations = max_iterations
        self.temperature_profile = np.zeros((self.stages, len(self.components)))



    def calculate_Pij(self, C1: np.array, C2: np.array, C3: np.array, C4: np.array, C5: np.array, T_j: np.array) -> np.array:
        """
        Calculate the saturation pressure of a component at a given temperature.

        Parameters:
        C1, C2, C3, C4, C5 : np.array
            Coefficients for the equation.
        T_j : np.array
            Temperature in Celsius.

        Returns:
        np.array
            Saturation pressure in atm.
        """
       # Convert temperature from Celsius to Kelvin
        T_j_K = T_j + 273.15

        # 𝑃𝑠 = exp[𝐶1 + 𝐶2𝑇 + 𝐶3ln𝑇 + 𝐶4𝑇**𝐶5]

        P = np.array([
            [
                np.exp(C1[i] + C2[i] / T_j_K[j] + C3[i] * np.log(T_j_K[j]) + C4[i] * T_j_K[j] ** C5[i])
                for j in range(len(T_j))
            ] 
            for i in range(len(self.components))
        ])
        P = P*1e-5
        return P
    
    def calculatate_Kij(self, Pij_sat: np.array)-> np.array:
        """
        Calculate the Ki values for the components.

        Parameters:
        Pi_sat : list
            List of saturation pressures for the components.

        Returns:
        list
            List of Ki values for the components.
        """
        return np.array([[Pij_sat[i][j]/self.pressure for i in range(len(Pij_sat))] for j in range(len(Pij_sat[0]))])
    
    def calculate_Lj(self, feed_tray: int , Vj: np.array)-> np.array:
        """
        Calculate the liquid flow rates for each component.

        Parameters:
        feed_tray : int
            NUmber of tray at which feed is given in the distillation column.
        Vj : list
            List of vapor flow rates.

        Returns:
        list
            List of liquid flow rates for each component.
        """
        Lj = np.zeros_like(Vj)

        # Lj Before feed tray
        for i in range(feed_tray-1):
            Lj[i] = Vj[i+1] - self.distillate_rate
            
        # Lj After feed tray
        for i in range(feed_tray-1, self.stages-1):
            Lj[i] = Vj[i+1] + self.feed_rate - self.distillate_rate
        
        for i in range(self.stages-1,self.stages):
            Lj[i] = Lj[i-1]

        return Lj
    
    def calculate_Aj(self, Lj: np.array)-> np.array:
        """
        Calculate the Aj values for each component.

        Parameters:
        Lj : list
            List of liquid flow rates for each component.

        Returns:
        list
            List of Aj values for each component.
        """
        # Aj = Lj-1
        Aj = np.zeros_like(Lj)
        Aj[0] = 0
        Aj[1:] = Lj[:-1]
        return Aj
    
    
    def calculate_Bij(self, Vj: np.array, Lj: np.array, Kij: np.array)-> np.array:
        """
        Calculate the Bij values for each component.

        Parameters:
        Lj : list
            List of liquid flow rates for each component.
        Kij : list
            List of Ki values for the components.

        Returns:
        list
            List of Bij values for each component.
        """
        # Bij[i][j] = Lj[i] + (Vj+1[i] +self.reflux) * Kij[i][j]
        
        return np.array([[self.distillate_rate - Lj[j] - (Vj[j]) * Kij[j][i] for j in range(len(Kij[0]))] for i in range(len(Kij))])
    

    def calculate_Cij(self, Vj: np.array, Kij: np.array)-> np.array:
        """
        Parameters:
        Vj : list
            List of vapor flow rates.
        Kij : list
            List of Ki values for the components.

        Returns:
        list
            List of Cj values for each component.
        """
        # Cij = - Vj * Kij
        # Ensure Vj is reshaped for broadcasting
        Vj_reshaped = Vj.reshape(1, -1)  # Shape (1, m)

        # Compute Cij using broadcasting
        Cij = Vj_reshaped * Kij  # (n, m) * (1, m) results in (n, m)
        Cij = np.transpose(Cij)
        # Shift left and set the last column to 0
        Cij_shifted = np.hstack((Cij[:, 1:], np.zeros((Cij.shape[0], 1))))
        return Cij_shifted
    
    def calculate_Dij(self,  feed_tray :int)-> np.array:
        """
        Calculate the Dj values for each component.

        Returns:
        list
            List of Dj values for each component.
        """
        # Dj = 0.0
        # Dij = 0
        Dij = np.zeros((self.num_components, self.num_components))
        Dij[feed_tray-1] = [-self.feed_rate* X for X in self.feed_composition]
        return np.transpose(Dij)

    # Thoams algorithm for bubble point method
    def thomas_algorithm(self, A: np.array, B: np.array, C: np.array, D: np.array)-> np.array:
        """
        Solve a tridiagonal system of equations using the Thomas algorithm.

        Parameters:
        A : list
            Lower diagonal elements.
        B : list
            Main diagonal elements.
        C : list
            Upper diagonal elements.
        D : list
            Right-hand side vector.

        Returns:
        list
            Solution vector.
        """
        n = len(B)
        c_prime = np.zeros(n)
        d_prime = np.zeros(n)
        x = np.zeros(n)
        
        c_prime[0] = C[0]/B[0]
        d_prime[0] = D[0]/B[0]
        
        for i in range(1, n):
            c_prime[i] = C[i]/(B[i] - A[i]*c_prime[i-1])
            d_prime[i] = (D[i] - A[i]*d_prime[i-1])/(B[i] - A[i]*c_prime[i-1])
        
        x[-1] = d_prime[-1]
        
        for i in range(n-2, -1, -1):
            x[i] = d_prime[i] - c_prime[i]*x[i+1]
        
        return x
    
    # Check summation of Xij = 1
    def check_summation(self, Xij: np.array) -> np.array:
        """
        Check if the summation of Xij is equal to 1.

        Parameters:
        Xij : list
            List of component mole fractions.

        Returns:
        bool
            True if the summation is equal to 1, False otherwise.
        """
        # NOrmalisation of Xij values to 1
        # Normalize each row to sum to 1
        Xij_nor = Xij / Xij.sum(axis=1, keepdims=True)  # Row-wise normalization

        return Xij_nor # Row-wise normalization
    
    def calculate_Tj(self, A: np.array, B: np.array, C: np.array, Xij: np.array) -> np.array:
        """
        Calculate the temperature profile (Tj) for each stage in the distillation column 
        using the Bubble Point method.

        Parameters:
        A, B, C : Antoine equation coefficients (arrays)
        Xij : Mole fractions of components in the liquid phase for each stage.

        Returns:
        T_j : np.array
            Temperature values for each stage.
        """
        T_j = np.zeros(self.stages)

        for stage in range(self.stages):
            T_guess = np.mean(self.boiling_points)  # Initial guess within boiling points range

            for _ in range(self.max_iterations):
                P_total = sum([10 ** (A[i] - (B[i] / (T_guess + C[i]))) * Xij[stage][i] 
                               for i in range(len(self.components))])

                if np.isclose(P_total, self.pressure, atol=self.convergence_tolerance):
                    break
                
                T_guess += (self.pressure - P_total) * 0.1  # Adjust the temperature guess

            T_j[stage] = T_guess

        return T_j

    def calculate_Vj_heat_balance(self, Lj: np.array, hL: np.array, hV: np.array, hF: float, Q: float, feed_tray: int) -> np.array:
        """
        Calculate vapor flow rates (Vj) for each stage using the Heat Balance Equation.

        Parameters:
        Lj : np.array
            Liquid flow rates per stage.
        hL : np.array
            Enthalpies of the liquid phase per stage.
        hV : np.array
            Enthalpies of the vapor phase per stage.
        hF : float
            Enthalpy of the feed.
        Q : float
            Heat added or removed from the system.
        feed_tray : int
            Index of the feed tray.

        Returns:
        Vj : np.array
            Vapor flow rates per stage.
        """
        Vj = np.zeros(self.stages)

        for j in range(self.stages):
            if j == 0:  # Condenser stage
                Vj[j] = (hL[j] * Lj[j] + Q) / hV[j]
            elif j == self.stages - 1:  # Reboiler stage
                Vj[j] = (hL[j-1] * Lj[j-1] + hF * self.feed_rate - hL[j] * Lj[j] + Q) / hV[j]
            else:  # Intermediate stages
                Vj[j] = (hL[j-1] * Lj[j-1] + hV[j+1] * Vj[j+1] + (hF * self.feed_rate if j == feed_tray else 0) - hL[j] * Lj[j] + Q) / hV[j]

        return Vj

    def calculate_enthalpies(self, Xij: np.array, Kij: np.array, hij: np.array) -> Tuple[np.array, np.array, float]:
        """
        Calculate the enthalpies of liquid (hL), vapor (hV), and feed (hF).

        Parameters:
        Xij : np.array
            Mole fractions of components in the liquid phase per stage.
        Kij : np.array
            Equilibrium constants per component per stage.
        hij : np.array
            Pure component enthalpies.

        Returns:
        Tuple containing:
        - hL : np.array (Liquid enthalpies per stage)
        - hV : np.array (Vapor enthalpies per stage)
        - hF : float (Feed enthalpy)
        """
        hL = np.zeros(self.stages)
        hV = np.zeros(self.stages)

        for j in range(self.stages):
            hL[j] = sum([Xij[j][i] * hij[i] for i in range(len(self.components))])
            hV[j] = sum([Kij[j][i] * Xij[j][i] * hij[i] for i in range(len(self.components))]) / \
                    sum([Kij[j][i] * Xij[j][i] for i in range(len(self.components))])

        hF = sum([self.feed_composition[i] * hij[i] for i in range(len(self.components))])

        return hL, hV, hF
    

    # Now comparing Values of the Vj and Tj for each stage and checking the convergence of the solution

    def mse_convergence(self, initial_T, initial_V, tol=1e-6, max_iters=1000):
        """
        Iterates the temperature and vapor flow update functions until convergence based on MSE.
        """
        T_old, V_old = np.array(initial_T, dtype=float), np.array(initial_V, dtype=float)
        
        for iteration in range(max_iters):
            T_new, V_new = process(T_old, V_old)  # Updated values from process function

            mse_T = np.mean((T_new - T_old) ** 2)
            mse_V = np.mean((V_new - V_old) ** 2)

            logging.debug("Iteration %d: MSE_T = %.10f, MSE_V = %.10f", iteration + 1, mse_T, mse_V)
            
            if mse_T < tol and mse_V < tol:
                logging.info("Converged in %d iterations", iteration + 1)
                return T_new, V_new
            
            T_old, V_old = T_new, V_new  # Update for next iteration

        logging.warning("Did not converge after %d iterations", max_iters)
        return T_old, V_old


def process(Tj_guess : np.array, Vj_guess : np.array):
    # Calculate saturation pressures
    Pij_sat = distillation.calculate_Pij(C1, C2, C3, C4, C5, Tj_guess)
    logging.debug("Saturation pressures: %s", Pij_sat)
    
    # Calculate equilibrium constants
    Kij = distillation.calculatate_Kij(Pij_sat)
    logging.debug("Equilibrium constants: %s", Kij)
    
    # Calculate liquid flow rates
    Lj = distillation.calculate_Lj(feed_tray, Vj_guess)
    logging.debug("Liquid flow rates: %s", Lj)
    
    # Calculate Aj, Bij, Cij, Dij
    Aj = distillation.calculate_Aj(Lj)
    Bij = distillation.calculate_Bij(Vj_guess,Lj, Kij)
    Cij = distillation.calculate_Cij(Vj_guess, Kij)
    Dij = distillation.calculate_Dij(feed_tray)
    logging.debug("Aj: %s", Aj)
    logging.debug("Bij: %s", Bij)
    logging.debug("Cij: %s", Cij)
    logging.debug("Dij: %s", Dij)
    
    # Solve for Xij using Thomas algorithm
    temp_values = []

    for i in range(3):
        temp = distillation.thomas_algorithm(Aj, Bij[i], Cij[i], Dij[i])
        temp_values.append(temp)
    Xij_= np.array(temp_values)
    logging.debug("Component mole fractions (Xij): %s", Xij_)
    
    Xij= distillation.check_summation(Xij_)
    logging.debug("Component mole fractions (Xij): %s", Xij)

    # Calculate new temperature profile
    Tj_new = distillation.calculate_Tj(A,B,C,Xij)
    logging.debug("New temperature profile: %s", Tj_new)
    
    # Calculate enthalpies
    hL, hV, hF = distillation.calculate_enthalpies(Xij, Kij, hij_l)
    logging.debug("Enthalpies of liquid: %s", hL)
    logging.debug("Enthalpies of vapor: %s", hV)
    logging.debug("Enthalpy of feed: %s", hF)
    
    # Calculate new vapor flow rates
    Vj_new = distillation.calculate_Vj_heat_balance(Lj, hL, hV, hF, Q, feed_tray)
    logging.debug("New vapor flow rates: %s", Vj_new)

    return Tj_new, Vj_new


import logging

# Configure logging
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

if __name__ == "__main__":

    # Antoine's coefficients for methanol, ethanol, and n-propanol
    A = [5.20409, 5.3229, 5.31384]
    B = [1581.341, 1670.409, 1690.864]
    C = [-33.5, -40.191, -51.804]

    # Coefficients for the Clasius claperyon equation for methanol, ethanol, and n-propanol
    C1 = [81.768, 74.475, 88.134]
    C2 = [-6876, -7164.3, -8438.6]
    C3 = [-8.7078, -7.327, -9.0766]
    C4 = [7.1926e-6, 3.1340e-6, 8.3303e-18]
    C5 = [2, 2, 6]

    # Enthalpies of liquid for methanol, ethanol, and n-propanol (kJ/kmol)
    hij_l = [-238.8, -277.6, -326.6]
    # Enthalpies of vapor for methanol, ethanol, and n-propanol (kJ/kmol)
    hij_v = [-35.3, -38.6, -47.6]


    # Feed tray number
    feed_tray = 2
    # Heat added or removed from the stage (kJ/h)
    Q = 0
    # Heat of feed (kJ/kmol)
    hF = 0
    # Number of stages
    stages = 3
    # Operating pressure (atm)
    P = 1
    # Feed composition
    feed_composition = [0.6, 0.2, 0.2]
    # Feed rate (kmol/h)
    feed_rate = 1000
    # Distillate rate (kmol/h)
    distillate_rate = 600
    # Reflux rate (kmol/h)
    reflux_rate = 2000
    # Convergence tolerance
    tolerance = 1e-3
    # Maximum number of iterations
    max_iterations = 10

    # Create an instance of the distillation column
    distillation = Distillation()

    # Initial guess for the temperature profile (Kelvin)
    Tj_guess = np.linspace(min(distillation.boiling_points) + 273.15, max(distillation.boiling_points) + 273.15, distillation.stages)

    # Initial guess for the vapor flow rates
    Vj_guess = np.array([2000, 2000, 2000])

    logging.debug("Initial temperature guess: %s", Tj_guess)
    logging.debug("Initial vapor flow rates: %s", Vj_guess)

    # Run the MSE convergence function to obtain the final values
    final_T, final_V = distillation.mse_convergence(Tj_guess, Vj_guess, max_iters=max_iterations)

    # Print the converged temperature profile and vapor flow rates
    logging.info("Converged Temperature Profile: %s", final_T)
    logging.info("Converged Vapor Flow Rates: %s", final_V)

    print("Final Converged Temperatures:", final_T)
    print("Final Converged Vapor Flow Rates:", final_V)