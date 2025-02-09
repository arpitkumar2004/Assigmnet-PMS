# Q3: BUBBLE POINT METHOD – MULTICOMPONENT MULTISTAGE DISTILLATION
""" 
1000 kmol/h of a saturated-liquid mixture of 60% methanol (normal boiling point 65 C), 20 mol% ethanol 
(normal boiling point 98 C), and 20 mol% n-propanol (normal boiling point 97 C) is fed to the middle 
stage of a distillation column having three equilibrium stages, a total condenser, a partial reboiler, all 
operated at 1 atm. The distillate rate is 600 kmol/h, and the external reflux rate is 2,000 kmol/h of 
saturated liquid.
 """
# Write your own code to obtain converged solutions of the temperature profiles.



# You can use the following code as a starting point:
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
        self.reflux_rate = reflux_rate
        self.stages = stages
        self.pressure = pressure
        self.components = components
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
        # 𝑃𝑠 = exp[𝐶1 + 𝐶2𝑇 + 𝐶3ln𝑇 + 𝐶4𝑇**𝐶5]
        return np.array([[np.exp(C1[i] + C2[i]*T_j[j] + C3[i]*np.log(T_j[j]) + C4[i]*T_j[j]**C5[i]) for j in range(len(T_j))] for i in range(len(self.components))])
    
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
            Lj[i] = Vj[i+1] + self.distillate_rate
            
        # Lj After feed tray
        for i in range(feed_tray-1, self.stages-1):
            Lj[i] = Vj[i+1] + self.feed_rate + self.distillate_rate
            
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
    
    def calculate_Bij(self, Lj: np.array, Kij: np.array)-> np.array:
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
        return np.array([[Lj[i] + (Lj[i+1] + self.reflux_rate) * Kij[i][j] for j in range(len(Kij[0]))] for i in range(len(Kij))])
    

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
        return np.array([[-Vj[j] * Kij[j] for j in range(len(Vj))] for i in range(len(Kij))])
    
    def calculate_Dij(self,  feed_tray :int)-> np.array:
        """
        Calculate the Dj values for each component.

        Returns:
        list
            List of Dj values for each component.
        """
        # Dij = 0
        Dij = np.zeros((self.num_components, self.num_components))
        Dij[feed_tray-1] = self.feed_rate * self.feed_composition
        return Dij

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
    def check_summation(self, Xij: np.array) -> bool:
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
        Xij = Xij / np.sum(Xij)

        return np.isclose(sum(Xij), 1, atol=self.convergence_tolerance)
    
    # calculate Tj from summation of Xij and Kij values equals to 1
    def calculate_Tj(self, A:np.array, B: np.array, C: np.array, Xij:np.array ) -> np.array:
        """
        Calculate the temperature profile for the distillation column.

        Parameters:
        A : list
            List of A coeff. of antoine equation.
        B : list
            List of B coeff. of antoine equation.
        C : list
            List of C coeff. of antoine equation.
        Xij : list
            List of Xij values for the components.

        Returns:
        list
            Temperature values for stages of the distillation column.
        """
        # returning that value of T_j for which summation of  

        # np.e(A[i] + (B[i]/T_j[j]+C[i])Xij[i][j])  

        # for all i and j values should be equal to self.pressure, 
        # for this we need to assume the values of T_j from the range of boiling points of the components for iteration

        T_j = np.zeros(self.stages)
        for stage in range(self.stages):
            T_guess = np.mean(self.boiling_points)
            for _ in range(self.max_iterations):
                P_total = sum([10**(A[i] - (B[i] / (T_guess + C[i]))) * Xij[stage][i] for i in range(len(self.components))])
                if np.isclose(P_total, self.pressure, atol=self.convergence_tolerance):
                    break
                T_guess += (self.pressure - P_total) * 0.1  # Adjust the temperature guess
            T_j[stage] = T_guess
        return T_j
    


    # Calculation Vj values from Heat balance eqution for each stage
    # hLj-1*Lj-1 + hVj+1*Vj+1 + hFj*Fj - hLj*(Lj+ self.distillate) - hVj*(Vj + self.reflux) + Q = 0

    def calculate_Vj_heat_balance(self, Lj: np.array, hL: np.array, hV: np.array, hF: float, Q: float, feed_tray: int) -> np.array:
        """
        Calculate the vapor flow rates for each component using the heat balance equation.

        Parameters:
        Lj : np.array
            List of liquid flow rates for each component.
        hL : np.array
            List of enthalpies of liquid for each stage.
        hV : np.array
            List of enthalpies of vapor for each stage.
        hF : float
            Enthalpy of the feed.
        Q : float
            Heat added or removed from the stage.
        feed_tray : int
            The feed tray number.

        Returns:
        np.array
            Vapor flow rates for each component.
        """
        Vj = np.zeros_like(Lj)

        # Heat balance equation for each stage
        for j in range(self.stages):
            if j == 0:
                Vj[j] = (hL[j] * Lj[j] + Q) / hV[j]
            elif j == self.stages - 1:
                Vj[j] = (hL[j-1] * Lj[j-1] + hF * self.feed_rate - hL[j] * Lj[j] + Q) / hV[j]
            else:
                Vj[j] = (hL[j-1] * Lj[j-1] + hV[j+1] * Vj[j+1] + (hF * self.feed_rate if j == feed_tray else 0) - hL[j] * Lj[j] + Q) / hV[j]

        return Vj
    
    # Calcualtion of enthalpy of liquid and vapor and feed for each stage
    def calculate_enthalpies(self, Xij: np.array, Kij: np.array, hij: np.array) -> Tuple[np.array, np.array, float]:
        """
        Calculate the enthalpies of liquid, vapor, and feed.

        Parameters:
        Xij : np.array
            Mole fractions of components in the liquid phase for each stage.
        Kij : np.array
            Equilibrium constants for the components.
        hij : np.array
            Enthalpies of components.

        Returns:
        hL : np.array
            Enthalpies of liquid for each stage.
        hV : np.array
            Enthalpies of vapor for each stage.
        hF : float
            Enthalpy of the feed.
        """
        hL = np.zeros(self.stages)
        hV = np.zeros(self.stages)
        
        for j in range(self.stages):
            hL[j] = sum([Xij[j][i] * hij[i] for i in range(len(self.components))])
            hV[j] = sum([Kij[j][i] * Xij[j][i] * hij[i] for i in range(len(self.components))]) / sum([Kij[j][i] * Xij[j][i] for i in range(len(self.components))])
        
        hF = sum([self.feed_composition[i] * hij[i] for i in range(len(self.components))])
        
        return hL, hV, hF
    

    # Now comparing Values of the Vj and Tj for each stage and checking the convergence of the solution
    def check_convergence(self, Vj: np.array, Tj: np.array, Vj_prev: np.array, Tj_prev: np.array) -> bool:
        """
        Check if the solution has converged.

        Parameters:
        Vj : np.array
            Vapor flow rates for each component.
        Tj : np.array
            Temperature values for stages of the distillation column.
        Vj_prev : np.array
            Previous values of vapor flow rates.
        Tj_prev : np.array
            Previous values of temperature values.

        Returns:
        bool
            True if the solution has converged, False otherwise.
        """
        return np.all(np.isclose(Vj, Vj_prev, atol=self.convergence_tolerance)) and np.all(np.isclose(Tj, Tj_prev, atol=self.convergence_tolerance))
    

    

if __name__ == "__main__":


    # Antoine's coefficients for methanol, ethanol, and n-propanol
    A = [5.20409, 5.3229, 5.31384]
    B = [1581.341, 1670.409, 1690.864]
    C = [-33.5, -40.191, -51.804]

    # Coefficients for the Clasius claperyon equation for methanol, ethanol, and n-propanol
    C1 = [81.768, 74.475,88.134]
    C2 = [-6876, -7164.3, -8438.6]
    C3 = [-8.7078, -7.327,-9.0766]
    C4 = [7.1926e-6, 3.1340e-6, 8.3303e-18]
    C5 = [2,2,6]

    # Enthalpies of liquid for methanol, ethanol, and n-propanol (kJ/kmol)
    hij_l = [-238.8, -277.6, -326.6]
    # Enthalpies of vapor for methanol, ethanol, and n-propanol (kJ
    hij_v = [-108.9, -167.8, -217.8]
    # Latent heat of vaporization for methanol, ethanol, and n-propanol (kJ/kmol)
    delta_H_vap = [35.3, 38.56, 44.0]
    # Critical temperatures for methanol, ethanol, and n-propanol (K)
    Tc = [513, 514, 536]
    # Critical pressures for methanol, ethanol, and n-propanol (bar)
    Pc = [80.9, 61.4, 61.4]
    # Acentric factors for methanol, ethanol, and n-propanol
    omega = [0.559, 0.644, 0.666]
    # Feed temperature (Celsius)
    T_feed = 45
    # Feed tray number
    feed_tray = 2
    # Heat added or removed from the stage (kJ/h)
    Q = 0
    # Heat of feed (kJ/kmol)
    hF = 0
    # Feed flow rate (kmol/h)
    F = 0
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
    tolerance = 1e-6
    # Maximum number of iterations
    max_iterations = 1000



    # Create an instance of the distillation column
    distillation = Distillation()
    #  Initial guess for the temperature profile (Celsius)
    Tj_guess = np.linspace(min(distillation.boiling_points),max(distillation.boiling_points), distillation.stages)
    # Initial guess for the vapor flow rates
    Vj = np.array([0, 500, 500])

    # write complete code here to obtain converged solutions of the temperature profiles by using above defined functions and classes
    # You can use the following code as a starting point:
    # distillation.solve(Tj_guess, Vj, feed_tray, feed_composition, A, B, C, C1, C2, C3, C4, C5, hij, hij_v, delta_H_vap, Tc, Pc, omega, T_feed, Q, hF, F, stages, P, feed_rate, distillate_rate, reflux_rate, tolerance, max_iterations)
    # print(distillation.temperature_profile)
    # print(distillation.Vj)
    # print(distillation.Xij)
    # print(distillation.Kij)
