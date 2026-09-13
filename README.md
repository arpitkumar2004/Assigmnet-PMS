# Process Modelling, Simulation & Computer-Aided Process Engineering (PMS & CAPE)

[![Python](https://img.shields.io/badge/Python-3.10%2B-blue.svg)](https://www.python.org/)
[![TensorFlow](https://img.shields.io/badge/TensorFlow-2.x-orange.svg)](https://tensorflow.org/)
[![Aspen HYSYS](https://img.shields.io/badge/Aspen-HYSYS%20%7C%20EDR-darkgreen.svg)](https://www.aspentech.com/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Institution](https://img.shields.io/badge/IIT%20Kharagpur-Chemical%20Engineering-red.svg)](http://www.iitkgp.ac.in/)

A comprehensive computational repository dedicated to **Chemical Process Modelling & Simulation (PMS)** and **Computer-Aided Process Engineering (CAPE)**. This repository houses mathematical modeling from first principles, iterative numerical algorithms, machine learning / deep learning for thermophysical property prediction, pinch analysis for energy integration, shortcut and rigorous distillation column design, and industrial heat exchanger rating using **Aspen HYSYS & EDR**.

---

## Table of Contents
1. [Portfolio Overview & CV Highlights](#portfolio-overview--cv-highlights)
2. [Important Project Reports (Ready for CV/Portfolio)](#important-project-reports-ready-for-cvportfolio)
   - [Project 1: Multicomponent Multistage Distillation Simulation (Wang-Henke BP Method)](#1-multicomponent-multistage-distillation-simulation-wang-henke-bp-method)
   - [Project 2: Machine Learning & Deep Learning for Thermophysical Property Estimation](#2-machine-learning--deep-learning-for-thermophysical-property-estimation)
   - [Project 3: Non-Linear Multicomponent Flash Separation & Recovery Optimization](#3-non-linear-multicomponent-flash-separation--recovery-optimization)
   - [Project 4: Debutanizer Column Design via Fenske-Underwood-Gilliland (FUG) Shortcut](#4-debutanizer-column-design-via-fenske-underwood-gilliland-fug-shortcut)
   - [Project 5: Pinch Analysis, Temperature Cascade & Heat Exchanger Network (HEN) Synthesis](#5-pinch-analysis-temperature-cascade--heat-exchanger-network-hen-synthesis)
   - [Project 6: Industrial Shell-and-Tube Heat Exchanger Rating & Retrofit (Aspen HYSYS & EDR)](#6-industrial-shell-and-tube-heat-exchanger-rating--retrofit-aspen-hysys--edr)
   - [Project 7: Multicomponent Reboiled Absorber Simulation (Simultaneous Correction Method)](#7-multicomponent-reboiled-absorber-simulation-simultaneous-correction-method)
3. [Resume / CV Ready Bullet Points](#resume--cv-ready-bullet-points)
4. [Mathematical Formulations & Algorithms](#mathematical-formulations--algorithms)
5. [Repository Structure](#repository-structure)
6. [Installation & Execution Guide](#installation--execution-guide)

---

## Portfolio Overview & CV Highlights

This repository contains academic and industrial-grade engineering projects completed as part of advanced chemical engineering coursework (**CH620003: Process Modelling and Simulation** and **CAPE Laboratory** at **IIT Kharagpur**).

### Core Competencies Demonstrated:
- **First-Principles Numerical Simulation**: Rigorous stage-by-stage MESH equation modeling, tridiagonal matrix solver (Thomas Algorithm), Newton-Raphson non-linear root finding, Simultaneous Correction (SC) methods.
- **Process Data Science & Scientific ML**: Closed-form Normal Equation regression, 3D loss surface optimization, deep neural network (MLP) architecture design in TensorFlow/Keras for thermodynamic property prediction with sparse training sets (10% data).
- **Process Energy Integration & Optimization**: Pinch Technology, composite curves ($T\text{-}H$ diagrams), Problem Table Algorithm (temperature cascade), minimum utility targeting, and grid diagram HEN synthesis saving **1.8 MW** of energy.
- **Commercial Process Simulators**: Industrial thermal design, rating, fouling analysis, and re-rating/retrofitting of shell-and-tube heat exchangers in **Aspen HYSYS** and **Aspen EDR**.

---

## Important Project Reports (Ready for CV/Portfolio)

Below is the structured breakdown of the primary engineering projects in this repository, including direct links to their detailed reports, mathematical models, code files, and key quantifiable outcomes.

---

### 1. Multicomponent Multistage Distillation Simulation (Wang-Henke BP Method)
* **Domain**: Rigorous Separation Processes & Numerical Modelling
* **Primary Report**: [22CH30008_PMS_A1_Report.pdf (Section 4)](file:///c:/Users/RDRL/Desktop/Assigmnet%20PMS/Assignment%201/22CH30008_PMS_A1_Report.pdf)
* **Hand Calculations & Verification**: [C_a hand iteration (1).pdf](file:///c:/Users/RDRL/Desktop/Assigmnet%20PMS/Assignment%201/Problem3/C_a%20hand%20iteration%20(1).pdf)
* **Source Code**: [C_b.py](file:///c:/Users/RDRL/Desktop/Assigmnet%20PMS/Assignment%201/Problem3/C_b.py), [C_c.py](file:///c:/Users/RDRL/Desktop/Assigmnet%20PMS/Assignment%201/Problem3/C_c.py), [test.py](file:///c:/Users/RDRL/Desktop/Assigmnet%20PMS/Assignment%201/test.py)

#### Problem Summary:
Simulation of a 3-stage (extendable to 10 stages) distillation column with a total condenser and partial reboiler at 1 atm separating a ternary alcohol feed: **60 mol% Methanol, 20 mol% Ethanol, 20 mol% n-Propanol** ($F = 1000\text{ kmol/h}$, $D = 600\text{ kmol/h}$, $L_0 = 2000\text{ kmol/h}$).

#### Key Technical Contributions:
- Formulated full **MESH** equations (Material, Equilibrium, Summation, Heat balances) for multicomponent vapor-liquid equilibrium (VLE).
- Implemented an object-oriented `Distillation` solver in Python utilizing the **Thomas Algorithm (TDMA)** to invert tridiagonal component matrices in $\mathcal{O}(N)$ computational complexity.
- Utilized the extended Clausius-Clapeyron equation ($P^s = \exp(C_1 + C_2/T + C_3\ln T + C_4 T^{C_5})$) for temperature-dependent vapor pressures and stage bubble-point updates.
- Extended the column model to 10 equilibrium stages and conducted feed tray location sensitivity analysis to minimize thermal fluctuations and optimize separation purity.
- Implemented rigorous unit test suites with `unittest` in [test.py](file:///c:/Users/RDRL/Desktop/Assigmnet%20PMS/Assignment%201/test.py) covering matrix tridiagonal solvers, enthalpy balances, and convergence routines.

---

### 2. Machine Learning & Deep Learning for Thermophysical Property Estimation
* **Domain**: Scientific Machine Learning & Molecular Thermodynamics
* **Primary Report**: [22CH30008_PMS_A1_Report.pdf (Section 2)](file:///c:/Users/RDRL/Desktop/Assigmnet%20PMS/Assignment%201/22CH30008_PMS_A1_Report.pdf)
* **Jupyter Notebook**: [Sol_A1.ipynb](file:///c:/Users/RDRL/Desktop/Assigmnet%20PMS/Assignment%201/Problem1/Sol_A1.ipynb)
* **Visualizations**: [fig1.png](file:///c:/Users/RDRL/Desktop/Assigmnet%20PMS/Assignment%201/Problem1/fig1.png), [fig2.png](file:///c:/Users/RDRL/Desktop/Assigmnet%20PMS/Assignment%201/Problem1/fig2.png), [fig4.png](file:///c:/Users/RDRL/Desktop/Assigmnet%20PMS/Assignment%201/Problem1/fig4.png)

#### Problem Summary:
Estimating critical thermophysical properties (normal boiling point $T_b$ and reduced boiling point $T_b/T_c$) of over 600 organic chemical compounds from molecular weight ($M_w$) and Pitzer acentric factor ($\omega$).

#### Key Technical Contributions:
- **Baseline Parametric Regression**: Derived simple and multivariable linear regressions. Formulated closed-form Ordinary Least Squares via the **Normal Equation**:
  $$\boldsymbol{\theta} = (\mathbf{X}^T \mathbf{X})^{-1} \mathbf{X}^T \mathbf{y}$$
  Achieved $R^2 = 0.7861$ and $MSE = 0.0004$ on reduced boiling point correlation ($\theta_0 = 0.5958$, $\theta_1 = 0.0002$, $\theta_2 = 0.1546$).
- **Cost Function Landscape Visualization**: Rendered 3D meshgrid contour plots of the Mean Squared Error loss surface against parameter vectors to evaluate convexity and gradient behavior.
- **Deep Neural Network (MLP)**: Designed a deep feedforward regression network in **TensorFlow/Keras** trained under extreme data scarcity (10% training, 90% test set).
- **Hyperparameter Optimization & TensorBoard**: Executed grid search across network depths (1–4 hidden layers) and widths (16–64 neurons/layer) with ReLU activations and Adam optimization.
- **Results**: Best architecture (3 hidden layers, 16 neurons/layer) attained a test **MSE of $0.00057$** and **MAE of $0.0166$** ($< 1.7\%$ error), demonstrating superior non-linear property representation over traditional empirical correlations.

---

### 3. Non-Linear Multicomponent Flash Separation & Recovery Optimization
* **Domain**: Non-Linear Thermodynamics & Equilibrium Optimization
* **Primary Report**: [22CH30008_PMS_A1_Report.pdf (Section 3)](file:///c:/Users/RDRL/Desktop/Assigmnet%20PMS/Assignment%201/22CH30008_PMS_A1_Report.pdf)
* **Source Code**: [B_a_b.py](file:///c:/Users/RDRL/Desktop/Assigmnet%20PMS/Assignment%201/Problem2/B_a_b.py), [B_c.py](file:///c:/Users/RDRL/Desktop/Assigmnet%20PMS/Assignment%201/Problem2/B_c.py)

#### Problem Summary:
Separation of a ternary aromatic mixture (100 kmol/h feed: **60 mol% Benzene, 25 mol% Toluene, 15 mol% o-Xylene**) in an isothermal flash drum operating at 1.01325 bar across temperatures from 100°C to 150°C.

#### Key Technical Contributions:
- Formulated the non-linear **Rachford-Rice objective function**:
  $$f(\psi) = \sum_{i=1}^C \frac{z_i (K_i - 1)}{1 + \psi (K_i - 1)} = 0$$
  and derived its exact analytical first derivative for second-order quadratic convergence via **Newton-Raphson iteration**.
- Dynamically integrated Antoine vapor pressure equations to update equilibrium ratios $K_i(T)$.
- Base Case ($100^\circ\text{C}$): Converged in 4 iterations to vapor fraction $\psi = 0.6840$ ($V = 68.4\text{ kmol/h}$, $L = 31.6\text{ kmol/h}$) with vapor phase composition $y = [0.6963, 0.2245, 0.0793]$.
- **Process Optimization**: Conducted fine-grid temperature sweeps ($0.0001^\circ\text{C}$ resolution) to maximize Benzene recovery in distillate, identifying the global optimum at **$T = 106.99^\circ\text{C}$** yielding a peak recovery of **$92.47\%$**.

---

### 4. Debutanizer Column Design via Fenske-Underwood-Gilliland (FUG) Shortcut
* **Domain**: Hydrocarbon Fractionation & Shortcut Column Design
* **Primary Report**: [22CH30008.pdf (Part 1)](file:///c:/Users/RDRL/Desktop/Assigmnet%20PMS/Assignment%202/22CH30008.pdf)
* **Source Code**: [1.py](file:///c:/Users/RDRL/Desktop/Assigmnet%20PMS/Assignment%202/Problem%201/1.py), [test.ipynb](file:///c:/Users/RDRL/Desktop/Assigmnet%20PMS/Assignment%202/test.ipynb)

#### Problem Summary:
Design and specification of an 8-component petroleum debutanizer column ($i\text{C}_4$, $n\text{C}_4$, $i\text{C}_5$, $n\text{C}_5$, $\text{C}_6$, $\text{C}_7$, $\text{C}_8$, $\text{C}_9$) operating at 80 psia (552 kPa) to fractionate Light Key ($n\text{C}_4$) and Heavy Key ($i\text{C}_5$).

#### Key Technical Contributions:
- **Fenske Equation**: Estimated minimum theoretical stages ($N_{min} = 14.93$) and calculated product splits for all 6 non-key components across top and bottoms streams.
- **Underwood Equation**: Numerically solved the non-linear root condition using `scipy.optimize.fsolve` for the active root $\theta$ between relative volatilities to evaluate the minimum reflux ratio ($R_{min} = 91.77$).
- **Gilliland Correlation**: Automated calculation of actual theoretical stages required at an operating reflux ratio of $1.3 R_{min}$, determining $N = 27.26$ equilibrium stages.
- Validated hydrocarbon K-values under elevated pressure and evaluated operational feasibility.

---

### 5. Pinch Analysis, Temperature Cascade & Heat Exchanger Network (HEN) Synthesis
* **Domain**: Process Integration, Energy Efficiency & Thermal Pinch Analysis
* **Primary Reports**: 
  - [Pinch Analysis and Heat Exchanger Network Integration.docx](file:///c:/Users/RDRL/Desktop/Assigmnet%20PMS/Assignment%202/Pinch%20Analysis%20and%20Heat%20Exchanger%20Network%20Integration.docx)
  - [22CH30008.pdf (Part 2)](file:///c:/Users/RDRL/Desktop/Assigmnet%20PMS/Assignment%202/22CH30008.pdf)
  - [PMS A2-1.pdf](file:///c:/Users/RDRL/Desktop/Assigmnet%20PMS/PMS%20A2-1.pdf), [A2 Q3 PMS.pdf](file:///c:/Users/RDRL/Desktop/Assigmnet%20PMS/Assignment%202/Problem%203/A2%20Q3%20PMS.pdf)
* **Visuals & Hand Synthesis**: [3 graph.jpeg](file:///c:/Users/RDRL/Desktop/Assigmnet%20PMS/Assignment%202/Problem%203/3%20graph.jpeg)

#### Problem Summary:
Energy integration of a multi-stream chemical process network consisting of 3 hot streams and 3 cold streams with a minimum temperature approach $\Delta T_{min} = 10^\circ\text{C}$.

#### Key Technical Contributions:
- Constructed **Composite Curves** (Hot and Cold curves on a $T\text{-}H$ diagram) to visualize maximum heat exchange potential.
- Formulated the **Problem Table Algorithm (Temperature Cascade Method)** using shifted interval temperatures ($T^* = T - \Delta T_{min}/2$ for hot streams, $T^* = T + \Delta T_{min}/2$ for cold streams) to identify enthalpy surpluses and deficits.
- **Pinch Point Determination**: Identified the exact pinch point at cold stream $90^\circ\text{C}$ / hot stream $100^\circ\text{C}$.
- **Utility Targets**: Evaluated minimum hot utility demand ($Q_{H,min} = 0.1\text{ MW}$), minimum cold utility demand ($Q_{C,min} = 0.05\text{ MW}$), and recovered **$1.8\text{ MW}$ of heat internally**.
- **HEN Synthesis & Flowsheet Integration**: Designed an optimal heat exchanger network respecting pinch heuristics (no heat transfer across pinch, no cold utilities above pinch, no hot utilities below pinch) and synthesized a full **Process Flow Diagram (PFD)** with labeled heat duties and split streams.

---

### 6. Industrial Shell-and-Tube Heat Exchanger Rating & Retrofit (Aspen HYSYS & EDR)
* **Domain**: Computer-Aided Process Engineering (CAPE) & Industrial Thermal Design
* **Primary Report**: [Assignment 2 CAPE LAb report.pdf](file:///c:/Users/RDRL/Desktop/Assigmnet%20PMS/CAPE%20Lab/Assignment%202%20CAPE%20LAb%20report.pdf)
* **Simulation Files**: [cape lab.hsc](file:///c:/Users/RDRL/Desktop/Assigmnet%20PMS/CAPE%20Lab/cape%20lab.hsc), [cape lab.bk0](file:///c:/Users/RDRL/Desktop/Assigmnet%20PMS/CAPE%20Lab/cape%20lab.bk0), [A2Q1.hsc](file:///c:/Users/RDRL/Desktop/Assigmnet%20PMS/CAPE%20Lab/A2Q1.hsc), [A2Q1.bk0](file:///c:/Users/RDRL/Desktop/Assigmnet%20PMS/CAPE%20Lab/A2Q1.bk0)

#### Problem Summary:
Thermal design, rating, and retrofitting of an industrial shell-and-tube heat exchanger in **Aspen HYSYS** and **Aspen EDR (Exchanger Design & Rating)**.

#### Key Technical Contributions:
- **Base Case Design**: Simulated a multi-pass (2-pass tube, 1-pass shell) exchanger heating 80,000 lb/hr Benzene (70°F to 140°F, 45 psia) using Toluene (235°F to 150°F, 40 psia) with allowable $\Delta P \le 5\text{ psia}$ and fouling resistance $R_f = 0.0015\text{ ft}^2\cdot\text{hr}\cdot^\circ\text{F/BTU}$. Calculated total heat duty of **$2,124,640\text{ BTU/hr}$ ($2.125\text{ MMBTU/hr}$)**.
- **Detailed EDR Rating**: Selected TEMA shell type, baffle cuts, tube bundle pitch, diameter, and metallurgy to satisfy pressure drop constraints and heat transfer coefficients.
- **Retrofit / Re-Rating Scenario**: Evaluated the existing physical unit for a severe new process duty: heating Methanol using high-temperature pressurized water at **$10.5\text{ MMBTU/hr}$** (a ~5x duty increase). Retuned operating pressures, stream allocations, and checked overdesign/vibration margins.

---

### 7. Multicomponent Reboiled Absorber Simulation (Simultaneous Correction Method)
* **Domain**: Complex Separation Processes & Simultaneous Non-Linear Solvers
* **Primary Code**: [2.py](file:///c:/Users/RDRL/Desktop/Assigmnet%20PMS/Assignment%202/2.py)
* **Problem Definition**: [A2.pdf (Question 2)](file:///c:/Users/RDRL/Desktop/Assigmnet%20PMS/Assignment%202/A2.pdf)

#### Problem Summary:
Rigorous simulation of a 13-stage deethanizing reboiled absorber column operating at 400 psia (2.76 MPa) with n-decane absorbent oil separating light hydrocarbons ($\text{C}_1, \text{C}_2, \text{C}_3, n\text{C}_4, n\text{C}_5$).

#### Key Technical Contributions:
- Implemented the **Simultaneous Correction (SC)** algorithm to simultaneously solve coupled stage material balances, phase equilibria, and non-linear enthalpy balances.
- Interfaced **CoolProp** thermodynamic library (`PropsSI`) for cryogenic and high-pressure liquid/vapor hydrocarbon enthalpies.
- Iteratively converged stage temperatures, interstage flow rates, and reboiler duty under high pressure.

---

## Mathematical Formulations & Algorithms

### 1. Rachford-Rice Equation (Flash VLE)
The objective function for vapor fraction $\psi = V/F$:
$$f(\psi) = \sum_{i=1}^C \frac{z_i (K_i - 1)}{1 + \psi(K_i - 1)} = 0$$
Analytical derivative for Newton-Raphson update:
$$f'(\psi) = - \sum_{i=1}^C \frac{z_i (K_i - 1)^2}{[1 + \psi(K_i - 1)]^2}$$
Newton-Raphson iteration step:
$$\psi^{(k+1)} = \psi^{(k)} - \frac{f(\psi^{(k)})}{f'(\psi^{(k)})}$$

### 2. Thomas Algorithm (Tridiagonal Matrix Algorithm - TDMA)
For stage-by-stage component mass balance tridiagonal systems:
$$A_j x_{j-1} + B_j x_j + C_j x_{j+1} = D_j$$
Forward elimination:
$$C'_j = \frac{C_j}{B_j - A_j C'_{j-1}}, \quad D'_j = \frac{D_j - A_j D'_{j-1}}{B_j - A_j C'_{j-1}}$$
Back substitution:
$$x_j = D'_j - C'_j x_{j+1}$$

### 3. Normal Equation (Linear Regression)
$$\boldsymbol{\theta} = (\mathbf{X}^T \mathbf{X})^{-1} \mathbf{X}^T \mathbf{y}$$
$$\text{Cost: } J(\boldsymbol{\theta}) = \frac{1}{2m} \sum_{i=1}^m (\mathbf{x}^{(i)\boldsymbol{\theta}} - y^{(i)})^2$$

### 4. Fenske-Underwood-Gilliland (FUG) Equations
- **Fenske Equation** (Minimum theoretical stages):
  $$N_{min} = \frac{\ln\left[\left(\frac{x_{LK,D}}{x_{HK,D}}\right)\left(\frac{x_{HK,B}}{x_{LK,B}}\right)\right]}{\ln \alpha_{LK/HK}}$$
- **Underwood Equation** (Minimum reflux ratio):
  $$\sum_{i=1}^C \frac{\alpha_i z_{i,F}}{\alpha_i - \theta} = 1 - q, \quad R_{min} = \sum_{i=1}^C \frac{\alpha_i x_{i,D}}{\alpha_i - \theta} - 1$$

---

## Repository Structure

```text
Assigmnet PMS/
├── README.md                                          # Comprehensive project documentation & CV guide
├── PMS A2-1.pdf                                       # Assignment 2 Question 3 handwritten solution & HEN
│
├── Assignment 1/                                      # ASSIGNMENT 1: ML, Flash, & Bubble Point Distillation
│   ├── A1.pdf                                         # Original assignment problem sheet
│   ├── 22CH30008_PMS_A1_Report.pdf                    # Complete compiled engineering report (12 pages)
│   ├── data_file.xlsx                                 # Thermophysical properties dataset (600+ compounds)
│   ├── test.py                                        # Unit testing suite for distillation module
│   ├── Sol_A!_q3.ipynb                                # Distillation prototyping notebook
│   ├── Problem1/                                      # Q1: Regression & Neural Networks
│   │   ├── Sol_A1.ipynb                               # End-to-end ML notebook (EDA, OLS, ANN, TensorBoard)
│   │   └── fig1.png, fig2.png, fig4.png               # Regression fits, 3D loss surface & residual plots
│   ├── Problem2/                                      # Q2: Multicomponent Flash Problem
│   │   ├── B_a_b.py                                   # Newton-Raphson Rachford-Rice solver & temperature sweep
│   │   └── B_c.py                                     # Benzene recovery optimization algorithm
│   └── Problem3/                                      # Q3: Multicomponent Distillation (Wang-Henke BP)
│       ├── C_a hand iteration (1).pdf                 # Step-by-step verified hand calculation
│       ├── C_b.py                                     # OOP Distillation class with Thomas Algorithm
│       └── C_c.py                                     # 10-stage distillation & feed location optimizer
│
├── Assignment 2/                                      # ASSIGNMENT 2: FUG, SC Method, & Pinch Analysis
│   ├── A2.pdf                                         # Original assignment problem sheet
│   ├── 22CH30008.pdf                                  # Debutanizer & Pinch Analysis engineering report
│   ├── Pinch Analysis and Heat Exchanger Network...docx # Comprehensive HEN design & cascade report
│   ├── 2.py                                           # Simultaneous Correction (SC) Reboiled Absorber
│   ├── test.ipynb                                     # FUG shortcut distillation calculations
│   ├── Problem 1/                                     # Q1: FUG Shortcut Distillation
│   │   └── 1.py                                       # Fenske, Underwood, & Gilliland Python solver
│   └── Problem 3/                                     # Q3: Pinch Analysis & Flowsheet Synthesis
│       ├── A2 Q3 PMS.pdf                              # Detailed handwritten cascade & HEN calculations
│       └── 3 graph.jpeg                               # Composite curves & heat exchanger grid diagram
│
└── CAPE Lab/                                          # COMPUTER-AIDED PROCESS ENGINEERING (ASPEN)
    ├── Assignment 2 CAPE LAb report.pdf               # Complete Aspen HYSYS/EDR rating report (10 pages)
    ├── cape lab.hsc, cape lab.bk0                     # Aspen HYSYS simulation case & backup files
    ├── A2Q1.hsc, A2Q1.bk0                             # Heat exchanger simulation files
    ├── Assignment 2/                                  # HYSYS output stream and block property reports
    │   ├── q1p1.txt, q1p2.txt, q1p3.txt, q2p1.txt
    └── arpit/                                         # Aspen Plus simulation models (.apwz)
```

---

## Installation & Execution Guide

### Prerequisites
- Python 3.8 to 3.10
- Jupyter Notebook / JupyterLab
- Commercial license for **Aspen HYSYS & EDR** (for `.hsc` and `.bk0` files)

### 1. Clone Repository & Setup Environment
```bash
git clone https://github.com/arpitkumar2004/Assigmnet-PMS.git
cd Assigmnet-PMS

# Create virtual environment
python -m venv venv
# Activate on Windows:
.\venv\Scripts\activate
# Activate on Linux/macOS:
source venv/bin/activate

# Install required packages
pip install numpy scipy pandas matplotlib seaborn tensorflow scikit-learn openpyxl coolprop
```

### 2. Running Simulations

#### Run Multicomponent Distillation Unit Tests:
```bash
cd "Assignment 1"
python test.py
```

#### Run Multicomponent Flash & Recovery Optimization:
```bash
cd "Assignment 1/Problem2"
python B_c.py
```

#### Run Debutanizer FUG Shortcut Calculations:
```bash
cd "Assignment 2/Problem 1"
python 1.py
```

#### Run Reboiled Absorber Simultaneous Correction Solver:
```bash
cd "Assignment 2"
python 2.py
```

#### Launch Jupyter Notebooks:
```bash
jupyter notebook
# Open Assignment 1/Problem1/Sol_A1.ipynb for ML/ANN workflows
# Open Assignment 2/test.ipynb for FUG distillation workflows
```

---

## Author
- **Arpit Kumar**
- Department of Chemical Engineering, Indian Institute of Technology (IIT) Kharagpur
- GitHub: [@arpitkumar2004](https://github.com/arpitkumar2004)
