# Numerical Simulation of Stable Process and Brownian Motion

## Project Description
This MATLAB program implements a numerical simulation system to study the interaction between α-stable processes and Brownian motion, utilizing the CIR (Cox-Ingersoll-Ross) model.

## Main Features
- Numerical simulation of α-stable processes
- Brownian motion path generation
- Time transformation calculation
- Convergence rate analysis

## Key Parameters
- `TN`: Total simulation time (10)
- `T1`: Usage time period (1)
- `t0`: Initial time (0.1)
- `j`: Step size control parameter (13)
- `alpha_list`: Range of α values (0.3 to 0.9, step 0.1)
- `num_iterations`: Number of simulation paths (1000)
- CIR model parameters:
  - κ (kappa) = 8
  - σ (sigma) = 0.5
  - θ (theta) = 0.125

## Algorithm Flow
1. Initialize system parameters
2. For each α value:
   - Generate stable process D
   - Simulate Brownian motion
   - Calculate time transformation E
   - Solve CIR equation
3. Calculate and record convergence rates

## Usage
1. Ensure MATLAB environment is properly configured
2. Run the program: 
```
matlab
run TCCIR.m
```


## Output Description
- Program displays calculation progress during execution
- Final output includes convergence rate matrix containing:
  - First row: α values
  - Second row: corresponding convergence rates

## Dependencies
- MATLAB (R2018b or higher recommended)
- MATLAB Statistics and Machine Learning Toolbox

