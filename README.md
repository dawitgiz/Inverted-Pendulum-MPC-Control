# Inverted Pendulum Control Using MPC and NMPC

## Overview

This repository contains MATLAB implementations of Model Predictive Control (MPC) and Nonlinear Model Predictive Control (NMPC) for an inverted pendulum system.

The project explores and compares three predictive control strategies for stabilizing an inverted pendulum:

1. Linear Model Predictive Control (Linear MPC)
2. Nonlinear Model Predictive Control (NMPC) without terminal cost
3. Nonlinear Model Predictive Control (NMPC) with LQR-based terminal cost

## System Description

The inverted pendulum system consists of a pendulum mounted on a horizontally moving cart. The control input is the horizontal force applied to the cart.

The system states are:

- Cart position (m)
- Cart velocity (m/s)
- Pendulum angle (rad)
- Pendulum angular velocity (rad/s)

The control objective is to stabilize the pendulum around its upright equilibrium while regulating cart position and respecting actuator constraints.

## Implemented Controllers

### 1. Linear MPC
Uses a linearized state-space model of the inverted pendulum for prediction and constrained optimization.

### 2. NMPC Without Terminal Cost
Uses the nonlinear system model for prediction, with numerical integration using the fourth-order Runge–Kutta (RK4) method.

### 3. NMPC With LQR-Based Terminal Cost
Extends the nonlinear MPC formulation by including an LQR-based terminal penalty to influence the predicted terminal state.

## MATLAB Files

| File | Description |
|---|---|
| `Linear_MPC_Comparison.m` | Linear MPC implementation |
| `NMPC_Without_Terminal_Cost.m` | Nonlinear MPC without terminal cost |
| `NMPC_With_Terminal_Cost.m` | Nonlinear MPC with LQR-based terminal cost |
| `Three_MPC_Comparison.m` | Comparative simulation of the three controllers |

## Software Requirements

- MATLAB
- Optimization Toolbox (`quadprog` and `fmincon`)
- Control System Toolbox (for state-space modeling and LQR calculations, where used)

## How to Run

1. Clone or download this repository.
2. Open MATLAB.
3. Add the downloaded project folder to the MATLAB path.
4. Open the desired `.m` file.
5. Run the script.

For the comparative study, run `Three_MPC_Comparison.m`.

## Performance Evaluation

The controllers can be evaluated using:

- Pendulum angular response
- Cart position response
- Settling time
- Integral of Absolute Error (IAE)
- Control energy
- Maximum control force
- Computational time

Simulation results and figures will be added to this repository as the project develops.

## Project Status

**Ongoing development.**

The repository currently provides the MATLAB controller implementations. Simulation results, additional tests, and documentation may be updated in future versions.

## Author

Dawit Gizachew Misikir

Electrical and Computer Engineering  
Control and Instrumentation Engineering

## License

No license has been assigned yet. All rights are reserved by the author unless a license is added to this repository.
