# ⚡ IEEE 14-Bus Power Flow Analysis | MATLAB Implementation

This project implements two classical load flow algorithms — **Gauss-Seidel (GS)** and **Newton-Raphson (NRF)** — to solve the steady-state power flow of the IEEE 14-bus test system, providing bus voltages, power injections, and network losses.

## 🎯 Objective
To compute:
- Bus voltage magnitudes and phase angles
- Active and reactive power flows
- Power mismatches and convergence analysis

## 🛠️ System Details
- **Buses:** 14 total (1 Slack, 2 PV, 11 PQ)
- **Transmission Lines:** 20 lines (with R, X, B data)
- **Transformers:** 3 transformers with off-nominal tap ratios (modeled in Ybus)
- **Shunt Compensation:** Bus 9 with 0.19 pu capacitor

## 🔍 Methodologies

### 1. Gauss-Seidel Method
- Formulated **Ybus** from line, transformer, and shunt data.
- Classified buses into Slack, PV, PQ.
- Used an acceleration factor (**α = 1.6**) to improve convergence.
- Iteratively updated voltages using power mismatch calculations.
- Reclassified PV buses as PQ when reactive limits were violated.
- Converged when voltage change < 1e-6 or after 100 iterations.

### 2. Newton-Raphson Method
- Modeled nonlinear power equations and formed **Jacobian matrix**.
- Solved for ΔP and ΔQ mismatches using matrix inversion.
- Updated voltage magnitudes and angles at each iteration.
- Handled PV bus Q-limit violations with reclassification.
- Achieved rapid convergence (~5–7 iterations) with high accuracy.

## 📈 Results
- Successfully computed the final voltages and power flows.
- Identified and managed PV bus reactive limits.
- Compared results between GS and NRF: NRF converged faster with fewer iterations.
- Calculated total power losses and validated solution accuracy.

## ✅ Tools & Implementation
- Developed in **MATLAB** using matrix operations and iterative solvers.
- Encapsulated Ybus formation, bus classification, and convergence logic in modular functions.
- Tested on the standard **IEEE 14-bus system**, matching reference results.

## 🔧 Future Improvements
- Extend to larger IEEE test cases (30-bus, 57-bus).
- Integrate optimal power flow (OPF) and contingency analysis.
- Simulate renewable energy sources and dynamic loads.

