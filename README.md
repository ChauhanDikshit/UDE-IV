# UDE-IV: Enhanced Unified Differential Evolution for Constrained Optimization

## Overview

UDE-IV is an advanced differential evolution algorithm designed for solving constrained real-parameter optimization problems (COPs). It builds upon the success of its predecessors, UDE-II and UDE-III, which achieved top rankings in the CEC 2018 and CEC 2024 competitions. UDE-IV introduces a dual-population framework, strategy adaptation, and enhanced exploration mechanisms for improved performance on complex constrained problems.

## Key Features

- **Four Trial Vector Generation Strategies**:
  - `DE/rand/1`
  - `DE/current-to-rand/1`
  - `DE/current-to-pbest/1`
  - `DE/rand-to-pbest/1`

- **Dual Population Mechanism**:
  - Top sub-population uses fixed strategies.
  - Bottom sub-population applies adaptive strategy selection.

- **Strategy Adaptation**:
  - Probabilistic selection based on historical success.

- **Improved Exploration**:
  - Enhanced use of `pbest`-guided strategies.

- **Parameter Adaptation**:
  - LSHADE44-based mechanism.

- **Constraint Handling**:
  - Combines feasibility rule and ε-constraint method (from C²oDE).

- **Stagnation Avoidance**:
  - Improved with optimal parameter tuning.

- **Simplified Parent Selection**:
  - Ranking-based selection removed for better efficiency.

## Benchmarking

UDE-IV has been rigorously tested on the CEC 2025 benchmark suite consisting of 28 constrained 30D problems and has shown superior performance compared to UDE-III.

## Getting Started

### Prerequisites

- MATLAB

# Run the main optimization script
# Example (MATLAB):
main_ude.m
