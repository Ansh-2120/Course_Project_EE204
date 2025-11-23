# Advanced Circuit Simulator

## Description

This is a comprehensive circuit analysis tool that implements conventional circuit analysis methodologies including matrix-based topology analysis, Kirchhoff's laws, and Laplace transform techniques. Built with Python, it provides both educational insight into circuit theory and powerful computational capabilities for analyzing RLC circuits.

Unlike simple circuit simulators,it explicitly shows the mathematical framework behind circuit analysis - from incidence matrices to KCL/KVL equations to Laplace domain solutions.

---

## Features

- **Matrix-Based Analysis**: Automatic generation of Incidence, Cutset, and Tieset matrices for circuit topology
- **Equation Visualization**: Step-by-step display of time-domain and Laplace-domain equations
- **Kirchhoff's Laws**: Explicit KCL (node equations) and KVL (loop equations) formulation
- **Laplace Transform Analysis**: Symbolic impedance calculations and s-domain solutions
- **Time-Domain Response**: Inverse Laplace transform for transient analysis
- **Frequency-Domain Response**: Bode plots (magnitude and phase) for AC analysis
- **Component Support**: Resistors, Inductors, Capacitors, DC/AC Voltage Sources, Current Sources
- **Interactive Visualization**: Dynamic circuit diagrams and real-time plotting

---

## Techniques Used

### Circuit Analysis Methods
- **Topology Analysis**: Graph theory for building incidence, cutset, and tieset matrices
- **Nodal Analysis**: KCL-based node voltage method
- **Loop Analysis**: KVL-based mesh current method
- **Laplace Transform**: Converting differential equations to algebraic equations in s-domain
- **Inverse Laplace Transform**: Converting s-domain solutions back to time-domain

### Computational Libraries
- **SymPy**: Symbolic mathematics, Laplace transforms, equation solving
- **NumPy**: Numerical computations and matrix operations
- **NetworkX**: Graph algorithms for circuit topology (spanning trees, cycle detection)
- **Streamlit**: Interactive web-based user interface
- **Plotly**: Interactive plotting for time and frequency domain
- **Matplotlib**: Circuit diagram generation
- **Pandas**: Data display and matrix visualization

---

## Project Structure

```
├── circuit_simulator.py    # Main UI and interface logic (5-tab layout)
├── calculations.py         # Core circuit analysis engine
└── README.md              # This file
```

### Key Components

**`calculations.py`** - CircuitSimulator class with:
- Matrix generation methods (incidence, cutset, tieset)
- Equation builders (time-domain, Laplace, KCL, KVL)
- Circuit solver using SymPy's symbolic solver
- Response calculators (time and frequency domain)

**`circuit_simulator.py`** - Streamlit interface with:
- Component input controls
- 5 analysis tabs: Circuit Diagram, Matrix Analysis, Equations, Time Domain, Frequency Domain
- Interactive plotting and equation display

---

## How to Run

### Installation

1. **Install Python 3.8+**

2. **Install required packages:**
```bash
pip install streamlit sympy numpy matplotlib plotly networkx pandas
```

3. **(Optional) For better circuit layouts:**
```bash
# Linux
sudo apt-get install graphviz graphviz-dev
pip install pygraphviz

# macOS
brew install graphviz
pip install pygraphviz
```

### Running the Simulator

```bash
streamlit run circuit_simulator.py
```

Then open your browser to `http://localhost:8501`

---

## Usage

1. **Add Components**: Select type (R/L/C/V/I), enter value, specify nodes
2. **Solve Circuit**: Click "Solve Circuit" to generate matrices and equations
3. **View Results**: 
   - **Matrix Analysis**: See topology matrices (A, Q, B)
   - **Equations**: View KCL/KVL and Laplace-domain equations
   - **Time Domain**: Transient response plots
   - **Frequency Domain**: Bode magnitude and phase plots

---

## Technical Approach

1. **Build Circuit Graph**: Use NetworkX to represent circuit topology
2. **Generate Matrices**: Create incidence (A), cutset (Q), and tieset (B) matrices
3. **Formulate Equations**: Build KCL and KVL equations symbolically
4. **Transform to s-domain**: Apply Laplace transform to get impedances
5. **Solve System**: Use SymPy to solve linear system symbolically
6. **Calculate Responses**:
   - Time: Inverse Laplace transform → numerical evaluation
   - Frequency: Substitute s = jω → magnitude and phase

---

## Contributors

- **Ansh Agarwal** - 240102120
- **Parate Aditya Nitin** - 240102123

**Course**: EE204 - Circuit Theory  
**Institution**: IIT Guwahati
