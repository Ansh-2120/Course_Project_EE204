# Course-Project-EE204 : Circuit Simulator

## Overview

A Python-based circuit Simulator tool designed to model and analyze electrical (RLC) circuits interactively. 
It supports both time-domain and frequency-domain analyses using numerical and symbolic computation.

### Key Features:
- Add and analyze RLC components (Resistors, Inductors, Capacitors) and sources (Voltage, Current).
- Solve circuit equations using Kirchhoff's Current Law (KCL).
- Perform Laplace transform-based frequency-domain analysis.
- Generate time-domain responses using inverse Laplace transforms.
- Visualize circuits as dynamic graphs.
- Interactive plots for time and frequency responses.

---

## Libraries Used

The following Python libraries are required to run the Circuit Simulator:

- **Streamlit**: For building the interactive web interface.
- **Sympy**: For symbolic mathematics and solving circuit equations.
- **NumPy**: For numerical computations and data handling.
- **Matplotlib**: For plotting and visualizing the circuit diagrams.
- **Plotly**: For interactive time-domain and frequency-domain plots.
- **NetworkX**: For graph-based visualization of circuit diagrams.

---

## Code Explanation

### `calculations.py`
This file contains the core class `CircuitSimulator`, which implements the computational logic. 
It includes methods for:

- **Component Management**:
  - `add_component()`: Adds RLC components and sources to the circuit.
- **Symbolic Equation Setup**:
  - `setup_node_variables()`: Creates symbolic variables for voltages and currents.
  - `build_equations()`: Constructs equations using KCL and component relations.
- **Solution and Analysis**:
  - `solve_circuit()`: Solves the equations for node voltages and branch currents.
  - `get_time_domain_response()`: Calculates time-domain responses via inverse Laplace transforms.
  - `get_frequency_response()`: Computes frequency response (magnitude and phase) over a given range.
- **Visualization**:
  - `create_circuit_visualization()`: Generates dynamic circuit diagrams.

### `circuit-simulator.py`
This file provides the interactive user interface built with Streamlit. 
Users can:
- Add circuit components via dropdowns and numeric inputs.
- Visualize the circuit dynamically.
- Analyze circuit responses (time and frequency domains) through interactive plots.
- Manage and clear components as needed.

---

## How to Run

### Prerequisites

1. Install Python 3.8+.
2. Install the required libraries using pip:
   ```bash
   pip install streamlit sympy numpy matplotlib plotly networkx
   ```

### Steps

1. Clone the repository:
   ```bash
   git clone <repository-url>
   cd <repository-folder>
   ```

2. Run the simulator using Streamlit:
   ```bash
   streamlit run circuit-simulator.py
   ```

3. Open the provided local URL (usually `http://localhost:8501`) in a web browser.

---

## Project Flow

1. **Add Components**:
   - Specify the type (Resistor, Capacitor, Inductor, Voltage Source, or Current Source).
   - Enter values and node connections.
   - For AC sources, specify additional parameters (amplitude, frequency, waveform type).
2. **Analyze Circuit**:
   - Solve the circuit equations to compute node voltages and branch currents.
   - View time-domain and frequency-domain plots.
3. **Visualize Circuit**:
   - Explore the dynamically generated circuit graph.
4. **Iterate**:
   - Modify components and analyze again.

---

## Example

### Input:
- Add components:
  - Resistor (100 Ω) between Node 0 and Node 1.
  - Voltage Source (10 V) between Node 1 and Node 0.

### Output:
- Voltage at Node 1: 10 V.
- Time-domain and frequency-domain responses plotted interactively.

---

## Future Improvements

- Support for additional component types (e.g., transformers, diodes).
- Enhanced visualization with 3D circuit layouts.
- Exportable reports for circuit analyses.

---

## Contributors

- **[Ansh Agarwal]** - 240102120
- **[Aditya Nitin Parate]** - 240102123
