import base64
from matplotlib import pyplot as plt
import streamlit as st
import sympy as sp
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from sympy.abc import t, s
from sympy.integrals import inverse_laplace_transform
import networkx as nx
from sympy import *
import io
from base64 import b64encode

class CircuitSimulator:
    def __init__(self):
        self.components = []
        self.nodes = set()
        self.equations = []
        self.node_vars = {}
        self.solutions = {}
        self.s = sp.Symbol('s')
        self.t = sp.Symbol('t')
        self.component_count = {'R': 0, 'L': 0, 'C': 0, 'V': 0, 'I': 0}
        
    def add_component(self, name, value, n1, n2, ac_params=None):
        """Add a component to the circuit with optional AC parameters"""
        if n1 == n2:
            raise ValueError("Start and end nodes cannot be the same")
            
        # Validate component values
        if float(value) <= 0:
            raise ValueError(f"Component value must be positive, got {value}")
            
        self.component_count[name] += 1
        component_id = f"{name}{self.component_count[name]}"
        component = {
            'id': component_id,
            'name': name,
            'value': float(value),
            'n1': int(n1),
            'n2': int(n2),
            'ac_params': ac_params
        }
        self.components.append(component)
        self.nodes.update([n1, n2])
        return True

    def setup_node_variables(self):
        """Setup variables for nodes and currents"""
        self.node_vars = {}
        for node in self.nodes:
            self.node_vars[node] = sp.Symbol(f'V{node}')
        
        self.current_vars = {}
        for comp in self.components:
            if comp['name'] in ['V', 'L']:
                self.current_vars[comp['id']] = sp.Symbol(f'I_{comp["id"]}')

    def build_equations(self):
        """Build circuit equations with support for AC sources"""
        self.equations = []
        
        # KCL equations for each node
        for node in self.nodes:
            if node == min(self.nodes):  # Skip reference node
                continue

            current_sum = 0

            for comp in self.components:
                if comp['n1'] == node or comp['n2'] == node:
                    v1 = self.node_vars[comp['n1']]
                    v2 = self.node_vars[comp['n2']]

                    # Define current direction: positive if entering node
                    if comp['name'] == 'R':
                        i = (v1 - v2) / comp['value']
                    elif comp['name'] == 'C':
                        i = self.s * comp['value'] * (v1 - v2)
                    elif comp['name'] == 'L':
                        i = (v1 - v2) / (comp['value'] * self.s)
                    elif comp['name'] == 'V':
                        i = self.current_vars[comp['id']]
                        if comp.get('ac_params'):
                            ac = comp['ac_params']
                            omega = 2 * sp.pi * ac['frequency']
                            v_source = ac['magnitude'] / (self.s + omega * 1j) if ac['type'] == 'sin' else \
                                       ac['magnitude'] / (self.s - omega * 1j)
                            self.equations.append(v1 - v2 - v_source)
                        else:
                            self.equations.append(v1 - v2 - comp['value'])
                    elif comp['name'] == 'I':
                        i = comp['value']

                    # Add current contribution based on direction
                    if comp['n1'] == node:
                        current_sum += i
                    else:
                        current_sum -= i

            if current_sum != 0:
                self.equations.append(current_sum)

        # Add component-specific equations
        for comp in self.components:
            if comp['name'] == 'L':
                v1 = self.node_vars[comp['n1']]
                v2 = self.node_vars[comp['n2']]
                i = self.current_vars[comp['id']]
                self.equations.append(v1 - v2 - comp['value'] * self.s * i)

        # Add reference node equation
        if self.nodes:
            ref_node = min(self.nodes)
            self.equations.append(self.node_vars[ref_node])

    def solve_circuit(self):
        """Solve the circuit equations"""
        try:
            if not self.components:
                raise ValueError("No components in circuit")
                
            self.setup_node_variables()
            self.build_equations()

            variables = list(self.node_vars.values()) + list(self.current_vars.values())

            if not variables:
                raise ValueError("No variables to solve for")

            sympy_equations = [sp.Eq(eq, 0) for eq in self.equations]
            solution = sp.solve(sympy_equations, variables, dict=True)

            if not solution:
                raise ValueError("No solution found")

            self.solutions = solution[0]
            return True
        except Exception as e:
            print(f"Error solving circuit: {str(e)}")
            return False

    def get_frequency_response(self, fmin=1, fmax=1000, points=100):
        """Calculate frequency response"""
        f = np.logspace(np.log10(fmin), np.log10(fmax), points)
        w = 2 * np.pi * f
        responses = {}

        for var, expr in self.solutions.items():
            magnitude = []
            phase = []
            for w_val in w:
                try:
                    v_complex = complex(expr.subs(self.s, 1j * w_val))
                    
                    # Convert to native float values to make them serializable
                    magnitude.append(float(abs(v_complex)))
                    phase.append(float(np.angle(v_complex, deg=True)))
                except Exception as e:
                    magnitude.append(0.0)  # Default to zero if there's an error
                    phase.append(0.0)

            responses[str(var)] = {'magnitude': magnitude, 'phase': phase}

        return f, responses


    from sympy.utilities.lambdify import lambdify

    def get_time_domain_response(self, tmax=0.01, points=1000):
        t_vals = np.linspace(0, tmax, points)
        responses = {}
        dt = t_vals[1] - t_vals[0]
        MAX_AMPLITUDE = 100  # Maximum amplitude limit

        if 'V0' in self.solutions:
            responses['V0'] = np.zeros(len(t_vals))

        for var, expr in self.solutions.items():
            if var == 'V0':
                continue

            try:
                print(f"Processing variable: {var}, Expression: {expr}")
                if isinstance(expr, (int, float)):
                    # If the expression is a constant, create a constant response
                    responses[str(var)] = np.full(len(t_vals), float(expr))
                else:
                    # Simplify the symbolic expression
                    expr = sp.simplify(expr)
                    print(f"Simplified expr for {var}: {expr}")

                    if expr.is_constant():
                        # Handle constant expressions
                        responses[str(var)] = np.full(len(t_vals), float(expr))
                    else:
                        # Perform inverse Laplace transform
                        time_expr = inverse_laplace_transform(expr, self.s, t)

                        # Handle special case for DiracDelta
                        if 'DiracDelta' in str(time_expr):
                            print(f"Handling DiracDelta for {var}")
                            coeff = 1.0
                            if isinstance(time_expr, sp.Mul):
                                for arg in time_expr.args:
                                    if not arg.has(sp.DiracDelta):
                                        coeff *= float(arg)

                            pulse_width = 5 * dt
                            sigma = dt
                            values = np.zeros(len(t_vals))
                            for i, t_val in enumerate(t_vals):
                                if t_val < pulse_width:
                                    values[i] = coeff * (1.0 / (sigma * np.sqrt(2 * np.pi))) * \
                                                np.exp(-0.5 * (t_val / sigma) ** 2)

                            pulse_area = np.trapz(values, t_vals)
                            if pulse_area != 0:
                                values *= coeff / pulse_area

                            if np.max(np.abs(values)) > MAX_AMPLITUDE:
                                values *= MAX_AMPLITUDE / np.max(np.abs(values))
                        else:
                            # Lambdify for numerical evaluation
                            try:
                                time_func = lambdify(t, time_expr, modules=["numpy"])
                                values = time_func(t_vals)
                            except Exception as e:
                                print(f"Error in lambdifying {var}: {e}")
                                values = np.zeros(len(t_vals))
                        values = np.real(values)
                        # Limit amplitude if necessary
                        if np.max(np.abs(values)) > MAX_AMPLITUDE:
                            values *= MAX_AMPLITUDE / np.max(np.abs(values))

                        # Handle large discontinuities
                        for i in range(1, len(values)):
                            if abs(values[i] - values[i - 1]) > MAX_AMPLITUDE:
                                values[i] = values[i - 1]

                        responses[str(var)] = values
            except Exception as e:
                print(f"Error processing {var}: {e}")
                responses[str(var)] = np.zeros(len(t_vals))

        return t_vals, responses


    @staticmethod
    def create_circuit_visualization(components, nodes):
        """
        Create an enhanced circuit visualization with distinct labels for AC sources and other components.
        Avoids overlapping wires and provides a professional circuit-like appearance.
        """
        # Create a graph for the circuit
        G = nx.Graph()
        edge_labels = {}  # Store edge labels for components

        # Add edges and labels for each component
        for comp in components:
            label = f"{comp['name']} ({comp['value']})"
            if comp['name'] == 'V' and comp.get('ac_params'):  # Handle AC source separately
                ac_params = comp['ac_params']
                label = f"AC\n{ac_params['magnitude']}V @ {ac_params['frequency']}Hz"
            G.add_edge(comp['n1'], comp['n2'], label=label)
            edge_labels[(comp['n1'], comp['n2'])] = label

        # Use Graphviz for circuit layout
        try:
            pos = nx.nx_agraph.graphviz_layout(G, prog="dot")  # Circuit layout using dot
        except ImportError:
            pos = nx.spring_layout(G)

        # Plot the graph
        plt.figure(figsize=(10, 10))
        nx.draw_networkx_nodes(G, pos, node_size=500, node_color="lightgray", edgecolors="black")
        nx.draw_networkx_edges(G, pos, edge_color="black", width=1.5)
        nx.draw_networkx_labels(G, pos, font_size=10, font_color="black", font_weight="bold")

        # Add edge labels
        nx.draw_networkx_edge_labels(
            G,
            pos,
            edge_labels=edge_labels,
            font_size=9,
            font_color="blue",
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.8),
        )
        for comp in components:
            if comp['name'] == 'V' and comp.get('ac_params'):  
                n1, n2 = comp['n1'], comp['n2']
                mid_x = (pos[n1][0] + pos[n2][0]) / 2
                mid_y = (pos[n1][1] + pos[n2][1]) / 2
                plt.text(
                    mid_x,
                    mid_y + 0.1,
                    "AC Source",
                    fontsize=10,
                    fontweight="bold",
                    color="red",
                    ha="center",
                    va="center",
                    bbox=dict(facecolor="yellow", edgecolor="red", boxstyle="round,pad=0.2", alpha=0.5),
                )

        # Add title
        plt.title("Circuit Diagram", fontsize=14)
        plt.axis("off")

        # Save the plot to a BytesIO object for embedding
        buffer = io.BytesIO()
        plt.savefig(buffer, format="png", bbox_inches="tight")
        buffer.seek(0)

        # Convert the image to base64
        encoded_image = base64.b64encode(buffer.read()).decode()
        buffer.close()

        # Return HTML for embedding in Streamlit
        return f'''
        <div style="display: flex; justify-content: center; align-items: center; margin: 20px 0;">
            <img src="data:image/png;base64,{encoded_image}" 
                style="max-width: 100%; height: auto; background-color: white; padding: 20px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1);"
            />
        </div>
        '''



    def _format_ac_params(self, component):
        """Format AC parameters for display"""
        if component.get('ac_params'):
            params = component['ac_params']
            return f" ({params['magnitude']}V, {params['frequency']}Hz {params['type']})"
        return ""
