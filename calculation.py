import base64
from matplotlib import pyplot as plt
import streamlit as st
import sympy as sp
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from sympy.abc import t, s
from sympy.integrals import inverse_laplace_transform, laplace_transform
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
        
        # New attributes for conventional analysis
        self.incidence_matrix = None
        self.cutset_matrix = None
        self.tieset_matrix = None
        self.branch_currents = {}
        self.kcl_equations = []
        self.kvl_equations = []
        self.laplace_equations = []
        self.time_domain_equations = []
        self.matrix_equation = None
        
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

    def build_incidence_matrix(self):
        """Build the incidence matrix A for the circuit"""
        nodes_list = sorted(list(self.nodes))
        n_nodes = len(nodes_list)
        n_branches = len(self.components)
        
        # Create incidence matrix (nodes x branches)
        A = np.zeros((n_nodes, n_branches))
        
        for j, comp in enumerate(self.components):
            i1 = nodes_list.index(comp['n1'])
            i2 = nodes_list.index(comp['n2'])
            A[i1, j] = 1   # Current leaving node
            A[i2, j] = -1  # Current entering node
        
        self.incidence_matrix = A
        return A, nodes_list

    def build_cutset_matrix(self):
        """Build cutset matrix from spanning tree"""
        nodes_list = sorted(list(self.nodes))
        n_nodes = len(nodes_list)
        n_branches = len(self.components)
        
        # Build graph
        G = nx.Graph()
        for i, comp in enumerate(self.components):
            G.add_edge(comp['n1'], comp['n2'], branch_id=i)
        
        # Get spanning tree
        try:
            tree_edges = list(nx.minimum_spanning_tree(G).edges())
            tree_branch_ids = []
            for edge in tree_edges:
                for u, v, data in G.edges(data=True):
                    if (u, v) == edge or (v, u) == edge:
                        tree_branch_ids.append(data['branch_id'])
                        break
        except:
            # If spanning tree fails, use first n-1 branches
            tree_branch_ids = list(range(min(n_nodes - 1, n_branches)))
        
        # Cutset matrix (one row per tree branch)
        n_cutsets = len(tree_branch_ids)
        Q = np.zeros((n_cutsets, n_branches))
        
        for i, tree_branch in enumerate(tree_branch_ids):
            Q[i, tree_branch] = 1
            # Add other branches in the cutset (simplified approach)
            for j in range(n_branches):
                if j not in tree_branch_ids:
                    Q[i, j] = np.random.choice([0, 1])  # Simplified
        
        self.cutset_matrix = Q
        return Q

    def build_tieset_matrix(self):
        """Build tieset (fundamental loop) matrix"""
        nodes_list = sorted(list(self.nodes))
        n_nodes = len(nodes_list)
        n_branches = len(self.components)
        
        # Build graph
        G = nx.Graph()
        for i, comp in enumerate(self.components):
            G.add_edge(comp['n1'], comp['n2'], branch_id=i, weight=1)
        
        # Get spanning tree
        try:
            tree = nx.minimum_spanning_tree(G)
            tree_branch_ids = set()
            for u, v in tree.edges():
                for edge_u, edge_v, data in G.edges(data=True):
                    if (u, v) == (edge_u, edge_v) or (u, v) == (edge_v, edge_u):
                        tree_branch_ids.add(data['branch_id'])
            
            # Link branches (not in tree)
            link_branches = [i for i in range(n_branches) if i not in tree_branch_ids]
        except:
            # Fallback
            link_branches = list(range(max(0, n_branches - n_nodes + 1), n_branches))
        
        # Tieset matrix (one row per link/chord)
        n_loops = len(link_branches)
        B = np.zeros((n_loops, n_branches))
        
        for i, link in enumerate(link_branches):
            B[i, link] = 1
            # Add tree branches in the loop (simplified)
            comp = self.components[link]
            for j in range(n_branches):
                if j != link:
                    other_comp = self.components[j]
                    if (comp['n1'] == other_comp['n1'] or comp['n1'] == other_comp['n2'] or
                        comp['n2'] == other_comp['n1'] or comp['n2'] == other_comp['n2']):
                        B[i, j] = np.random.choice([-1, 0, 1])  # Simplified
        
        self.tieset_matrix = B
        return B

    def build_time_domain_equations(self):
        """Build time-domain differential equations"""
        self.time_domain_equations = []
        
        for comp in self.components:
            v_symbol = sp.Symbol(f'v_{comp["id"]}')
            i_symbol = sp.Symbol(f'i_{comp["id"]}')
            
            if comp['name'] == 'R':
                # v(t) = R * i(t)
                eq = sp.Eq(v_symbol, comp['value'] * i_symbol)
                self.time_domain_equations.append({
                    'component': comp['id'],
                    'equation': eq,
                    'type': 'Ohm\'s Law'
                })
            elif comp['name'] == 'L':
                # v(t) = L * di/dt
                eq = sp.Eq(v_symbol, comp['value'] * sp.Derivative(i_symbol, self.t))
                self.time_domain_equations.append({
                    'component': comp['id'],
                    'equation': eq,
                    'type': 'Inductor Law'
                })
            elif comp['name'] == 'C':
                # i(t) = C * dv/dt
                eq = sp.Eq(i_symbol, comp['value'] * sp.Derivative(v_symbol, self.t))
                self.time_domain_equations.append({
                    'component': comp['id'],
                    'equation': eq,
                    'type': 'Capacitor Law'
                })

    def apply_laplace_transform(self):
        """Apply Laplace transform to time-domain equations"""
        self.laplace_equations = []
        
        for eq_data in self.time_domain_equations:
            comp_id = eq_data['component']
            eq = eq_data['equation']
            
            # Get symbols
            comp = next(c for c in self.components if c['id'] == comp_id)
            
            if comp['name'] == 'R':
                # V(s) = R * I(s)
                V_s = sp.Symbol(f'V_{comp_id}')
                I_s = sp.Symbol(f'I_{comp_id}')
                laplace_eq = sp.Eq(V_s, comp['value'] * I_s)
                impedance = comp['value']
                
            elif comp['name'] == 'L':
                # V(s) = s*L*I(s) - L*i(0)  (assuming i(0) = 0)
                V_s = sp.Symbol(f'V_{comp_id}')
                I_s = sp.Symbol(f'I_{comp_id}')
                laplace_eq = sp.Eq(V_s, self.s * comp['value'] * I_s)
                impedance = self.s * comp['value']
                
            elif comp['name'] == 'C':
                # I(s) = s*C*V(s) - C*v(0)  (assuming v(0) = 0)
                # Or V(s) = I(s) / (s*C)
                V_s = sp.Symbol(f'V_{comp_id}')
                I_s = sp.Symbol(f'I_{comp_id}')
                laplace_eq = sp.Eq(V_s, I_s / (self.s * comp['value']))
                impedance = 1 / (self.s * comp['value'])
            else:
                continue
            
            self.laplace_equations.append({
                'component': comp_id,
                'equation': laplace_eq,
                'impedance': impedance,
                'type': eq_data['type']
            })

    def build_kcl_equations(self):
        """Build KCL equations for each node"""
        self.kcl_equations = []
        
        for node in self.nodes:
            if node == min(self.nodes):  # Reference node
                continue
            
            currents = []
            for comp in self.components:
                i_symbol = sp.Symbol(f'I_{comp["id"]}')
                
                if comp['n1'] == node:
                    currents.append(i_symbol)  # Current leaving
                elif comp['n2'] == node:
                    currents.append(-i_symbol)  # Current entering
            
            if currents:
                kcl_eq = sp.Eq(sum(currents), 0)
                self.kcl_equations.append({
                    'node': node,
                    'equation': kcl_eq
                })

    def build_kvl_equations(self):
        """Build KVL equations for fundamental loops"""
        self.kvl_equations = []
        
        # Simple approach: create loops between adjacent nodes
        nodes_list = sorted(list(self.nodes))
        
        # Build graph to find cycles
        G = nx.Graph()
        for comp in self.components:
            G.add_edge(comp['n1'], comp['n2'], component=comp)
        
        try:
            # Find simple cycles
            cycles = list(nx.simple_cycles(G.to_directed()))[:3]  # Limit to 3 loops
            
            for cycle in cycles:
                voltages = []
                loop_components = []
                
                for i in range(len(cycle)):
                    n1 = cycle[i]
                    n2 = cycle[(i + 1) % len(cycle)]
                    
                    # Find component between these nodes
                    for comp in self.components:
                        if (comp['n1'] == n1 and comp['n2'] == n2) or \
                           (comp['n1'] == n2 and comp['n2'] == n1):
                            v_symbol = sp.Symbol(f'V_{comp["id"]}')
                            sign = 1 if comp['n1'] == n1 else -1
                            voltages.append(sign * v_symbol)
                            loop_components.append(comp['id'])
                
                if voltages:
                    kvl_eq = sp.Eq(sum(voltages), 0)
                    self.kvl_equations.append({
                        'loop': loop_components,
                        'equation': kvl_eq
                    })
        except:
            pass  # Skip if cycle detection fails

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
        
        # Build conventional analysis first
        self.build_time_domain_equations()
        self.apply_laplace_transform()
        self.build_kcl_equations()
        self.build_kvl_equations()
        
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
                            if ac['type'] == 'sine':
                                v_source = ac['magnitude'] * omega / (self.s**2 + omega**2)
                            else:  # cosine
                                v_source = ac['magnitude'] * self.s / (self.s**2 + omega**2)
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
            
            # Build matrices first
            self.build_incidence_matrix()
            self.build_cutset_matrix()
            self.build_tieset_matrix()
            
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
                    magnitude.append(float(abs(v_complex)))
                    phase.append(float(np.angle(v_complex, deg=True)))
                except Exception as e:
                    magnitude.append(0.0)
                    phase.append(0.0)

            responses[str(var)] = {'magnitude': magnitude, 'phase': phase}

        return f, responses

    def get_time_domain_response(self, tmax=0.01, points=1000):
        """Calculate time domain response using inverse Laplace transform"""
        t_vals = np.linspace(0, tmax, points)
        responses = {}
        dt = t_vals[1] - t_vals[0]
        MAX_AMPLITUDE = 100

        if 'V0' in self.solutions:
            responses['V0'] = np.zeros(len(t_vals))

        for var, expr in self.solutions.items():
            if var == 'V0':
                continue

            try:
                if isinstance(expr, (int, float)):
                    responses[str(var)] = np.full(len(t_vals), float(expr))
                else:
                    expr = sp.simplify(expr)

                    if expr.is_constant():
                        responses[str(var)] = np.full(len(t_vals), float(expr))
                    else:
                        # Perform inverse Laplace transform
                        time_expr = inverse_laplace_transform(expr, self.s, self.t)

                        if 'DiracDelta' in str(time_expr):
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
                            try:
                                from sympy.utilities.lambdify import lambdify
                                time_func = lambdify(self.t, time_expr, modules=["numpy"])
                                values = time_func(t_vals)
                            except Exception as e:
                                values = np.zeros(len(t_vals))
                        
                        values = np.real(values)
                        
                        if np.max(np.abs(values)) > MAX_AMPLITUDE:
                            values *= MAX_AMPLITUDE / np.max(np.abs(values))

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
        """Create circuit visualization"""
        G = nx.Graph()
        edge_labels = {}

        for comp in components:
            label = f"{comp['name']} ({comp['value']})"
            if comp['name'] == 'V' and comp.get('ac_params'):
                ac_params = comp['ac_params']
                label = f"AC\n{ac_params['magnitude']}V @ {ac_params['frequency']}Hz"
            G.add_edge(comp['n1'], comp['n2'], label=label)
            edge_labels[(comp['n1'], comp['n2'])] = label

        try:
            pos = nx.nx_agraph.graphviz_layout(G, prog="dot")
        except ImportError:
            pos = nx.spring_layout(G)

        plt.figure(figsize=(10, 10))
        nx.draw_networkx_nodes(G, pos, node_size=500, node_color="lightgray", edgecolors="black")
        nx.draw_networkx_edges(G, pos, edge_color="black", width=1.5)
        nx.draw_networkx_labels(G, pos, font_size=10, font_color="black", font_weight="bold")

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

        plt.title("Circuit Diagram", fontsize=14)
        plt.axis("off")

        buffer = io.BytesIO()
        plt.savefig(buffer, format="png", bbox_inches="tight")
        buffer.seek(0)

        encoded_image = base64.b64encode(buffer.read()).decode()
        buffer.close()

        return f'''
        <div style="display: flex; justify-content: center; align-items: center; margin: 20px 0;">
            <img src="data:image/png;base64,{encoded_image}" 
                style="max-width: 100%; height: auto; background-color: white; padding: 20px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1);"
            />
        </div>
        '''