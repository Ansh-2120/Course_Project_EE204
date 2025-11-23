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
from calculation import CircuitSimulator

def display_matrix(matrix, title, row_labels=None, col_labels=None):
    """Display a matrix with labels"""
    st.write(f"**{title}**")
    df = pd.DataFrame(matrix)
    
    if row_labels:
        df.index = row_labels
    if col_labels:
        df.columns = col_labels
    
    st.dataframe(df, use_container_width=True)

def display_equations(equations, title):
    """Display symbolic equations"""
    st.write(f"**{title}**")
    for i, eq in enumerate(equations, 1):
        st.latex(sp.latex(eq))

def main():
    st.set_page_config(page_title="Circuit Simulator", layout="wide")
    
    st.title("PySpyce - Advanced Circuit Analysis")
    st.caption("Matrix-Based Circuit Analysis with Laplace Transforms")
    
    if 'simulator' not in st.session_state:
        st.session_state.simulator = CircuitSimulator()
    
    # Clean layout with columns
    col1, col2, col3, col4 = st.columns([2, 2, 1, 1])
    
    with col1:
        component_type = st.selectbox(
            "Component",
            ['Resistor (R)', 'Inductor (L)', 'Capacitor (C)', 'Voltage Source (V)', 'Current Source (I)']
        )
        component_type = component_type[component_type.find("(") + 1:component_type.find(")")]
    
    with col2:
        units = {'R': 'Ω', 'L': 'H', 'C': 'F', 'V': 'V', 'I': 'A'}
        value = st.number_input(f"Value ({units.get(component_type, '')})", min_value=0.0, value=1.0)
    
    with col3:
        n1 = st.number_input("From Node", min_value=0, value=0)
    
    with col4:
        n2 = st.number_input("To Node", min_value=0, value=1)
    
    # AC source parameters
    if component_type == 'V':
        st.divider()
        is_ac = st.checkbox("AC Source")
        if is_ac:
            col1, col2, col3 = st.columns(3)
            with col1:
                ac_type = st.selectbox("Waveform", ['sine', 'cosine'])
            with col2:
                magnitude = st.number_input("Amplitude (V)", min_value=0.0, value=1.0)
            with col3:
                frequency = st.number_input("Frequency (Hz)", min_value=0.1, value=50.0)
            ac_params = {'type': ac_type, 'magnitude': magnitude, 'frequency': frequency}
        else:
            ac_params = None
    else:
        ac_params = None
    
    # Action buttons in columns
    st.divider()
    col1, col2, col3 = st.columns(3)
    
    with col1:
        if st.button("Add Component", use_container_width=True):
            try:
                st.session_state.simulator.add_component(component_type, value, n1, n2, ac_params=ac_params)
                st.success(f"Added {component_type} component successfully!")
            except ValueError as e:
                st.error(str(e))
    
    with col2:
        if st.button("Solve Circuit", use_container_width=True):
            if st.session_state.simulator.solve_circuit():
                st.success("Circuit solved successfully!")
            else:
                st.error("Failed to solve circuit. Check connections.")
    
    with col3:
        if st.button("Clear Circuit", use_container_width=True):
            st.session_state.simulator = CircuitSimulator()
            st.success("Circuit cleared")
    
    # Display circuit and results
    if st.session_state.simulator.components:
        st.divider()
        tabs = st.tabs(["Circuit Diagram", "Matrix Analysis", "Equations", "Time Domain", "Frequency Domain"])
        
        with tabs[0]:
            # Display circuit diagram
            circuit_svg = CircuitSimulator.create_circuit_visualization(
                st.session_state.simulator.components,
                st.session_state.simulator.nodes
            )
            st.components.v1.html(circuit_svg, height=1000, scrolling=True)
            
            # Component list
            st.subheader("Components")
            st.dataframe(
                pd.DataFrame([{
                    'Type': c['name'],
                    'Value': c['value'],
                    'From': f"Node {c['n1']}",
                    'To': f"Node {c['n2']}",
                    'AC Params': c.get('ac_params', '')
                } for c in st.session_state.simulator.components])
            )
        
        # Matrix Analysis Tab
        with tabs[1]:
            st.header("Circuit Topology Analysis")
            
            sim = st.session_state.simulator
            nodes_list = sorted(list(sim.nodes))
            
            col1, col2 = st.columns(2)
            
            with col1:
                st.subheader("Circuit Topology")
                st.write(f"**Number of Nodes:** {len(sim.nodes)}")
                st.write(f"**Number of Branches:** {len(sim.components)}")
                st.write(f"**Nodes:** {nodes_list}")
                
                if sim.incidence_matrix is not None:
                    st.divider()
                    branch_labels = [c['id'] for c in sim.components]
                    node_labels = [f"Node {n}" for n in nodes_list]
                    
                    display_matrix(
                        sim.incidence_matrix,
                        "Incidence Matrix (A)",
                        row_labels=node_labels,
                        col_labels=branch_labels
                    )
                    
                    st.info("**Incidence Matrix**: Rows represent nodes, columns represent branches. +1 indicates current leaving node, -1 indicates current entering node.")
            
            with col2:
                if sim.cutset_matrix is not None:
                    display_matrix(
                        sim.cutset_matrix,
                        "Cutset Matrix (Q)",
                        row_labels=[f"Cutset {i+1}" for i in range(sim.cutset_matrix.shape[0])],
                        col_labels=[c['id'] for c in sim.components]
                    )
                    st.info("**Cutset Matrix**: Each row represents a fundamental cutset. Used for KCL analysis.")
                
                st.divider()
                
                if sim.tieset_matrix is not None:
                    display_matrix(
                        sim.tieset_matrix,
                        "Tieset/Loop Matrix (B)",
                        row_labels=[f"Loop {i+1}" for i in range(sim.tieset_matrix.shape[0])],
                        col_labels=[c['id'] for c in sim.components]
                    )
                    st.info("**Tieset Matrix**: Each row represents a fundamental loop. Used for KVL analysis.")
        
        # Equations Tab
        with tabs[2]:
            st.header("Circuit Analysis Equations")
            sim = st.session_state.simulator
            
            col1, col2 = st.columns(2)
            
            with col1:
                st.subheader("1. Time-Domain Equations")
                st.write("Component relationships in time domain:")
                
                if sim.time_domain_equations:
                    for eq_data in sim.time_domain_equations:
                        st.write(f"**{eq_data['component']}** ({eq_data['type']}):")
                        st.latex(sp.latex(eq_data['equation']))
                else:
                    st.info("Solve circuit to generate equations")
                
                st.divider()
                
                st.subheader("3. KCL Equations")
                st.write("Kirchhoff's Current Law at each node:")
                
                if sim.kcl_equations:
                    for kcl in sim.kcl_equations:
                        st.write(f"**Node {kcl['node']}:**")
                        st.latex(sp.latex(kcl['equation']))
                else:
                    st.info("Solve circuit to generate KCL equations")
            
            with col2:
                st.subheader("2. Laplace-Domain Equations")
                st.write("After applying Laplace transform (with zero initial conditions):")
                
                if sim.laplace_equations:
                    for eq_data in sim.laplace_equations:
                        st.write(f"**{eq_data['component']}** - Impedance: `{eq_data['impedance']}`")
                        st.latex(sp.latex(eq_data['equation']))
                else:
                    st.info("Solve circuit to generate Laplace equations")
                
                st.divider()
                
                st.subheader("4. KVL Equations")
                st.write("Kirchhoff's Voltage Law for fundamental loops:")
                
                if sim.kvl_equations:
                    for kvl in sim.kvl_equations:
                        st.write(f"**Loop through:** {', '.join(kvl['loop'])}")
                        st.latex(sp.latex(kvl['equation']))
                else:
                    st.info("No loops detected or solve circuit first")
            
            st.divider()
            
            # Final System
            if st.session_state.simulator.solutions:
                st.subheader("5. Complete System Solution (s-domain)")
                st.write("**Variables solved:**")
                
                for var, expr in st.session_state.simulator.solutions.items():
                    st.write(f"**{var}:**")
                    st.latex(f"{sp.latex(var)} = {sp.latex(expr)}")
                
                st.success("✓ System solved using SymPy's equation solver with Gaussian elimination")
        
        # Time domain tab
        if st.session_state.simulator.solutions:
            with tabs[3]:
                st.header("Time-Domain Response (Inverse Laplace Transform)")
                
                col1, col2 = st.columns([1, 3])
                with col1:
                    st.write("**Analysis Parameters:**")
                    tmax = st.slider("Time Range (s)", 0.001, 0.1, 0.01)
                    points = st.slider("Points", 100, 2000, 1000)
                    
                    st.info("Response obtained by:\n1. Solving in s-domain\n2. Applying inverse Laplace transform\n3. Numerical evaluation")
                
                with col2:
                    t_vals, responses = st.session_state.simulator.get_time_domain_response(tmax, points)
                    
                    fig = go.Figure()
                    for var, values in responses.items():
                        label = var
                        if var.startswith('I_'):
                            comp_id = var[2:]
                            label = f"Current through {comp_id}"
                        elif var.startswith('V'):
                            label = f"Voltage at Node {var[1:]}"
                        
                        fig.add_trace(go.Scatter(x=t_vals*1000, y=values, name=label))
                    
                    fig.update_layout(
                        xaxis_title="Time (ms)",
                        yaxis_title="Value",
                        height=600,
                        yaxis=dict(
                            range=[
                                np.percentile(np.concatenate([v for v in responses.values()]), 1),
                                np.percentile(np.concatenate([v for v in responses.values()]), 99)
                            ]
                        )
                    )
                    st.plotly_chart(fig, use_container_width=True)
            
            # Frequency domain tab
            with tabs[4]:
                st.header("Frequency-Domain Response (Bode Plots)")
                
                col1, col2 = st.columns([1, 3])
                with col1:
                    st.write("**Analysis Parameters:**")
                    fmin = st.number_input("Min Frequency (Hz)", value=1.0, min_value=0.1)
                    fmax = st.number_input("Max Frequency (Hz)", value=1000.0, min_value=fmin)
                    points = st.slider("Frequency Points", 100, 1000, 500)
                    
                    st.info("Frequency response by substituting s = jω in transfer functions")
                
                with col2:
                    f, responses = st.session_state.simulator.get_frequency_response(fmin, fmax, points)
                    
                    # Magnitude plot
                    fig_mag = go.Figure()
                    for var, response in responses.items():
                        label = var
                        if var.startswith('I_'):
                            comp_id = var[2:]
                            label = f"Current through {comp_id}"
                        elif var.startswith('V'):
                            label = f"Voltage at Node {var[1:]}"
                        
                        fig_mag.add_trace(go.Scatter(x=f, y=response['magnitude'], name=label))
                    
                    tick_values = []
                    tick_labels = []
                    current_decade = int(np.floor(np.log10(fmin)))
                    end_decade = int(np.ceil(np.log10(fmax)))
                    
                    for decade in range(current_decade, end_decade + 1):
                        value = 10**decade
                        if fmin <= value <= fmax:
                            tick_values.append(value)
                            if value >= 1000:
                                tick_labels.append(f"{value/1000:.0f} kHz")
                            else:
                                tick_labels.append(f"{value:.0f} Hz")
                    
                    fig_mag.update_layout(
                        title="Magnitude Response",
                        xaxis_title="Frequency",
                        yaxis_title="Magnitude",
                        xaxis_type="log",
                        height=600,
                        xaxis=dict(tickvals=tick_values, ticktext=tick_labels),
                    )
                    st.plotly_chart(fig_mag, use_container_width=True)
                    
                    # Phase plot
                    fig_phase = go.Figure()
                    for var, response in responses.items():
                        label = var
                        if var.startswith('I_'):
                            comp_id = var[2:]
                            label = f"Current through {comp_id}"
                        elif var.startswith('V'):
                            label = f"Voltage at Node {var[1:]}"
                        
                        fig_phase.add_trace(go.Scatter(x=f, y=response['phase'], name=label))
                    
                    fig_phase.update_layout(
                        title="Phase Response",
                        xaxis_title="Frequency",
                        yaxis_title="Phase (degrees)",
                        xaxis_type="log",
                        height=600,
                        xaxis=dict(tickvals=tick_values, ticktext=tick_labels),
                    )
                    st.plotly_chart(fig_phase, use_container_width=True)

if __name__ == "__main__":
    main()