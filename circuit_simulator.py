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
from calculations import CircuitSimulator

def main():
    st.set_page_config(page_title="Circuit Simulator", layout="wide")
    
    st.title("PySpyce")
    
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
        tabs = st.tabs(["Circuit", "Time Domain", "Frequency Domain"])
        
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
        
        # Time domain tab
        if st.session_state.simulator.solutions:
            with tabs[1]:
                col1, col2 = st.columns([1, 3])
                with col1:
                    tmax = st.slider("Time Range (s)", 0.001, 0.1, 0.01)
                    points = st.slider("Points", 100, 2000, 1000)
                
                t_vals, responses = st.session_state.simulator.get_time_domain_response(tmax, points)
                
                fig = go.Figure()
                for var, values in responses.items():
                    # Convert variable names to descriptive labels
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
                    # Set y-axis range to exclude extreme spikes
                    yaxis=dict(
                        range=[
                            np.percentile(np.concatenate([v for v in responses.values()]), 1),
                            np.percentile(np.concatenate([v for v in responses.values()]), 99)
                        ]
                    )
                )
                st.plotly_chart(fig, use_container_width=True)
            
            # Frequency domain tab
            with tabs[2]:
                col1, col2 = st.columns([1, 3])
                with col1:
                    fmin = st.number_input("Min Frequency (Hz)", value=1.0, min_value=0.1)
                    fmax = st.number_input("Max Frequency (Hz)", value=1000.0, min_value=fmin)
                    points = st.slider("Frequency Points", 100, 1000, 500)
                
                f, responses = st.session_state.simulator.get_frequency_response(fmin, fmax, points)
                
                # Magnitude plot
                fig_mag = go.Figure()
                for var, response in responses.items():
                    # Convert variable names to descriptive labels
                    label = var
                    if var.startswith('I_'):
                        comp_id = var[2:]
                        label = f"Current through {comp_id}"
                    elif var.startswith('V'):
                        label = f"Voltage at Node {var[1:]}"
                    
                    fig_mag.add_trace(go.Scatter(x=f, y=response['magnitude'], name=label))
                
                # Calculate appropriate tick values based on the frequency range
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
                    xaxis_title="Frequency",
                    yaxis_title="Magnitude",
                    xaxis_type="log",
                    height=600,
                    xaxis=dict(tickvals=tick_values, ticktext=tick_labels),
                )
                st.plotly_chart(fig_mag, use_container_width=True)
                
                # Phase plot with same x-axis formatting
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
                    xaxis_title="Frequency",
                    yaxis_title="Phase (degrees)",
                    xaxis_type="log",
                    height=600,
                    xaxis=dict(tickvals=tick_values, ticktext=tick_labels),
                )
                st.plotly_chart(fig_phase, use_container_width=True)

if __name__ == "__main__":
    main()
