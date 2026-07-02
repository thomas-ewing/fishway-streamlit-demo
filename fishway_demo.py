import streamlit as st
import matplotlib.pyplot as plt
import plotly.express as px
import pandas as pd
import numpy as np
from verticalslot import FishwayDefinition as fd

def plot_fishway_slot_loss(f):

    number_of_slots = len(f.slot_inverts)
    number_of_pools = number_of_slots-1

    fig = px.line(x=np.arange(number_of_slots)+1, y=f.slot_differential, title='Slot Loss',
                     labels={'x': 'Slot Number', 'y': 'Slot Loss [m]'},
                     range_x=[1-0.2, number_of_slots+0.2],
                     height=300, template='plotly_dark',
                     markers=True)
        
    fig.update_traces(marker=dict(size=10, symbol='x', color='red'),
                      hovertemplate='Slot %{x}: %{y:.3f}')
    fig.update_yaxes(range=[0, None], rangemode='tozero')
    fig.update_xaxes(tickmode='array', tickvals=np.arange(number_of_slots)+1)

    return fig

def plot_cell_turb(f):

    number_of_slots = len(f.slot_inverts)
    number_of_pools = number_of_slots-1

    fig = px.line(x=np.arange(number_of_pools)+1, y=f.pool_diss_rate, title='Cell Turbulence',
                     labels={'x': 'Cell Number', 'y': 'Cell Turbulence [W/m^3]'},
                     range_x=[1-0.2, number_of_pools+0.2],
                     height=300, template='plotly_dark',
                     markers=True)
    
    fig.update_traces(marker=dict(size=10, symbol='x', color='green'),
                      hovertemplate='Cell %{x}: %{y:.3f}')
    fig.update_yaxes(range=[0, None], rangemode='tozero')
    fig.update_xaxes(tickmode='array', tickvals=np.arange(number_of_pools)+1)

    return fig

def plot_fishway_result(f):

    plt.style.use('dark_background')

    number_of_slots = len(f.slot_inverts)
    number_of_pools = number_of_slots-1

    # Set up figure

    fig, ax = plt.subplots(figsize=[12,4])

    # Plot slot inverts

    ax.plot(slot_invert_array, 'o', color='yellow', label='Slot Inverts')

    # Plot pool floor levels

    floors_x = np.ravel((np.arange(number_of_pools+1),np.arange(number_of_pools+1),), order='F')[1:]
    floors_y = np.append(np.ravel((f.pool_inverts,f.pool_inverts,), order='F'), f.slot_inverts[-1])
    ax.plot(floors_x, floors_y, color='grey', label='Cell Floor Levels', linewidth=3)

    # Plot water levels

    water_level_x = np.append(np.append(-1, np.ravel((np.arange(number_of_pools+1),np.arange(number_of_pools+1),), order='F')), number_of_pools+1)
    water_level_y = np.pad(np.append(f.downstream_water_level, np.append(np.ravel((f.pool_levels,f.pool_levels,), order='F'), f.upstream_water_level)),1,mode='edge')
    ax.plot(water_level_x, water_level_y, color='deepskyblue', label='Cell Water Levels', linewidth=3)

    # Water level labels

    ax.annotate('DS: {:.3f}m'.format(f.downstream_water_level), 
            xy=(-0.5, f.downstream_water_level,), 
            xytext=(0, -5), 
            textcoords='offset points',
            fontweight='bold',
            rotation_mode='anchor',
            rotation=90,
            size='small',
            ha='right', va='center')
    
    ax.annotate('US: {:.3f}m'.format(f.upstream_water_level), 
            xy=(number_of_slots-0.5, f.upstream_water_level,), 
            xytext=(0, -5), 
            textcoords='offset points',
            fontweight='bold',
            rotation_mode='anchor',
            rotation=90,
            size='small',
            ha='right', va='center')

    for p in range(number_of_pools):

        ax.annotate('P{:.0f}: {:.3f}m'.format(p+1, f.pool_levels[p]), 
                xy=(p+0.5, f.pool_levels[p]), 
                xytext=(0, -5), 
                textcoords='offset points',
                rotation_mode='anchor',
                rotation=90,
                size='small',
                fontweight='bold',
                ha='right', va='center')


    # Decorate Plots

    ax.set_xlim([-2, number_of_slots+1])
    ax.set_xticks(np.arange(0,number_of_slots))
    ax.set_xticklabels(np.arange(0,number_of_slots)+1)
    ax.set_xlabel('Slot Number')
    ax.set_ylabel('Elevation [m]')
    plt.grid(axis='x', alpha=0.2)
    plt.legend()
    
    return fig

st.set_page_config(layout="wide")

st.title("Interactive Fishway Demo!")
st.write("Welcome to a simple fishway calculation demo!")

col1, col2 = st.columns([1,3])

downstream_water_level = col1.number_input("Downstream Water Level: ", value=1.0)
upstream_water_level = col1.number_input("Upstream Water Level: ", value=2.0)
entry_slot_level = col1.number_input("Entry Slot Level: ", value=0.0)
slot_increment = col1.number_input("Vertical Slot Increment: ", value=0.1, step=0.001, format="%0.3f")
number_of_slots = col1.number_input("Number of Slots: ", value=10)
slot_width = col1.number_input("Slot Width: ", value=0.1, step=0.001, format="%0.3f")
pool_area = col1.number_input("Pool Area: ", value=6.0)

slot_invert_array = np.linspace(entry_slot_level, entry_slot_level+((number_of_slots-1)*slot_increment), number_of_slots)
number_of_pools = number_of_slots-1
pool_invert_array = slot_invert_array[0:number_of_pools]
pool_area_array = pool_area*np.ones(number_of_pools)

temp_fw = fd(slot_invert_array, pool_invert_array, pool_area_array, upstream_water_level, downstream_water_level)
temp_fw.add_slot_definition(np.array([0.0,5.0]), np.array([0,5])*slot_width)
temp_fw.solve_hydraulics()

col2.pyplot(plot_fishway_result(temp_fw))

col2a, col2b = col2.columns([1,1])

col2a.plotly_chart(plot_fishway_slot_loss(temp_fw))

col2b.plotly_chart(plot_cell_turb(temp_fw))


