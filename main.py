# -----------------------------------------------------------------------------
# Script: main
# Author: Ludwig Horvath
# Email: ludhor@kth.se
# Date: 2024-10-05
# Description: This script runs the main simulation, orchestrating the interactions 
#              between different classes and managing the overall workflow of the 
#              orbital simulation.
# -----------------------------------------------------------------------------

# Link with orbital parameters: https://ssd.jpl.nasa.gov/orbits.html 
from matplotlib.patches import Circle
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from Classes import *

def plot_orbital_elements(orbit, detail='low'):
    """
    Plots the orbital elements for a given orbit object.

    Parameters:
        orbit: The orbit object containing the orbital parameters.
        detail (str): The detail level of the plot ('low', 'medium', or 'high').

    Returns:
        ax: The axes with the plotted orbital elements.
    """
    e = orbit

    # Plot the orbital curve
    ax.plot(e.curve[0,:], e.curve[1,:], e.curve[2,:], color=e.color)

    if detail in ('medium', 'high'):
        # Plot the key points: Center (C), Focus (F), Aphelion (A), Perihelion (P)
        ax.scatter(e.C[0], e.C[1], e.C[2], color=e.color, label='Z') 
        ax.scatter(e.F[0], e.F[1], e.F[2], color=e.color, label='F')
        ax.scatter(e.AP[0], e.AP[1], e.AP[2], color=e.color, label='A')
        ax.scatter(e.PE[0], e.PE[1], e.PE[2], color=e.color, label='P')

        # Plot the lines connecting key points
        ax.plot(e.line_AP[:,0], e.line_AP[:,1], e.line_AP[:,2], color=e.color, linewidth=0.5)
        ax.plot(e.line_Vert[:,0], e.line_Vert[:,1], e.line_Vert[:,2], color=e.color, linewidth=0.5)

    if detail in ('high'):
        # Add arrows indicating the orbital direction
        ax.add_artist(e.arrow_x)
        ax.add_artist(e.arrow_y)
        ax.add_artist(e.arrow_z)

    return ax

# Constants
G = 6.67430 * 10**(-11)  # Universal gravitational constant
AU = 149597871 * 10**3    # [m/1 AU]
errtol = 0.01
t0 = 0
sec_per_tunit = 3600 * 24  # time unit is [day]
dt = sec_per_tunit         # time step [1 day]
earthyear = 365
close = False

while not close:
    # Input for simulation time
    sim_time = float(input('⏳ Enter simulation time in Earth-years: '))
    tend = sim_time * earthyear * sec_per_tunit
    time = np.arange(t0, tend, dt)

    # Create solar system simulation
    solarsystem = System('SolarSystem')
    solarsystem.create()

    print('🔄 Solar system successfully created! Press Enter to continue...')

    # Initialize position array
    positions = np.zeros((len(solarsystem.orbits) * 3, len(time)))
    print('🚀 Running simulation... Please wait...\n')

    # Run simulation for the specified time
    for index, t in enumerate(time):
        for nr, orbit in enumerate(solarsystem.orbits):
            orbit.time_step()
            positions[nr * 3:3 + 3 * nr, index] = orbit.r_vec

    print('✅ Simulation complete! Now entering visualization phase...\n')

    # Define planet colors
    planet_colors = {
        'Sun': '#ffd700',
        'Mercury': '#b0b0b0',
        'Venus': '#e1ce7a',
        'Earth': '#1f77b4',
        'Mars': '#ff5733',
        'Jupiter': '#ffad33',
        'Saturn': '#eedc82',
        'Uranus': '#7fdbff',
        'Neptune': '#5073b8',
        'Pluto': '#be93c5'
    }

    solarsystem.colordictionary = planet_colors
    markers = {}

    print('📊 Entering visualization phase...\n')

    rerun = True

    while rerun:
        # Define default visualization settings
        visualized_planet = {
            'Mercury': 1,
            'Venus': 1,
            'Earth': 1,
            'Mars': 1,
            'Jupiter': 1,
            'Saturn': 1,
            'Uranus': 1,
            'Neptune': 1,
            'Pluto': 1
        }

        # Input for planet selection
        selection_choice = input('All planets ("Enter") Selected planets ("S"): ')

        if selection_choice == "S":
            for planet in visualized_planet.keys():
                planet_choice = input('Include ("Enter") ' + str(planet) + ', Exclude ("D") ' + str(planet) + ': ')
                if planet_choice == "D":
                    visualized_planet[planet] = 0

        planet_array = [solarsystem.planets[planet] for planet in visualized_planet]

        # Input for plot choice
        plot_choice = input('2D (1) or 3D (2) animation?: ').strip()

        if plot_choice == '1':
            # 2D Animation
            fig, ax = plt.subplots()

            ax.set_facecolor('black')
            ax.set_aspect('equal')
            fig.patch.set_facecolor('black')
            ax.set_axis_off()

            for nr, planet_obj in enumerate(planet_array):
                if visualized_planet[planet_obj.name]:
                    point, = ax.plot([], [], 'o', color=solarsystem.colordictionary[planet_obj.name], linewidth=0.4)
                    markers[planet_obj.name] = point
                    ax.plot(positions[3 * nr, :], positions[3 * nr + 1, :],
                            color=solarsystem.colordictionary[planet_obj.name], linewidth=0.4)

            # Plot the Sun as a circle
            gui_sun = Circle((0, 0), solarsystem.centerbody.radius,
                             color=solarsystem.colordictionary[solarsystem.centerbody.name],
                             label=solarsystem.centerbody.name)
            ax.add_patch(gui_sun)

            # Update function for 2D animation
            def update(frame):
                for nr, planet_obj in enumerate(planet_array):
                    if visualized_planet[planet_obj.name]:
                        x, y = positions[3 * nr, frame], positions[3 * nr + 1, frame]
                        markers[planet_obj.name].set_data(x, y)
                return list(markers.values())

            ani = FuncAnimation(fig, update, frames=range(1, len(time), 1), interval=10, repeat=False)
            plt.show()

        elif plot_choice == '2':
            # 3D Animation
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')

            # Hide grid lines and ticks
            ax.grid(False)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_zticks([])
            ax.set_facecolor('black')
            fig.patch.set_facecolor('black')
            ax.set_axis_off()
            ax.xaxis.set_pane_color((0.0, 0.0, 0.0, 1.0))
            ax.yaxis.set_pane_color((0.0, 0.0, 0.0, 1.0))
            ax.zaxis.set_pane_color((0.0, 0.0, 0.0, 1.0))

            # Set axes limits based on positions
            max_range = np.array([positions.max() - positions.min() for positions in [positions[i::3, :] for i in range(3)]]).max() / 2.0
            mid_x = (positions[0, :].max() + positions[0, :].min()) * 0.5
            mid_y = (positions[1, :].max() + positions[1, :].min()) * 0.5
            mid_z = (positions[2, :].max() + positions[2, :].min()) * 0.5

            ax.set_xlim(mid_x - max_range, mid_x + max_range)
            ax.set_ylim(mid_y - max_range, mid_y + max_range)
            ax.set_zlim(mid_z - max_range, mid_z + max_range)
            ax.set_aspect('equal')

            # Static graphics for 3D animation
            for nr, planet_obj in enumerate(planet_array):
                if visualized_planet[planet_obj.name]:
                    ax.plot(positions[3 * nr, :], positions[3 * nr + 1, :], positions[3 * nr + 2, :],
                            color=solarsystem.colordictionary[planet_obj.name], linewidth=0.4)
                    point, = ax.plot([], [], [], 'o', color=solarsystem.colordictionary[planet_obj.name], markersize=5)
                    markers[planet_obj.name] = point
            
            # Plot orbital elements for each planet
            for orbit in solarsystem.orbits:
                if visualized_planet[orbit.body2.name]:
                    plot_orbital_elements(orbit)

            # Plot the Sun in 3D
            ax.scatter([0], [0], [0], color=solarsystem.colordictionary[solarsystem.centerbody.name],
                       s=500, marker='o', label=solarsystem.centerbody.name)

            # Update function for 3D animation
            def update(frame):
                for nr, planet_obj in enumerate(planet_array):
                    if visualized_planet[planet_obj.name]:
                        x, y, z = positions[3 * nr, frame], positions[3 * nr + 1, frame], positions[3 * nr + 2, frame]
                        markers[planet_obj.name].set_data(x, y)
                        markers[planet_obj.name].set_3d_properties(z)
                return list(markers.values())

            ani = FuncAnimation(fig, update, frames=range(1, len(time), 1), interval=10, repeat=False)
            plt.show()
        
        # Prompt to rerun the visualization
        rerun_choice = input('Rerun simulation and visualization? (y/n): ')
        rerun = rerun_choice.lower() == 'y'

print("Simulation and visualization complete!")
