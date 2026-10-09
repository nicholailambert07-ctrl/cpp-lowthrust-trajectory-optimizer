import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# 1. Load the C++ simulation output
try:
    df = pd.read_csv('optimal_trajectory.csv')
except FileNotFoundError:
    print("Error: 'optimal_trajectory.csv' not found. Run your compiled C++ binary first!")
    exit(1)

# 2. Setup Figure and Layout
fig = plt.figure(figsize=(14, 10))
fig.suptitle('C++ Low-Thrust Optimal Orbit Transfer Simulation', fontsize=16, fontweight='bold')

# --- Graph 1: 3D Trajectory Spiral ---
ax1 = fig.add_subplot(2, 2, 1, projection='3d')
ax1.plot(df['x_km'], df['y_km'], df['z_km'], color='cyan', label='Spacecraft Path', linewidth=0.8)

# Render Earth sphere for scale
u = np.linspace(0, 2 * np.pi, 30)
v = np.linspace(0, np.pi, 30)
r_earth_km = 6378.137
x_e = r_earth_km * np.outer(np.cos(u), np.sin(v))
y_e = r_earth_km * np.outer(np.sin(u), np.sin(v))
z_e = r_earth_km * np.outer(np.ones(np.size(u)), np.cos(v))
ax1.plot_surface(x_e, y_e, z_e, color='blue', alpha=0.3)

ax1.set_title('3D Trajectory Spiral')
ax1.set_xlabel('X (km)')
ax1.set_ylabel('Y (km)')
ax1.set_zlabel('Z (km)')

# --- Graph 2: Altitude Profile ---
ax2 = fig.add_subplot(2, 2, 2)
ax2.plot(df['time_days'], df['altitude_km'], color='crimson', linewidth=2)
ax2.set_title('Altitude Profile')
ax2.set_xlabel('Time (days)')
ax2.set_ylabel('Altitude (km)')
ax2.grid(True)

# --- Graph 3: Orbital Inclination ---
ax3 = fig.add_subplot(2, 2, 3)
ax3.plot(df['time_days'], df['inclination_deg'], color='purple', linewidth=2)
ax3.set_title('Orbital Inclination Change')
ax3.set_xlabel('Time (days)')
ax3.set_ylabel('Inclination (deg)')
ax3.grid(True)

# --- Graph 4: Spacecraft Mass Depletion ---
ax4 = fig.add_subplot(2, 2, 4)
ax4.plot(df['time_days'], df['mass_kg'], color='forestgreen', linewidth=2)
ax4.set_title('Spacecraft Mass Depletion')
ax4.set_xlabel('Time (days)')
ax4.set_ylabel('Mass (kg)')
ax4.grid(True)

plt.tight_layout()
plt.savefig('orbital_transfer_plots.png', dpi=300) # Saves high-res image
plt.show() # Displays interactive window