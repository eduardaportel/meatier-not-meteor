import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from sgp4.api import Satrec, jday

PI = np.pi

class Rentry_analysis:
    def __init__(self):
        self.e2 = 0.00669437999014
        self.a = 6378137.0     # WGS84 semi-major axis [m]
        self.b = 6356752.3142  # WGS84 semi-minor axis [m]

    def Earth_WGS84_ECEF(self):
        # Mesh in latitude/longitude
        phi = np.linspace(-PI/2, PI/2, 100)   # latitude
        lam = np.linspace(-PI, PI, 200)       # longitude
        phi, lam = np.meshgrid(phi, lam)

        # Constant altitude
        h = np.full_like(phi, 0.0)  # ellipsoidal height

        # Prime vertical radius
        N = self.a / np.sqrt(1 - self.e2 * np.sin(phi)**2)

        # Coordinates in km
        X = (N + h) * np.cos(phi) * np.cos(lam) / 1e3
        Y = (N + h) * np.cos(phi) * np.sin(lam) / 1e3
        Z = ((1 - self.e2) * N + h) * np.sin(phi) / 1e3

        return X, Y, Z

    def propagate_TLE(self, tle_line1, tle_line2, minutes=3600, step=0.5):
        """Propagate a TLE using SGP4 (works for LEO & deep-space)."""
        
        satellite = Satrec.twoline2rv(tle_line1, tle_line2)
        
        # Choose start epoch (JDAY)
        jd, fr = jday(2025, 10, 4, 0, 0, 0.0)  # Year, month, day, h,m,s

        times = np.arange(0, minutes, step)  # propagation times in minutes
        sat_positions = []

        for t in times:
            e, r, v = satellite.sgp4(jd, fr + t/1440.0)  # propagate
            if e != 0:
                raise RuntimeError(f"SGP4 error code {e}")

            sat_positions.append(r)

        return np.array(sat_positions)

# ------------------------------
# Usage
# ------------------------------
class1 = Rentry_analysis()
X, Y, Z = class1.Earth_WGS84_ECEF()

# Example ISS TLE (2025-09-30 epoch)
tle1 = "1 XXXXU XXXXX   2461000.5  0.00016912  00000-0  30874-3 0  9992"
tle2 = "2 XXXXX  12.58732  306.50 067  195.6472  194.37 1.7496"

sat_pos = class1.propagate_TLE(tle1, tle2, minutes=3600, step=1)

# ------------------------------
# Plot
# ------------------------------
fig, ax = plt.subplots(subplot_kw={"projection": "3d"}, figsize=(8,8))

# Earth ellipsoid
ax.plot_surface(X, Y, Z, cmap=cm.Blues, alpha=0.6)

# Orbit trajectory
ax.plot(sat_pos[:,0], sat_pos[:,1], sat_pos[:,2], 'r', label="Orbit")

ax.set_xlabel("X [km]")
ax.set_ylabel("Y [km]")
ax.set_zlabel("Z [km]")
ax.set_title("SGP4 Orbit around WGS-84 Earth")
ax.legend()
plt.show()