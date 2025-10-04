import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from sgp4.api import Satrec, jday
from astropy import units as u
from astropy.time import Time
from poliastro.bodies import Earth, Sun
from poliastro.ephem import Ephem
from poliastro.twobody import Orbit
import json

PI = np.pi

class Rentry_analysis:
    
    def __init__(self):
        self.e2 = 0.00669437999014
        self.a = 6378137.0     # WGS84 semi-major axis [m]
        self.b = 6356752.3142  # WGS84 semi-minor axis [m]

    def readData_from_json_NEO(self, path):
        with open(path, "r") as file:
            data = json.load(file)  
        orbital_dict = data["orbital_data"]
        diameter_dict = data["estimated_diameter"]["meters"]
        return orbital_dict, diameter_dict
    
    def earth_WGS84_ECEF(self):
        phi = np.linspace(-PI/2, PI/2, 200)   # latitude
        lam = np.linspace(-PI, PI, 200)       # longitude
        phi, lam = np.meshgrid(phi, lam)
        h = np.full_like(phi, 0.0)  
        N = self.a / np.sqrt(1 - self.e2 * np.sin(phi)**2)
        X = (N + h) * np.cos(phi) * np.cos(lam) / 1e3
        Y = (N + h) * np.cos(phi) * np.sin(lam) / 1e3
        Z = ((1 - self.e2) * N + h) * np.sin(phi) / 1e3
        return X, Y, Z

    def propagate_TLE(self, tle_line1, tle_line2, minutes=3600, step=0.5):
        satellite = Satrec.twoline2rv(tle_line1, tle_line2)
        jd, fr = jday(2025, 10, 4, 0, 0, 0.0)  
        times = np.arange(0, minutes, step)  
        sat_positions = []
        for t in times:
            e, r, v = satellite.sgp4(jd, fr + t/1440.0)
            if e != 0:
                raise RuntimeError(f"SGP4 error code {e}")
            sat_positions.append(r)
        return np.array(sat_positions)
    
    def keplerian_to_RV(self, body, orbital_dict):
        epoch = Time(orbital_dict["orbit_determination_date"], scale="tdb")
        ecc   = float(orbital_dict["eccentricity"])
        sma   = float(orbital_dict["semi_major_axis"]) 
        inc   = float(orbital_dict["inclination"])
        LAN   = float(orbital_dict["ascending_node_longitude"])
        argp  = float(orbital_dict["perihelion_argument"])
        M0    = float(orbital_dict["mean_anomaly"])
        orb = Orbit.from_classical(
            body,
            sma * u.AU,           
            ecc * u.one,
            inc * u.deg,
            LAN * u.deg,
            argp * u.deg,
            M0 * u.deg,
            epoch=epoch
        )
        r, v = orb.rv()
        return epoch, r.to_value(u.km), v.to_value(u.km / u.s)
    
    def RV_to_keplerian(self, body, epoch, r, v):
        orb = Orbit.from_vectors(body, r * u.km, v * u.km/u.s, epoch=epoch)
        a    = orb.a.to(u.AU)
        ecc  = orb.ecc
        inc  = orb.inc.to(u.deg)
        raan = orb.raan.to(u.deg)
        argp = orb.argp.to(u.deg)
        nu   = orb.nu.to(u.deg)
        M    = orb.M.to(u.deg)
        n    = orb.n.to(u.deg/u.day)
        T    = orb.period.to(u.day)

        orbit_dict = {
            "epoch": epoch.iso,   # ISO string
            "eccentricity": float(ecc),
            "semi_major_axis": float(a.value),
            "inclination": float(inc.value),
            "ascending_node_longitude": float(raan.value),
            "argument_of_periapsis": float(argp.value),
            "true_anomaly": float(nu.value),
            "mean_anomaly": float(M.value),
            "mean_motion": float(n.value),
            "orbital_period": float(T.value)
        }

        return orbit_dict
    
    def relative_reference_frame(self, body1, body2, epoch):

        eph1 = Ephem.from_body(body1, epoch)
        
        eph2 = Ephem.from_body(body2, epoch)
        
        orb1 = Orbit.from_ephem(Sun, eph1, epoch)
        
        orb2 = Orbit.from_ephem(Sun, eph2, epoch)
        
        r1, v1 = orb1.rv()
        
        r2, v2 = orb2.rv()
        
        r_rel = (r1 - r2).to(u.km).value
        
        v_rel = (v1 - v2).to(u.km/u.s).value
        
        return r_rel, v_rel
    
    def change_reference_frame(self, rel_state_vector, state_vector2):
        
        r_rel, v_rel = rel_state_vector
        
        r2, v2 = state_vector2
        
        r1 = -np.array(r_rel) + np.array(r2)
        
        v1 = -np.array(v_rel) + np.array(v2)
        
        return r1, v1
    

# ------------------------------
# Usage
# ------------------------------
class1 = Rentry_analysis()

path = "C:\\Users\\Claudio Manuel\\Desktop\\prog\\linguagens\\python\\satellite\\NASA\\meatier-not-meteor\\simulations\\orbits.json"

X, Y, Z = class1.earth_WGS84_ECEF()
orbit_dict, dia_dict = class1.readData_from_json_NEO(path)
epoch, rs, vs = class1.keplerian_to_RV(Sun, orbit_dict)
r_rel, v_rel = class1.relative_reference_frame(Sun, Earth, epoch)
re, ve = class1.change_reference_frame((r_rel, v_rel), (rs, vs))

print("Transformed state vector wrt Earth:")
print("r [km] =", re)
print("v [km/s] =", ve)

# Example ISS TLE (dummy data)
tle1 = "1 XXXXU XXXXX   2461000.5  0.00016912  00000-0  30874-3 0  9992"
tle2 = "2 XXXXX  12.58732  306.50 067  195.6472  194.37 1.7496"
sat_pos = class1.propagate_TLE(tle1, tle2, minutes=60, step=1)

fig, ax = plt.subplots(subplot_kw={"projection": "3d"}, figsize=(8,8))
ax.plot_surface(X, Y, Z, cmap=cm.Blues, alpha=0.6)
ax.plot(sat_pos[:,0], sat_pos[:,1], sat_pos[:,2], 'r', label="Orbit")
ax.set_xlabel("X [km]"); ax.set_ylabel("Y [km]"); ax.set_zlabel("Z [km]")
ax.set_title("SGP4 Orbit around WGS-84 Earth")
ax.legend()
plt.show()