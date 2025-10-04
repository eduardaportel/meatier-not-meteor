import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from sgp4.api import Satrec, jday
from astropy import units as u
from datetime import datetime, timedelta
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

    def propagate_TLE(self, tle_line1, tle_line2, minutes, step):

        # Get the epoch string from TLE line 1 (positions 18–32)
        tle_epoch_str = tle_line1[18:32].strip()  # e.g. "25001.00000000"
        year_short = int(tle_epoch_str[:2])
        doy = float(tle_epoch_str[2:])

        year_full = 2000 + year_short if year_short < 57 else 1900 + year_short

        epoch_datetime = datetime(year_full, 1, 1) + timedelta(days=doy - 1)
        jd, fr = jday(epoch_datetime.year, epoch_datetime.month, epoch_datetime.day,
                    epoch_datetime.hour, epoch_datetime.minute, epoch_datetime.second)

        satellite = Satrec.twoline2rv(tle_line1, tle_line2)
        
        times = np.arange(0, minutes, step)
        
        sat_positions = []
        
        altitudes = []

        for t in times:

            e, r, v = satellite.sgp4(jd, fr + t / 1440.0)

            if e == 0:
                sat_positions.append(r)
                r_norm = np.linalg.norm(r)
                altitudes.append(r_norm - 6378)
            else:   
                reentry_time = t
                print(f"SGP4 error code {e} at t = {reentry_time} min → Reentrada detectada!")
                break

        return np.array(sat_positions), altitudes
    
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

        # Extraímos os valores e convertemos para float logo aqui
        a    = orb.a.to(u.AU).value
        
        ecc  = float(orb.ecc)
        
        inc  = orb.inc.to(u.deg).value
        
        raan = orb.raan.to(u.deg).value
        
        argp = orb.argp.to(u.deg).value
        
        nu   = orb.nu.to(u.rad).value
        
        n    = orb.n.to(u.deg/u.day).value
        
        T    = orb.period.to(u.day).value

        # --- Mean anomaly ---
        if ecc < 1.0:
            E = 2 * np.arctan2(
                np.tan(nu / 2.0),
                np.sqrt((1 + ecc) / (1 - ecc))
            )
            if E < 0:
                E += 2*np.pi
            M = np.degrees(E - ecc*np.sin(E))
        else:
            H = 2 * np.arctanh(
                np.sqrt((ecc-1)/(ecc+1)) * np.tan(nu/2.0)
            )
            M = np.degrees(ecc*np.sinh(H) - H)

        orbit_dict = {
            "epoch": epoch.iso,       # string ISO
            "eccentricity": ecc,
            "semi_major_axis": a,     # AU
            "inclination": inc,       # deg
            "ascending_node_longitude": raan,  # deg
            "argument_of_periapsis": argp,     # deg
            "true_anomaly": np.degrees(nu),    # deg
            "mean_anomaly": M,                 # deg
            "mean_motion": n,                  # deg/day
            "orbital_period": T                # days
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
        
    def write_tle_from_elements(self, satnum, elements):
        
        mean_motion_revs_per_day = elements["mean_motion"] / 360.0  # deg/day -> rev/day

        # Build epoch string in YYDDD.DDDDD format
        epoch_time = Time(elements["epoch"])
        
        year = epoch_time.datetime.year % 100
        
        day_of_year = epoch_time.datetime.timetuple().tm_yday
        
        frac_of_day = (epoch_time.datetime.hour*3600 + 
                    epoch_time.datetime.minute*60 + 
                    epoch_time.datetime.second) / 86400.0
        
        epoch_str = f"{year:02d}{day_of_year + frac_of_day:012.8f}"

        # Line 1 (dummy drag terms etc.)
        tle_line1 = f"1 {satnum:05d}U 25001A   {epoch_str}  .00000000  00000-0  00000-0 0  9991"

        # Line 2 (eccentricity is without decimal point in TLE!)
        ecc_str = f"{elements['eccentricity']:.7f}".split(".")[1][:7]  

        tle_line2 = (f"2 {satnum:05d} "
                    f"{elements['inclination']:8.4f} "
                    f"{elements['ascending_node_longitude']:8.4f} "
                    f"{ecc_str:7s} "
                    f"{elements['argument_of_periapsis']:8.4f} "
                    f"{elements['mean_anomaly']:8.4f} "
                    f"{mean_motion_revs_per_day:11.8f}00001")

        return tle_line1, tle_line2


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

#re = np.array([2.76118303e+04, -1.35428444e+04, -1.29941920e+04])
#ve = np.array([1.59638675, 2.08307728, 2.46499025])

orbit_dict_earth = class1.RV_to_keplerian(Earth, epoch, re, ve)

#tle1 = "1 99999U 25001A   25278.00000000  .00500000  00000-0  50000-3 0  9991"
#tle2 = "2 99999  97.5000  0.0000 0001000  0.0000  0.0000 15.22000000    01"

tle1, tle2 = class1.write_tle_from_elements(99999, orbit_dict_earth)

step=100
duration = 525599.42184*5
sat_pos, altitudes = class1.propagate_TLE(tle1, tle2, minutes=duration, step=step)

fig, ax = plt.subplots(subplot_kw={"projection": "3d"}, figsize=(8,8))

ax.scatter(X, Y, Z, color="blue", s=200, label="Earth")

ax.plot(sat_pos[:,0], sat_pos[:,1], sat_pos[:,2], 'r', label="Orbit")
ax.set_xlabel("X [m]"); ax.set_ylabel("Y [m]"); ax.set_zlabel("Z [m]")
ax.set_title("SGP4 Orbit around WGS-84 Earth")
ax.legend()
plt.show()

plt.figure(figsize=(10, 5))

times = np.arange(0, len(altitudes))

plt.plot(times, altitudes)
plt.axhline(80, color='red', linestyle='--', label="Reentry limit")
plt.title("Altitude decay over 5 months")
plt.xlabel("Time [min]")
plt.ylabel("Altitude [km]")
plt.grid(True)
plt.legend()
plt.show()