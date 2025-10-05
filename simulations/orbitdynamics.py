import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from sgp4.api import Satrec, jday
from astropy import units as u
from scipy.integrate import solve_ivp
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
        self.Mu = 398600.4418  # [km^3/s^2]
        self.H = 8500
        self.rho0 = 1.225
        self.rho_body = 3000
        self.diameter = 0
        self.mass = 0
        self.volume = 0

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

    def propagate_TLE(self, minutes, orbit_dict, step):
       
        sat_positions = []
        sat_velocity = []
        altitudes = []
        
        reentry_flag = False

        times = np.arange(0, minutes, step)

        diameter = self.diameter

        for t in times:
            tle_line1, tle_line2 = self.write_tle_from_elements(99999, orbit_dict, diameter)

            # Get the epoch string from TLE line 1 (positions 18–32)
            tle_epoch_str = tle_line1[18:32].strip()  
            
            year_short = int(tle_epoch_str[:2])
            
            doy = float(tle_epoch_str[2:])
            
            year_full = 2000 + year_short if year_short < 57 else 1900 + year_short
            
            epoch_datetime = datetime(year_full, 1, 1) + timedelta(days=doy - 1)

            jd, fr = jday(epoch_datetime.year, epoch_datetime.month, epoch_datetime.day,
                          epoch_datetime.hour, epoch_datetime.minute, epoch_datetime.second)

            satellite = Satrec.twoline2rv(tle_line1, tle_line2)
            
            e, r, v = satellite.sgp4(jd, fr + t / 1440.0)

            r_norm = np.linalg.norm(r)
            
            print(f"Altitude {r_norm - 6378}")

            if r_norm - 6378 > diameter/2 and reentry_flag == True:

                r_reentry, v_reentry = self.reentry_calculation(r, v, diameter, times)
        
                sat_positions.append(r_reentry.flatten())
        
                sat_velocity.append(v_reentry.flatten())
                break

            if e == 0 and r_norm - 6378 > 250:
                sat_positions.append(r)
                
                sat_velocity.append(v)
                
                altitudes.append(r_norm - 6378)
                
                reentry_flag = False
            else:
                print(f"SGP4 error code {e} at t = {t} min")
                print("Switching to reentry integration (solve_ivp)...")

                # Perform reentry trajectory integration using IVP
                t_reentry, y_reentry = self.reentry_calculation(r, v, t_span=(t, minutes))

                for x, y, z, vx, vy, vz in y_reentry:
                    sat_positions.append([x, y, z])
                    sat_velocity.append([vx, vy, vz])
                    altitudes.append(np.linalg.norm([x, y, z]) - 6378)               
                break

        return np.array(sat_positions), np.array(sat_velocity), altitudes
                
    def reentry_ivp(self, t, y, diameter):
       
        r = y[:3]   # km
        v = y[3:]   # km/s

        r = np.array(r)
        v = np.array(v)
        r_norm = np.linalg.norm(r)

        altitude = r_norm - 6378
        a_g = -self.Mu * r / (r_norm**3)

        Bstar_drag, Ballistic_coef = self.atmospheric_geo_drag(
            diameter=diameter, rho_body=self.rho_body, altitude=altitude
        )

        wEarth = np.array([0, 0, 7.2921159*1e-5])
        vrel = np.array(v - np.cross(wEarth, r))
        vrel_norm = np.linalg.norm(vrel)

        a_drag = -0.5 * Bstar_drag * vrel_norm * vrel      
        a_reentry = a_g + a_drag
        v_reentry = v

        return np.hstack((v_reentry, a_reentry))

    def reentry_trigger(self, t, y):
            r_norm = np.linalg.norm(y[:3])
            return r_norm - 6378 - self.diameter/2
    
    def get_reentry_event(self):
        
        def event(t, y):
            return self.reentry_trigger(t, y)

        event.terminal = True

        event.direction = -1

        return event

    def reentry_calculation(self, r0, v0, t_span):
        
        t_span = (t_span[0], t_span[-1])

        y0 = np.hstack((r0, v0))

        sol = solve_ivp(
            fun=lambda t, y: self.reentry_ivp(t, y, self.diameter),
            t_span=t_span,
            y0=y0,
            method="RK45",
            max_step=10,
            events=self.get_reentry_event(),
            rtol=1e-4,
            atol=1e-4
        )
        
        return sol.t, sol.y.T
    
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
            "epoch": epoch.iso,       
            "eccentricity": ecc,
            "semi_major_axis": a,     
            "inclination": inc,       
            "ascending_node_longitude": raan,  
            "argument_of_periapsis": argp,     
            "true_anomaly": np.degrees(nu),    
            "mean_anomaly": M,                 
            "mean_motion": n,                  
            "orbital_period": T                
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
    
    def write_tle_from_elements(self, satnum, elements, diameter):
        
        mean_motion_revs_per_day = elements["mean_motion"] / 360.0
    
        epoch_time = Time(elements["epoch"])
    
        year = epoch_time.datetime.year % 100
    
        day_of_year = epoch_time.datetime.timetuple().tm_yday
    
        frac_of_day = (epoch_time.datetime.hour * 3600 +
                       epoch_time.datetime.minute * 60 +
                       epoch_time.datetime.second +
                       epoch_time.datetime.microsecond / 1e6) / 86400.0

        epoch_str = f"{year:02d}{day_of_year + frac_of_day:012.8f}"

        altitude_km = elements['semi_major_axis'] - 6378.0  

        Bstar, Ballistic = self.atmospheric_geo_drag(diameter, self.rho_body, altitude_km)

        bstar_str = f"{Bstar:.5e}".replace("e", "").replace("+", "").zfill(8)

        bstar_tle = f"{bstar_str[0:5]}{bstar_str[5]}{bstar_str[6]}"

        ecc_str = f"{elements['eccentricity']:.7f}".split(".")[1][:7]

        tle_line1 = (
            f"1 {satnum:05d}U 25001A   {epoch_str}  "
            f".00000000  00000-0 {bstar_tle} 0  9991"
        )

        tle_line2 = (
            f"2 {satnum:05d} "
            f"{elements['inclination']:8.4f} "
            f"{elements['ascending_node_longitude']:8.4f} "
            f"{ecc_str:7s} "
            f"{elements['argument_of_periapsis']:8.4f} "
            f"{elements['mean_anomaly']:8.4f} "
            f"{mean_motion_revs_per_day:11.8f}00001"
        )

        return tle_line1, tle_line2

    def atmospheric_geo_drag(self, diameter, rho_body, altitude):
        
        radius = diameter/2

        Volume = 4*PI/3 * radius**3

        Area = PI*radius**2

        mass = rho_body * Volume

        self.mass = mass
        self.volume = Volume

        Cd = 0.7

        Balllistic_coef = mass / (Cd*Area)

        rho = self.rho0 * np.exp(-altitude/self.H)

        Bstar_drag = rho * Cd * Area / mass

        return Bstar_drag, Balllistic_coef

# ------------------------------
# Usage
# ------------------------------
class1 = Rentry_analysis()

path = "C:\\Users\\Claudio Manuel\\Desktop\\prog\\linguagens\\python\\satellite\\NASA\\meatier-not-meteor\\simulations\\orbits.json"

X, Y, Z = class1.earth_WGS84_ECEF()
orbit_dict, dia_dict = class1.readData_from_json_NEO(path)

class1.diameter = dia_dict["estimated_diameter_max"]

epoch, rs, vs = class1.keplerian_to_RV(Sun, orbit_dict)
r_rel, v_rel = class1.relative_reference_frame(Sun, Earth, epoch)
re, ve = class1.change_reference_frame((r_rel, v_rel), (rs, vs))

print("Transformed state vector wrt Earth:")
print("r [km] =", re)
print("v [km/s] =", ve)

orbit_dict_earth = class1.RV_to_keplerian(Earth, epoch, re, ve)

tle1 = "1 13343U 82092A   82348.50000000  .00035000  00000-0  12000-3 0  9991"
tle2 = "2 13343  65.0000 150.0000 0005000   0.0000  90.0000 16.00000000    01"

satellite = Satrec.twoline2rv(tle1, tle2)
tle_epoch_str = tle1[18:32].strip()
year_short = int(tle_epoch_str[:2])
doy = float(tle_epoch_str[2:])
year_full = 2000 + year_short if year_short < 57 else 1900 + year_short
epoch_datetime = datetime(year_full, 1, 1) + timedelta(days=doy - 1)

jd, fr = jday(epoch_datetime.year, epoch_datetime.month, epoch_datetime.day,
              epoch_datetime.hour, epoch_datetime.minute, epoch_datetime.second)

e, r, v = satellite.sgp4(jd, fr)
orbit_dict_earth = class1.RV_to_keplerian(Earth, epoch, r, v)

step = 0.1
duration = 63933.0 * 10
sat_pos, sat_vel, altitudes = class1.propagate_TLE(minutes=duration, 
                                            orbit_dict=orbit_dict_earth, 
                                            step=step)

'''
fig, ax = plt.subplots(subplot_kw={"projection": "3d"}, figsize=(8,8))
ax.plot_surface(X, Y, Z, rstride=5, cstride=5, color='blue', alpha=0.4, linewidth=0)
ax.scatter(sat_pos[0,0], sat_pos[0,1], sat_pos[0,2], 'g', label="Starting point")
ax.scatter(sat_pos[-1,0], sat_pos[-1,1], sat_pos[-1,2], 'o', label="Ending point")
ax.plot(sat_pos[:,0], sat_pos[:,1], sat_pos[:,2], 'r', label="Orbit")
ax.set_xlabel("X [m]"); ax.set_ylabel("Y [m]"); ax.set_zlabel("Z [m]")
ax.set_title("SGP4 Orbit around WGS-84 Earth")
ax.legend()
plt.show()
'''

plt.figure(figsize=(10, 5))
times = np.arange(0, len(altitudes))
plt.plot(times, altitudes)
plt.axhline(80, color='red', linestyle='--', label="Reentry limit")
plt.title(f"Altitude decay over {times[-1]/43800} months")
plt.xlabel("Time [min]")
plt.ylabel("Altitude [km]")
plt.grid(True)
plt.legend()
plt.show()

plt.figure(figsize=(10, 5))

times = np.arange(0, len(altitudes))

plt.plot( np.arange(0, len(sat_vel)) , class1.mass*np.linalg.norm(sat_vel)**2/2)
plt.axhline(80, color='red', linestyle='--', label="Reentry limit")
plt.title(f"Altitude decay over {times[-1]/43800} months")
plt.xlabel("Time [min]")
plt.ylabel("Altitude [km]")
plt.grid(True)
plt.legend()
plt.show()