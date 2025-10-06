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

class Reentry_analysis:
    """
    Class to simulate reentry trajectories of NEOs or satellites
    using simplified SGP4 propagation + atmospheric drag modeling.
    """
    
    def __init__(self):
        """
        Initialize Earth/WGS-84 parameters and default object properties.
        """
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
        """
        Load NEO data from JSON file.
        Returns: name, orbital elements dict, diameter dict.
        """
        with open(path, "r") as file:
       
            data = json.load(file)  
       
        orbital_dict = data["orbital_data"]
       
        diameter_dict = data["estimated_diameter"]["meters"]
       
        name = data["name"]

        return name, orbital_dict, diameter_dict
    
    def earth_WGS84_ECEF(self):
        """
        Compute Earth ellipsoid mesh (ECEF coordinates in km)
        for plotting purposes using WGS-84.
        Returns: X, Y, Z meshgrids.
        """
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
        """
        Propagate satellite orbit using TLE format + SGP4.
        Stops when altitude < 80 km (reentry).
        Args:
            minutes: total propagation duration [min]
            orbit_dict: orbital elements
            step: propagation step size [min]
        Returns: positions [km], velocities [km/s], altitudes [km]
        """
        sat_positions = []
        sat_velocity = []
        altitudes = []
        
        reentry_flag = False

        times = np.arange(0, minutes, step)

        diameter = self.diameter

        velocity = 0.0

        for t in times:
            tle_line1, tle_line2 = self.write_tle_from_elements(99999, orbit_dict, diameter, velocity)

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

            if e == 0 and r_norm - 6378 > 80:
                sat_positions.append(r)
                
                sat_velocity.append(v)
                
                altitudes.append(r_norm - 6378)
                
                velocity = np.linalg.norm(v)
            else:
                print(f"SGP4 error code {e} at t = {t} min")
                '''
                print("Switching to reentry integration (solve_ivp)...")

                # Perform reentry trajectory integration using IVP
                t_reentry, y_reentry = self.reentry_calculation(r, v)

                for x, y, z, vx, vy, vz in y_reentry:
                    altitude = (np.linalg.norm([x, y, z]) - 6378) >= 0
                    if altitude <= 0:
                        break  
                    sat_positions.append([x, y, z])
                    sat_velocity.append([vx, vy, vz])
                    altitudes.append(altitude)
                ''' 
                break

        altitudes = np.array(altitudes)
        sat_positions = np.array(sat_positions)
        sat_velocity = np.array(sat_velocity)

        return sat_positions, sat_velocity, altitudes
                
    
    def keplerian_to_RV(self, body, orbital_dict):
        """
        Convert classical Keplerian orbital elements → position/velocity vectors.
        Args:
            body: reference central body (Earth or Sun)
            orbital_dict: dict with orbital elements
        Returns: epoch [Time], r [km], v [km/s]
        """
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
        """
        Convert position/velocity vectors → Keplerian elements.
        Args:
            body: central body
            epoch: epoch of state
            r: position vector [km]
            v: velocity vector [km/s]
        Returns: orbital elements dict
        """

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
        """
        Compute relative state (r,v) of body1 w.r.t. body2 at given epoch.
        Useful for switching from heliocentric to geocentric frame.
        Returns: r_rel [km], v_rel [km/s]
        """
        
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
        """
        Change reference frame using relative and absolute vectors.
        Args:
            rel_state_vector: (r_rel, v_rel)
            state_vector2: (r2, v2)
        Returns: new (r1, v1) in km and km/s
        """
        
        r_rel, v_rel = rel_state_vector
    
        r2, v2 = state_vector2
    
        r1 = -np.array(r_rel) + np.array(r2)
    
        v1 = -np.array(v_rel) + np.array(v2)
    
        return r1, v1
    
    def write_tle_from_elements(self, satnum, elements, diameter, velocity):
        """
        Construct TLE lines from orbital elements + drag term (Bstar).
        Args:
            satnum: satellite number
            elements: orbital dict
            diameter: body diameter [m]
            velocity: body velocity magnitude [km/s]
        Returns: tle_line1, tle_line2
        """
        
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

        Bstar, Ballistic = self.atmospheric_geo_drag(diameter, self.rho_body, altitude_km, velocity)

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

    def atmospheric_geo_drag(self, diameter, rho_body, altitude, velocity):
        """
        Estimate ballistic coefficient and drag factor Bstar
        using exponential atmosphere model.
        Args:
            diameter: body diameter [m]
            rho_body: density of body [kg/m^3]
            altitude: altitude above Earth [km]
            velocity: velocity magnitude [km/s]
        Returns: Bstar_drag, Ballistic_coef
        """
        
        radius = diameter/2

        volume = 4*PI/3 * radius**3

        Area = PI*radius**2

        mass = rho_body * volume

        self.mass = mass
        self.volume = volume

        M = 0

        if velocity == 0:
            Cd = 0.7
        else:
            M = self.compute_mach(altitude, velocity)

            if M <= 0.722:

                Cd = 0.45*M**2 + 0.424
            
            else:
                Cd = 2.1*np.exp(-1.2*(M+0.35)) - 8.9*np.exp(-2.2*(M+0.35)) + 0.92

        Balllistic_coef = mass / (Cd*Area)

        rho = self.rho0 * np.exp(-altitude/self.H)

        Bstar_drag = rho * Cd * Area / mass

        return Bstar_drag, Balllistic_coef
    
    def compute_mach(self, altitude_km, v):
        """
        Compute Mach number at given altitude using simplified
        piecewise sound speed model.
        Args:
            altitude_km: altitude [km]
            v: velocity [km/s]
        Returns: Mach number
        """

        # Example: simplified lookup for speed of sound (m/s)
        # Replace with interpolation from table

        if altitude_km < 11:
            a = 340.0   # sea level
        elif altitude_km < 20:
            a = 295.0
        elif altitude_km < 32:
            a = 303.0
        elif altitude_km < 47:
            a = 330.0
        elif altitude_km < 51:
            a = 340.0
        elif altitude_km < 71:
            a = 355.0
        else:
            a = 270.0   # ~80 km

        M = v * 1000.0
        return M / a
    
    def ecef_to_geodetic(self, x, y, z, a=6378.137, b=6356.7523142):
        """
        Convert ECEF (km) to geodetic latitude, longitude, altitude (WGS-84).
        Returns: lat [deg], lon [deg], alt [km]
        """

        f = (a - b) / a
        e2 = f * (2 - f)

        lon = np.arctan2(y, x)
        r = np.sqrt(x**2 + y**2)

        # initial guess for latitude
        lat = np.arctan2(z, r * (1 - e2))
        alt = 0.0

        # iterate to improve latitude/altitude
        for _ in range(5):
            N = a / np.sqrt(1 - e2 * np.sin(lat)**2)
            alt = r / np.cos(lat) - N
            lat = np.arctan2(z + e2 * N * np.sin(lat), r)

        lat = np.degrees(lat)
        lon = np.degrees(lon)

        return lat, lon, alt



    '''
    def reentry_calculation(self, r0, v0):

        t_span = (0, 2)
        
        y0 = np.hstack((r0, v0))

        sol = solve_ivp(
            fun=lambda t, y: self.reentry_ivp(t, y, self.diameter),
            t_span=t_span,
            y0=y0,
            method="RK45",
            first_step=0.01,
            max_step=0.1,
            events=self.get_reentry_event(),
            rtol=1e-4,
            atol=1e-4
        )
        
        return sol.t, sol.y.T

    def reentry_ivp(self, t, y, diameter):
       
        r = y[:3]   # km
        v = y[3:]   # km/s

        r = np.array(r)
        v = np.array(v)
        r_norm = np.linalg.norm(r)

        altitude = r_norm - 6378
        a_g = -self.Mu * r / (r_norm**3)

        Bstar_drag, Ballistic_coef = self.atmospheric_geo_drag(
            diameter=diameter,
            rho_body=self.rho_body,
            altitude=altitude,
            velocity=np.linalg.norm(v)   # km/s
        )

        wEarth = np.array([0, 0, 7.2921159*1e-5])
        vrel = np.array(v - np.cross(wEarth, r))
        vrel_norm = np.linalg.norm(vrel)

        a_drag = -0.5 * Bstar_drag * vrel_norm * vrel      
        a_reentry = a_g + a_drag
        v_reentry = v

        return np.hstack((v_reentry, a_reentry))

    def reentry_trigger(self, t, y):

        r_norm = np.linalg.norm(y[:3])  # km

        altitude = r_norm - 6378.0  # km above Earth's mean radius

        if altitude > 0:

            return altitude   # stop at 0 km altitude
        
        else:
            return -1
    
    def get_reentry_event(self):
        
        def event(t, y):
            return self.reentry_trigger(t, y)

        event.terminal = True

        event.direction = -1

        return event
    '''

if __name__ == "__main__":
    
    EnvSimulation = Reentry_analysis()

    path = ".\\orbits.json"

    name, orbit_dict, dia_dict = EnvSimulation.readData_from_json_NEO(path)

    EnvSimulation.diameter = dia_dict["estimated_diameter_max"]

    epoch, rs, vs = EnvSimulation.keplerian_to_RV(Sun, orbit_dict)

    r_rel, v_rel = EnvSimulation.relative_reference_frame(Sun, Earth, epoch)

    re, ve = EnvSimulation.change_reference_frame((r_rel, v_rel), (rs, vs))

    orbit_dict_earth = EnvSimulation.RV_to_keplerian(Earth, epoch, re, ve)

    step = 1000
    duration = 52560000
    sat_pos, sat_vel, altitudes = EnvSimulation.propagate_TLE(minutes=duration, 
                                                orbit_dict=orbit_dict_earth, 
                                                step=step)

    lat, lon, alt = EnvSimulation.ecef_to_geodetic(sat_pos[-1, 0], sat_pos[-1, 1], sat_pos[-1, 2])

    vx, vy, vz = sat_vel[-1, 0], sat_vel[-1, 1], sat_vel[-1, 2]
    print(f"Final positon: lat, lon, alt: {[lat, lon, alt]}")
    print(f"Final velocity: vx, vy, vz: {[vx, vy, vz]}")

    results = {
        'name':name,
        '2Dresults' : {
            'lat': lat,
            'lon': lon,
            'alt': alt,
            'vx': vx,
            'vy': vy,
            'vz': vz,
            'v_norm': np.linalg.norm( [vx, vy, vz] ),
            'mass': EnvSimulation.mass,
            'volume': EnvSimulation.volume,
            'diameter': dia_dict["estimated_diameter_max"]
        },
        '3Dresults' : {
            'rs_x': rs[0],
            'rs_y': rs[1],
            'rs_z': rs[2],
            'vs_x': vs[0],
            'vs_y': vs[1],
            'vs_z': vs[2]
        } 
    }
    # Write JSON file
    with open("results_2D_and_3D.json", "w") as f:
        json.dump(results, f, indent=4)