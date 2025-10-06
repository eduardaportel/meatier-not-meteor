import math
from clients.usgs_client import USGSClient
from enum import Enum

class AsteroidDensity(Enum):
    """Esta enumeração representa as estações do ano."""
    AVERAGE_ICE = 1000
    AVERAGE_POROUS_ROCK = 1500
    AVERAGE_DENSE_ROCK = 3000
    AVERAGE_IRON = 8000

class AsteroidEstimations:
    """
    Estimativas relacionadas ao impacto de asteroides.
    """
# física do impacto 
    def estimate_crater_size(self, energy_megatons: float):
        """
        Estima o diâmetro e profundidade da cratera com base na energia.
        Fórmula aproximada derivada do Earth Impact Effects Program.
        """
        diameter_km = 1.161 * (energy_megatons ** 0.294)
        depth_km = 0.2 * diameter_km
        return diameter_km, depth_km

    def estimate_earthquake_magnitude(self, energy_megatons: float):
        """
        Converte energia de impacto em magnitude sísmica equivalente (Mw).
        """
        energy_joules = energy_megatons * 4.184e15
        seismic_energy = energy_joules * 10**-4
        magnitude = (2/3) * (math.log10(seismic_energy) - 4.8)
        return magnitude

    def estimate_impact_energy(self, diameter_m: float, velocity_kms, density: float = float(AsteroidDensity.AVERAGE_DENSE_ROCK.value)):
        """
        Calcula a energia de impacto de um asteroide (em megatons de TNT).
        - diameter_m: diâmetro médio do asteroide [m]
        - velocity_kms: velocidade relativa [km/s]
        - density: densidade média (rocha ~3000 kg/m³)
        """
        radius = diameter_m / 2
        volume = (4/3) * math.pi * radius**3
        mass = volume * density  # kg
        velocity_ms = velocity_kms * 1_000
        energy_joules = 0.5 * mass * velocity_ms**2
        energy_megatons = energy_joules / 4.184e15  # 1 megaton TNT = 4.184e15 J

        return {
            "mass_kg": mass, "energy_megatons": energy_megatons,"energy_joules": energy_joules
        }
    
    def estimate_crater_area(self, energy_megatons: float):
        """
        Estima o diâmetro e profundidade da cratera com base na energia.
        Fórmula aproximada derivada do Earth Impact Effects Program.
        """
        diameter_km = 1.161 * (energy_megatons ** 0.294)
        depth_km = 0.2 * diameter_km
        cratera_area_km2 = diameter_km * depth_km
        return cratera_area_km2

# checar oceano ou terra
    def is_ocean_or_land(self, lat: float, lon: float):
        """
        Uses USGS elevation data to determine whether the point is on land or ocean.
        Returns:
            (is_ocean: bool, elevation_m: float)
        """
        usgs = USGSClient()
        elevation_m = usgs.get_elevation(lat, lon)
        return (True, elevation_m) if elevation_m < 0 else (False, elevation_m)

#tsunami
    def estimate_transient_crater(self, diameter_m, velocity_ms, density):
        """Modelo baseado em Ward & Asphaug para cavidade transitória."""
        transient_crater_diameter_m = 1.161 * (diameter_m ** 0.78) * (velocity_ms ** 0.44) * (density / 1000) ** 0.26 * 9.81 ** (-0.22)
        transient_crater_radius_m = transient_crater_diameter_m / 2
        average_crater_depth_m = 0.25 * transient_crater_diameter_m
        return transient_crater_diameter_m, transient_crater_radius_m, average_crater_depth_m

    def estimate_initial_tsunami_height(self, average_crater_depth_m):
        """Altura inicial da onda logo após o impacto."""
        initial_height_m = 0.25 * average_crater_depth_m
        return initial_height_m

    def variable_tsunami_height(self, initial_height_m, transient_crater_radius_m, distance_km):
        """Altura da tsunami em função da distância (km)."""
        tsunami_height_m = initial_height_m * (transient_crater_radius_m / (distance_km * 1000)) ** 0.5
        return tsunami_height_m
    
    def shore_tsunami_height(self, tsunami_height_m, average_oceanic_depth, shore_depth):
        """Altura da tsunami na costa."""
        shore_tsunami_height_m = tsunami_height_m * (average_oceanic_depth / shore_depth) ** 0.25
        return shore_tsunami_height_m
