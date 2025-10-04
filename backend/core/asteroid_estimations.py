# core/crater_model.py
import math

class AsteroidEstimations:
    """
    Estimativas relacionadas ao impacto de asteroides.
    """
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

    def estimate_impact_energy(self, diameter_m, velocity_kms, density=3000):
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
            "mass_kg": mass,
            "energy_megatons": energy_megatons,
            "energy_joules": energy_joules
        }

    def estimate_tsunami_height(self, energy_megatons, distance_km):
        """
        Modelo simplificado de altura de tsunami em função da energia e distância.
        """
        height_m = 100 * (energy_megatons ** 0.35) / (distance_km ** 0.5)
        return height_m
