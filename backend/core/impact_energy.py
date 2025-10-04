# core/impact_energy.py
# import math

# def estimate_impact_energy(diameter_m, velocity_kms, density=3000):
#     """
#     Calcula a energia de impacto de um asteroide (em megatons de TNT).
#     - diameter_m: diâmetro médio do asteroide [m]
#     - velocity_kms: velocidade relativa [km/s]
#     - density: densidade média (rocha ~3000 kg/m³)
#     """
#     radius = diameter_m / 2
#     volume = (4/3) * math.pi * radius**3
#     mass = volume * density  # kg
#     velocity_ms = velocity_kms * 1_000
#     energy_joules = 0.5 * mass * velocity_ms**2
#     energy_megatons = energy_joules / 4.184e15  # 1 megaton TNT = 4.184e15 J

#     return {
#         "mass_kg": mass,
#         "energy_megatons": energy_megatons,
#         "energy_joules": energy_joules
#     }
