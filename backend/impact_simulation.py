#codigo antigo meu
import os
import sys

# Ensure local package imports work when running this file directly
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
if CURRENT_DIR not in sys.path:
    sys.path.insert(0, CURRENT_DIR)

from core.asteroid_estimations import AsteroidEstimations
from clients.nasa_client import NASAClient

def run_impact_simulation(neo_id: str, lat: float, lon: float, distance_km: float = 500):
    nasa = NASAClient()
    asteroid = nasa.get_neo_lookup(neo_id)
    estimations = AsteroidEstimations()

# parâmetros principais
    diameter = (
            float(asteroid["estimated_diameter"]["meters"]["estimated_diameter_min"]) +
            float(asteroid["estimated_diameter"]["meters"]["estimated_diameter_max"])
        ) / 2
    velocity_kms = float(asteroid["close_approach_data"][0]["relative_velocity"]["kilometers_per_second"])
    velocity_ms = velocity_kms * 1000

# impact = estimate_impact_energy(asteroid["diameter_m"], asteroid["velocity_kms"])
    impact = estimations.estimate_impact_energy(diameter_m=diameter, velocity_kms=velocity_kms)

    # localização: Terra ou Oceano
    try:
        is_ocean, elevation = estimations.is_ocean_or_land(lat, lon)
        surface_type = "Ocean" if is_ocean else "Land"
    except Exception as e:
        surface_type = f"Unknown ({e})"
        elevation = None
        is_ocean = False

   # cratera e magnitude sísmica 
    crater_area_km2 = estimations.estimate_crater_size(impact["energy_megatons"])
    magnitude = estimations.estimate_earthquake_magnitude(impact["energy_megatons"])
 
    # tsunami (apenas se oceânico)
    tsunami = None
    if is_ocean:
        transient_d, transient_r, avg_depth = estimations.estimate_transient_crater(diameter, velocity_ms, 3000)
        initial_wave = estimations.estimate_initial_tsunami_height(avg_depth)
        wave_at_distance = estimations.variable_tsunami_height(initial_wave, transient_r, distance_km)
        wave_at_shore = estimations.shore_tsunami_height(
            wave_at_distance, average_oceanic_depth=4000, shore_depth=200
        )

        tsunami = {
            "transient_crater_diameter_m": transient_d,
            "initial_height_m": initial_wave,
            "height_500km_m": wave_at_distance,
            "shore_height_m": wave_at_shore
        }

    return {
        "asteroid_name": asteroid.get("name"),
        "neo_id": neo_id,
        "diameter_m": diameter,
        "velocity_kms": velocity_kms,
        "energy_megatons": impact["energy_megatons"],
        "surface_type": surface_type,
        "elevation_m": elevation,
        "crater_area_km2": crater_area_km2,
        "earthquake_magnitude": magnitude,
        "tsunami": tsunami
    }

def main():
    lat, lon = -20.0, -30.0  # Exemplo no Atlântico Sul
    results = run_impact_simulation("3542519", lat, lon)
    

    print(f"\nSimulação — {results['asteroid_name']} ({results['neo_id']})")
    print(f"Localização: {lat:.2f}, {lon:.2f}")
    print(f"Diâmetro: {results['diameter_m']:.1f} m | Velocidade: {results['velocity_kms']:.2f} km/s")
    print(f"Energia: {results['energy_megatons']:.2f} Mt TNT")
    print(f"Cratera: {results['crater_area_km2']:.2f} km²")
    print(f"Terremoto equivalente: Mw {results['earthquake_magnitude']:.2f}")

    if results["tsunami"]:
        print("\n Tsunami estimada:")
        print(f"- Cavidade transitória: {results['tsunami']['transient_crater_diameter_m']:.1f} m")
        print(f"- Altura inicial: {results['tsunami']['initial_height_m']:.2f} m")
        print(f"- Altura a 500 km: {results['tsunami']['height_500km_m']:.2f} m")
        print(f"- Altura na costa: {results['tsunami']['shore_height_m']:.2f} m")
    else:
        print("\nSem geração significativa de tsunami (impacto terrestre).")

if __name__ == "__main__":
    main()
