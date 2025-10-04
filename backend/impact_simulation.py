# from core.asteroid_data import get_asteroid_data
# from core.impact_energy import estimate_impact_energy
# from core.crater_model import estimate_crater_size, estimate_earthquake_magnitude
# from core.tsunami_model import estimate_tsunami_height
from core.asteroid_estimations import AsteroidEstimations
from clients.nasa_client import NASAClient


def main():
    # asteroid = get_asteroid_data("3542519")  # (2010 PK9)
    asteroid = NASAClient().get_neo_lookup("3542519")

    # print(asteroid)
    # print(asteroid["diameter_min"])
    # print(asteroid["velocity_kms"])
    estimations = AsteroidEstimations()

    diameter = (
            int(asteroid["estimated_diameter"]["meters"]["estimated_diameter_min"]) +
            int(asteroid["estimated_diameter"]["meters"]["estimated_diameter_max"])
        ) / 2
    velocity = float(asteroid["close_approach_data"][0]["relative_velocity"]["kilometers_per_second"])

    # impact = estimate_impact_energy(asteroid["diameter_m"], asteroid["velocity_kms"])
    impact = estimations.estimate_impact_energy(diameter_m=diameter, velocity_kms=velocity)

    crater_d, crater_h = estimations.estimate_crater_size(impact["energy_megatons"])
    magnitude = estimations.estimate_earthquake_magnitude(impact["energy_megatons"])
    tsunami_height = estimations.estimate_tsunami_height(impact["energy_megatons"], distance_km=500)

    print(f"Asteroide: {asteroid['name']}")
    print(f"Diâmetro: {diameter:.1f} m")
    print(f"Velocidade: {velocity} km/s")
    print(f"Energia: {impact['energy_megatons']:.2f} Mt TNT")
    print(f"Cratera: {crater_d:.2f} km de diâmetro, {crater_h:.2f} km de profundidade")
    print(f"Terremoto equivalente: Mw {magnitude:.2f}")
    print(f"Tsunami estimado: {tsunami_height:.2f} m de altura a 500 km")

if __name__ == "__main__":
    main()
