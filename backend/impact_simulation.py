from core.asteroid_data import get_asteroid_data
from core.impact_energy import estimate_impact_energy
from core.crater_model import estimate_crater_size, estimate_earthquake_magnitude
from core.tsunami_model import estimate_tsunami_height

def main():
    asteroid = get_asteroid_data("3542519")  # (2010 PK9)
    impact = estimate_impact_energy(asteroid["diameter_m"], asteroid["velocity_kms"])
    
    crater_d, crater_h = estimate_crater_size(impact["energy_megatons"])
    magnitude = estimate_earthquake_magnitude(impact["energy_megatons"])
    tsunami_height = estimate_tsunami_height(impact["energy_megatons"], distance_km=500)
    
    print(f"Asteroide: {asteroid['name']}")
    print(f"Diâmetro: {asteroid['diameter_m']:.1f} m")
    print(f"Velocidade: {asteroid['velocity_kms']} km/s")
    print(f"Energia: {impact['energy_megatons']:.2f} Mt TNT")
    print(f"Cratera: {crater_d:.2f} km de diâmetro, {crater_h:.2f} km de profundidade")
    print(f"Terremoto equivalente: Mw {magnitude:.2f}")
    print(f"Tsunami estimado: {tsunami_height:.2f} m de altura a 500 km")

if __name__ == "__main__":
    main()
