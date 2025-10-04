# core/asteroid_data.py
# import requests

# def get_asteroid_data(neo_id: str, api_key="DEMO_KEY"):
#     """
#     Obtém dados de um asteroide específico pela API da NASA NEO.
#     """
#     url = f"https://api.nasa.gov/neo/rest/v1/neo/{neo_id}?api_key={api_key}"
#     response = requests.get(url)
#     if response.status_code != 200:
#         raise ValueError("Erro ao acessar NASA NEO API.")
    
#     data = response.json()
#     return {
#         "name": data["name"],
#         "diameter_m": (
#             data["estimated_diameter"]["meters"]["estimated_diameter_min"] +
#             data["estimated_diameter"]["meters"]["estimated_diameter_max"]
#         ) / 2,
#         "velocity_kms": float(data["close_approach_data"][0]["relative_velocity"]["kilometers_per_second"]),
#         "mass_kg": None,  # será calculado depois
#         "orbit": data["orbital_data"]
#     }
