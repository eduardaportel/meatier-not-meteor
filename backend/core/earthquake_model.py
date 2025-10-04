# core/earthquake_model.py
import requests

def get_recent_earthquakes(region="world", min_magnitude=5):
    """
    Obtém dados do catálogo NEIC (USGS) para comparação.
    """
    url = f"https://earthquake.usgs.gov/fdsnws/event/1/query?format=geojson&minmagnitude={min_magnitude}"
    response = requests.get(url)
    if response.status_code != 200:
        raise ValueError("Erro ao consultar USGS Earthquake Catalog.")
    return response.json()
