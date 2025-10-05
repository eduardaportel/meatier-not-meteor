# utils/geo_utils.py
import requests

def get_elevation(lat, lon):
    url = f"https://nationalmap.gov/epqs/pqs.php?x={lon}&y={lat}&units=Meters&output=json"
    response = requests.get(url)
    return response.json()["USGS_Elevation_Point_Query_Service"]["Elevation_Query"]["Elevation"]
