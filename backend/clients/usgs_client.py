"""Client for USGS Earthquake Catalog API."""
import config
from .base_client import BaseAPIClient

class USGSClient(BaseAPIClient):
    """Client for USGS Earthquake Catalog API."""

    def __init__(self):
        super().__init__(base_url=config.USGS_BASE_URL)

    def get_earthquakes(self, start_time: str, end_time: str, min_magnitude: float):
        """
        Fetches earthquake events within a given time range and minimum magnitude.
        API Docs: https://earthquake.usgs.gov/fdsnws/event/1/
        """
        endpoint = "/query"
        params = {
            "format": "geojson",
            "starttime": start_time,
            "endtime": end_time,
            "minmagnitude": min_magnitude,
        }
        return self.get(endpoint, params=params)

    def get_elevation(self, lat: float, lon: float, units: str = "Meters") -> float:
            """
            Fetches elevation for a given latitude and longitude.
            Returns elevation in meters (can be negative for ocean floor).
            """
            endpoint = "https://nationalmap.gov/epqs/pqs.php"
            params = {
                "x": lon,
                "y": lat,
                "units": units,
                "output": "json",
            }
            data = self.get(endpoint, params=params)
            try:
                service = data.get("USGS_Elevation_Point_Query_Service", {})
                query = service.get("Elevation_Query") or service.get("Elevation_Queries", [{}])[0]
                elevation = query.get("Elevation")
                if elevation is None:
                    raise ValueError
                return float(elevation)
            except Exception:
                raise ValueError(f"Unexpected response format from USGS: {data}")