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
