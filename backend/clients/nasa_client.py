"""Client for NASA's APIs."""
import config # Import our configuration
from clients.base_client import BaseAPIClient

class NASAClient(BaseAPIClient):
    """Client for NASA's APIs."""

    def __init__(self):
        super().__init__(base_url=config.NASA_BASE_URL)
        self.api_key = config.NASA_API_KEY


    def get_neos_feed(self, start_date: str, end_date: str):
        """
        Fetches a list of Near-Earth Objects based on their closest approach date to Earth.
        API Docs: https://api.nasa.gov/ (Neo-Feed)
        """
        endpoint = "/neo/rest/v1/feed"
        params = {
            "start_date": start_date,
            "end_date": end_date,
            "api_key": self.api_key,
        }
        return self.get(endpoint, params=params)

    def get_neo_lookup(self, neo_id: str):
        """
        Fetches a specific Near-Earth Object by its NASA JPL small-body ID.
        API Docs: https://api.nasa.gov/ (Neo-Lookup)
        """
        endpoint = f"/neo/rest/v1/neo/{neo_id}"
        params = {
            "api_key": self.api_key,
        }
        return self.get(endpoint, params=params)

    def get_neo_lookup_by_url(self, url):
        """
        Fetches a specific Near-Earth Object by its entire URL provided by Neo-Feed or Neo-Browse. E.g.:
        {
            "links": {
                "self": "http://api.nasa.gov/neo/rest/v1/neo/2000433?api_key=DEMO_KEY"
            },
            "id": "2000433",
            ...
        }
        API Docs: https://api.nasa.gov/ (Neo-Lookup)
        """
        return self.get_by_url(url)
    

    def get_neo_browse(self, page: int = 0, size: int = 20):
        """
        Fetches a paginated list of Near-Earth Objects.
        API Docs: https://api.nasa.gov/ (Neo-Browse)
        """
        endpoint = "/neo/rest/v1/neo/browse"
        params = {
            "page": page,
            "size": size,
            "api_key": self.api_key,
        }
        return self.get(endpoint, params=params)
