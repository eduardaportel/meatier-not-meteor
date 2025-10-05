"""A base client for making API requests with error handling."""

from typing import Optional, Dict, Any
import requests

class APIError(Exception):
    """Custom exception for API errors."""
    def __init__(self, status_code: int, message: str):
        self.status_code = status_code
        self.message = f"API Error {status_code}: {message}"
        super().__init__(self.message)

class BaseAPIClient:
    """A base client for making API requests."""

    def __init__(self, base_url: str):
        self.base_url = base_url

    def _request(
        self,
        method: str,
        endpoint: str,
        params: Optional[Dict[str, Any]] = None,
        data: Optional[Dict[str, Any]] = None,
        headers: Optional[Dict[str, str]] = None,
    ) -> Any:
        """
        Makes a request to the API.
        Args:
            method: The HTTP method (e.g., 'GET', 'POST').
            endpoint: The API endpoint path.
            params: URL parameters.
            data: Request body for POST/PUT requests.
            headers: Request headers.
        Returns:
            The JSON response from the API.
        Raises:
            APIError: If the API returns a non-2xx status code.
        """
        url = f"{self.base_url}{endpoint}"

        try:
            response = requests.request(
                method=method,
                url=url,
                params=params,
                json=data,
                headers=headers,
                timeout=10 # seconds
            )
            response.raise_for_status()  # Raises HTTPError for bad responses (4xx or 5xx)

        except requests.exceptions.HTTPError as e:
            raise APIError(
                status_code=e.response.status_code,
                message=e.response.text
            ) from e

        except requests.exceptions.RequestException as e:
            raise APIError(status_code=503, message=str(e)) from e
        
        return response.json()

    def _request_by_url(
            self,
            method: str,
            url: str,
            headers: Optional[Dict[str, str]] = None,
        ) -> Any:
        """
        Makes a request to the API through provided URL.
        Args:
            method: The HTTP method (e.g., 'GET', 'POST').
            url: The full URL to request.
            headers: Request headers.
        Returns:
            The JSON response from the API.
        Raises:
            APIError: If the API returns a non-2xx status code.
        """

        try:
            response = requests.request(
                method=method,
                url=url,
                headers=headers,
                timeout=10 # seconds
            )
            # print(response.request.url)
            response.raise_for_status()  # Raises HTTPError for bad responses (4xx or 5xx)

        except requests.exceptions.HTTPError as e:
            raise APIError(
                status_code=e.response.status_code,
                message=e.response.text
            ) from e

        except requests.exceptions.RequestException as e:
            # For network-related errors (e.g., DNS failure, connection refused)
            raise APIError(status_code=503, message=str(e)) from e
        # print(response.json())
        return response.json()

    def get(self, endpoint: str, params: Optional[Dict[str, Any]] = None) -> Any:
        """Performs a GET request."""
        r = self._request("GET", endpoint, params=params)
        return r
    
    def get_by_url(self, url: str) -> Any:
        """Performs a GET request."""
        r = self._request_by_url("GET", url)
        return r

    def post(self, endpoint: str, data: Optional[Dict[str, Any]] = None) -> Any:
        """Performs a POST request."""
        return self._request("POST", endpoint, data=data)
