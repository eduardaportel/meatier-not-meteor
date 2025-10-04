"""Configuration module for API keys and base URLs."""
import os
from dotenv import load_dotenv

# Load environment variables from .env file
load_dotenv()

# API Keys
NASA_API_KEY = os.getenv("NASA_API_KEY", "DEMO_KEY")
duda_api_key = 'o9eX4otm7mdz9oSxXmlE5pTJG63QS4uy2cGCpyo2'

# Base URLs
NASA_BASE_URL = "https://api.nasa.gov"
USGS_BASE_URL = "https://earthquake.usgs.gov/fdsnws/event/1"
