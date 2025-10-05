"""Controller module for handling API requests related to space data."""

import os
import sys
from flask import Flask, jsonify
import json
from datetime import datetime

from operator import itemgetter

from clients.nasa_client import NASAClient
from core.asteroid_estimations import AsteroidEstimations

app = Flask(__name__)

@app.route("/")
def index():
    """Health check endpoint."""
    return jsonify({
        "status": "ok",
        "message": "Space Apps API is running.",
        "endpoints": ["/api/v1/neos", "/api/v1/earthquakes"]
    })

@app.route("/api/v1/neos")
def get_neos():
    """
    Endpoint to retrieve Near-Earth Objects (NEOs) from NASA.
    Query Parameters:
        - start_date (str, YYYY-MM-DD): Defaults to today.
        - end_date (str, YYYY-MM-DD): Defaults to 7 days from start_date.
    """
    nasa_client = NASAClient()
    asteroid_estimations = AsteroidEstimations()

    result = []
    count = 0
    neo_data = nasa_client.get_neo_browse(page=0, size=20)

    total_pages = int(neo_data["page"]["total_pages"])

    for i in range(total_pages):
        for data in neo_data["near_earth_objects"]:
            print(f"NEO ID: {data['id']} | Name: {data['name']} | {count}")
            count += 1

            diameter_m = sum([
                float(data["estimated_diameter"]["meters"]["estimated_diameter_min"]),
                float(data["estimated_diameter"]["meters"]["estimated_diameter_max"])
            ]) / 2
            
            velocity_kms = 0
            now = datetime.now()
            
            date_format = '%Y-%m-%d'
            
            impact_estimations = {}
            mass_kg = ""
            
            for d in data["close_approach_data"]:
                
                date_obj = datetime.strptime(d["close_approach_date"], date_format)
                
                if date_obj >= now:
                    velocity_kms = float(d["relative_velocity"]["kilometers_per_second"])
                    impact_estimations = asteroid_estimations.estimate_impact_energy(diameter_m, velocity_kms)
                    break

            # mass_kg, energy_megatons, energy_joules = itemgetter("mass_kg", "energy_megatons", "energy_joules")
            # (asteroid_estimations.estimate_impact_energy(diameter_m, velocity_kms))

            print(impact_estimations)
            if 'mass_kg' in impact_estimations:
                mass_kg = impact_estimations["mass_kg"]

            dict = {
                "name": data["name"],
                "diameter_m": diameter_m,
                "velocity_kms": velocity_kms,
                "mass_kg": mass_kg,
                "orbit": data["orbital_data"]
            }
            result.append(dict)

    # print(result)

    # try:

    #   if nasa_client is None:
    #     # Open and load the JSON file
    #     with open('backend/controller/MOCK_NEO_DATA.json', 'r', encoding="utf-8") as json_file:
    #         data = json.load(json_file)

    return jsonify(result), 200

    # except FileNotFoundError:
    #     return jsonify({"error": "File not found"}), 404
    # except json.JSONDecodeError:
    #     return jsonify({"error": "Invalid JSON format in file"}), 500

if __name__ == "__main__":
    app.run(debug=True, port=5001)
