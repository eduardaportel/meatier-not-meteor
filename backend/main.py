"""Main script to demonstrate usage of NASA and USGS API clients."""
from clients.nasa_client import NASAClient
from clients.usgs_client import USGSClient
from clients.base_client import APIError

def fetch_data():
    """Example function to fetch data from NASA and USGS APIs."""
    # Initialize the clients
    nasa_client = NASAClient()
    usgs_client = USGSClient()

    print("Fetching Near-Earth Object Data from NASA ---")
    try:
        result = []

        page = 0
        count = 0
        # Usage for get_neo_browse
        neo_data = nasa_client.get_neo_browse(page=page, size=20)
        total_pages = int(neo_data["page"]["total_pages"])

        for i in range(total_pages):
            for data in neo_data["near_earth_objects"]:
                print(f"NEO ID: {data['id']} | Name: {data['name']} | {count}")
                count += 1

        print(result)

    except APIError as e:
        print(f"Error fetching NASA data: {e}")

    print("\nFetching Earthquake Data from USGS ---")
    try:
        # Example: Get significant earthquakes from the last month
        earthquakes = usgs_client.get_earthquakes(
            start_time="2025-09-03T00:00:00",
            end_time="2025-10-03T23:59:59",
            min_magnitude=5.0
        )

        feature_count = len(earthquakes.get("features", []))
        print(f"Found {feature_count} earthquakes with magnitude 5.0+ in the last month.")

        # Print info about the first earthquake found
        if feature_count > 0:
            first_quake = earthquakes["features"][0]["properties"]
            print(f"Most recent quake: {first_quake['mag']} magnitude at {first_quake['place']}")

    except APIError as e:
        print(f"Error fetching USGS data: {e}")


if __name__ == "__main__":
    fetch_data()
