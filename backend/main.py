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
        # Usage for get_neos_feed
        # start_date = "2025-10-03"
        # end_date = "2025-10-10"
        # neo_data = nasa_client.get_neos_feed(start_date=start_date, end_date=end_date)

        # element_count = neo_data.get("element_count", 0)
        # first_day_neos = neo_data.get("near_earth_objects", {}).get(start_date, [])

        # print(f"Found {element_count} Near-Earth Objects between {start_date} and {end_date}.")
        # if first_day_neos:
        #     print(f"First NEO: {first_day_neos[0]['name']}")
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

                # dict = {
                #     "name": data["name"],
                #     "diameter_m": (
                #         data["estimated_diameter"]["meters"]["estimated_diameter_min"] +
                #         data["estimated_diameter"]["meters"]["estimated_diameter_max"]
                #     ) / 2,
                #     "velocity_kms": float(data["close_approach_data"][0]["relative_velocity"]["kilometers_per_second"]),
                #     "mass_kg": None,  # será calculado depois
                #     "orbit": data["orbital_data"]
                # }
                # result.append(dict)

        print(result)

        # size = neo_data.get("page", {}).get("size", [])
        # pages = neo_data.get("page", {}).get("total_pages", [])
        # print(f"Size {size} | Pages: {pages}")


        # Use for get_neo_lookup
        # if first_day_neos:
        #     neo_id = first_day_neos[0]['id']
        #     neo_details = nasa_client.get_neo_lookup(neo_id=neo_id)
        #     print(f"Details for NEO ID {neo_id}: Name {neo_details['name']}")

        #     # or use get_neo_lookup_by_url
        #     neo_url = first_day_neos[0]['links']['self'] # Attention: URL exposes the API key
        #     neo_details = nasa_client.get_neo_lookup_by_url(neo_url)
        #     print(f"Details for NEO {neo_id} by URL: Name - {neo_details['name']}")

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
