# core/tsunami_model.py
import math

def estimate_transient_crater(diameter_m, velocity_ms, density):
    """ 
    Modelo matemático validado experimentalmente por Ward & Asphaug 
    """
    transient_crater_diameter_m = 1.161 * (diameter_m ** 0.78) * (velocity_ms ** 0.44) * (density/1000) ** 0.26 * 9.81 ** (-0.22)
    transient_crater_radius_m = transient_crater_diameter_m / 2
    average_crater_depth_m = 0.25 * transient_crater_diameter_m
    
    return transient_crater_diameter_m, transient_crater_radius_m, average_crater_depth_m

def estimate_initial_tsunami_height(average_crater_depth_m):
    """
    Modelo simplificado de altura de tsunami em função da energia e distância.
    """
    initial_height_m = 0.25 * average_crater_depth_m
    return initial_height_m

def variable_tsunami_height(initial_height_m, transient_crater_radius_m, distance_km):
    """
    Altura da tsunami em função da distância em km do ponto inicial do impacto
    """
    tsunami_height_m = initial_height_m * (transient_crater_radius_m / (distance_km * 1000)) ** 0.5
    return tsunami_height_m

def shore_tsunami_height(tsunami_height_m, average_oceanic_depth, shore_depth):
    """
    Altura da tsunami na costa 
    """
    shore_tsunami_height_m = tsunami_height_m * (average_oceanic_depth / shore_depth) ** 0.25
    return shore_tsunami_height_m
