# core/crater_model.py
import math

def estimate_crater_size(energy_megatons):
    """
    Estima o diâmetro e profundidade da cratera com base na energia.
    Fórmula aproximada derivada do Earth Impact Effects Program.
    """
    diameter_km = 1.161 * (energy_megatons ** 0.294) 
    depth_km = 0.2 * diameter_km
    return diameter_km, depth_km

def estimate_earthquake_magnitude(energy_megatons):
    """
    Converte energia de impacto em magnitude sísmica equivalente (Mw).
    """
    energy_joules = energy_megatons * 4.184e15
    magnitude = (2/3) * (math.log10(energy_joules) - 4.8)
    return magnitude
