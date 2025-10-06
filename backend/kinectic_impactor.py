import math
from impact_simulation import ImpactSimulation
from core.asteroid_estimations import AsteroidEstimations


# ---- Classe base: Asteroide ----
class Asteroid:
    def __init__(self, mass_kg: float, velocity_ms: float):
        self.mass_kg = mass_kg
        self.velocity_ms = velocity_ms

    @property
    def diameter_m(self):
        """Estimativa do diâmetro (densidade média 3000 kg/m³)"""
        return (6 * self.mass_kg / (math.pi * 3000)) ** (1 / 3)


# ---- Classe principal: Kinetic Impactor ----
class KineticImpactor(Asteroid):
    def __init__(
        self,
        mass_kg: float,
        velocity_ms: float,
        impactor_mass_kg: float,
        impactor_velocity_ms: float,
        target_lat: float,
        target_lon: float
    ):
        super().__init__(mass_kg, velocity_ms)
        self.impactor_mass_kg = impactor_mass_kg
        self.impactor_velocity_ms = impactor_velocity_ms
        self.target_lat = target_lat
        self.target_lon = target_lon
        self.simulation = ImpactSimulation()
        self.estimations = AsteroidEstimations()

    def intercept_asteroid(self):
        m1 = self.mass_kg
        v1 = self.velocity_ms
        m2 = self.impactor_mass_kg
        v2 = self.impactor_velocity_ms

        # Conservação de momento linear
        v_final = (m1 * v1 - m2 * v2) / (m1 + m2)
        slowdown_pct = (1 - abs(v_final) / v1) * 100

        # Avaliação do sucesso
        if slowdown_pct >= 99.9:
            result = "🌎 Sucesso! O asteroide foi completamente desviado."
            tip = "Excelente cálculo! O impacto anulou quase toda a velocidade do asteroide."
        elif slowdown_pct >= 80:
            result = "⚠️ Parcialmente bem-sucedido. O asteroide foi desacelerado, mas ainda representa risco."
            tip = "Tente aumentar a massa ou a velocidade do impactor."
        elif slowdown_pct >= 50:
            result = "❗ Impacto fraco. O asteroide ainda mantém boa parte da velocidade."
            tip = "Aumente consideravelmente a velocidade ou massa do impactor."
        else:
            result = "💥 Falha. O impacto foi insuficiente para alterar significativamente a trajetória."
            tip = "Use um impactor mais massivo ou muito mais rápido."

        # Simulação opcional (para integração futura)
        try:
            impact_result = self.simulation.run_custom_impact(
                diameter_m=self.diameter_m,
                velocity_kms=v_final / 1000,
                lat=self.target_lat,
                lon=self.target_lon,
                distance_km=500
            )
        except Exception:
            impact_result = None

        # ---- Relatório ----
        print("\n🛰️  KINETIC IMPACTOR MISSION REPORT")
        print("=========================================")
        print(f"🌍 Local do impacto: lat {self.target_lat:.2f}, lon {self.target_lon:.2f}")
        print("-----------------------------------------")
        print("☄️ ASTEROID")
        print(f"  Massa: {m1:,.2f} kg")
        print(f"  Velocidade inicial: {v1:,.2f} m/s")
        print("-----------------------------------------")
        print("🛰️ IMPACTOR (dados inseridos)")
        print(f"  Massa: {m2:,.2f} kg")
        print(f"  Velocidade: {v2:,.2f} m/s")
        print("-----------------------------------------")
        print("💥 COLISÃO")
        print(f"  Velocidade após impacto: {v_final:,.2f} m/s")
        print(f"  Redução de velocidade: {slowdown_pct:.2f}%")
        print("\n" + result)
        print(f"Dica: {tip}")
        print("=========================================\n")

        return {
            "asteroid_mass_kg": m1,
            "asteroid_velocity_ms": v1,
            "asteroid_final_velocity_ms": v_final,
            "asteroid_slowdown_pct": slowdown_pct,
            "impactor_mass_kg": m2,
            "impactor_velocity_ms": v2,
            "impact_result": impact_result,
        }


# ---- Execução interativa ----
if __name__ == "__main__":
    print("\n=== 🌠 SIMULADOR DE IMPACTO CINÉTICO ===\n")
    asteroid_mass = float(input("Massa do asteroide (kg): "))
    asteroid_velocity = float(input("Velocidade do asteroide (m/s): "))
    impactor_mass = float(input("Massa do impactor (kg): "))
    impactor_velocity = float(input("Velocidade do impactor (m/s): "))
    target_lat = float(input("Latitude do impacto (°): "))
    target_lon = float(input("Longitude do impacto (°): "))

    mission = KineticImpactor(
        mass_kg=asteroid_mass,
        velocity_ms=asteroid_velocity,
        impactor_mass_kg=impactor_mass,
        impactor_velocity_ms=impactor_velocity,
        target_lat=target_lat,
        target_lon=target_lon,
    )

    mission.intercept_asteroid()
