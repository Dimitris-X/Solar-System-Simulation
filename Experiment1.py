from PlanetarySystem import PlanetarySystem
from Planet import Planet


def main():
    system = PlanetarySystem("parameters-solar (1).txt", integrator="beeman")
    system.run_simulation()
    system.print_years()


if __name__ == "__main__":
    main()
