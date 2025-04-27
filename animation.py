from PlanetarySystem import PlanetarySystem


def main():
    system = PlanetarySystem("parameters-solar (animation).txt", integrator="beeman")
    system.run_animation()


if __name__ == "__main__":
    main()
