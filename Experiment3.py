from MarsMission import MarsMission


def main():
    # mass 2180 kg = 3.65025e-22 earth masses similar to centaur stage on perseverance mission
    # starting height is 0.001 AU to both escape the earth's gravity and be close to it at the start
    # velocity range similar to perseverance mission
    # velocity 36,900 km/h, available range for search ±10% meaning about 2.0 au/year to 2.6 au/year
    # mars radius 2.26608e-5 AU
    # height of areosynchronous orbit is the goal distance: 0.00013588429 AU

    mission = MarsMission(3.65025e-22, 2.0, 2.6, 1.e-3, 0.00013588429, "parameters-solar (3).txt")
    # Note, this may take multiple minutes to run
    mission.calculate_trajectory()


if __name__ == "__main__":
    main()
