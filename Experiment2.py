from PlanetarySystem import PlanetarySystem
import matplotlib.pyplot as plt
import numpy as np


def main():
    # creates two identical systems, one using beeman integration and one using euler integration
    system_beeman = PlanetarySystem("parameters-solar (2).txt", integrator="beeman")
    system_euler = PlanetarySystem("parameters-solar (2).txt", integrator="euler")
    system_beeman.run_simulation()
    system_euler.run_simulation()

    # creates the axes of the graph
    plt.style.use("default")
    x_values = np.linspace(0, 1000, num=len(system_beeman.energy_history))
    plt.ylabel("Total energy (joules)")
    plt.xlabel("Time elapsed (earth years)")
    plt.title("Energy vs Time, Beeman and Euler (timestep=0.01)")
    system_beeman.print_energy_stats()

    # displays the graph
    plt.scatter(x_values, system_beeman.energy_history, s=1, label="Beeman integration")
    plt.scatter(x_values, system_euler.energy_history, s=1, label="Euler integration")
    plt.legend()
    plt.show()


if __name__ == "__main__":
    main()
