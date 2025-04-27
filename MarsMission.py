import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from Planet import Planet
from PlanetarySystem import PlanetarySystem
from Rocket import Rocket
import copy


class MarsMission:
    """
    this class handles the calculations for the optimal mars mission
    """
    def __init__(self, rocket_mass, min_vel, max_vel, distance_earth, goal_distance, filename_read):
        self.system_plain = PlanetarySystem(filename_read)  # planetary system without rocket
        self.step = self.system_plain.step
        self.rocket_mass = rocket_mass
        self.min_vel = min_vel
        self.max_vel = max_vel
        self.distance = distance_earth  # the initial distance from the earth
        self.goal_distance = goal_distance  # the desired minimum distance from the earth

    def create_simulation(self, angle, velocity, max_years):
        """
        creates a simulation where the rocket has the specified
        """
        system = copy.deepcopy(self.system_plain)  # copies the prototype system
        system.limit = max_years/self.step
        for planet in system.planets:
            if planet.name == "earth":
                earth = planet
            elif planet.name == "mars":
                mars = planet
        rocket = Rocket(earth.pos, earth.vel, self.rocket_mass, self.distance, angle, velocity, mars)
        system.planets.insert(2, rocket)  # inserts the rocket into the list of planets of the system
        system.run_simulation()
        print(f"closest approach = {rocket.closest_dist} AU at time {rocket.closest_time} years with angle {angle} radians and speed {velocity} AU/year")
        return rocket.closest_dist, rocket.closest_time

    def run_animation(self, angle, velocity):
        """
        runs an animation with the specified initial conditions
        """
        system = copy.deepcopy(self.system_plain)  # creates a copy of the prototype system
        for planet in system.planets:
            if planet.name == "earth":
                earth = planet
            elif planet.name == "mars":
                mars = planet
        rocket = Rocket(earth.pos, earth.vel, self.rocket_mass, self.distance, angle, velocity, mars)
        system.planets.insert(1, rocket)  # inserts the rocket into the list of planets of the system
        system.run_animation()

        print(f"closest approach = {rocket.closest_dist} AU")
        print(f"time of closest approach = {rocket.closest_time} years")

    def search_range(self, angle_nums, vel_nums, max_years):
        """
        divides the range of available speeds and angles into a specified number of values,
        checks every combination of speed and angle,
        and finds which one of those results in the smallest closest approach to mars
        """
        outcome_dict = {}  # dictonary of all the possibilities checked
        for i in range(angle_nums):
            for j in range(vel_nums):
                angle = 2*np.pi*i/angle_nums
                velocity = self.min_vel + j*(self.max_vel-self.min_vel)/(vel_nums - 1)
                outcome = self.create_simulation(angle, velocity, max_years)

                # the stating conditions work as a key
                # a tuple containing the closest approach and corresponding time are the value
                outcome_dict.update({(angle, velocity): outcome})

        min_dist = (0, self.min_vel)  # default value for best starting conditions
        # the optimal starting conditions are calculated
        for condition, outcome in outcome_dict.items():
            if outcome_dict[condition][0] < outcome_dict[min_dist][0]:
                min_dist = condition

        # return a tuple containing the starting conditions and
        # a tuple containing the closest approach and corresponding time
        return min_dist, outcome_dict[min_dist]

    def hill_climb(self, initial_distance, initial_angle, initial_vel, angle_step, vel_step, max_years, max_depth):
        """
        this method start with a set of initial conditions and tries altering them by the specified "steps" to improve
        the closest distance to mars. It then recursively calls itself with the improved conditions. If none of the
        possible changes are improvements, it calls itself again, this time halving the size of the angle and velocity
        increments. This process repeats until the desired closest distance is achieved,
        or the maximum recursion depth is reached
        """
        # checks if the desired distance has been reached
        if initial_distance < self.goal_distance:
            print("Goal reached!")
            print(f"the optimal angle is {initial_angle} radians")
            print(f"the optimal speed is {initial_vel} AU/year")
            dist, time = self.create_simulation(initial_angle, initial_vel, max_years)
            print(f"the final closest distance is {dist} AU")
            print(f"the final travel time is {time} years")
            print(f"The final increments were: vel={vel_step} and angle={angle_step}")
            print(f"Max depth remaining = {max_depth}")
            return initial_angle, initial_vel

        # if the maximum depth is reached the recursive chain is terminated
        if max_depth == 0:
            print("maximum recursion depth reached")
            print(f"the optimal angle is {initial_angle} radians")
            print(f"the optimal speed is {initial_vel} AU/year")
            dist, time = self.create_simulation(initial_angle, initial_vel, max_years)
            print(f"the final closest distance is {dist} AU")
            print(f"the final travel time is {time} years")
            print(f"The final increments were: vel={vel_step} and angle={angle_step}")
            print(f"The final increments were: vel={vel_step} and angle={angle_step}")
            return

        # checks if increasing the angle by an increment improves the distance
        potential_improvement = self.create_simulation(initial_angle + angle_step, initial_vel, max_years)
        if potential_improvement[0] < initial_distance:
            return self.hill_climb(potential_improvement[0], initial_angle + angle_step, initial_vel, angle_step,
                                   vel_step, max_years, max_depth - 1)

        # checks if decreasing the angle by an increment improves the distance
        potential_improvement = self.create_simulation(initial_angle - angle_step, initial_vel, max_years)
        if potential_improvement[0] < initial_distance:
            return self.hill_climb(potential_improvement[0], initial_angle - angle_step, initial_vel, angle_step,
                                   vel_step, max_years, max_depth - 1)

        # checks if increasing the velocity by an increment improves the distance
        # as long this doesn't violate the maximum velocity
        if initial_vel + vel_step < self.max_vel:
            potential_improvement = self.create_simulation(initial_angle, initial_vel + vel_step, max_years)
            if potential_improvement[0] < initial_distance:
                return self.hill_climb(potential_improvement[0], initial_angle, initial_vel + vel_step, angle_step,
                                       vel_step, max_years, max_depth - 1)

        # checks if decreases the velocity by an increment improves the distance
        # as long this doesn't violate the minimum velocity
        if initial_vel - vel_step > self.min_vel:
            potential_improvement = self.create_simulation(initial_angle, initial_vel - vel_step, max_years)
            if potential_improvement[0] < initial_distance:
                return self.hill_climb(potential_improvement[0], initial_angle, initial_vel - vel_step, angle_step,
                                       vel_step, max_years, max_depth - 1)

        # if none of the changes are improvements, the size of the increments is decreased
        return self.hill_climb(initial_distance, initial_angle, initial_vel, angle_step / 2, vel_step / 2,
                               max_years, max_depth - 1)

    def calculate_trajectory(self, angle_nums=15, vel_nums=7, angle_step=0.1, vel_step=0.1, max_years=1, max_depth=200):
        """
        this method combines search_rage and hill_climb to give calculate the optimal launch conditions
        """
        # the search_range method is used to produce an initial guess for the hill_climb method
        timestep = 0.001
        self.change_time_step(timestep)
        initial_search = self.search_range(angle_nums, vel_nums, max_years)

        # the answer from the search_range method becomes the initial guess fo the hill_climb method
        initial_distance = initial_search[1][0]
        initial_angle = initial_search[0][0]
        initial_vel = initial_search[0][1]

        temp_solution = self.hill_climb(initial_distance, initial_angle, initial_vel, angle_step, vel_step, max_years, max_depth)

        timestep = timestep / 2
        angle_step = angle_step * 1 / 3
        vel_step = vel_step * 1 / 3
        self.change_time_step(timestep)
        initial_distance = self.create_simulation(temp_solution[0], temp_solution[1], max_years)[0]
        temp_solution_old = temp_solution
        temp_solution = self.hill_climb(initial_distance, temp_solution[0], temp_solution[1], angle_step, vel_step,
                                        max_years, max_depth)
        while temp_solution != temp_solution_old:
            timestep = timestep/2
            angle_step = angle_step * 1/4
            vel_step = vel_step * 1/4
            self.change_time_step(timestep)

            initial_distance = self.create_simulation(temp_solution[0], temp_solution[1], max_years)[0]
            temp_solution_old = temp_solution
            temp_solution = self.hill_climb(initial_distance, temp_solution[0], temp_solution[1], angle_step, vel_step,
                                            max_years, max_depth)

        # the final solutions are printed
        print(f"Final solution: {temp_solution}")
        print(f"Final timestep: {timestep}")

    def change_time_step(self, new_time_step):
        """
        changes the timestep used in the class to a specified value
        """
        self.step = new_time_step
        self.system_plain.step = new_time_step
