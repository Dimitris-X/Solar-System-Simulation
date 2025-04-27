import numpy as np


class Planet:
    """
    represents a planet that affects other with its gravity and is affected by the gravity of others
    """
    def __init__(self, mass, position, velocity, name="unnamed", colour=(0, 0, 0)):
        self.mass = mass
        self.pos = position
        self.pos_old = position
        self.vel = velocity
        self.name = name
        self.acc = 0
        self.acc_old = 0
        self.force = 0
        self.colour = colour
        self.potential = 0
        self.new_years_list = [0]

    def update_euler(self, time_step):
        """
        makes all the changes necessary for a single timestep using euler integration
        """
        self.pos_old = self.pos

        self.acc_old = self.acc

        self.acc = self.force / self.mass

        self.pos = self.pos + self.vel * time_step

        self.vel += self.acc * time_step

    def update_position_beeman(self, time_step):
        """
        updates the position of the planet according to beeman integration and stores the previous position
        """
        self.pos_old = self.pos
        self.pos = self.pos + self.vel * time_step + (self.acc / 2 + (self.acc - self.acc_old) / 6) * time_step ** 2

    def update_velocity_beeman(self, time_step):
        """
        updates the velocity of the planet according to beeman integration and stores the previous position
        """
        new_acc = self.force / self.mass
        self.vel = self.vel + (2*new_acc+5*self.acc-self.acc_old)*time_step / 6
        self.acc_old = self.acc
        self.acc = new_acc

    def check_new_year(self, total_time, time_step):
        # this code checks if the planet has crossed the +x axis (if it's orbiting counterclockwise)
        # or the -x axis (if it's orbiting clockwise)
        if self.pos_old[1] < 0 < self.pos[1]:
            # assumes that the planet moved linearly between the two positions and calculates when the planet crossed
            # the x axis as a percentage of the last timestep
            timestep_percentage = self.pos_old[1]/(self.pos_old[1] - self.pos[1])

            self.new_years_list.append(total_time + time_step * timestep_percentage)

    def get_year_stats(self):
        """
        returns the planets average orbital period and the associated deviation by integrating
        through the list of new years
        """

        # since the list is initialized as [0] the number of years passed in the planet
        # is equal to the length of the list minus 1

        if len(self.new_years_list) == 1:
            # when the length of the list is one, the planet never crossed the +x axis so flot('inf')
            # is given as default value
            return float('inf'), float('inf')
        else:
            average_period = self.new_years_list[-1] / (len(self.new_years_list) - 1)
            variance = 0
            for i in range(len(self.new_years_list)-1):
                period = self.new_years_list[i+1] - self.new_years_list[i]
                variance += (period - average_period)**2
            variance = variance / (len(self.new_years_list)-1)
            standard_deviation = np.sqrt(variance)
            return average_period, standard_deviation

    def get_kinetic_energy(self):
        """
        returns the kinetic energy of the planet
        """
        return self.mass * np.dot(self.vel, self.vel) / 2

    def __str__(self):
        # return f"{self.name} at {self.pos}, with speed {self.vel} and velocity {self.vel}"
        return self.name

    def add_force(self, force):
        """
        adds a force to the planets total force
        """
        self.force += force

    def add_potential(self, potential):
        """
        adds a potential to the planets total gravitational potential
        """
        self.potential += potential
