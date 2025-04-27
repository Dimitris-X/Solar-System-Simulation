from Planet import Planet
import numpy as np


class Rocket(Planet):
    """
    represents the rocket that is sent to mars
    expands the Planet class
    """
    def __init__(self, earth_position, earth_velocity, mass, distance, angle, velocity, mars):
        self.initial_velocity = velocity
        self.angle = angle
        unit_direction = np.array([np.cos(angle), np.sin(angle)])
        self.time = 0

        # the initial conditions for the planet are translated to positions and velocity vectors in the Planet class
        super().__init__(mass, earth_position + unit_direction * distance, earth_velocity + unit_direction * velocity,
                         name="rocket", colour="w")
        self.mars = mars

        # the initial distance to mars is determined
        distance_vector = self.mars.pos - self.pos
        self.closest_dist = np.sqrt(distance_vector @ distance_vector)
        self.closest_time = 0

    def update_position_beeman(self, time_step):
        """
        Overrides the corresponding method in the Planet class so the distance to mars is also considered
        """
        self.pos_old = self.pos
        self.pos = self.pos + self.vel * time_step + (self.acc / 2 + (self.acc - self.acc_old) / 6) * time_step ** 2
        self.time += time_step
        self.check_mars_distance()  # checks if this is the closest distance to mars reached yet

    def check_mars_distance(self):
        """
        Checks the distance of the rocket with mars and updates the closest distance and correpsonding time
        """
        distance_vector = self.mars.pos - self.pos
        distance = np.sqrt(distance_vector @ distance_vector)
        if distance < self.closest_dist:
            self.closest_dist = distance
            self.closest_time = self.time
