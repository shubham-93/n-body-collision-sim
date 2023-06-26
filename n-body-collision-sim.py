"""
Simulation of n particles of random masses, positions and velocities colliding
elastically inside a 2D box.

Author: Shubham Maheshwari
email: shubham.93@gmail.com
GitHub: https://github.com/shubham-93
"""

import random
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


class Particle:
    """
    Particle class holding particle attributes of label, mass, position and velocity.
    In addition, it has 2 methods: one to move the particle ahead in time,
    and another to check if, after moving, it hit a boundary.
    """

    def __init__(self, label, mass, x_start, y_start, vel_x, vel_y):
        """
        Initialize the label, mass, starting position coordinates and
        starting velocities of a particle.
        """
        self.label = label
        self.mass = mass
        self.x = x_start
        self.y = y_start
        self.vel_x = vel_x
        self.vel_y = vel_y

        # Set minimum resolutions for checking collisions with a boundary of the box
        self.minimum_boundary_res_x = abs(self.vel_x * dt)
        self.minimum_boundary_res_y = abs(self.vel_y * dt)


    def collision_checker(self):
        """
        Check if a particle collides with a boundary.
        Then check if it collides with another particle.
        If collision is true, then update the particle velocity (or velocities of all
        colliding particles).
        """

        # Check if a particle collides with either the bottom or top boundary,
        # and then if it collides with either the left or right boundary
        if (
            ((self.y - bottom_boundary < self.minimum_boundary_res_y) and (self.vel_y < 0))
            or
            ((top_boundary - self.y < self.minimum_boundary_res_y) and (self.vel_y > 0))
        ):
            self.vel_y = - self.vel_y

        if (
            ((self.x - left_boundary < self.minimum_boundary_res_x) and (self.vel_x < 0))
            or
            ((right_boundary - self.x < self.minimum_boundary_res_x) and (self.vel_x > 0))
        ):
            self.vel_x = - self.vel_x

        # After checking for collision with a boundary, check
        # if there are collisions between the particles themselves.
        # Note that in the main body of the program, we have created
        # a list of instances of Particle class with
        # random masses, positions and velocities initialized

        for particle in particles:
            if self.label != particle.label:
                # position vectors of self and other particle (1=self, 2=other)
                r_1_vec = np.array([self.x, self.y])
                r_2_vec = np.array([particle.x, particle.y])

                v_1_vec = np.array([self.vel_x, self.vel_y])
                v_2_vec = np.array([particle.vel_x, particle.vel_y])

                # relative displacement and velocity vectors (of particle 1 wrt 2)
                r_12_vec = r_1_vec - r_2_vec
                v_12_vec = v_1_vec - v_2_vec

                r_21_vec = -r_12_vec
                v_21_vec = -v_12_vec

                # Check if distance between particles is less than the
                # threshold distance where collision is guaranteed to occur
                # and also check if the relative displacement is opposite to
                # the relative velociy (which would imply that the 
                # particles are approaching to collide)
                if (np.linalg.norm(r_12_vec) < np.linalg.norm(v_12_vec)*dt) and (np.dot(r_12_vec, v_12_vec) < 0):
                    v_1_new_vec = v_1_vec - (2*particle.mass/(self.mass + particle.mass))*(np.dot(v_12_vec, r_12_vec)/np.power(np.linalg.norm(r_12_vec), 2))*(r_12_vec)
                    v_2_new_vec = v_2_vec - (2*self.mass/(self.mass + particle.mass))*(np.dot(v_21_vec, r_21_vec)/np.power(np.linalg.norm(r_12_vec), 2))*(r_21_vec)

                    self.vel_x = v_1_new_vec[0]
                    self.vel_y = v_1_new_vec[1]

                    particle.vel_x = v_2_new_vec[0]
                    particle.vel_y = v_2_new_vec[1]


    def update_time_step(self, dt):
        """
        Update the particle position in one time step dt, and then check for collisions
        using self.collision_checker().
        """
        self.x = self.x + self.vel_x * dt
        self.y = self.y + self.vel_y * dt

        self.collision_checker()


# Main program body starts here

# Time step and range
dt = 0.02
t = np.arange(0, 100, dt)
total_frames = len(t)


# Boundaries of the box
left_boundary = -60
right_boundary = 60
bottom_boundary = 0
top_boundary = 60


# Number of particles
num_particles = 15


# Initializing instances of Particle class with random masses, positions and velocities
particles = [Particle(label = label,
                     mass = random.randint(5,20),
                     x_start = random.randint(left_boundary, right_boundary),
                     y_start = random.randint(bottom_boundary, top_boundary),
                     vel_x = random.randint(-100,100),
                     vel_y = random.randint(-100,100)
                     )for label in range(num_particles)]


# Create the figure
fig, ax = plt.subplots()
xdata, ydata = [], []
ax.grid()
ln, = ax.plot([], [], marker='o', markersize=5)
ln.set_markerfacecolor('r')
ln.set_markeredgecolor('k')
ln.set_linestyle('')


# Figure initialization function
def init_fig():
    ax.set_xlim(-100, 100)
    ax.set_ylim(-100, 100)

    ax.plot([left_boundary, right_boundary],[bottom_boundary, bottom_boundary], linewidth=2, color='k')
    ax.plot([left_boundary, right_boundary],[top_boundary, top_boundary], linewidth=2, color='k')
    ax.plot([left_boundary, left_boundary],[bottom_boundary, top_boundary], linewidth=2, color='k')
    ax.plot([right_boundary, right_boundary],[bottom_boundary, top_boundary], linewidth=2, color='k')
    
    return ln,

# Evolve the system of particles in time
def evolution(frame):
    [particles[label].update_time_step(dt) for label in range(num_particles)]

    x_data = [particles[label].x for label in range(num_particles)]
    y_data = [particles[label].y for label in range(num_particles)]


    ln.set_data(x_data, y_data)
    return ln,


# Animation
anim = FuncAnimation(fig, evolution, init_func=init_fig, blit=True, interval=50, frames=total_frames, repeat=False)
plt.show()