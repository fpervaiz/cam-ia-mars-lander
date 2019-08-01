import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

# mass, initial position and velocity
m = 1
position0 = [0, 5000000, 0]
velocity0 = [3200, 0, 0]

# mars attributes
mars_mass = 6.42e23
mars_radius = 3386000.0

G = 6.673e-11
M = mars_mass

# simulation time, timestep and time
t_max = 20000

# critical dt for Verlet stability appears to be 1.95
dt = 1

t_array = np.arange(0, t_max, dt)

# initialise empty arrays to record trajectories
position = np.array([position0])
velocity = np.array([velocity0])

# dynamics calculation

for t in t_array:
    # calculate r and angles from last position
    last_position = position[-1]
    last_r = np.abs(np.linalg.norm(last_position))
    last_x = last_position[0]
    last_y = last_position[1]
    last_z = last_position[2]
    phi = np.arctan(last_y/last_x)
    theta = np.arccos(last_z/last_r)

    a = (-G*M)/(last_r*last_r)    
    # hack
    if last_x < 0:
        a = -a        
    
    # calculate new position and velocity    
    aX = a*np.sin(theta)*np.cos(phi)
    aY = a*np.sin(theta)*np.sin(phi)
    aZ = a*np.cos(theta)

    a_vector = np.array([aX, aY, aZ])

    if t == 0:
        # Euler integration
        new_pos = [last_position + dt * velocity[-1]]
        new_vel = [(last_position + dt * a_vector)]
        position = np.append(position, new_pos, axis=0)
        velocity = np.append(velocity, new_vel, axis=0)
    else:
        # Verlet integration
        new_pos = [2*last_position - position[-2] + dt*dt*a_vector]
        new_vel =[(last_position - position[-2])/dt]
        position = np.append(position, new_pos, axis=0)
        velocity = np.append(velocity, new_vel, axis=0)

    new_r = np.abs(np.linalg.norm(new_pos))
    if new_r < mars_radius:
        print("Collided with Mars at time {}s".format(t))
        break


x_pos = position[:,0]
y_pos = position[:,1]
z_pos = position[:,2]

fig = plt.figure()
ax = plt.axes(projection='3d')
#ax._axis3don = False

# draw Mars sphere
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x = mars_radius * np.outer(np.cos(u), np.sin(v))
y = mars_radius * np.outer(np.sin(u), np.sin(v))
z = mars_radius * np.outer(np.ones(np.size(u)), np.cos(v))
ax.plot_surface(x, y, z,  rstride=4, cstride=4, color='r')

# plot orbit
ax.plot3D(x_pos, y_pos, z_pos)
ax.set_aspect('equal')

# start set equal axes
max_range = np.array([x_pos.max()-x_pos.min(), y_pos.max()-y_pos.min(), z_pos.max()-z_pos.min()]).max() / 2.0

mid_x = (x_pos.max()+x_pos.min()) * 0.5
mid_y = (y_pos.max()+y_pos.min()) * 0.5
mid_z = (z_pos.max()+z_pos.min()) * 0.5
ax.set_xlim(mid_x - max_range, mid_x + max_range)
ax.set_ylim(mid_y - max_range, mid_y + max_range)
ax.set_zlim(mid_z - max_range, mid_z + max_range)
# end set equal axes
plt.show()



