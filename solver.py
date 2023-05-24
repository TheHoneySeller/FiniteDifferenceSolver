import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import math
from functools import partial

def solve(u):
    coefficient = float(delta_t) / (delta_s ** 2)

    for n in range(max_time - 1):
        for r in range(1, space_shape[0] - 1):
            for s in range(1, space_shape[1] - 1):
                u[n+1, r, s] = u[n, r, s] + coefficient * (second_order_difference(u, n, r, s, x_axis_index) + \
                                                           second_order_difference(u, n, r, s, y_axis_index)+ k2*u[n, r, s] ) 
                cauchy_boundary(u, n+1, center, boundary_radius, r, s)
                
def solve2(u, p):
    coefficient = float(delta_t) / (delta_s ** 2)

    for n in range(max_time - 1):
        for r in range(1, space_shape[0] - 1):
            for s in range(1, space_shape[1] - 1):
                p[n+1,r,s]=p[n,r,s]-k3*delta_t*p[n,r,s]*u[n,r,s]
                u[n+1, r, s] = u[n, r, s] + coefficient * (second_order_difference(u, n, r, s, x_axis_index) + \
                                                           second_order_difference(u, n, r, s, y_axis_index)+ k2*u[n, r, s] ) \
                                                           *p[n,r,s]
                cauchy_boundary(u, n+1, center, boundary_radius, r, s)

def second_order_difference(u, n, r, s, axis_index):

    index1 = [n, r, s]
    index1[axis_index + 1] += 1
    index1 = tuple(index1)

    index2 = [n, r, s]
    index2 = tuple(index2)

    index3 = [n, r, s]
    index3[axis_index + 1] -= 1
    index3 = tuple(index3)

    return u[index1] - 2*u[index2] + u[index3]

def plot(n, u):
    plt.clf()

    plt.title(f"time: {n*delta_t}")
    
    plt.pcolormesh(u[n], cmap=plt.cm.jet, vmin=-1000, vmax=1000)
    plt.colorbar()

def neumann_boundary(u, n, center, radius, r, s):
    pass

def cauchy_boundary(u, n, center, radius, r, s):
    if is_on_boundary(r, s, center, radius):
        u[n, r, s] = cauchy_boundary_value
    elif distance((r,s), center) > radius:
        u[n, r, s] = cauchy_boundary_value

def is_on_boundary(r, s, center, radius):
    corners = ((r-0.5, s+0.5), (r-0.5, s-0.5), (r+0.5,s+0.5), (r+0.5, s-0.5))
    is_in = False
    is_out = False

    for corner in corners:
        if distance(corner, center) <= radius:
            is_in = True
        if distance(corner, center) > radius:
            is_out = True
    
    return is_in and is_out

def distance(x1, x2):
    return math.sqrt((x1[0] - x2[0])**2 + (x1[1] - x2[1])**2)        

def normal(center, r, s):
    n = (r - center[0], s - center[1])
    norm = norm(n)
    n = ( i/norm for i in n)
    return n

def norm(x):
    return math.sqrt(x[0]**2 + x[1]**2)

def set_boundary_conditions(u, u_initial):
    u.fill(u_initial)

    for n in range(max_time):
        for r in range(space_shape[0]):
            for s in range(space_shape[1]):
                if distance(center, (r, s)) > boundary_radius:
                    u[n, r, s] = boundary_value

k2 = 0.1
k3 = 0.01*k2/2.4

cauchy_boundary_value = 0
boundary_radius = 24

boundary_value = cauchy_boundary_value

x_axis_index = 0
y_axis_index = 1

space_shape = (50,50) #2d for now todo: make it arbitrary dimensional
max_time = 1000

delta_s = 1 #space delta (delta x, y etc)
delta_t = 0.1 #time delta

u_initial = 200 #ambient initial tempreture
p_initial = 1

center = (space_shape[0]/2,space_shape[1]/2)

u = np.empty( (max_time,) + space_shape)
v = np.empty( (max_time,) + space_shape)

p = np.empty( (max_time,) + space_shape)

p.fill(p_initial)

# set boundry conditions
set_boundary_conditions(u, u_initial)
set_boundary_conditions(v, u_initial)

solve2(u, p)
solve(v)


anim = animation.FuncAnimation(plt.figure(), partial(plot, u=np.subtract(v,u)), interval=1, frames=max_time, repeat=False)
anim.save("solution.gif")