import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation

def solve(u):
    coefficient = delta_t / (delta_s ** 2)

    for n in range(max_time - 1):
        for r in range(space_shape[0] - 1):
            for s in range(space_shape[1] - 1):
                u[n+1, r, s] = u[n, r, s] + coefficient * (second_order_difference(u, n, r, s, x_axis_index, delta_s) + \
                                                           second_order_difference(u, n, r, s, y_axis_index, delta_s) ) 
                

def second_order_difference(u, n, r, s, axis_index, delta):

    index1 = [n, r, s]
    index1[axis_index + 1] += 1
    index1 = tuple(index1)

    index2 = [n, r, s]
    index2 = tuple(index2)

    index3 = [n, r, s]
    index3[axis_index + 1] -= 1
    index3 = tuple(index3)

    return u[index1] - 2*u[index2] + u[index3]

x_axis_index = 0
y_axis_index = 1

space_shape = (10,10) #2d for now todo: make it arbitrary dimensional
max_time = 100

delta_s = 1 #space delta (delta x, y etc)
delta_t = 1 #time delta

u_initial = 0 #ambient initial tempreture
boundry_top = 100
boundry_bottom = 0
boundry_left = 0
boundry_right = 0

u = np.empty( (max_time,) + space_shape)

u.fill(u_initial)

# set boundry conditions
u[:, :1, :] = boundry_top
u[:, space_shape[0]-1:, :] = boundry_bottom
# no need to set others they are all zero anyways for now

solve(u)