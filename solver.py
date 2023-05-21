import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def solve(u):
    coefficient = float(delta_t) / (delta_s ** 2)

    for n in range(max_time - 1):
        for r in range(1, space_shape[0] - 1):
            for s in range(1, space_shape[1] - 1):
                u[n+1, r, s] = u[n, r, s] + coefficient * (second_order_difference(u, n, r, s, x_axis_index) + \
                                                           second_order_difference(u, n, r, s, y_axis_index) ) 
                

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

def plot(n):
    plt.clf()

    plt.title(f"time: {n*delta_t}")
    
    plt.pcolormesh(u[n], cmap=plt.cm.jet, vmin=0, vmax=100)
    plt.colorbar()




x_axis_index = 0
y_axis_index = 1

space_shape = (50,50) #2d for now todo: make it arbitrary dimensional
max_time = 500

delta_s = 1 #space delta (delta x, y etc)
delta_t = 0.1 #time delta

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

anim = animation.FuncAnimation(plt.figure(), plot, interval=1, frames=max_time, repeat=False)
anim.save("solution.gif")