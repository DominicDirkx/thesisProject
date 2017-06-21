from matplotlib import pyplot as plt
from matplotlib import animation
import numpy as np
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d.axes3d as p3
import seaborn as sns
import json
import pylab
from load_data import load_orbit, load_manifold, load_bodies_location, load_lagrange_points_location


def init():
    for line in lines:
        line.set_data([], [])
        line.set_3d_properties([])
    return lines


def animate(i):
    if i > 0:
        for j, line in enumerate(lines):
            if j < numberOfOrbitsPerManifolds:
                x = manifold_S_plus.xs(j + 1)['x'].tolist()
                y = manifold_S_plus.xs(j + 1)['y'].tolist()
                z = manifold_S_plus.xs(j + 1)['z'].tolist()

            if numberOfOrbitsPerManifolds <= j < numberOfOrbitsPerManifolds * 2:
                x = manifold_S_min.xs(j - numberOfOrbitsPerManifolds + 1)['x'].tolist()
                y = manifold_S_min.xs(j - numberOfOrbitsPerManifolds + 1)['y'].tolist()
                z = manifold_S_min.xs(j - numberOfOrbitsPerManifolds + 1)['z'].tolist()

            if numberOfOrbitsPerManifolds * 2 <= j < numberOfOrbitsPerManifolds * 3:
                x = manifold_U_plus.xs(j - numberOfOrbitsPerManifolds * 2 + 1)['x'].tolist()
                y = manifold_U_plus.xs(j - numberOfOrbitsPerManifolds * 2 + 1)['y'].tolist()
                z = manifold_U_plus.xs(j - numberOfOrbitsPerManifolds * 2 + 1)['z'].tolist()

            if numberOfOrbitsPerManifolds * 3 <= j:
                x = manifold_U_min.xs(j - numberOfOrbitsPerManifolds * 3 + 1)['x'].tolist()
                y = manifold_U_min.xs(j - numberOfOrbitsPerManifolds * 3 + 1)['y'].tolist()
                z = manifold_U_min.xs(j - numberOfOrbitsPerManifolds * 3 + 1)['z'].tolist()

            for k in range(len(x)):
                if x[k] < 0.5 or x[k] > 1.5 or (y[k] or z[k]) < -0.5 or (y[k] or z[k]) > 0.5:
                    x[k] = pylab.NaN
                    y[k] = pylab.NaN
                    z[k] = pylab.NaN
                else:
                    pass

            line.set_data(x[:i], y[:i])
            line.set_3d_properties(z[:i])

    # ax.view_init(elev=10., azim=i)

    try:
        t = manifold_U_min.xs(1).index.values[i]
        time_text.set_text('t = {:.2f}'.format(round(abs(t), 2)))
    except IndexError:
        pass
    return lines


def cr3bp_velocity(x_loc, y_loc, c):
    r_1 = np.sqrt((x_loc + massParameter) ** 2 + y_loc ** 2)
    r_2 = np.sqrt((x_loc - 1 + massParameter) ** 2 + y_loc ** 2)
    v = x_loc ** 2 + y_loc ** 2 + 2 * (1 - massParameter) / r_1 + 2 * massParameter / r_2 - c
    return v


with open("../config/config.json") as data_file:
    config = json.load(data_file)


for orbit_type in config.keys():
    for orbit_name in config[orbit_type].keys():
        print(orbit_name)

        numberOfOrbitsPerManifolds = 100

        fig = plt.figure(figsize=(20, 20))
        ax = fig.add_subplot(111, projection='3d')
        # plt.rcParams['animation.ffmpeg_path'] = 'ffmpeg-git-20170607-64bit-static/ffmpeg'

        ax.set_xlim3d([0.5, 1.5])
        ax.set_xlabel('x')
        ax.set_ylim3d([-0.5, 0.5])
        ax.set_ylabel('y')
        ax.set_zlim3d([-0.5, 0.5])
        ax.set_zlabel('z')

        time_text = ax.text(1, 1, 1, s='', transform=ax.transAxes, size=22)

        orbit = load_orbit('../data/raw/' + orbit_name + '_final_orbit.txt')
        manifold_S_plus = load_manifold('../data/raw/' + orbit_name + '_W_S_plus.txt')
        manifold_S_min = load_manifold('../data/raw/' + orbit_name + '_W_S_min.txt')
        manifold_U_plus = load_manifold('../data/raw/' + orbit_name + '_W_U_plus.txt')
        manifold_U_min = load_manifold('../data/raw/' + orbit_name + '_W_U_min.txt')

        plt.plot(orbit['x'], orbit['y'], orbit['z'], color='blue')

        EARTH_GRAVITATIONAL_PARAMETER = 3.986004418E14
        SUN_GRAVITATIONAL_PARAMETER = 1.32712440018e20
        MOON_GRAVITATIONAL_PARAMETER = SUN_GRAVITATIONAL_PARAMETER / (328900.56 * (1.0 + 81.30059))
        massParameter = MOON_GRAVITATIONAL_PARAMETER / (MOON_GRAVITATIONAL_PARAMETER + EARTH_GRAVITATIONAL_PARAMETER)

        C = float(config[orbit_type][orbit_name]['C'])
        x_range = np.arange(0.5, 1.5, 0.001)
        y_range = np.arange(-0.5, 0.5, 0.001)
        X, Y = np.meshgrid(x_range, y_range)
        Z = cr3bp_velocity(X, Y, C)
        if Z.min() < 0:
            plt.contourf(X, Y, Z, 0, colors='black', alpha=0.05)

        title = 'C = ' + str(round(C, 3)) + \
                ', T = ' + str(round(float(config[orbit_type][orbit_name]['T']), 3))
        plt.title(title, size=30)

        phi = np.linspace(0, 2 * np.pi, 100)
        theta = np.linspace(0, np.pi, 100)
        bodies = load_bodies_location()
        body = 'Moon'
        x_body = bodies[body]['r'] * np.outer(np.cos(phi), np.sin(theta)) + bodies[body]['x']
        y_body = bodies[body]['r'] * np.outer(np.sin(phi), np.sin(theta)) + bodies[body]['y']
        z_body = bodies[body]['r'] * np.cos(theta) + bodies[body]['z']
        ax.plot_surface(x_body, y_body, z_body, color='black')

        lagrange_points = load_lagrange_points_location()
        for lagrange_point in ['L1', 'L2']:
            ax.scatter3D(lagrange_points[lagrange_point]['x'],
                         lagrange_points[lagrange_point]['y'],
                         lagrange_points[lagrange_point]['z'], color='grey', marker='d', alpha=0.75)
            ax.text(lagrange_points[lagrange_point]['x'],
                    lagrange_points[lagrange_point]['y'],
                    lagrange_points[lagrange_point]['z'], lagrange_point, size=16)

        numberOfOrbits = numberOfOrbitsPerManifolds * 4
        color_palette_green = sns.dark_palette('green', n_colors=numberOfOrbitsPerManifolds)
        color_palette_red = sns.dark_palette('red', n_colors=numberOfOrbitsPerManifolds)
        lines = [ax.plot(manifold_S_plus.xs(idx + 1)['x'].tolist()[:1], manifold_S_plus.xs(idx + 1)['y'].tolist()[:1], manifold_S_plus.xs(idx + 1)['z'].tolist()[:1], color=color_palette_green[idx])[0] for idx in range(numberOfOrbitsPerManifolds)]
        lines.extend([ax.plot(manifold_S_min.xs(idx + 1)['x'].tolist()[:1], manifold_S_min.xs(idx + 1)['y'].tolist()[:1], manifold_S_min.xs(idx + 1)['z'].tolist()[:1], color=color_palette_green[idx])[0] for idx in range(numberOfOrbitsPerManifolds)])
        lines.extend([ax.plot(manifold_U_plus.xs(idx + 1)['x'].tolist()[:1], manifold_U_plus.xs(idx + 1)['y'].tolist()[:1], manifold_U_plus.xs(idx + 1)['z'].tolist()[:1], color=color_palette_red[idx])[0] for idx in range(numberOfOrbitsPerManifolds)])
        lines.extend([ax.plot(manifold_U_min.xs(idx + 1)['x'].tolist()[:1], manifold_U_min.xs(idx + 1)['y'].tolist()[:1], manifold_U_min.xs(idx + 1)['z'].tolist()[:1], color=color_palette_red[idx])[0] for idx in range(numberOfOrbitsPerManifolds)])

        numberOfFrames = 0
        for index in range(1, numberOfOrbitsPerManifolds + 1):
            numberOfFrames = max(numberOfFrames, len(manifold_S_plus.xs(index)['x']))
            numberOfFrames = max(numberOfFrames, len(manifold_S_min.xs(index)['x']))
            numberOfFrames = max(numberOfFrames, len(manifold_U_plus.xs(index)['x']))
            numberOfFrames = max(numberOfFrames, len(manifold_U_min.xs(index)['x']))

        anim = animation.FuncAnimation(fig, animate, init_func=init,
                                       frames=int(numberOfFrames), interval=1, blit=True)

        Writer = animation.writers['ffmpeg']
        writer = Writer(fps=30, metadata=dict(artist='Koen Langemeijer'))
        anim.save(('../data/animations/' + orbit_name + '_3d_zoom.mp4'), writer=writer)

# plt.show()
