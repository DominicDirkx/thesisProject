import pandas as pd
import numpy as np


def load_manifold(file_path):
    pd.options.mode.chained_assignment = None  # Turn off SettingWithCopyWarning

    input_data = pd.read_table(file_path, delim_whitespace=True, header=None).filter(list(range(7)))
    input_data.columns = ['time', 'x', 'y', 'z', 'xdot', 'ydot', 'zdot']
    split_on_index = list(input_data[input_data['time'] == 0].index)

    output_data = []
    for idx, start_index in enumerate(split_on_index):

        if idx != len(split_on_index)-1:
            data_per_orbit = input_data[start_index:split_on_index[idx+1]]
        else:
            data_per_orbit = input_data[start_index:]

        data_per_orbit['orbitNumber'] = idx
        # data_per_orbit['orbitNumber'] = idx + 1
        output_data.append(data_per_orbit)
        pass

    output_data = pd.concat(output_data).reset_index(drop=True).set_index(['orbitNumber', 'time'])
    return output_data

def load_manifold_refactored(file_path):
    pd.options.mode.chained_assignment = None  # Turn off SettingWithCopyWarning

    input_data = pd.read_table(file_path, delim_whitespace=True, header=None).filter(list(range(7)))
    input_data.columns = ['time', 'x', 'y', 'z', 'xdot', 'ydot', 'zdot']
    split_on_index = list(input_data[input_data['time'] == 0].index)

    output_data = []

    if input_data['time'].mean() > 0:
        # For stable manifolds
        for idx, start_index in enumerate(split_on_index):

            if idx != len(split_on_index)-1:
                data_per_orbit = input_data[start_index:split_on_index[idx+1]]
            else:
                data_per_orbit = input_data[start_index:]

            data_per_orbit['orbitNumber'] = idx
            output_data.append(data_per_orbit)
    else:
        # For unstable manifolds
        for idx, end_index in enumerate(split_on_index):

            if idx != 0:
                data_per_orbit = input_data[split_on_index[idx - 1] + 1:end_index + 1]
            else:
                data_per_orbit = input_data[:end_index + 1]

            data_per_orbit['orbitNumber'] = idx
            output_data.append(data_per_orbit)

    output_data = pd.concat(output_data).reset_index(drop=True).set_index(['orbitNumber', 'time'])
    return output_data


def load_manifold_incl_stm(file_path):
    pd.options.mode.chained_assignment = None  # Turn off SettingWithCopyWarning

    input_data = pd.read_table(file_path, delim_whitespace=True, header=None).rename(columns={0: 'time'})
    split_on_index = list(input_data[input_data['time'] == 0].index)

    output_data = []
    for idx, start_index in enumerate(split_on_index):

        if idx != len(split_on_index)-1:
            data_per_orbit = input_data[start_index:split_on_index[idx+1]]
        else:
            data_per_orbit = input_data[start_index:]

        data_per_orbit['orbitNumber'] = idx

        output_data.append(data_per_orbit)
        pass

    output_data = pd.concat(output_data).reset_index(drop=True).set_index(['orbitNumber', 'time'])
    return output_data


def load_orbit(file_path):
    data = pd.read_table(file_path, delim_whitespace=True, header=None).filter(list(range(7)))
    data.columns = ['time', 'x', 'y', 'z', 'xdot', 'ydot', 'zdot']
    return data


def load_initial_conditions(file_path):
    data = pd.read_table(file_path, delim_whitespace=True, header=None)
    data.columns = ['orbitId', 'C', 'T', 'x', 'y', 'z', 'xdot', 'ydot', 'zdot',
                    'lambda1real', 'lambda1imag', 'lambda2real', 'lambda2imag', 'lambda3real', 'lambda3imag',
                    'lambda4real', 'lambda4imag', 'lambda5real', 'lambda5imag', 'lambda6real', 'lambda6imag']
    return data


def load_initial_conditions_incl_M(file_path):
    data = pd.read_table(file_path, delim_whitespace=True, header=None)
    # data.columns = ['orbitId', 'C', 'T', 'x', 'y', 'z', 'xdot', 'ydot', 'zdot']
    return data


def load_differential_corrections(file_path):
    data = pd.read_table(file_path, delim_whitespace=True, header=None)
    # data.columns = ['numberOfIterations', 'C', 'T', 'x', 'y', 'z', 'xdot', 'ydot', 'zdot']
    return data


def load_lagrange_points_location():
    location_lagrange_points = {'L1': [0.8369151483688, 0, 0],
                                'L2': [1.1556821477825, 0, 0],
                                'L3': [-1.003037581609, 0, 0],
                                'L4': [0.4878494189827, 0.8660254037844, 0],
                                'L5': [0.4878494189827, -0.8660254037844, 0]}

    location_lagrange_points = pd.DataFrame.from_dict(location_lagrange_points)
    location_lagrange_points.index = ['x', 'y', 'z']
    return location_lagrange_points


def load_bodies_location():
    EARTH_GRAVITATIONAL_PARAMETER = 3.986004418E14
    SUN_GRAVITATIONAL_PARAMETER = 1.32712440018e20
    MOON_GRAVITATIONAL_PARAMETER = SUN_GRAVITATIONAL_PARAMETER / (328900.56 * (1.0 + 81.30059))
    massParameter = MOON_GRAVITATIONAL_PARAMETER / (MOON_GRAVITATIONAL_PARAMETER + EARTH_GRAVITATIONAL_PARAMETER)

    location_bodies = {'Earth': [-massParameter, 0, 0, 6371 / 384400],
                       'Moon': [1-massParameter, 0, 0, 1737 / 384400]}
    location_bodies = pd.DataFrame.from_dict(location_bodies)
    location_bodies.index = ['x', 'y', 'z', 'r']
    return location_bodies


def cr3bp_velocity(x_loc, y_loc, c):
    EARTH_GRAVITATIONAL_PARAMETER = 3.986004418E14
    SUN_GRAVITATIONAL_PARAMETER = 1.32712440018e20
    MOON_GRAVITATIONAL_PARAMETER = SUN_GRAVITATIONAL_PARAMETER / (328900.56 * (1.0 + 81.30059))
    massParameter = MOON_GRAVITATIONAL_PARAMETER / (MOON_GRAVITATIONAL_PARAMETER + EARTH_GRAVITATIONAL_PARAMETER)

    r_1 = np.sqrt((massParameter + x_loc) ** 2 + y_loc ** 2)
    r_2 = np.sqrt((1 - massParameter - x_loc) ** 2 + y_loc ** 2)
    v = x_loc ** 2 + y_loc ** 2 + 2 * (1 - massParameter) / r_1 + 2 * massParameter / r_2 - c
    return v


def computeJacobiEnergy(x, y, z, xdot, ydot, zdot):
    EARTH_GRAVITATIONAL_PARAMETER = 3.986004418E14
    SUN_GRAVITATIONAL_PARAMETER = 1.32712440018e20
    MOON_GRAVITATIONAL_PARAMETER = SUN_GRAVITATIONAL_PARAMETER / (328900.56 * (1.0 + 81.30059))
    mass_parameter = MOON_GRAVITATIONAL_PARAMETER / (MOON_GRAVITATIONAL_PARAMETER + EARTH_GRAVITATIONAL_PARAMETER)
    r_1 = np.sqrt((mass_parameter + x) ** 2 + y ** 2 + z ** 2)
    r_2 = np.sqrt((1 - mass_parameter - x) ** 2 + y ** 2 + z ** 2)
    v = np.sqrt(xdot**2 + ydot**2 + zdot**2)
    c = x**2 + y**2 + 2*(1-mass_parameter)/r_1 + 2*mass_parameter/r_2 - v**2
    return c


if __name__ == '__main__':
    # manifold_file_path = '../data/near_vertical_1_W_S_min.txt'
    # manifold_df = load_manifold(manifold_file_path)
    # print(manifold_df)

    # orbit_file_path = '../data/near_vertical_1_final_orbit.txt'
    # orbit_df = load_orbit(orbit_file_path)
    # print(orbit_df)

    # initial_conditions_file_path = '../data/raw/horizontal_L2_initial_conditions.txt'
    # initial_conditions_df = load_initial_conditions(initial_conditions_file_path)
    # print(initial_conditions_df)

    # initial_conditions_file_path = '../data/raw/horizontal_L1_initial_conditions.txt'
    # initial_conditions_incl_M_df = load_initial_conditions_incl_M(initial_conditions_file_path)
    # print(initial_conditions_incl_M_df)

    initial_conditions_file_path = '../data/raw/L1_halo_initial_conditions.txt'
    initial_conditions_incl_M_df = load_initial_conditions_incl_M(initial_conditions_file_path)
    print(initial_conditions_incl_M_df)


    # differential_correction_file_path = '../data/raw/horizontal_L1_differential_correction.txt'
    # differential_correction_df = load_differential_corrections(differential_correction_file_path)
    # print(differential_correction_df)

    import math
    t_scale = 27.3217/(2*math.pi)
    x_scale = 0.3844e6
    v_scale = x_scale/t_scale/(24*3600)