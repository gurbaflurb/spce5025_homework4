from keplarianElements import KeplerianElements
import keHelperFunctions

import math

import tabulate

def main():

    vectors_file = 'vectors.yaml'
    vector_data = keHelperFunctions.read_in_yaml(vectors_file)

    ke1 = KeplerianElements(vector_data['vectors'][f'vector1']['x_pos'],
                               vector_data['vectors'][f'vector1']['y_pos'],
                               vector_data['vectors'][f'vector1']['z_pos'],
                               vector_data['vectors'][f'vector1']['x_velocity'],
                               vector_data['vectors'][f'vector1']['y_velocity'],
                               vector_data['vectors'][f'vector1']['z_velocity'])

    inital_parameters = [[ke1.r_vector[0], ke1.r_vector[1], ke1.r_vector[2], ke1.r_dot_vector[0], ke1.r_dot_vector[1], ke1.r_dot_vector[2]]]

    print('Initial Parameters for Problems 1 and 2')
    print(tabulate.tabulate(inital_parameters, headers=['X', 'Y', 'Z', 'XD', 'YD', 'ZD'], floatfmt=".6f"))
    print()

    print(f'----- Problem 1 -----')

    new_pos = ke1.r_vector
    new_vel = ke1.r_dot_vector
    p1_table = []

    # Use RK4 to integrate the vector 10 steps using a_arrow = -(mu/r^3)*r_vector
    for i in range(0,10):
        new_pos, new_vel = keHelperFunctions.keplarian_rk4(new_pos, new_vel, vector_data['vectors'][f'vector1']['step_size'], ke1.mu)
        current_step = vector_data['vectors'][f'vector1']['step_size']*(i+1)
        table_entry = [current_step, new_pos[0], new_pos[1], new_pos[2], new_vel[0], new_vel[1], new_vel[2]]
        p1_table.append(table_entry)

    print('RK4 Estimated Positions and Velocities')
    print(tabulate.tabulate(p1_table, headers=['Step (Seconds)', 'X', 'Y', 'Z', 'XD', 'YD', 'ZD'], floatfmt=".8f"))
    print()

    # Check your work: compare RK4 results using f and g
    # We should get ~same~ results
    # Use delta_E approach as in Exam 1
    
    cur_nu = ke1.nu
    cur_pos = ke1.r_vector
    cur_vel = ke1.r_dot_vector

    fg_table = []
    diff_table = []

    for i in range(0, 10):
        loc, perigee_passes = ke1.determine_location_after_n_seconds(vector_data['vectors'][f'vector1']['step_size'], math.degrees(cur_nu))
        f = ke1.determine_f(math.degrees(cur_nu), math.degrees(loc))
        g = ke1.determine_g(math.degrees(cur_nu), math.degrees(loc))
        g_dot = ke1.determine_g_dot(math.degrees(loc))
        f_dot = ke1.determine_f_dot_new(f, g, g_dot)

        new_pos = ke1.determine_arbitrary_new_position(f, g, cur_pos, cur_vel)
        new_vel = ke1.determine_arbitrary_new_velocity(f_dot, g_dot, cur_pos, cur_vel)

        current_step = vector_data['vectors'][f'vector1']['step_size']*(i+1)
        table_entry = [current_step, new_pos[0], new_pos[1], new_pos[2], new_vel[0], new_vel[1], new_vel[2]]

        
        
        
        
        
        

        diff_entry = [current_step,
                      new_pos[0] - p1_table[i][1],
                      new_pos[1] - p1_table[i][2],
                      new_pos[2] - p1_table[i][3],
                      new_vel[0] - p1_table[i][4],
                      new_vel[1] - p1_table[i][5],
                      new_vel[2] - p1_table[i][6]]

        fg_table.append(table_entry)
        diff_table.append(diff_entry)

        cur_nu = loc
        cur_pos = new_pos
        cur_vel = new_vel


    print('f and g Positions and Velocities')
    print(tabulate.tabulate(fg_table, headers=['Step (Seconds)', 'X', 'Y', 'Z', 'XD', 'YD', 'ZD'], floatfmt=".8f"))
    print()

    print('Differences table')
    print(tabulate.tabulate(diff_table, headers=['Step (Seconds)', 'X', 'Y', 'Z', 'XD', 'YD', 'ZD'], floatfmt=".8f"))
    print()

    
    # Redo RK4 integration using units of Earth Radii and Hours
    # Dont forget to change the units of r_vector and r_dot_vector and mu
    # Radius of Earth is given as 6378137 meters
    # Seconds to hours -> 1s * 1m/60s * 1hr/60m -> 3600 seconds in an hour

    # Mu is in meters cubed per second squared, converting to Earth radii per hour
    converted_mu = ke1.mu * (math.pow(3600, 2)/math.pow(6378137, 3))

    # Each position coordinate is in meters, so we need to convert to Earth Radii
    er_pos = ke1.r_vector / 6378137

    # Each velocity component is in meters per second, so we need to convert to Earth Radii per hour
    er_vel = ke1.r_dot_vector * (3600/6378137)

    # Step is in seconds, so we need to convert it to hours
    converted_step = vector_data['vectors'][f'vector1']['step_size']/3600

    p2_table = []

    for i in range(0,10):
        er_pos, er_vel = keHelperFunctions.keplarian_rk4(er_pos, er_vel, converted_step, converted_mu)
        current_step = vector_data['vectors'][f'vector1']['step_size']*(i+1)
        table_entry = [current_step, er_pos[0], er_pos[1], er_pos[2], er_vel[0], er_vel[1], er_vel[2]]
        p2_table.append(table_entry)

    print('Earth Radii Converted position and velocities')
    print(tabulate.tabulate(p2_table, headers=['Step (Seconds)', 'X', 'Y', 'Z', 'XD', 'YD', 'ZD'], floatfmt=".8f"))
    
    print()

    print(f'----- Problem 2 -----')
    # Use initial ECI vector and same step size as in problem 1
    # Use RK4 to integrate 10 steps using the given equation
    new_pos = ke1.r_vector
    new_vel = ke1.r_dot_vector
    table = []

    for i in range(0, 10):
        new_pos, new_vel = keHelperFunctions.keplarian_rk4_oblate_earth(new_pos, new_vel, vector_data['vectors'][f'vector1']['step_size'])
        current_step = vector_data['vectors'][f'vector1']['step_size']*(i+1)
        table_entry = [current_step, new_pos[0], new_pos[1], new_pos[2], new_vel[0], new_vel[1], new_vel[2]]
        table.append(table_entry)

    print('RK4 Estimated Positions and Velocities accounting for Earths Oblateness')
    print(tabulate.tabulate(table, headers=['Step (Seconds)', 'X', 'Y', 'Z', 'XD', 'YD', 'ZD'], floatfmt=".8f"))



if __name__ == '__main__':
    main()
