from keplarianElements import KeplerianElements

import keHelperFunctions

import numpy as np

def main():

    vectors_file = 'vectors.yaml'
    vector_data = keHelperFunctions.read_in_yaml(vectors_file)

    print(f'----- Problem 1 -----')
    ke1 = KeplerianElements(vector_data['vectors'][f'vector1']['x_pos'],
                               vector_data['vectors'][f'vector1']['y_pos'],
                               vector_data['vectors'][f'vector1']['z_pos'],
                               vector_data['vectors'][f'vector1']['x_velocity'],
                               vector_data['vectors'][f'vector1']['y_velocity'],
                               vector_data['vectors'][f'vector1']['z_velocity'])
    
    print(f'Initial Position: {ke1.r_vector}')
    print(f'Initial Velocity: {ke1.r_dot_vector}')
    print()

    new_pos = ke1.r_vector
    new_vel = ke1.r_dot_vector

    # Use RK4 to integrate the vector 10 steps using a_arrow = -(mu/r^3)*r_vector
    for i in range(0,10):
        new_pos, new_vel = keHelperFunctions.keplarian_rk4(new_pos, new_vel, vector_data['vectors'][f'vector1']['step_size'], ke1.mu)

        print(f'Position after step {i+1}: {new_pos}')
        print(f'Velocity after step {i+1}: {new_vel}')
        print()


    # Check your work: compare RK4 results using f and g
    # We should get ~same results
    # Use delta_E approach as in Exam 1
    f = ke1.determine_f()
    
    # Redo RK4 integration using units of Eearth Radii and Hours
    # Dont forget to change the units of r and mu
    

    print(f'----- Problem 2 -----')
    # Use initial ECI vector and same step size as in problem 1
    # Use RK4 to integrate 10 steps using the given equation




if __name__ == '__main__':
    main()
