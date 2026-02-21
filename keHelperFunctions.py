import yaml
import math

import numpy as np


# Read in a yaml that has all the initial vectors for position and velocity
def read_in_yaml(file_name):
    with open(file_name, 'r') as f:
        data = yaml.load(f.read(), Loader=yaml.SafeLoader)
        return data



def convert_arbitrary_perifocal_to_eci(a, e, inclination, raan, aop, nu) -> tuple:
        '''Converts manually provided perifocal values to ECI coordinates. Uses radians and not degrees. Passing in degrees will mess up the calculation'''

        # Hard coded because a nice solution for this exists (Pulled from slide 74)
        x = [ math.cos(raan)*math.cos(aop) - math.sin(raan)*math.sin(aop)*math.cos(inclination), -math.cos(raan)*math.sin(aop) - math.sin(raan)*math.cos(aop)*math.cos(inclination), math.sin(raan)*math.sin(inclination)]
        y = [ math.sin(raan)*math.cos(aop)+math.cos(raan)*math.sin(aop)*math.cos(inclination), -math.sin(raan)*math.sin(aop)+math.cos(raan)*math.cos(aop)*math.cos(inclination), -math.cos(raan)*math.sin(inclination)]
        z = [ math.sin(inclination)*math.sin(aop), math.sin(inclination)*math.cos(aop), math.cos(inclination)]

        return (x, y, z)



def find_arbitrary_position_and_velocity_vector(a: float, eccentricity: float, nu: float) -> tuple:
    '''Returns a tuple in the form position_vector, velocity_vector'''
    mu = 398600441800000.0 # From WGS84
    
    perifocal = a*(1.0-math.pow(eccentricity, 2))
    radius = perifocal/(1.0 + eccentricity*math.cos(nu))

    return ([radius*math.cos(nu), radius*math.sin(nu), 0.0],
            [-math.sqrt(mu/perifocal)*math.sin(nu), math.sqrt(mu/perifocal)*(eccentricity+math.cos(nu)), 0.0])



def keplarian_rk4(r_vector, r_dot_vector, step, mu):
     '''Takes in the position vector (r_vector), velocity vector (r_dot_vector), r (norm of r_vector), step (Size of step between each round), and mu (Probably from WGS 84).
     Mu is typically provided as meters cubed over seconds squared.'''

     # Commented out sections can help to bug things 

     # Given on slide 22
     r0norm = np.linalg.norm(r_vector)
     rd_a_pt = (-mu/math.pow(r0norm, 3))*r_vector

     r_a = r_vector + (step/2)*r_dot_vector
     rd_a = r_dot_vector + (step/2)*rd_a_pt
     k1 = [r_dot_vector, rd_a_pt]

#      print(f'k1  : {k1[0]}')
#      print(f'k1 a: {k1[1]}')
#      print(f'ra: {r_a}')
#      print(f'va: {rd_a}')
#      print(f'Norm of r: {r0norm}')
#      print()

     # Given on slide 23
     ranorm = np.linalg.norm(r_a)
     rd_b_pt = (-mu/math.pow(ranorm, 3))*r_a

     r_b = r_vector+(step/2)*rd_a
     rd_b = r_dot_vector+(step/2)*rd_b_pt
     k2 = [rd_a, rd_b_pt]

#      print(f'k2  : {k2[0]}')
#      print(f'k2 a: {k2[1]}')
#      print(f'ra: {r_b}')
#      print(f'va: {rd_b}')
#      print(f'Norm of r: {ranorm}')
#      print()

     # Given on slide 24
     rbnorm = np.linalg.norm(r_b)
     rd_c_pt = (-mu/math.pow(rbnorm, 3))*r_b

     r_c = r_vector+(step)*rd_b # These are wrong on slide 24, we multiply by the step instead of step / 2
     rd_c = r_dot_vector+(step)*rd_c_pt # These are wrong on slide 24, we multiply by the step instead of step / 2
     k3 = [rd_b, rd_c_pt]

#      print(f'k3  : {k3[0]}')
#      print(f'k3 a: {k3[1]}')
#      print(f'rc: {r_c}')
#      print(f'vc: {rd_c}')
#      print(f'Norm of r: {rbnorm}')
#      print()


     # Given on slide 25
     k4 = [rd_c, (-mu/math.pow(np.linalg.norm(r_c), 3))*r_c]
#      print(f'k4  : {k4[0]}')
#      print(f'k4 a: {k4[1]}')
#      print()

     step_position_solution = r_vector + step*( ((k1[0])/6) + ((k2[0])/3) + ((k3[0])/3) + ((k4[0])/6) )

     step_velocity_solution = r_dot_vector + step*( ((k1[1])/6) + ((k2[1])/3) + ((k3[1])/3) + ((k4[1])/6) )

     return (step_position_solution, step_velocity_solution)



def keplarian_rk4_oblate_earth(r_vector, r_dot_vector, step):
     '''Takes into account the Earth is Oblate. Takes in the position vector (r_vector), velocity vector (r_dot_vector), r (norm of r_vector), and step (Size of step between each round)
     Mu in this case is not provided by the user but instead defined by the WGS 84 standard for this specific implementation.'''

     # Defined by WGS 84
     earth_radius = 6378137 # This value is in meters
     earth_j2 = 0.00108262998905194
     mu = 398600441800000 # This value is in meters cubed per second squared

     # Commented out sections can help to bug things 

     # K1
     # Given on slide 52
     r = math.sqrt(math.pow(r_vector[0], 2) + math.pow(r_vector[1], 2) + math.pow(r_vector[2], 2))

     a_pt_s = r_vector[2]/r
     a_pt_1 = -(mu/math.pow(r, 3))
     a_pt_2 = (1+((3*earth_j2)/2)*math.pow(earth_radius/r, 2)*(1-5*math.pow(a_pt_s, 2)))*r_vector
     
     a = a_pt_1 * a_pt_2

     r_a = r_vector + (step/2)*r_dot_vector
     rd_a = r_dot_vector + (step/2)*a
     k1 = [r_dot_vector, a]

#      print(f'k1  : {k1[0]}')
#      print(f'k1 a: {k1[1]}')
#      print(f'ra: {r_a}')
#      print(f'va: {rd_a}')
#      print(f'Norm of r: {r0norm}')
#      print()

     # K2
     r = math.sqrt(math.pow(r_a[0], 2) + math.pow(r_a[1], 2) + math.pow(r_a[2], 2))

     b_pt_s = r_a[2]/r
     b_pt_1 = -(mu/math.pow(r, 3))
     b_pt_2 = (1+((3*earth_j2)/2)*math.pow(earth_radius/r, 2)*(1-5*math.pow(b_pt_s, 2)))*r_a
     
     a = b_pt_1 * b_pt_2

     r_b = r_vector+(step/2)*rd_a
     rd_b = r_dot_vector+(step/2)*a
     k2 = [rd_a, a]

#      print(f'k2  : {k2[0]}')
#      print(f'k2 a: {k2[1]}')
#      print(f'ra: {r_b}')
#      print(f'va: {rd_b}')
#      print(f'Norm of r: {ranorm}')
#      print()

     # K3
     r = math.sqrt(math.pow(r_b[0], 2) + math.pow(r_b[1], 2) + math.pow(r_b[2], 2))

     c_pt_s = r_b[2]/r
     c_pt_1 = -(mu/math.pow(r, 3))
     c_pt_2 = (1+((3*earth_j2)/2)*math.pow(earth_radius/r, 2)*(1-5*math.pow(c_pt_s, 2)))*r_b
     
     a = c_pt_1 * c_pt_2

     r_c = r_vector+(step)*rd_b # These are wrong on slide 24, we multiply by the step instead of step / 2
     rd_c = r_dot_vector+(step)*a # These are wrong on slide 24, we multiply by the step instead of step / 2
     k3 = [rd_b, a]

#      print(f'k3  : {k3[0]}')
#      print(f'k3 a: {k3[1]}')
#      print(f'rc: {r_c}')
#      print(f'vc: {rd_c}')
#      print(f'Norm of r: {rbnorm}')
#      print()

     # K4
     # Given on slide 25
     r = math.sqrt(math.pow(r_c[0], 2) + math.pow(r_c[1], 2) + math.pow(r_c[2], 2))

     d_pt_s = r_c[2]/r
     d_pt_1 = -(mu/math.pow(r, 3))
     d_pt_2 = (1+((3*earth_j2)/2)*math.pow(earth_radius/r, 2)*(1-5*math.pow(d_pt_s, 2)))*r_c
     
     a = d_pt_1 * d_pt_2

     k4 = [rd_c, a]
#      print(f'k4  : {k4[0]}')
#      print(f'k4 a: {k4[1]}')
#      print()

     step_position_solution = r_vector + step*( ((k1[0])/6) + ((k2[0])/3) + ((k3[0])/3) + ((k4[0])/6) )

     step_velocity_solution = r_dot_vector + step*( ((k1[1])/6) + ((k2[1])/3) + ((k3[1])/3) + ((k4[1])/6) )

     return (step_position_solution, step_velocity_solution)
