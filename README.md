SPCE 5025 - Homework 4

Using previous homeworks as a base, and extending it to meet homework 4 requirements

Takes in a provided set of vectors, one for the position, and one for the velocity of a satellite, and now a step size. This then calculates the following:

# Problem 1
- Applies RK4 to estimate the change in position and velocity 10 times with a step size of 60
- Applies the f and g functions to the same initial position and velocity 10 times and compares to RK4 results
- Changes the units to Earth Radii and Hours and then re-performs RK4

# Problem 2
- Applies RK4 to estimate the change in position and velocity 10 times with a step size of 60
- Instead of using meters per second (Like in problem 1), instead converts units to Earth Radii per hour

# Output
Outputs tables with the step in seconds, the X, Y, and Z coordinates, and the Ẋ (aka XD), ẏ (aka YD), and Ż (aka ZD)
 