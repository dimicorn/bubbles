Run the main.py file.

The program will ask you: "Data or plot? "

If you input "data" at the start of the program, 
the program will create a Excel table (make sure to change the file_path at the constants.py file)
with some necessary data for all the cases with an intersection (the odd case (*) doesn't have an intersection).
Currently most the values that could be calculated are commented out, but that could be changed.
Be mind reading comments, not all values/functions are working correctly.

If you input "plot" at the start of the program, the program will ask you: "Velocity or pressure? ".

If you input "vel", the program will create and automatically save figures of V(lambda_c) and velocity slopes
for a specific value of gamma, k_rho and all the possible values of n_int.

If you input "pres", the program will create and automatically save scaled figures of P(lambda_c) for all the cases
that have an intersection. Intersection is determined mathematically. If the distance between the last point of
the ODE solution and the velocity slope is less than 0.001, we consider it as an intersection. There is a odd
case (*), where the ODE curve intersects with the velocity slope and goes even further, it is considered separately.

If there will be any updates, I will write them here.

Date: 13.09