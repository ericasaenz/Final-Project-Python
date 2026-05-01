# Final-Project-Python
Plant Disease Spread Simulation

This project simulates how a disease spreads through a population of plants using a logistic growth model. It compares Euler’s Method and RK4 to the exact solution to evaluate their accuracy.

The user inputs values such as the initial number of infected plants, infection rate, total time, step size, and carrying capacity. The program then calculates how the infection changes over time, displays a results table, computes the final error, and generates a graph comparing all three methods.

To run the program, open a terminal and navigate the file. Then run the command:

python filename.py

After running the program, you will be prompted to enter the required input values. For example:

Enter the initial number of infected plants: 5  
Enter the infection rate: 0.4  
Enter the total simulation time: 20  
Enter the step size: 1  
Enter the total number of plants (maximum infected possible): 100

Once the inputs are provided, the program will output a table of results, display the final error comparison, and show a graph of the simulation. This example produces a smooth S-shaped growth curve and clearly shows the difference in accuracy between Euler’s Method and RK4.

This project demonstrates how numerical methods can be used to approximate solutions to differential equations and highlights the differences in accuracy between Euler’s Method and RK4.
