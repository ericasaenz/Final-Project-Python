import math
import matplotlib.pyplot as plt


# This function gets a float input from the user and checks if it is valid
def get_float(prompt, min_value=None, allow_zero=True):
    # keep asking until the user enters a valid number
    while True:
        try:
            value = float(input(prompt))

            # check if value is below the allowed minimum
            if min_value is not None and value < min_value:
                print(f"Value must be at least {min_value}.")
                continue

            # check if zero is allowed
            if not allow_zero and value == 0:
                print("Value cannot be zero.")
                continue

            return value

        except ValueError:
            print("Invalid input. Please enter a number.")


# This is the logistic growth equation (rate of change of infection)
def disease_model(t, infected, rate, carrying_capacity):
    # returns dI/dt = rI(1 - I/K)
    return rate * infected * (1 - infected / carrying_capacity)


# Euler's Method: uses one slope per step to approximate the solution
def euler_method(rate, carrying_capacity, initial_infected, total_time, step_size):
    times = [0.0]                  # stores time values
    infected_vals = [initial_infected]  # stores infected values

    t = 0.0
    infected = initial_infected

    # loop until total simulation time is reached
    while t < total_time:
        # make sure we do not go past the final time
        step = min(step_size, total_time - t)

        # update infected value using Euler formula
        infected = infected + step * disease_model(t, infected, rate, carrying_capacity)

        # move time forward
        t += step

        # keep values within realistic bounds (0 to max plants)
        infected = max(0.0, min(infected, carrying_capacity))

        # store results
        times.append(t)
        infected_vals.append(infected)

    return times, infected_vals


# Runge-Kutta 4th Order Method (RK4): more accurate than Euler
def rk4_method(rate, carrying_capacity, initial_infected, total_time, step_size):
    times = [0.0]
    infected_vals = [initial_infected]

    t = 0.0
    infected = initial_infected

    while t < total_time:
        step = min(step_size, total_time - t)

        # compute 4 slopes (k1, k2, k3, k4)
        k1 = disease_model(t, infected, rate, carrying_capacity)
        k2 = disease_model(t + step / 2, infected + step * k1 / 2, rate, carrying_capacity)
        k3 = disease_model(t + step / 2, infected + step * k2 / 2, rate, carrying_capacity)
        k4 = disease_model(t + step, infected + step * k3, rate, carrying_capacity)

        # combine slopes to update infected value
        infected = infected + (step / 6) * (k1 + 2 * k2 + 2 * k3 + k4)

        t += step

        # keep within realistic limits
        infected = max(0.0, min(infected, carrying_capacity))

        times.append(t)
        infected_vals.append(infected)

    return times, infected_vals


# Exact solution of the logistic equation 
# I(t) = K / (1 + A * e^(-rt))
def exact_solution(rate, carrying_capacity, initial_infected, times):
    values = []

    # if no infection at start, solution stays zero
    if initial_infected == 0:
        return [0.0 for _ in times]

    # compute constant A
    a_value = (carrying_capacity - initial_infected) / initial_infected

    # compute exact value at each time
    for t in times:
        infected = carrying_capacity / (1 + a_value * math.exp(-rate * t))
        values.append(infected)

    return values


# Prints a table comparing Euler, RK4, and exact solution
def print_results_table(times, euler_vals, rk4_vals, exact_vals):
    print("\nResults Table")
    print(f"{'Time':<12}{'Euler':<16}{'RK4':<16}{'Exact':<16}")

    for t, e_val, r_val, ex_val in zip(times, euler_vals, rk4_vals, exact_vals):
        print(f"{t:<12.2f}{e_val:<16.4f}{r_val:<16.4f}{ex_val:<16.4f}")


# Computes final error of each method compared to exact solution
def compute_final_errors(euler_vals, rk4_vals, exact_vals):
    euler_error = abs(euler_vals[-1] - exact_vals[-1])
    rk4_error = abs(rk4_vals[-1] - exact_vals[-1])
    return euler_error, rk4_error


# Plots all three solutions on one graph
def plot_results(times, euler_vals, rk4_vals, exact_vals):
    plt.figure(figsize=(10, 6))

    plt.plot(times, euler_vals, marker='o', label="Euler Method")
    plt.plot(times, rk4_vals, marker='s', label="RK4 Method")
    plt.plot(times, exact_vals, linestyle='--', label="Exact Solution")

    plt.title("Plant Disease Spread Simulation")
    plt.xlabel("Time")
    plt.ylabel("Number of Infected Plants")

    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()


# Main program that runs everything
def main():
    print("Plant Disease Spread Simulation")
    print("This program compares Euler's Method and RK4 for a disease spread model.\n")

    # get user inputs
    initial_infected = get_float("Enter the initial number of infected plants: ", min_value=0)
    infection_rate = get_float("Enter the infection rate: ", min_value=0)
    total_time = get_float("Enter the total simulation time: ", min_value=0, allow_zero=False)
    step_size = get_float("Enter the step size: ", min_value=0, allow_zero=False)
    carrying_capacity = get_float("Enter the total number of plants (maximum infected possible): ", min_value=0, allow_zero=False)

    # check for invalid case
    if initial_infected > carrying_capacity:
        print("Initial infected plants cannot be greater than total number of plants.")
        return

    # run Euler method
    times_euler, euler_vals = euler_method(
        infection_rate,
        carrying_capacity,
        initial_infected,
        total_time,
        step_size,
    )

    # run RK4 method
    times_rk4, rk4_vals = rk4_method(
        infection_rate,
        carrying_capacity,
        initial_infected,
        total_time,
        step_size,
    )

    # compute exact solution
    exact_vals = exact_solution(infection_rate, carrying_capacity, initial_infected, times_euler)

    # print results table
    print_results_table(times_euler, euler_vals, rk4_vals, exact_vals)

    # compute and display errors
    euler_error, rk4_error = compute_final_errors(euler_vals, rk4_vals, exact_vals)

    print("\nFinal Error Comparison")
    print(f"Euler final absolute error: {euler_error:.6f}")
    print(f"RK4 final absolute error:   {rk4_error:.6f}")

    # compare accuracy
    if rk4_error < euler_error:
        print("RK4 was more accurate for this simulation.")
    elif euler_error < rk4_error:
        print("Euler was more accurate for this simulation.")
    else:
        print("Both methods had the same final error.")

    # plot graph
    plot_results(times_euler, euler_vals, rk4_vals, exact_vals)


# run the program
if __name__ == "__main__":
    main()