from bayes_opt import BayesianOptimization, JSONLogger, Events
import subprocess
import json
import os
import pickle  # Import pickle to save the model
import matplotlib.pyplot as plt

# take case from the sim as a test case

# use simular geometry to ensure convergence

def cup_function(blank_radius):

    # target geometry
    ideal_height = 20
    ideal_radius = 33/2
    ideal_profile_radius = 5

    punch_profile_radius = 5
    punch_min_radius = 33/2
    die_profile_radius = 5
    die_punch_gap = 1.15
    punch_depth = 35

    abaqus_command = [
        r"C:\SIMULIA\Commands\abaqus.bat", "cae", "noGUI=stampy.py", "--",
        str(ideal_height), str(ideal_radius), str(ideal_profile_radius),
        str(punch_depth), str(punch_profile_radius), str(punch_min_radius),
        str(die_profile_radius), str(die_punch_gap), str(blank_radius)
    ]

    # Run the subprocess and capture height_diff from stdout
    result = subprocess.run(abaqus_command, capture_output=True, text=True, shell=True)

    # If Abaqus failed, return None
    if result.returncode != 0:
        print("Abaqus execution failed!")
        return None

    try:
        with open(r"C:\Users\kam97\OneDrive - University of Bath\Documents\build\stamping_op_build3\sim_output.json", "r") as f:
            data_out = json.load(f)

        height_diff = data_out["height difference"]

        if height_diff == None:
            height_diff = -100
        print(f"Simulation Complete, height difference: {height_diff}")
        return height_diff

    except Exception as exc:
        print("Could not parse sim_output.json:", exc)
        print("STDOUT:\n", result.stdout)
        print("STDERR:\n", result.stderr)
        return None

# Define a callback function to save the model after each optimization step.
def save_model_callback(optimizer):
    model_path = r"C:\Users\kam97\OneDrive - University of Bath\Documents\build\stamping_op_build3\bayes_model.pkl"
    try:
        with open(model_path, "wb") as f:
            pickle.dump(optimizer, f)
        print(f"Model saved to {model_path}")
    except Exception as e:
        print("Error saving model:", e)

def main():
    # Bounded region of parameter space
    pbounds = {
        "blank_radius": (40, 45), 
    }

    optimizer = BayesianOptimization(
        f=cup_function,
        pbounds=pbounds,    
        random_state=1,
    )

    # Subscribe the JSON logger
    logger = JSONLogger(path=r"C:\Users\kam97\OneDrive - University of Bath\Documents\build\stamping_op_build3\op_logs.log")

    optimizer.subscribe(Events.OPTIMIZATION_STEP, logger)

    optimizer.maximize(init_points=10, n_iter=20)

    save_model_callback(optimizer)

    op_path = r"C:\Users\kam97\OneDrive - University of Bath\Documents\build\stamping_op_build3"

    # optimisation convergence plot
    # plot_convergence(op_path)



def mesh_conv_final_max():

    import pandas as pd

    sim_numbers = [1]
    data_list = ["thickness_results.csv","strain_eq_results.csv", "springback_displacement_results.csv"]
    var_list = ["SectionThickness","PEEQ", "U_mag"]

    main_sim_output =  r"C:\Users\kam97\OneDrive - University of Bath\Documents\build\build_temp295"

    for sim in sim_numbers:

        print(sim)

        sim_path = rf"{main_sim_output}\sim_{sim}"
        # print(sim_path)

        for (data, var) in zip(data_list, var_list):

            data_path = rf"{sim_path}\{data}"

            # print(data_path)

            data_temp = pd.read_csv(data_path)

            max_val = data_temp[var].max()

            print(var, max_val)

def mesh_conv_final_min():

    import pandas as pd

    sim_numbers = [1]
    data_list = ["thickness_results.csv", "reaction_force_results.csv"]
    var_list = ["SectionThickness", "rf"]

    main_sim_output =  r"C:\Users\kam97\OneDrive - University of Bath\Documents\build\build_temp295"

    for sim in sim_numbers:

        print(sim)

        sim_path = rf"{main_sim_output}\sim_{sim}"
        # print(sim_path)

        for (data, var) in zip(data_list, var_list):

            data_path = rf"{sim_path}\{data}"

            # print(data_path)

            data_temp = pd.read_csv(data_path)

            max_val = data_temp[var].min()

            print(var, max_val)


def plot_exp():

    import numpy as np
    import matplotlib.pyplot as plt
    import pandas as pd

    # Load data from file (adjust the filename as needed)
    data = np.loadtxt('punch_rf_experimental.txt', delimiter=',')

    # Split data into x and y arrays
    x = data[:, 0]
    y = data[:, 1]

    # Create the plot
    plt.plot(x, y, marker='o', linestyle='--', color='b', label='Experimental Data')


    cwd = r"C:\Users\kam97\OneDrive - University of Bath\Documents\build\build_temp263\sim_16"

    rf_file_path = f"{cwd}/reaction_force_results.csv"

    rf_df = pd.read_csv(rf_file_path)
    
    rf_df["displacement"] = rf_df["displacement"] * -1
    rf_df["rf"] = rf_df["rf"] * -1e-3

    plt.plot(rf_df["displacement"], rf_df["rf"] * 4, color="r", label="Simulation Data") # for quarter model
    plt.xlabel("Displacement (mm)")
    plt.ylabel("Reaction Force (KN)")
    plt.legend()
    plt.grid(True)

    plt.savefig("rf_plot_exp_sim.png")




def plot_thickness_exp():

    import numpy as np
    import matplotlib.pyplot as plt
    import pandas as pd

    exp_data = np.loadtxt('exp_thickness_profile.txt', delimiter=',')
    sim_data = np.loadtxt('thickness_sim_data.txt', delimiter=',')

    plt.plot(exp_data[:, 0], exp_data[:, 1], marker="^", linestyle='--', color='b', label='Experimental Data')
    plt.plot(sim_data[:, 0], sim_data[:, 1], linestyle='-', color='r', label='Sim Data')

    plt.xlabel("Curvilinear distance from cup's center (mm)")
    plt.xlim([9, 35])
    plt.ylabel("Thickness (mm)")
    plt.legend()
    plt.grid(True)
    
    plt.savefig("thickness_exp_sim.png")


def plot_convergence():

    # Path to your log file containing one JSON record per line
    log_file = "op_logs.log"

    # Lists to store iteration number, target values, and elapsed times
    iterations = []
    targets = []
    elapsed_times = []

    # Read and parse the log file
    with open(log_file, 'r') as f:
        for i, line in enumerate(f):
            data = json.loads(line)
            iterations.append(i + 1)
            targets.append(data["target"])
            elapsed_times.append(data["datetime"]["elapsed"])

    # Compute best-so-far values: since your target is negative absolute difference,
    # the best value is the one closest to zero (i.e. the maximum)
    best_so_far = []
    current_best = float('-inf')
    for t in targets:
        if t > current_best:
            current_best = t
        best_so_far.append(current_best)

    # Plot best-so-far target value vs. iteration number
    plt.figure(figsize=(8, 6))
    plt.plot(iterations, best_so_far, marker='o', linestyle='-', color='blue')
    plt.xlabel('Iteration')
    plt.ylabel('Best Cup height error (mm)')
    # plt.yscale("log")
    plt.title('Convergence Plot (Iterations)')
    plt.grid(True)
    plt.savefig("convergence_iteratons.png")

    # Plot best-so-far target value vs. elapsed time (in seconds)
    plt.figure(figsize=(8, 6))
    plt.plot(elapsed_times, best_so_far, marker='o', linestyle='-', color='red')
    plt.xlabel('Elapsed Time (seconds)')
    plt.ylabel('Best cup height error (mm)')
    # plt.yscale("log")
    plt.title('Convergence Plot (Elapsed Time)')
    plt.grid(True)
    plt.savefig("convergence_time.png")

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec

def posterior(optimizer, grid):
    # Compute GP predictions on a grid: mean and standard deviation.
    mu, sigma = optimizer._gp.predict(grid, return_std=True)
    return mu, sigma

def plot_gp_default():
    # Define a dense grid over the search space for blank_radius (from 25 to 33)
    grid = np.linspace(25, 33, 10000).reshape(-1, 1)
    
    # Change to your simulation output directory
    os.chdir(r"C:\Users\kam97\OneDrive - University of Bath\Documents\build\stamping_op_build2")
    op_path = "bayes_model.pkl"
    
    # Load the optimizer object
    with open(op_path, "rb") as f:
        optimizer = pickle.load(f)
    
    independant_var = "blank_radius"
    dependant_var = "cup_height_error"  # or "cup_height", adjust as needed
    
    # Create the figure with two subplots
    fig = plt.figure(figsize=(16, 10))
    steps = len(optimizer.res)
    fig.suptitle('GP Prediction and Acquisition Utility After {} Steps'.format(steps),
                 fontsize=15)
    
    gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1])
    ax_gp = plt.subplot(gs[0])
    ax_acq = plt.subplot(gs[1])
    
    # Extract observed x and target values from the optimizer results
    x_obs = np.array([[res["params"][independant_var]] for res in optimizer.res])
    y_obs = np.array([res["target"] for res in optimizer.res])
    
    # Get GP predictions (mean and standard deviation) on the grid
    mu, sigma = posterior(optimizer, grid)
    
    # Compute the target (analytical cup height) on the grid using your analytical function
    target_vals = np.array([-abs(cup_height_analyt(x) - 20) for x in grid.flatten()])
    
    # Plot the target function (analytical) as a solid blue line
    ax_gp.plot(grid, target_vals, linewidth=3, color='blue', label='Target Function (Analytical)')
    
    # Plot the GP surrogate prediction, observations, and 95% confidence interval
    ax_gp.plot(grid, mu, '--', color='k', label='GP Prediction')
    ax_gp.scatter(x_obs.flatten(), y_obs, marker='D', s=80, color='r', label='Observations')
    ax_gp.fill_between(grid.flatten(),
                       mu.flatten() - 1.96 * sigma.flatten(),
                       mu.flatten() + 1.96 * sigma.flatten(),
                       alpha=0.6, fc='c', ec='none', label='95% Confidence Interval')
    ax_gp.set_xlim(grid.min(), grid.max())
    ax_gp.set_xlabel(independant_var, fontsize=15)
    ax_gp.set_ylabel(dependant_var, fontsize=15)
    
    # Plot the acquisition function using the optimizer's default acquisition function
    utility = optimizer.acquisition_function._get_acq(gp=optimizer._gp)(grid)
    ax_acq.plot(grid.flatten(), utility, label='Acquisition Utility', color='purple')
    
    # Mark the point with the highest acquisition utility as the next best guess
    next_index = np.argmax(utility)
    next_x = grid.flatten()[next_index]
    ax_acq.plot(next_x, utility[next_index], '*', markersize=15,
                label='Next Best Guess', markerfacecolor='gold',
                markeredgecolor='k', markeredgewidth=1)
    ax_acq.set_xlim(grid.min(), grid.max())
    ax_acq.set_xlabel(independant_var, fontsize=15)
    ax_acq.set_ylabel('Utility', fontsize=15)
    
    ax_gp.legend(loc='upper left')
    ax_acq.legend(loc='upper left')
    
    plt.savefig("gp_check.png")

def plot_sampling_history():

    independant_var = "blank_radius"

    # Change to your simulation output directory
    os.chdir(r"C:\Users\kam97\OneDrive - University of Bath\Documents\build\stamping_op_build2")
    op_path = "bayes_model.pkl"
    
    # Load the optimizer object
    with open(op_path, "rb") as f:
        optimizer = pickle.load(f)


    # Scatter plot of sample locations vs. iteration number
    x_obs = [res["params"][independant_var] for res in optimizer.res]
    iterations = np.arange(1, len(x_obs) + 1)
    plt.figure(figsize=(10, 4))
    plt.scatter(iterations, x_obs, color='green')
    plt.xlabel('Iteration')
    plt.ylabel(independant_var)
    plt.title('Sampling History')
    plt.grid(True)
    plt.savefig("sampling_history.png")

def cup_height_analyt(blank_radius):
    """
    Calculate the cup height based on area conservation, assuming constant thickness.
    
    The blank area is partitioned into:
      - The bottom area: a circle with radius (cup_radius - cup_fillet)
      - The fillet area: approximated as a quarter circle with radius cup_fillet
      - The wall area: the remaining area, which forms the vertical wall
    
    The lateral wall area is given by: A_wall = 2π * cup_radius * cup_height.
    Conserving area:
      π * blank_radius^2 = A_bottom + A_fillet + 2π * cup_radius * cup_height
    Solving for cup_height:
      cup_height = (π * blank_radius^2 - A_bottom - A_fillet) / (2π * cup_radius)
    
    Parameters:
        blank_radius (float): The radius of the initial blank.
    
    Returns:
        cup_height (float): The calculated height of the cup.
    """
    # Fixed cup design parameters
    cup_radius = 33 / 2     # Overall cup radius (e.g., 16.5)
    cup_fillet = 5          # Fillet width (or effective thickness for the fillet region)
    
    # Bottom: flat circular region with effective radius = cup_radius - cup_fillet
    bottom_radius = cup_radius - cup_fillet
    A_bottom = np.pi * bottom_radius**2
    
    # Fillet: assume the fillet spans a fixed fraction of the cup's circumference
    fillet_fraction = 0.25  # 25% of the cup circumference
    L_fillet = fillet_fraction * (2 * np.pi * cup_radius)  # Arc length of the fillet region
    A_fillet = cup_fillet * L_fillet  # Approximated as a rectangle: width * arc length
    
    # Blank area (conserved material area, assuming constant thickness)
    A_blank = np.pi * blank_radius**2
    
    # Remaining area forms the vertical wall
    A_wall = A_blank - A_bottom - A_fillet
    
    # The lateral wall area of a cylinder is: A_wall = 2π * cup_radius * cup_height
    cup_height = A_wall / (2 * np.pi * cup_radius)
    return cup_height


if __name__ == "__main__":
    

    # plot_gp_default()
    # plot_sampling_history()
    main()

    # mesh_conv_final_max()
    # print("##############")
    # mesh_conv_final_min()
    # plot_exp()

    # plot_thickness_exp()

    # os.chdir(r"C:\Users\kam97\OneDrive - University of Bath\Documents\build\stamping_op_build2")
    # plot_convergence()
