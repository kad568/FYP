from bayes_opt import BayesianOptimization, JSONLogger, Events, acquisition
import subprocess
import json
import os
import pickle  # Import pickle to save the model
import matplotlib.pyplot as plt
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern, RBF, RationalQuadratic, WhiteKernel, ConstantKernel
from scipy.optimize import NonlinearConstraint
from skopt.space import Real, Space
from skopt.sampler import Lhs

# take case from the sim as a test case

# use simular geometry to ensure convergence

def cup_function(BHF, die_profile_radius):

    # target geometry
    ideal_height = 20
    ideal_radius = 33/2
    ideal_profile_radius = 5
    blank_radius = 30.321158122674316

    punch_die_gap = 1.15


    punch_depth = 35
    friction = 0.09

    abaqus_command = [
        r"C:\SIMULIA\Commands\abaqus.bat", "cae", "noGUI=stampy_process.py", "--",
        str(ideal_height), str(ideal_radius), str(ideal_profile_radius),
        str(punch_depth), str(BHF), str(friction),
        str(die_profile_radius), str(punch_die_gap), str(blank_radius)
    ]

    # Run the subprocess and capture height_diff from stdout
    result = subprocess.run(abaqus_command, capture_output=True, text=True, shell=True)

    # If Abaqus failed, return None
    if result.returncode != 0:
        print("Abaqus execution failed!")
        return None

    try:
        with open(r"C:\Users\kam97\OneDrive - University of Bath\Documents\build\process_optimisation2\sim_output.json", "r") as f:
            data_out = json.load(f)

        min_thickness = data_out["min_thickness"]
        max_thickness = data_out["max_thickness"]

        if (min_thickness == None) or (max_thickness > 1.2):
            min_thickness = 0
        print(f"Simulation Complete, min_thickness: {min_thickness}")
        return min_thickness

    except Exception as exc:
        print("Could not parse sim_output.json:", exc)
        print("STDOUT:\n", result.stdout)
        print("STDERR:\n", result.stderr)
        return None

# Define a callback function to save the model after each optimization step.
def save_model_callback(optimizer):
    model_path = r"C:\Users\kam97\OneDrive - University of Bath\Documents\build\process_optimisation2\bayes_model.pkl"
    try:
        with open(model_path, "wb") as f:
            pickle.dump(optimizer, f)
        print(f"Model saved to {model_path}")
    except Exception as e:
        print("Error saving model:", e)

def main():
    # Bounded region of parameter space
    pbounds = {
        "BHF": (0, 1),
        # "punch_die_gap": (0, 1),
        "die_profile_radius": (0, 1),
    }

    acquisition_function = acquisition.ExpectedImprovement(xi=0.0)

    kernel = ConstantKernel(1.0, (0.1, 10)) * Matern(length_scale=[1.0, 1.0], nu=2.5, length_scale_bounds=[(1e-2, 10.0), (1e-2, 10.0)]) \
                   + WhiteKernel(noise_level=1e-5, noise_level_bounds=(1e-8, 1e-3))

    
    gp = GaussianProcessRegressor(alpha=1e-06,
                                    kernel=kernel,
                                    n_restarts_optimizer=20,
                                    normalize_y=False,
                                    random_state=0)  # set a random state for reproducibility

    optimizer = BayesianOptimization(
        f=cup_function,
        acquisition_function=acquisition_function,
        pbounds=pbounds,    
        random_state=0,
    )

    optimizer._gp = gp

    # Subscribe the JSON logger
    logger = JSONLogger(path=r"C:\Users\kam97\OneDrive - University of Bath\Documents\build\process_optimisation2\op_logs.log")

    optimizer.subscribe(Events.OPTIMIZATION_STEP, logger)

    # 5) use scikit-optimize’s LHS to generate n_initial “init” points
    space = Space([
        Real(*pbounds["BHF"], name="BHF"),
        Real(*pbounds["die_profile_radius"],    name="die_profile_radius"),
    ])

    lhs = Lhs(lhs_type="classic", criterion="maximin")
    n_initial = 20
    lhs_points = lhs.generate(space.dimensions, n_samples=n_initial, random_state=0)

    # 6) register each LHS point with the BO engine
    for x in lhs_points:
        params = {"BHF": x[0], "die_profile_radius": x[1]}
        y = cup_function(**params)
        optimizer.register(params=params, target=y)


    optimizer.maximize(init_points=0, n_iter=60)

    save_model_callback(optimizer)



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

    os.chdir(r"C:\Users\kam97\OneDrive - University of Bath\Documents\build\stamping_op_build3")

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
    os.chdir(r"C:\Users\kam97\OneDrive - University of Bath\Documents\build\stamping_op_build3")
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
    os.chdir(r"C:\Users\kam97\OneDrive - University of Bath\Documents\build\stamping_op_build3")
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


import os
import json
import pickle
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec

def plot_convergence_2d():
    """
    Reads the optimization log file (one JSON record per line) where each record contains:
      - "target": negative RMSE (closer to zero is better),
      - "params": with keys "punch_min_radius" and "punch_profile_radius",
      - "datetime": with "elapsed" time.
      
    Since the target is negative RMSE, the best value is the maximum (i.e. the one closest to zero).
    This function plots:
      - Best-so-far target value vs. iteration number.
      - Best-so-far target value vs. elapsed time.
    """
    os.chdir(r"C:\Users\kam97\OneDrive - University of Bath\Documents\build\process_optimisation2")
    log_file = "op_logs.log"

    iterations = []
    targets = []
    elapsed_times = []

    with open(log_file, 'r') as f:
        for i, line in enumerate(f):
            data = json.loads(line)
            iterations.append(i + 1)
            targets.append(data["target"])  # negative RMSE values
            elapsed_times.append(data["datetime"]["elapsed"])

    # Since target is negative RMSE (closer to 0 is better), best-so-far is the maximum value.
    best_so_far = []
    current_best = float('-inf')
    for t in targets:
        if t > current_best:
            current_best = t
        best_so_far.append(current_best)

    # Plot best-so-far target vs. iteration number.
    plt.figure(figsize=(8, 6))
    plt.plot(iterations, best_so_far, marker='o', linestyle='-', color='blue')
    plt.xlabel('Iteration')
    plt.ylabel('Best Target (-RMSE)')
    plt.title('Convergence Plot (Iterations)')
    plt.grid(True)
    plt.savefig("convergence_iterations_2d.png")

    # Plot best-so-far target vs. elapsed time.
    plt.figure(figsize=(8, 6))
    plt.plot(elapsed_times, best_so_far, marker='o', linestyle='-', color='red')
    plt.xlabel('Elapsed Time (seconds)')
    plt.ylabel('Best Target (-RMSE)')
    plt.title('Convergence Plot (Elapsed Time)')
    plt.grid(True)
    plt.savefig("convergence_time_2d.png")


def plot_gp_2d():
    """
    Creates two contour plots over a 2D grid spanning the design space:
      - GP Mean Prediction,
      - GP Prediction Uncertainty (standard deviation).
      
    The design space is defined by:
      - "punch_profile_radius": ranging from 3 to 6,
      - "punch_min_radius": ranging from 14 to 19.
      
    Observations are overlaid on both plots.
    """
    # Define grid ranges based on the given bounds
    r1 = np.linspace(0, 1, 100)   # punch_profile_radius
    r2 = np.linspace(0, 1, 100)  # punch_min_radius
    R1, R2 = np.meshgrid(r1, r2)
    grid = np.vstack([R1.ravel(), R2.ravel()]).T

    os.chdir(r"C:\Users\kam97\OneDrive - University of Bath\Documents\build\process_optimisation2")
    op_path = "bayes_model_norm.pkl"
    with open(op_path, "rb") as f:
        optimizer = pickle.load(f)

    # Extract observed sample locations and targets
    # X-axis: punch_profile_radius, Y-axis: punch_min_radius.
    x_obs = np.array([[res["params"]["punch_profile_radius"], res["params"]["punch_min_radius"]] 
                      for res in optimizer.res])
    y_obs = np.array([res["target"] for res in optimizer.res])

    # Get GP predictions: mean and standard deviation on the grid.
    mu, sigma = optimizer._gp.predict(grid, return_std=True)
    MU = mu.reshape(R1.shape)
    SIGMA = sigma.reshape(R1.shape)

    # Create side-by-side contour plots for GP mean and uncertainty.
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

    cp1 = ax1.contourf(R1, R2, MU, levels=50, cmap='viridis')
    cb = fig.colorbar(cp1, ax=ax1)
    cb.set_label('Normalized RMSE', labelpad=15)
    ax1.scatter(x_obs[:, 0], x_obs[:, 1], color='red', marker='D', s=50, label='Observations')
    ax1.set_title('GP Mean Prediction (RMSE)')
    ax1.set_xlabel('Punch Profile Radius')
    ax1.set_ylabel('Punch Min Radius')
    ax1.legend()

    cp2 = ax2.contourf(R1, R2, SIGMA, levels=50, cmap='inferno')
    cb2 = fig.colorbar(cp2, ax=ax2)
    cb2.set_label('Std Dev', labelpad=15)
    ax2.scatter(x_obs[:, 0], x_obs[:, 1], color='red', marker='D', s=50, label='Observations')
    ax2.set_title('GP Prediction Uncertainty (Std Dev)')
    ax2.set_xlabel('Punch Profile Radius')
    ax2.set_ylabel('Punch Min Radius')
    ax2.legend()

    plt.suptitle('GP Prediction and Uncertainty After {} Steps'.format(len(optimizer.res)), fontsize=16)
    plt.savefig("gp_2d.png")


def plot_acq_2d():
    """
    Plots the 2D acquisition function as a contour plot over the design space.
    The next best guess (i.e. the point with the highest acquisition utility) is marked,
    and the observed sample points are overlaid.
    """
    r1 = np.linspace(0, 1, 100)   # punch_profile_radius
    r2 = np.linspace(0, 1, 100)  # punch_min_radius
    R1, R2 = np.meshgrid(r1, r2)
    grid = np.vstack([R1.ravel(), R2.ravel()]).T

    os.chdir(r"C:\Users\kam97\OneDrive - University of Bath\Documents\build\process_optimisation2")
    op_path = "bayes_model_norm.pkl"
    with open(op_path, "rb") as f:
        optimizer = pickle.load(f)

    # Compute the acquisition function values on the grid.
    utility = optimizer.acquisition_function._get_acq(gp=optimizer._gp)(grid)
    UTILITY = utility.reshape(R1.shape)

    fig, ax = plt.subplots(figsize=(8, 6))
    cp = ax.contourf(R1, R2, UTILITY, levels=50, cmap='plasma')
    cb = fig.colorbar(cp, ax=ax)
    cb.set_label('Utility', labelpad=15)

    # Mark the next best guess (maximum acquisition utility).
    next_index = np.argmax(utility)
    next_x = grid[next_index, 0]
    next_y = grid[next_index, 1]
    ax.plot(next_x, next_y, '*', markersize=15, markerfacecolor='gold',
            markeredgecolor='k', markeredgewidth=1, label='Next Best Guess')

    # Overlay observed sample points.
    x_obs = np.array([[res["params"]["punch_profile_radius"], res["params"]["punch_min_radius"]] 
                      for res in optimizer.res])
    ax.scatter(x_obs[:, 0], x_obs[:, 1], color='red', marker='D', s=50, label='Observations')

    ax.set_xlabel('Punch Profile Radius')
    ax.set_ylabel('Punch Min Radius')
    ax.set_title('Acquisition Function Utility')
    ax.legend()

    plt.savefig("acquisition_2d.png")


def plot_sampling_history_2d():
    """
    Creates a 2D scatter plot of the sample locations over the design space.
    Points are color-coded by iteration number to visualize the sampling history.
    """
    os.chdir(r"C:\Users\kam97\OneDrive - University of Bath\Documents\build\process_optimisation2")
    op_path = "bayes_model_norm.pkl"
    with open(op_path, "rb") as f:
        optimizer = pickle.load(f)

    x_obs = np.array([[res["params"]["punch_profile_radius"], res["params"]["punch_min_radius"]] 
                      for res in optimizer.res])
    iterations = np.arange(1, len(x_obs) + 1)

    plt.figure(figsize=(8, 6))
    sc = plt.scatter(x_obs[:, 0], x_obs[:, 1], c=iterations, cmap='viridis', s=80)
    plt.xlabel('Punch Profile Radius')
    plt.ylabel('Punch Min Radius')
    plt.title('Sampling History (Color-coded by Iteration)')
    plt.colorbar(sc, label='Iteration')
    plt.grid(True)
    plt.savefig("sampling_history_2d.png")
    
def plot_sampling_history_rmse_2d():
    """
    Creates a 2D scatter plot of the sample locations over the design space.
    Points are color-coded by the RMSE value (obtained as the absolute value of the target)
    to visualize the performance.
    """
    import os
    import pickle
    import numpy as np
    import matplotlib.pyplot as plt

    os.chdir(r"C:\Users\kam97\OneDrive - University of Bath\Documents\build\process_optimisation2")
    op_path = "bayes_model_norm.pkl"
    with open(op_path, "rb") as f:
        optimizer = pickle.load(f)

    # Extract the design parameters from the optimizer results
    x_obs = np.array([[res["params"]["punch_profile_radius"], res["params"]["punch_min_radius"]] 
                      for res in optimizer.res])
    
    # Since target is defined as -RMSE, we take the absolute value to get the actual RMSE
    rmse_values = np.array([abs(res["target"]) for res in optimizer.res])
    
    plt.figure(figsize=(8, 6))
    sc = plt.scatter(x_obs[:, 0], x_obs[:, 1], c=rmse_values, cmap='viridis', s=80)
    plt.xlabel('Punch Profile Radius')
    plt.ylabel('Punch Min Radius')
    plt.title('Sampling History (Color-coded by RMSE)')
    plt.colorbar(sc, label='RMSE')
    plt.grid(True)
    plt.savefig("sampling_history_rmse_2d.png")


def op_gp():

    os.chdir(r"C:\Users\kam97\OneDrive - University of Bath\Documents\build\process_optimisation2")
    op_path = "bayes_model.pkl"

    with open(op_path, "rb") as f:
        optimizer = pickle.load(f)

    gp = optimizer._gp

    # print(gp)

    # print("Optimized kernel:", gp.kernel_)
    # print("Log marginal likelihood:", gp.log_marginal_likelihood())
    # print("##################### normalized #################")
    # print("Training inputs (X_train):", gp.X_train_)

    x_data = gp.X_train_

    punch_min_radius = [data[0] for data in x_data]
    punch_fillet_radius = [data[1] for data in x_data]

    max_punch_min_radius = max(punch_min_radius)
    min_punch_min_radius = min(punch_min_radius)

    max_punch_fillet_radius = max(punch_fillet_radius)
    min_punch_fillet_radius = min(punch_fillet_radius)

    x_data = []

    for radius, fillet in zip(punch_min_radius, punch_fillet_radius):

        punch_min_radius_norm = (radius - min_punch_min_radius) / (max_punch_min_radius - min_punch_min_radius)

        punch_filet_radius_norm = (fillet - min_punch_fillet_radius) / (max_punch_fillet_radius - min_punch_fillet_radius)

        data = [punch_filet_radius_norm, punch_min_radius_norm]

        x_data.append(data)

    x_train_norm = np.array(x_data)
    y_train = gp.y_train_

    # Now create a new GaussianProcessRegressor with your desired kernel and settings
    # kernel = Matern(length_scale=1, nu=2.5)
    # kernel = RBF(length_scale=1)
    # kernel = RationalQuadratic(length_scale=1, alpha=1)
    # kernel = ConstantKernel(1.0, (0.1, 10)) * Matern(length_scale=1, nu=2.5) \
    #                + WhiteKernel(noise_level=1e-5, noise_level_bounds=(1e-8, 1e-3))
    kernel = ConstantKernel(1.0, (0.1, 10)) * RationalQuadratic(length_scale=1.0, alpha=1.0, 
                                        length_scale_bounds=(1e-2, 10), alpha_bounds=(1e-2, 10)) \
                                        + WhiteKernel(noise_level=1e-8, noise_level_bounds=(1e-15, 1e-3))
    
    gp_norm = GaussianProcessRegressor(alpha=1e-06,
                                    kernel=kernel,
                                    n_restarts_optimizer=20,
                                    normalize_y=False,
                                    random_state=0)  # set a random state for reproducibility

    # Retrain the GP on the normalized inputs
    gp_norm.fit(x_train_norm, y_train)

    optimizer._gp = gp_norm

    # Create a new results list with normalized parameter values
    new_res = []
    for res in optimizer.res:
        orig_params = res["params"]
        new_params = {
            "punch_min_radius": (orig_params["punch_min_radius"] - min_punch_min_radius) / (max_punch_min_radius - min_punch_min_radius),
            "punch_profile_radius": (orig_params["punch_profile_radius"] - min_punch_fillet_radius) / (max_punch_fillet_radius - min_punch_fillet_radius)
        }
        # Keep the target value the same; if you need to update that too, apply similar normalization.
        new_res.append({
            "target": res["target"],
            "params": new_params
        })

    # Replace the optimizer's results with the new normalized results
    optimizer._space.res = new_res

    model_path = r"C:\Users\kam97\OneDrive - University of Bath\Documents\build\process_optimisation2\bayes_model_norm.pkl"
    
    with open(model_path, "wb") as f:
        pickle.dump(optimizer, f)
    print(f"Model saved to {model_path}")

    # print(optimizer.res)



    # Optionally, print the optimized kernel and log marginal likelihood
    # print(gp_norm)
    # print(str(2.5))
    print("Optimized kernel:", gp_norm.kernel_)
    print("Log marginal likelihood:", gp_norm.log_marginal_likelihood())

        


if __name__ == "__main__":
    

    main()



    # plot_gp_default()
    # plot_sampling_history()
    # plot_convergence()

    # op_gp()

    # # BO post pro
    # plot_convergence_2d()
    # plot_gp_2d()
    # plot_acq_2d()
    # plot_sampling_history_2d()
    # plot_sampling_history_rmse_2d()

    # Example usage:
    # Test UCB with kappa=1.0
    # plot_acq_2d(acq_kind="ucb", kappa=1.0)

    # # Test EI with xi=0.01
    # plot_acq_2d(acq_kind="ei", xi=0.01)

    # # Test PI with xi=0.01 (or another value if you prefer)
    # plot_acq_2d(acq_kind="poi", xi=0.01)



    # mesh_conv_final_max()
    # print("##############")
    # mesh_conv_final_min()
    # plot_exp()

    # plot_thickness_exp()

    # os.chdir(r"C:\Users\kam97\OneDrive - University of Bath\Documents\build\stamping_op_build2")
    # plot_convergence()
