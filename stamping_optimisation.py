from bayes_opt import BayesianOptimization, JSONLogger, Events
import subprocess
import json
import os


def cup_function(punch_depth, punch_profile_radius, punch_min_radius,
    die_profile_radius, die_punch_gap):

    # target geometry
    ideal_height = 10
    ideal_radius = 15
    ideal_profile_radius = 6.5

    abaqus_command = [
        r"C:\SIMULIA\Commands\abaqus.bat", "cae", "noGUI=stampy.py", "--",
        str(ideal_height), str(ideal_radius), str(ideal_profile_radius),
        str(punch_depth), str(punch_profile_radius), str(punch_min_radius),
        str(die_profile_radius), str(die_punch_gap)
    ]

    # Run the subprocess and capture RMSE from stdout
    result = subprocess.run(abaqus_command, capture_output=True, text=True, shell=True)

    # # If Abaqus failed, return None
    # if result.returncode != 0:
    #     print("Abaqus execution failed!")
    #     return None

    try:
        with open(r"C:\Users\kam97\OneDrive - University of Bath\Documents\build\sim_output.json", "r") as f:
            data_out = json.load(f)
        rmse = data_out["rmse"]

        print(f"Simulation Complete, RMSE: {rmse}")

        return rmse
    
    except Exception as exc:
        print("Could not parse sim_output.json:", exc)
        print("STDOUT:\n", result.stdout)
        print("STDERR:\n", result.stderr)
        return None


def main():

    # Bounded region of parameter space
    pbounds = {
        "punch_depth": (15, 25), 
        "punch_profile_radius": (4, 7), 
        "punch_min_radius": (10, 20),
        "die_profile_radius": (6, 9), 
        "die_punch_gap": (1,5)
        }

    optimizer = BayesianOptimization(
        f=cup_function,
        pbounds=pbounds,
        random_state=1,
    )


    logger = JSONLogger(path=r"C:\Users\kam97\OneDrive - University of Bath\Documents\build\op_logs.log")
    optimizer.subscribe(Events.OPTIMIZATION_STEP, logger)

    optimizer.maximize(init_points=2, n_iter=3)

    print(optimizer.max)


if __name__ == "__main__":
    main()

