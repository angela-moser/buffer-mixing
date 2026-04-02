# cstr equation solver for buffer mixing
import os, time
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from . import buffers, helpers, equilibration


class Mixer:
    def __init__(self, volume, flow, init_conc):
        self.volume = volume
        self.flow = flow
        self.init_conc = init_conc


def cstr_model(c, p, t):
    inlet_time, inlet_data, flow_rate, volume = p
    # Interpolate the inlet concentration based on current time t
    itp = interp1d(inlet_time, inlet_data, fill_value="extrapolate")
    c_inlet = float(itp(t))
    # Calculate the rate of change of concentration in the reactor
    dc_dt = (flow_rate / volume) * (c_inlet - c)
    return dc_dt


def solve_cstr_model(mixer, inlet_time, inlet_data):
    # Initial concentration in the mixer
    initial_concentration = mixer.init_conc

    # Ensure inlet time is properly formatted
    if not (isinstance(inlet_time, (list, tuple)) or hasattr(inlet_time, "__array__")) or len(inlet_time) == 0:
        raise ValueError("Inlet time must be a non-empty array.")
    if not (isinstance(inlet_data, (list, tuple)) or hasattr(inlet_data, "__array__")) or len(inlet_data) == 0:
        raise ValueError("Inlet data must be a non-empty array.")

    # Define the ODE system
    def ode_func(t, c, p):
        return cstr_model(c, p, t)

    t_span = (inlet_time[0], inlet_time[-1])
    t_eval = inlet_time
    params = (inlet_time, inlet_data, mixer.flow, mixer.volume)

    sol = solve_ivp(ode_func, t_span, [initial_concentration], args=(params,), t_eval=t_eval, rtol=2.3e-14, atol=1e-16, vectorized=False)

    output_time = sol.t
    output_array = sol.y[0]

    return output_time, output_array


def mix_buffers(buffers, fractions):
    # Validate inputs
    if len(buffers) != len(fractions):
        raise ValueError("Number of buffers and fractions must match.")

    # Threshold to avoid very small values
    fractions = [max(f, 1e-12) for f in fractions]

    # Initialize a dictionary to store the concentrations
    mixed_concentrations = {}

    # Compute concentrations for Na and Cl
    mixed_concentrations["Na"] = sum(buffers[i].C_Na * fractions[i] for i in range(len(buffers)))
    mixed_concentrations["Cl"] = sum(buffers[i].C_Cl * fractions[i] for i in range(len(buffers)))

    for i in range(len(buffers)):
        fraction = fractions[i]
        C_tot_list = buffers[i].C_tot_list
        types = buffers[i].types
        # Iterate over buffer types and total concentrations
        for type_name, conc in zip(types, C_tot_list):
            if type_name in mixed_concentrations:
                # Add to the existing concentration if the type is already in the dictionary
                mixed_concentrations[type_name] += conc * fraction
            else:
                # Initialize the concentration for this buffer type
                mixed_concentrations[type_name] = conc * fraction

    return mixed_concentrations


def get_outlet_concentrations(mixer, buffers, inlet_time, inlet_data):
    # Get the cstr output fraction profile
    output_time, cstr_fraction_array = solve_cstr_model(mixer, inlet_time, inlet_data)

    # Mix the buffers 50/50 and determine unique buffer types present
    mixed_concentrations = mix_buffers(buffers, [0.5, 0.5])

    # Collect the buffer types and remove Na and Cl from the list
    buffer_keys = [k for k in mixed_concentrations.keys() if k not in ("Na", "Cl")]

    # Preallocate solution array [time, Na, Cl, buffer comps]
    n_time = len(cstr_fraction_array)
    n_buffers = len(mixed_concentrations)
    solution_array = np.zeros((n_time, 1 + n_buffers))

    C_Na = np.zeros(n_time)
    C_Cl = np.zeros(n_time)
    C_a_array = np.zeros((n_time, len(buffer_keys)))

    # Mix the buffers according to the CSTR output at each time point and update concentrations in solution array
    for i in range(n_time):
        fractions = [cstr_fraction_array[i], 1 - cstr_fraction_array[i]]
        mixed_conc = mix_buffers(buffers, fractions)
        C_Na[i] = mixed_conc["Na"]
        C_Cl[i] = mixed_conc["Cl"]
        for j, key in enumerate(buffer_keys):
            C_a_array[i, j] = mixed_conc[key]

    # Fill the solution array
    solution_array[:, 0] = output_time
    solution_array[:, 1] = C_Na
    solution_array[:, 2] = C_Cl
    solution_array[:, 3:] = C_a_array

    return solution_array, buffer_keys


def equilibrate_cstr(solution_array: np.ndarray, buffer_keys: list):
    # Solve first row to determine species array size
    first_row = solution_array[0, :]

    # Convert buffer keys (symbols) back to strings
    types = [str(key) for key in buffer_keys]

    C_Na_1 = first_row[1]
    C_Cl_1 = first_row[2]
    C_tot_1 = first_row[3:]

    pKa_list, zA_list, A_list, b_list = buffers.get_buffer_info(types=types)

    # Solve for equilibrium concentrations
    species_1 = equilibration.equilibrate(C_Na_1, C_Cl_1, C_tot_1, pKa_list, zA_list, A_list, b_list)
    C_H_1, C_OH_1, C_Na_1, C_Cl_1, C_a_list_1 = species_1
    combined_C_a_1 = np.concatenate(C_a_list_1)

    ionic_strength_1 = helpers.calculate_ionic_strength(C_H_1, C_OH_1, C_Na_1, C_Cl_1, C_a_list_1, zA_list)
    activity_coeff_1 = helpers.calculate_activity_coeff(ionic_strength_1)
    pH_1 = helpers.calculate_pH(activity_coeff_1 * C_H_1 * 1e-3)  # Convert C_H from mM to M

    # Collect names for buffer species
    column_names = ["time", "ionic strength", "pH", "H+", "OH-", "Na+", "Cl-"]
    for i, species_group in enumerate(C_a_list_1):
        for j in range(len(species_group)):
            column_names.append(f"{types[i]} {j+1}")

    # Preallocate solution array [time, IS, pH, H, OH, Na, Cl, species_comps]
    equilibrated_solution_array = np.zeros((solution_array.shape[0], 7 + len(combined_C_a_1)))

    # Assign time column
    equilibrated_solution_array[:, 0] = solution_array[:, 0]

    # Assign first row solution
    equilibrated_solution_array[0, 1:7] = [ionic_strength_1, pH_1, C_H_1, C_OH_1, C_Na_1, C_Cl_1]
    equilibrated_solution_array[0, 7:] = list(combined_C_a_1)

    # Loop over remaining rows
    for i in range(1, solution_array.shape[0]):
        row = solution_array[i, :]
        C_Na = row[1]
        C_Cl = row[2]
        C_tot = row[3:]

        species = equilibration.equilibrate(C_Na, C_Cl, C_tot, pKa_list, zA_list, A_list, b_list)
        C_H, C_OH, C_Na, C_Cl, C_a_list = species
        combined_C_a = np.concatenate(C_a_list)

        ionic_strength = helpers.calculate_ionic_strength(C_H, C_OH, C_Na, C_Cl, C_a_list, zA_list)
        activity_coeff = helpers.calculate_activity_coeff(ionic_strength)
        pH = helpers.calculate_pH(activity_coeff * C_H * 1e-3)

        equilibrated_solution_array[i, 1:7] = [ionic_strength, pH, C_H, C_OH, C_Na, C_Cl]
        equilibrated_solution_array[i, 7:] = list(combined_C_a)

    # Create DataFrame
    equilibrated_solution_df = pd.DataFrame(equilibrated_solution_array, columns=column_names)

    return equilibrated_solution_df


def make_csv_cstr(df, save_dir):
    # Define the file path with .csv extension
    path = os.path.join(save_dir, "cstr_out.csv")

    # Write the DataFrame to a CSV file
    df.to_csv(path, index=False)

    print(f"CSV saved to {path}")


def plot_cstr(df: pd.DataFrame, save_dir: str):
    # Settings
    font_size = 12
    RPI_cherry = (0.816, 0, 0.086)  # RGB values for the plot color

    time = df["time"]
    pH = df["pH"]
    ionic_strength = df["ionic strength"]

    # Create subplots
    fig, axes = plt.subplots(1, 3, figsize=(13, 4), dpi=300)
    ax_pH, ax_IS, ax_species = axes

    # Plot pH
    ax_pH.plot(time, pH, label="pH", color=RPI_cherry, linewidth=2)
    ax_pH.set_xlim([0, time.max()])
    ax_pH.set_ylim([pH.min() * 0.9, pH.max() * 1.1])
    ax_pH.set_xlabel("Time [sec]", fontsize=font_size)
    ax_pH.set_ylabel("pH", fontsize=font_size)
    ax_pH.set_title("CSTR output pH", fontsize=font_size)

    # Plot ionic strength
    max_c = ionic_strength.max()
    ax_IS.plot(time, ionic_strength, label="ionic strength", color=RPI_cherry, linewidth=2)
    ax_IS.set_xlim([0, time.max()])
    ax_IS.set_ylim([0, max_c * 1.1])
    ax_IS.set_xlabel("Time [sec]", fontsize=font_size)
    ax_IS.set_ylabel("Ionic strength [mM]", fontsize=font_size)
    ax_IS.set_title("CSTR output ionic strength", fontsize=font_size)

    # Plot species concentrations
    exclude_columns = ["time", "ionic strength", "pH", "H+", "OH-"]
    for col in df.columns:
        if col not in exclude_columns:
            col_max = df[col].max()
            max_c = max(max_c, col_max)
            ax_species.plot(time, df[col], label=str(col), linewidth=2)

    ax_species.set_xlim([0, time.max()])
    ax_species.set_ylim([0, max_c * 1.1])
    ax_species.set_xlabel("Time [sec]", fontsize=font_size)
    ax_species.set_ylabel("Concentration [mM]", fontsize=font_size)
    ax_species.set_title("CSTR output concentrations", fontsize=font_size)
    ax_species.legend(loc="center left", bbox_to_anchor=(1, 0.5), fontsize=font_size * 0.8)

    # Adjust layout and save figure
    plt.tight_layout(pad=3.0)
    output_path = os.path.join(save_dir, "cstr_out.png")
    plt.savefig(output_path, bbox_inches="tight")

    print(f"Plot saved to {output_path}")


def run_cstr_sim(mixer, buffers, inlet_time, inlet_data, save_dir):
    print("Running CSTR simulation...")
    cstr_start_time = time.time()

    # Compute outlet concentrations
    solution_array, buffer_keys = get_outlet_concentrations(mixer, buffers, inlet_time, inlet_data)

    # Equilibrate the output concentrations
    equilibrated_solution_df = equilibrate_cstr(solution_array, buffer_keys)

    # Report elapsed time
    cstr_elapsed_time = time.time() - cstr_start_time
    print(f"Time for CSTR simulation: {cstr_elapsed_time:.2f} seconds")

    # Save results and plot
    make_csv_cstr(equilibrated_solution_df, save_dir)
    plot_cstr(equilibrated_solution_df, save_dir)

    return solution_array, buffer_keys
