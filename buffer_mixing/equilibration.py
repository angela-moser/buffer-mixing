"""
Functions for buffer equilibration calculations
"""
import numpy as np
from scipy.optimize import least_squares
from . import buffers, helpers


# =============================================================================
# Buffer equations using no nested lists
def buffer_equations(variables, C_Na, C_Cl, C_tot_list, Ka_all, zA_all, lengths, num_equations):
    eqs = np.zeros(num_equations)

    C_H = variables[0]
    C_OH = 1e-8 / C_H

    # Electroneutrality equation
    offset = 0
    idx = 1
    sum_charge = 0.0
    for b in range(len(lengths)):
        n = lengths[b]
        zA_list = zA_all[offset : offset + n]
        for i in range(n):
            sum_charge += variables[idx] * zA_list[i]
            idx += 1
        offset += n
    eqs[0] = C_H + C_Na - C_OH - C_Cl + sum_charge

    # Total buffer equations
    idx = 1
    eq_i = 1
    for b in range(len(lengths)):
        total = 0.0
        for i in range(lengths[b]):
            total += variables[idx]
            idx += 1
        eqs[eq_i] = C_tot_list[b] - total
        eq_i += 1

    # Dissociation equations
    ka_offset = 0
    c_offset = 0
    for b in range(len(lengths)):
        n = lengths[b]  # number of species in this buffer
        for i in range(n - 1): # One less Ka equation than buffer species
            C_i = variables[1 + c_offset + i]
            C_ip1 = variables[1 + c_offset + i + 1]
            eqs[eq_i] = Ka_all[ka_offset + i] * C_i - C_ip1 * C_H
            eq_i += 1
        ka_offset += n - 1
        c_offset += n

    return eqs

# =============================================================================
def flatten(l):
    # Flatten a list of lists
    return [item for sublist in l for item in sublist]


def unpack_array(array, lengths):
    # Unpack a flat array into a list of arrays with given lengths
    result = []
    idx = 0
    for l in lengths:
        result.append(array[idx:idx+l])
        idx += l
    return result


def solve_acid_base_multi(C_Na, C_Cl, C_tot_list, Ka_list_list, zA_list_list):
    # Convert from M to mM for Ka values
    Ka_list_list = [np.array(Ka) * 1e3 for Ka in Ka_list_list]

    # Ensure all concentrations are positive
    C_tot_list = [max(c, 0) for c in C_tot_list]

    # Calculate number of equations (Electroneutrality + total + dissociation equations)
    num_equations = 1 + len(C_tot_list) + sum(len(i) for i in Ka_list_list)

    # Flatten nested lists and record lengths of sublists
    Ka_all = np.concatenate(Ka_list_list)
    zA_all = np.concatenate(zA_list_list)
    lengths = np.array([len(l) for l in zA_list_list]) # number of species for each buffer

    def wrapped_equations(x):
        return buffer_equations(x, C_Na, C_Cl, C_tot_list, Ka_all, zA_all, lengths, num_equations)

    # Initial guesses for C_H and concentrations of acids
    highest_C_tot_index = np.argmax(C_tot_list)
    mean_Kas = [np.mean(Ka_list) for Ka_list in Ka_list_list]
    C_H_guess = mean_Kas[highest_C_tot_index]
    C_a_guess = [c/n for c, n in zip(C_tot_list, lengths) for _ in range(n)]

    # Create a list of possible initial H+ guesses to try
    H_guesses = [C_H_guess, C_H_guess*1e2, C_H_guess*1e-2, C_H_guess*1e5,
                 C_H_guess*1e-5, C_H_guess*1e7, C_H_guess*1e-7]
    bounds = (0, 1000)
    H_guesses = np.clip(H_guesses, bounds[0], bounds[1]).tolist()

    # Define a function to check for negative concentrations
    def is_valid_solution(sol):
        return sol.success and np.all(sol.x >= 0)

    # Solve the equations trying multiple guesses if needed
    for guess in H_guesses:
        initial_guess = np.concatenate([[guess], C_a_guess])
        sol = least_squares(wrapped_equations, initial_guess, bounds=bounds)
        # sol = root(buffer_equations, initial_guess, method='hybr') # Faster but doesn't always work
        if is_valid_solution(sol):
            break
        else:
            print(f'No valid solution. C = {sol.x}')
            # Return NaNs if no solution found
            return np.nan, np.nan, C_Na, C_Cl, [np.nan] * len(C_a_guess)

    # Extract concentrations
    C_H = sol.x[0]
    C_a_list = unpack_array(sol.x[1:], lengths)
    C_OH = 1e-8 / C_H  # Kw = 1e-8 mmol^2/L^2

    return C_H, C_OH, C_Na, C_Cl, C_a_list


def equilibrate(C_Na, C_Cl, C_tot_list, pKa_list_list, zA_list_list, A, b):
    # Initial pKa values
    Ka_list_list = helpers.pK_to_K(pKa_list_list)  # Precompute Ka values once
    tolerance = 1e-3
    difference = 1e2  # Start with a large difference to enter the loop

    # Initialize values for the loop
    C_H, C_OH, C_Na, C_Cl, C_a_list = solve_acid_base_multi(C_Na, C_Cl, C_tot_list, Ka_list_list, zA_list_list)
    species = [C_H, C_OH, C_Na, C_Cl, C_a_list]
    ionic_strength = helpers.calculate_ionic_strength(C_H, C_OH, C_Na, C_Cl, C_a_list, zA_list_list)

    # Precompute corrected pKa based on initial ionic strength
    pKa_corrected_list = helpers.corrected_pKa(ionic_strength, A, b, zA_list_list, pKa_list_list)
    Ka_corrected_list = helpers.pK_to_K(pKa_corrected_list)

    previous_pKa_list = pKa_corrected_list

    # Iterative loop to converge the pKa correction
    while difference > tolerance:
        # Solve the acid-base system with the corrected Ka values
        C_H, C_OH, C_Na, C_Cl, C_a_list = solve_acid_base_multi(C_Na, C_Cl, C_tot_list, Ka_corrected_list, zA_list_list)
        species = [C_H, C_OH, C_Na, C_Cl, C_a_list]
        ionic_strength = helpers.calculate_ionic_strength(C_H, C_OH, C_Na, C_Cl, C_a_list, zA_list_list)

        # Get corrected pKa values based on the updated ionic strength
        pKa_corrected_list = helpers.corrected_pKa(ionic_strength, A, b, zA_list_list, pKa_list_list)
        Ka_corrected_list = helpers.pK_to_K(pKa_corrected_list)

        # Calculate the differences for convergence check
        flat_pKa = np.array(flatten(pKa_corrected_list))
        flat_prev = np.array(flatten(previous_pKa_list))
        difference = np.max(np.abs(flat_pKa - flat_prev))

        previous_pKa_list = pKa_corrected_list  # Update for next iteration

    return species


def solve_buffer_system(types, C_Na, C_Cl, C_tot_list, do_print=True):
    # Get buffer info
    pKa_list, zA_list, A_list, b_list = buffers.get_buffer_info(types)

    # Solve for equilibrium concentrations
    species = equilibrate(C_Na, C_Cl, C_tot_list, pKa_list, zA_list, A_list, b_list)
    C_H, C_OH, C_Na, C_Cl, C_a_list = species

    ionic_strength = helpers.calculate_ionic_strength(C_H, C_OH, C_Na, C_Cl, C_a_list, zA_list)
    pKa_corrected = helpers.corrected_pKa(ionic_strength, A_list, b_list, zA_list, pKa_list)
    activity_coeff = helpers.calculate_activity_coeff(ionic_strength)
    pH = helpers.calculate_pH(activity_coeff * C_H * 1e-3)  # Convert C_H from mM to M first

    # Print results if necessary
    if do_print:
        helpers.print_results(types, species, ionic_strength, pKa_corrected, activity_coeff, pH, pKa_list)

    return species, ionic_strength, pKa_corrected, activity_coeff, pH


def solve_for_Na(target_pH, Na_low, Na_high, types, C_Cl, C_tot_list):
    tol = 1e-6
    max_iter = 1000

    # Unpack return values but keep only pH
    pH_low = solve_buffer_system(types=types, C_Na=Na_low, C_Cl=C_Cl, C_tot_list=C_tot_list, do_print=False)[4]
    pH_high = solve_buffer_system(types=types, C_Na=Na_high, C_Cl=C_Cl, C_tot_list=C_tot_list, do_print=False)[4]

    diff_low = pH_low - target_pH
    diff_high = pH_high - target_pH

    if diff_low * diff_high > 0:
        raise ValueError("Function does not change sign in the given interval; no root guaranteed.")

    iter_count = 0
    while (Na_high - Na_low) > tol and iter_count < max_iter:
        Na_mid = (Na_low + Na_high) / 2
        pH_mid = solve_buffer_system(types=types, C_Na=Na_mid, C_Cl=C_Cl, C_tot_list=C_tot_list, do_print=False)[4]

        diff = pH_mid - target_pH

        if abs(diff) < tol:
            C_Na = Na_mid
            break
        elif diff * diff_low < 0:
            Na_high = Na_mid
            diff_high = diff
        else:
            Na_low = Na_mid
            diff_low = diff

        iter_count += 1

    if iter_count == max_iter:
        raise RuntimeError("Maximum iterations reached without convergence.")

    C_Na = (Na_low + Na_high) / 2

    # Print result
    print("\n-------------------------------------")
    print(f"\nSolution for target pH = {target_pH:.3f} : CNa = {C_Na:.3f}")

    species, ionic_strength, pKa_corrected, activity_coeff, pH = solve_buffer_system(
        types=types, C_Na=C_Na, C_Cl=C_Cl, C_tot_list=C_tot_list, do_print=True
    )

    return species, ionic_strength, pKa_corrected, activity_coeff, pH
