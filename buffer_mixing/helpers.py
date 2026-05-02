"""
This module contains helper functions for the buffer calculator.
=========================
activity_coeff: equation 23 from Pabst and Carta 2007
corrected pKa: equation 24 from Pabst and Carta 2007
"""
import math
import os
from datetime import datetime


def calculate_ionic_strength(C_H, C_OH, C_Na, C_Cl, C_a, zA_list):
    # Sum contributions from C_H, C_OH, C_Na, C_Cl (all charge 1)
    charge_ions = sum([C_H, C_OH, C_Na, C_Cl])

    # Handle single buffer (1D) or multiple buffers (2D)
    if isinstance(zA_list[0], int):
        buffer_charge = sum([c * z**2 for c, z in zip(C_a, zA_list)])
    else:
        buffer_charge = sum(
            sum([c * z**2 for c, z in zip(ca, za)])
            for ca, za in zip(C_a, zA_list)
        )

    ionic_strength = 0.5 * (charge_ions + buffer_charge)
    return ionic_strength


def calculate_activity_coeff(IS, A=0.5114, b=0.1, z=1):
    IS = IS * 1e-3  # change to molar units from mM
    g = 10.0 ** (-z**2 * (A * math.sqrt(IS) / (1 + math.sqrt(IS)) - b * IS))
    return g


def corrected_pKa(IS, A, b, zA, pKa):
    IS = IS * 1e-3  # change to molar units from mM

    # Handle the case for an individual buffer (1D)
    if isinstance(zA[0], int):
        zA = zA[1:]  # skip first element (not starting at 1 because we don't subtract 1)
        pKa_cor = [pk + 2 * z * (A * math.sqrt(IS) / (1 + math.sqrt(IS)) - b * IS)
                   for pk, z in zip(pKa, zA)]
    else:
        # Handle the case for multiple buffers (2D)
        zA = [z[1:] for z in zA]
        pKa_cor = [
            [pk + 2 * zz * (Ai * math.sqrt(IS) / (1 + math.sqrt(IS)) - bi * IS)
             for pk, zz in zip(pKa_i, zA_i)]
            for Ai, bi, zA_i, pKa_i in zip(A, b, zA, pKa)
        ]

    return pKa_cor


def pK_to_K(pK_array):
    result = []
    for x in pK_array:
        if isinstance(x, (list, tuple)):
            result.append(pK_to_K(x))
        else:
            result.append(10.0 ** (-x))
    return result


def calculate_pH(C_H):
    if C_H <= 0:
        print(f"Negative C_H encountered (C_H = {C_H:.3e}), pH will be set to 0")
        return 0
    return -math.log10(C_H)


def make_results_folder(run_type, output_folder_path):
    # Check if the main output folder exists; if not, create it
    if not os.path.isdir(output_folder_path):
        os.makedirs(output_folder_path)

    # Create the "output" folder path
    output_path = os.path.join(output_folder_path, "output")

    # Create the "output" folder if it doesn't exist
    if not os.path.isdir(output_path):
        os.makedirs(output_path)

    # Create the new folder name with timestamp
    timestamp = datetime.now().strftime("%Y%m%d%H%M%S")
    folder_name = f"{run_type}_results_{timestamp}"
    new_folder_path = os.path.join(output_path, folder_name)

    # Create the new folder
    os.makedirs(new_folder_path, exist_ok=True)

    return new_folder_path


def print_results(types, species, ionic_strength, pKa_corrected, activity_coeff, pH, pKa_list):
    C_H, C_OH, C_Na, C_Cl, C_a = species

    print("\n-------------------------------------")

    # Print corrected pKa values
    if isinstance(pKa_corrected[0], (float, int)):
        # Individual buffer
        for i, (pka_cor, pka_orig) in enumerate(zip(pKa_corrected, pKa_list)):
            print(f"\ncorrected {types} pKa[{i}] = {pka_cor:.3f} (original pKa = {pka_orig:.3f})")
    else:
        # Multiple buffers
        for i, (type_name, pka_cor_list, pka_orig_list) in enumerate(zip(types, pKa_corrected, pKa_list)):
            for j, (pka_cor, pka_orig) in enumerate(zip(pka_cor_list, pka_orig_list)):
                print(f"\ncorrected {type_name} pKa[{j}] = {pka_cor:.3f} (original pKa = {pka_orig:.3f})")

    print("\n\nEquilibrium concentrations (mM):")
    print(f"C_H = {C_H:.3e}")
    print(f"C_OH = {C_OH:.3e}")
    print(f"C_Na = {C_Na:.3f}")
    print(f"C_Cl = {C_Cl:.3f}")

    # Print buffer concentrations
    if isinstance(C_a[0], (float, int)):
        # Individual buffer
        for i, ca in enumerate(C_a):
            print(f"{types} C_a[{i}] = {ca:.2f}")
    else:
        # Multiple buffers
        for i, (type_name, ca_list) in enumerate(zip(types, C_a)):
            for j, ca in enumerate(ca_list):
                print(f"{type_name} C_a[{j}] = {ca:.3f}")

    print(f"\nactivity coefficient = {activity_coeff:.3f}")
    print(f"\nionic strength = {ionic_strength:.2f} (mM)")
    print(f"\npH = {pH:.3f}")
