"""
Example for CSTR model usage.

1) Set buffer 1 and buffer 2
2) Provide inlet time and data as numpy arrays
3) Create a mixer with defined volume, flow rate, and initial concentration
4) Create the save directory using the helper function
5) Run the mixing simulation with cstr.run_cstr_sim
"""
import os, sys
import numpy as np
from pathlib import Path

root = Path(__file__).resolve().parents[2]
if str(root) not in sys.path:
    sys.path.insert(0, str(root))

from buffer_mixing import buffers, helpers, cstr


# ======================================================================================================
# Buffers for the methyl orange cstr experiments

b1 = buffers.Buffer(["tris"], 0.0, 42.35, [50.0]) # 50mM tris pH 7.40
# b1 = buffers.Buffer(["tris"], 150.0, 193.12, [50.0]) # 50mM tris + 150mM NaCl pH 7.40
# b1 = buffers.Buffer(["tris"], 0.0, 128.75, [150.0]) # 150mM tris pH 7.40

# b1 = buffers.Buffer(["phosphate"], 79.295, 0.0, [50.0]) # 50mM phos pH 7.00
# b1 = buffers.Buffer(["phosphate"], 232.362, 150.0, [50.0]) # 50mM phos + 150mM NaCl pH 7.00
# b1 = buffers.Buffer(["phosphate"], 249.699, 0.0, [150.0]) # 150mM phos pH 7.00

b2 = buffers.Buffer(["acetate"], 2.415, 0.0, [50.0]) # 50mM acetate pH 3.50
# b2 = buffers.Buffer(["acetate"], 52.736, 50.0, [50.0]) # 50mM acetate + 50mM NaCl pH 3.50
# b2 = buffers.Buffer(["acetate"], 8.20, 0.0, [150.0]) # 150mM acetate pH 3.50

# b2 = buffers.Buffer(["citrate"], 40.95, 0.0, [50.0]) # 50mM citrate pH 3.50
# b2 = buffers.Buffer(["citrate"], 130.40, 0.0, [150.0]) # 150mM citrate pH 3.50

# b2 = buffers.Buffer(["glycine"], 0.0, 3.78, [50.0]) # 50mM glycine pH 3.50
# b2 = buffers.Buffer(["glycine"], 0.0, 11.10, [150.0]) # 150mM glycine pH 3.50

# ======================================================================================================
### CSTR simulation

inlet_time = np.linspace(0, 1320, 1321)
inlet_data = array = np.concatenate([np.ones(121), np.zeros(600), np.ones(600)])
mixer = cstr.Mixer(6e-7, 1.67e-8, 1.0)
buffers = [b1, b2]

# Create results folder
save_dir = helpers.make_results_folder(run_type="cstr", output_folder_path=os.getcwd())

# Run CSTR simulation
cstr.run_cstr_sim(mixer, buffers, inlet_time, inlet_data, save_dir)
