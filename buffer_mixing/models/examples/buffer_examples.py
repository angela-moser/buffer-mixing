"""
Examples for buffer equilibration model usage.
"""
import sys
from pathlib import Path

current = Path(__file__).resolve()

for parent in current.parents:
    if (parent / "buffer_mixing").exists():
        sys.path.insert(0, str(parent))
        break

from buffer_mixing import equilibration


# ======================================================================================================
### Single buffer example
equilibration.solve_buffer_system(types=["formate"], C_Na=69.0, C_Cl=0.0, C_tot_list=[100.0]) # 4.001

# ======================================================================================================
### Mixture of two buffers example

# equilibration.solve_buffer_system(types=["acetate", "phosphate"], C_Na=104.4, C_Cl=0.0, C_tot_list=[50.0, 50.0]) # 6.02

# ======================================================================================================
### Getting Na from pH example

# equilibration.solve_for_Na(target_pH=3.5, Na_low=0.0, Na_high=5.0, types=["acetate"], C_Cl=0.0, C_tot_list=[50.0])

# ======================================================================================================
### Methyl orange experiment buffers

# equilibration.solve_buffer_system(types=["tris"], C_Na=0.0, C_Cl=42.35, C_tot_list=[50.0]) # 7.400
# equilibration.solve_buffer_system(types=["tris"], C_Na=150.0, C_Cl=193.12, C_tot_list=[50.0]) # 7.400
# equilibration.solve_buffer_system(types=["tris"], C_Na=0.0, C_Cl=128.75, C_tot_list=[150.0]) # 7.400

# equilibration.solve_buffer_system(types=["phosphate"], C_Na=79.295, C_Cl=0.0, C_tot_list=[50.0]) # 7.000
# equilibration.solve_buffer_system(types=["phosphate"], C_Na=249.699, C_Cl=0.0, C_tot_list=[150.0]) # 7.000
# equilibration.solve_buffer_system(types=["phosphate"], C_Na=232.362, C_Cl=150.0, C_tot_list=[50.0]) # 7.000

# equilibration.solve_buffer_system(types=["acetate"], C_Na=2.415, C_Cl=0.0, C_tot_list=[50.0]) # 3.500
# equilibration.solve_buffer_system(types=["acetate"], C_Na=8.20, C_Cl=0.0, C_tot_list=[150.0]) # 3.500
# equilibration.solve_buffer_system(types=["acetate"],  C_Na=52.736, C_Cl=50.0, C_tot_list=[50.0]) # 3.500

# equilibration.solve_buffer_system(types=["citrate"], C_Na=40.95, C_Cl=0.0, C_tot_list=[50.0]) # 3.500
# equilibration.solve_buffer_system(types=["citrate"], C_Na=130.40, C_Cl=0.0, C_tot_list=[150.0]) # 3.500

# equilibration.solve_buffer_system(types=["glycine"], C_Na=0.0, C_Cl=3.78, C_tot_list=[50.0]) # 3.500
# equilibration.solve_buffer_system(types=["glycine"], C_Na=0.0, C_Cl=11.10, C_tot_list=[150.0]) # 3.500

# ======================================================================================================
### Methyl orange calibration curve buffers

# No added salt
# equilibration.solve_buffer_system(types=["acetate"], C_Na=0.738, C_Cl=0.0, C_tot_list=[50.0]) # 3.200, IS 1.40
# equilibration.solve_buffer_system(types=["acetate"], C_Na=2.415, C_Cl=0.0, C_tot_list=[50.0]) # 3.500, IS 2.750
# equilibration.solve_buffer_system(types=["acetate"], C_Na=5.12, C_Cl=0.0, C_tot_list=[50.0]) # 3.800, IS 5.29
# equilibration.solve_buffer_system(types=["acetate"], C_Na=9.65, C_Cl=0.0, C_tot_list=[50.0]) # 4.100, IS 9.74
# equilibration.solve_buffer_system(types=["acetate"], C_Na=16.5, C_Cl=0.0, C_tot_list=[50.0]) # 4.400, IS 16.55

# 150 mM added salt
# equilibration.solve_for_Na(target_pH=3.2, Na_low=150.0, Na_high=200.0, types=["acetate"], C_Cl=150.0, C_tot_list=[50.0]) # IS 151.71
# equilibration.solve_for_Na(target_pH=3.5, Na_low=150.0, Na_high=200.0, types=["acetate"], C_Cl=150.0, C_tot_list=[50.0]) # IS 153.31
# equilibration.solve_for_Na(target_pH=3.8, Na_low=150.0, Na_high=200.0, types=["acetate"], C_Cl=150.0, C_tot_list=[50.0]) # IS 156.19
# equilibration.solve_for_Na(target_pH=4.1, Na_low=150.0, Na_high=200.0, types=["acetate"], C_Cl=150.0, C_tot_list=[50.0]) # IS 161.01
# equilibration.solve_for_Na(target_pH=4.4, Na_low=150.0, Na_high=200.0, types=["acetate"], C_Cl=150.0, C_tot_list=[50.0]) # IS 168.04

# 350 mM added salt
# equilibration.solve_for_Na(target_pH=3.2, Na_low=350.0, Na_high=400.0, types=["acetate"], C_Cl=350.0, C_tot_list=[50.0]) # IS 351.72
# equilibration.solve_for_Na(target_pH=3.5, Na_low=350.0, Na_high=400.0, types=["acetate"], C_Cl=350.0, C_tot_list=[50.0]) # IS 353.33
# equilibration.solve_for_Na(target_pH=3.8, Na_low=350.0, Na_high=400.0, types=["acetate"], C_Cl=350.0, C_tot_list=[50.0]) # IS 356.22
# equilibration.solve_for_Na(target_pH=4.1, Na_low=350.0, Na_high=400.0, types=["acetate"], C_Cl=350.0, C_tot_list=[50.0]) # IS 361.04
# equilibration.solve_for_Na(target_pH=4.4, Na_low=350.0, Na_high=400.0, types=["acetate"], C_Cl=350.0, C_tot_list=[50.0]) # IS 368.05

# 1000 mM added salt
# equilibration.solve_buffer_system(types=["acetate"], C_Na=1000.55, C_Cl=1000.0, C_tot_list=[50.0]) # 3.200, IS 1001.45
# equilibration.solve_buffer_system(types=["acetate"], C_Na=1002.36, C_Cl=1000.0, C_tot_list=[50.0]) # 3.500, IS 1002.81
# equilibration.solve_buffer_system(types=["acetate"], C_Na=1005.083, C_Cl=1000.0, C_tot_list=[50.0]) # 3.800, IS 1005.31
# equilibration.solve_buffer_system(types=["acetate"], C_Na=1009.456, C_Cl=1000.0, C_tot_list=[50.0]) # 4.100, IS 1009.57
# equilibration.solve_buffer_system(types=["acetate"], C_Na=1015.957, C_Cl=1000.0, C_tot_list=[50.0]) # 4.400, IS 1016.01

# Citrate
# equilibration.solve_for_Na(target_pH=3.2, Na_low=0.0, Na_high=200.0, types=["citrate"], C_Cl=0.0, C_tot_list=[50.0]) # IS 32.3
# equilibration.solve_for_Na(target_pH=3.5, Na_low=0.0, Na_high=200.0, types=["citrate"], C_Cl=0.0, C_tot_list=[50.0]) # IS 44.77
# equilibration.solve_for_Na(target_pH=3.8, Na_low=0.0, Na_high=200.0, types=["citrate"], C_Cl=0.0, C_tot_list=[50.0]) # IS 59.25
# equilibration.solve_for_Na(target_pH=4.1, Na_low=0.0, Na_high=200.0, types=["citrate"], C_Cl=0.0, C_tot_list=[50.0]) # IS 77.54
# equilibration.solve_for_Na(target_pH=4.4, Na_low=0.0, Na_high=200.0, types=["citrate"], C_Cl=0.0, C_tot_list=[50.0]) # IS 100.23

# Glycine
# equilibration.solve_buffer_system(types=["glycine"], C_Na=0.0, C_Cl=7.25, C_tot_list=[50.0]) # pH 3.200, IS 7.25
# equilibration.solve_buffer_system(types=["glycine"], C_Na=0.0, C_Cl=3.78, C_tot_list=[50.0]) # pH 3.500, IS 3.78
# equilibration.solve_buffer_system(types=["glycine"], C_Na=0.0, C_Cl=1.925, C_tot_list=[50.0]) # pH 3.800, IS 1.925
# equilibration.solve_buffer_system(types=["glycine"], C_Na=0.0, C_Cl=0.967, C_tot_list=[50.0]) # pH 4.100, IS 0.97
# equilibration.solve_buffer_system(types=["glycine"], C_Na=0.0, C_Cl=0.483, C_tot_list=[50.0]) # pH 4.400, IS 0.48

# Acetate + Tris
# equilibration.solve_buffer_system(types=["acetate", "tris"], C_Na=0.0, C_Cl=49.16, C_tot_list=[50.0, 50.0]) # 3.200, IS 50.77
# equilibration.solve_buffer_system(types=["acetate", "tris"], C_Na=0.0, C_Cl=47.27, C_tot_list=[50.0, 50.0]) # 3.500, IS 50.39
# equilibration.solve_buffer_system(types=["acetate", "tris"], C_Na=0.0, C_Cl=44.34, C_tot_list=[50.0, 50.0]) # 3.800, IS 50.19
# equilibration.solve_buffer_system(types=["acetate", "tris"], C_Na=0.0, C_Cl=39.63, C_tot_list=[50.0, 50.0]) # 4.100, IS 50.09
# equilibration.solve_buffer_system(types=["acetate", "tris"], C_Na=0.0, C_Cl=32.78, C_tot_list=[50.0, 50.0]) # 4.400, IS 50.04

# Acetate + Phosphate
# equilibration.solve_for_Na(target_pH=3.2, Na_low=0.0, Na_high=500.0, types=["acetate", "phosphate"], C_Cl=0.0, C_tot_list=[50.0, 50.0]) # IS 48.25
# equilibration.solve_for_Na(target_pH=3.5, Na_low=0.0, Na_high=500.0, types=["acetate", "phosphate"], C_Cl=0.0, C_tot_list=[50.0, 50.0]) # IS 51.41
# equilibration.solve_for_Na(target_pH=3.8, Na_low=0.0, Na_high=500.0, types=["acetate", "phosphate"], C_Cl=0.0, C_tot_list=[50.0, 50.0]) # IS 55.07
# equilibration.solve_for_Na(target_pH=4.1, Na_low=0.0, Na_high=500.0, types=["acetate", "phosphate"], C_Cl=0.0, C_tot_list=[50.0, 50.0]) # IS 60.26
# equilibration.solve_for_Na(target_pH=4.4, Na_low=0.0, Na_high=500.0, types=["acetate", "phosphate"], C_Cl=0.0, C_tot_list=[50.0, 50.0]) # IS 67.57
