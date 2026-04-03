import MDAnalysis as mda
from MDAnalysis.lib.distances import capped_distance
import numpy as np
from collections import defaultdict

# ================= Configuration =================
TOPOLOGY = "md.gro"  
TRAJECTORY = "md.xtc"
CUTOFF = 3.0  # Cutoff radius in A
# =================================================

print("Loading trajectory...")
u = mda.Universe(TOPOLOGY, TRAJECTORY)

# Atom selections
li_atoms = u.select_atoms("resname Li")
o_atoms = u.select_atoms("name O*")
bf_f_atoms = u.select_atoms("resname BF and name F*")

if len(li_atoms) == 0:
    raise ValueError("No atoms with resname Li found. Please check the topology.")

# Variables for statistics
o_coord_counts = defaultdict(int)
f_coord_count = 0

n_frames = len(u.trajectory)
n_li = len(li_atoms)
print(f"Found {n_li} Li atoms. Total frames: {n_frames}.")
print(f"Analyzing coordination within {CUTOFF} A, please wait...")

# Iterate through the trajectory
for ts in u.trajectory:
    # 1. Find O* around Li
    pairs_o = capped_distance(li_atoms.positions, o_atoms.positions, 
                              max_cutoff=CUTOFF, return_distances=False)
    
    if len(pairs_o) > 0:
        o_indices = pairs_o[:, 1]
        found_resnames = o_atoms.resnames[o_indices]
        unique_res, counts = np.unique(found_resnames, return_counts=True)
        for res, count in zip(unique_res, counts):
            o_coord_counts[res] += count

    # 2. Find F* (from BF) around Li
    pairs_f = capped_distance(li_atoms.positions, bf_f_atoms.positions, 
                              max_cutoff=CUTOFF, return_distances=False)
    f_coord_count += len(pairs_f)

# ================= Data Processing and Output =================
total_li_observations = n_li * n_frames

print("\n" + "="*50)
print("Average Coordination Number (per Li atom)")
print("="*50)

total_o_coord = 0
for res_name, count in o_coord_counts.items():
    avg_coord = count / total_li_observations
    total_o_coord += avg_coord
    print(f"O* from [{res_name}]: {avg_coord:.4f}")

avg_f_coord = f_coord_count / total_li_observations
print(f"F* from [BF] : {avg_f_coord:.4f}")

total_coord = total_o_coord + avg_f_coord
print("-" * 50)
print(f"Total average coordination number of Li: {total_coord:.4f}")

print("\n" + "="*50)
print("First Coordination Shell Composition (Total = 100%)")
print("="*50)
if total_coord > 0:
    for res_name, count in o_coord_counts.items():
        avg_coord = count / total_li_observations
        pct = (avg_coord / total_coord) * 100
        print(f"O* from [{res_name}] ratio: {pct:.2f} %")
    
    f_pct = (avg_f_coord / total_coord) * 100
    print(f"F* from [BF] ratio: {f_pct:.2f} %")
else:
    print("No coordinating atoms found within the specified cutoff.")
print("="*50)