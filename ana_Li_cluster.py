# -*- coding: utf-8 -*-
import MDAnalysis as mda
from MDAnalysis.analysis import distances
import os
from collections import Counter
import numpy as np

# === CONFIGURATION ===
GRO_FILE = "md.gro" 
OUTPUT_DIR = "coordination_clusters" 
CUTOFF = 3.0 

LI_RES = "Li"
TARGET_RESNAMES = ["18C", "BF", "EtO", "THF"]

# === MAIN SCRIPT ===
if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

u = mda.Universe(GRO_FILE)

# 确保读取到了盒子信息，否则 PBC 无法工作
if u.dimensions is None:
    print("WARNING: No box dimensions found! PBC will be disabled.")

li_atoms = u.select_atoms(f"resname {LI_RES}")
all_target_atoms = u.select_atoms("resname " + " ".join(TARGET_RESNAMES))

print(f"Total Li: {len(li_atoms)}")
print(f"Total Target Atoms: {len(all_target_atoms)}")

# === DIAGNOSIS STEP ===
# 计算第 1 个 Li 离子到所有目标原子的最小距离
if len(li_atoms) > 0 and len(all_target_atoms) > 0:
    # 强制计算第一对距离
    dist_single = distances.distance_array(li_atoms[0:1].positions, 
                                          all_target_atoms.positions, 
                                          box=u.dimensions)
    print(f"DEBUG: Min distance for the first Li is: {np.min(dist_single):.4f} Angstrom")
# ======================

cluster_results = []
snapshot_counts = Counter()

for li in li_atoms:
    # 这里的 selection 逻辑如果还是不行，说明 select_atoms 内部的 around 没跑通
    # 我们换一种更稳健的写法：直接用 distance_array 过滤
    dists = distances.distance_array(li.position.reshape(1,3), 
                                    all_target_atoms.positions, 
                                    box=u.dimensions)
    
    # 找到距离小于 CUTOFF 的原子索引
    nearby_indices = np.where(dists[0] < CUTOFF)[0]
    
    if len(nearby_indices) == 0:
        tag = "Uncoordinated"
    else:
        # 获取这些原子的残基名
        resnames = all_target_atoms[nearby_indices].resnames
        tag = "-".join(sorted(list(set(resnames))))
    
    cluster_results.append(tag)
    
    # 保存逻辑
    if tag != "Uncoordinated" and snapshot_counts[tag] < 3:
        nearby_res = all_target_atoms[nearby_indices].residues
        (nearby_res.atoms + li).write(os.path.join(OUTPUT_DIR, f"type_{tag}_Li{li.id}.pdb"))
        snapshot_counts[tag] += 1

# 输出结果
stats = Counter(cluster_results)
for tag, count in stats.most_common():
    print(f"{tag:<30} | {count:<6} | {(count/len(li_atoms)*100):>7.2f}%")