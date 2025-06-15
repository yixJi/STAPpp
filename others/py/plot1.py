import re
import numpy as np
import matplotlib.pyplot as plt
import os
from pathlib import Path

user_file = "patch_test"  
scale = 60        
filename = Path(fr'C:\Users\jiangyx\Desktop\hw\finite\STAPpp-master\src\data\{user_file}.out')

# === 读取文件内容 ===
with open(filename, 'r', encoding='utf-8', errors='ignore') as f:
    lines = f.readlines()

# === 定位各段开始行 ===
node_start = elem_start = disp_start = None
for i, line in enumerate(lines):
    if 'N O D A L   P O I N T   D A T A' in line:
        node_start = i + 3
    if 'E L E M E N T   I N F O R M A T I O N' in line:
        elem_start = i + 2
    if 'D I S P L A C E M E N T S' in line:
        disp_start = i + 2

# === 读取节点坐标（Q4） ===
def read_nodes_q4(lines, node_start):
    node_coords = []
    for line in lines[node_start:]:
        if re.match(r'\s*\d+\s+\d+\s+\d+', line):
            parts = line.split()
            x, y = map(float, parts[-2:])
            node_coords.append([x, y])
        elif line.strip() == '':
            break
    return np.array(node_coords)

# === 读取单元连通性（Q4） ===
def read_elem_connect_q4(lines, elem_start):
    elem_connect = []
    for line in lines[elem_start:]:
        if re.match(r'\s*\d+\s+\d+\s+\d+\s+\d+\s+\d+\s+\d+', line):
            parts = line.split()
            elem_connect.append([int(parts[1]), int(parts[2]), int(parts[3]), int(parts[4])])
        elif line.strip() == '':
            break
    return elem_connect

# === 读取位移（Q4） ===
def read_disp_q4(lines, disp_start):
    displacements = []
    for line in lines[disp_start:]:
        if re.match(r'\s*\d+\s+[-\deE.+]+\s+[-\deE.+]+', line):
            parts = line.split()
            dx = float(parts[1])
            dy = float(parts[2])
            displacements.append([dx, dy])
        elif line.strip() == '':
            break
    return np.array(displacements)

# === 读取数据 ===
node_coords = read_nodes_q4(lines, node_start)
elem_connect = read_elem_connect_q4(lines, elem_start)
displacements = read_disp_q4(lines, disp_start)

# === 绘图函数 ===
def plot_mesh(node_coords, elem_connect, title, color='b', elem_color='g', node_num_color='k'):
    for elem in elem_connect:
        pts = node_coords[np.array(elem) - 1]
        pts = np.vstack([pts, pts[0]])  # 封闭多边形
        plt.plot(pts[:, 0], pts[:, 1], color + '-')
    plt.scatter(node_coords[:, 0], node_coords[:, 1], c=color)
    for i, (x, y) in enumerate(node_coords):
        plt.text(x, y, str(i + 1), color=node_num_color, fontsize=10)
    plt.title(title)
    plt.axis('equal')
    plt.xlabel("x")
    plt.ylabel("y")
    plt.grid(True)

# === 生成图像并保存 ===
output_dir = filename.parent
img_base = output_dir / user_file

plt.figure(figsize=(8, 6))
plot_mesh(node_coords, elem_connect, 'Q4 Element Original Mesh', color='b')
plt.savefig(f"{img_base}_original.png", dpi=300)
plt.show()

plt.figure(figsize=(8, 6))
deformed_coords = node_coords + scale * displacements
plot_mesh(deformed_coords, elem_connect, f'Deformed Mesh (scale={scale})', color='r')
plt.savefig(f"{img_base}_deformed.png", dpi=300)
plt.show()
