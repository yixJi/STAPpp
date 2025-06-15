import numpy as np
from pathlib import Path

user_file = "test3"
file_path = Path(fr"C:\Users\jiangyx\Desktop\hw\finite\STAPpp-master\src\data\{user_file}.out")

L = 10      
h = 2       
b = 1      
E = 10000   
P = 2      
I = b * h**3 / 12  
def u_exact(x, y):
    y_c = y - h / 2
    return P * y_c * x / (E * I)

def v_exact(x):
    return -P * h * x**2 / (2 * E * I)

lines = file_path.read_text(encoding="utf-8", errors="ignore").splitlines()

node_coords = {}
for idx, line in enumerate(lines):
    if "N O D A L   P O I N T   D A T A" in line:
        node_start = idx + 3  
        break

for line in lines[node_start:]:
    parts = line.strip().split()
    if len(parts) >= 5 and parts[0].isdigit():
        try:
            node_id = int(parts[0])
            x = float(parts[-2])
            y = float(parts[-1])
            node_coords[node_id] = (x, y)
        except:
            continue
    elif line.strip() == "":
        break

displacements = {}
reading_disp = False
for idx, line in enumerate(lines):
    if "X-DISPLACEMENT" in line and "Y-DISPLACEMENT" in line:
        reading_disp = True
        continue
    if reading_disp:
        parts = line.strip().split()
        if len(parts) >= 3 and parts[0].isdigit():
            try:
                node_id = int(parts[0])
                dx = float(parts[1])
                dy = float(parts[2])
                displacements[node_id] = (dx, dy)
            except:
                continue
        elif line.strip() == "":
            break

print("节点编号 |    x    |    y    |  u_calc  |  u_exact | u误差 |  v_calc  |  v_exact | v误差")
print("-------------------------------------------------------------------------------------")

u_errors = []
v_errors = []

for node in sorted(node_coords):
    x, y = node_coords[node]
    dx, dy = displacements.get(node, (0.0, 0.0))
    u_ex = u_exact(x, y)
    v_ex = v_exact(x)
    e_u = dx - u_ex
    e_v = dy - v_ex
    u_errors.append(e_u**2)
    v_errors.append(e_v**2)
    print(f"{node:9d} | {x:6.2f} | {y:6.2f} | {dx:8.5f} | {u_ex:8.5f} | {e_u:7.2e} | {dy:8.5f} | {v_ex:8.5f} | {e_v:7.2e}")

N = len(node_coords)
l2_total = np.sqrt((np.sum(u_errors) + np.sum(v_errors)) / N)
print("-------------------------------------------------------------------------------------")
print(f"节点总数: {N}")
print(f"L2误差 (总位移): {l2_total:.6e}")
