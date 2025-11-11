import numpy as np

nodes = {}  # id : (x, y, z)

mesh_file = "bin_A\\old_channel.geo"

with open(mesh_file, "r") as f:
    lines = f.readlines()

reading_nodes = False
reading_cells = False

for line in lines:
    line = line.strip()

    if line.startswith("#"):
        if line.startswith("# Nodes"):
            reading_nodes = True
            reading_cells = False
        elif line.startswith("# Cells"):
            reading_nodes = False
            reading_cells = True
    elif reading_nodes and line:
        parts = line.split()        
        node_id = int(parts[0])
        x = float(parts[1])
        y = float(parts[2])
        z = float(parts[3])
        nodes[node_id] = (x, y, z)


new_lines = []
reading_cells = False

g=9.81
q0 = 2
z0 = 0

for line in lines:
    line_strip = line.strip()
    if line.startswith("#"):
        if line_strip.startswith("# Cells"):
            reading_cells = True
        else :
            reading_cells = False

    if reading_cells :
        if not line_strip or line_strip.startswith("#"):
            new_lines.append(line)
            continue

        parts = line_strip.split()
        cell_id = int(parts[0])
        first_node_id = int(parts[1])
        x_first_node = nodes[first_node_id][0]  # récupérer x du premier noeud
        print(x_first_node)
        h = (4/g)**(1/3)*(1+0.5*np.exp(-16*x_first_node/1000 -1))**2

        if first_node_id == 1:
            u0 = q0 / h
            z0 = 0
            new_alt = z0
            c0 = 0.5*u0**2 + g*(h + z0)
        else:
            new_alt = c0/g - q0**2/(2*g*h**2) - h

        # Remplacer la dernière colonne
        parts[-1] = f"{new_alt:.7E}"

        # Reconstituer la ligne
        new_line = "   ".join(parts)
        new_lines.append(new_line)
    else:
        new_lines.append(line)

with open("channel.geo", "w") as f:
    for line in new_lines:
        f.write(line.rstrip() + "\n")