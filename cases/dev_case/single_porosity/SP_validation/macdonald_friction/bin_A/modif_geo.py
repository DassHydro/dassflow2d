import numpy as np

###########################
## PARAMETERS DEFINITION ##
###########################

def h_ex(x, lx):
    return (1 + 0.5 * np.exp(-16 * (x / lx - 0.5) ** 2))

def int_friction(x, lx):
    N = 10000
    t = np.linspace(0, x, N)
    integrand = h_ex(t, lx)**(-10/3)
    result = np.trapz(integrand, t)
    return result

g = 9.81
q0 = 2
z0 = 0
lx = 100
n = 0.05

###########################
######## MODIF GEO ########
###########################

nodes = {}  # id : (x, y, z)

mesh_file = "channel_base.geo"

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

for line in lines:
    line_strip = line.strip()
    if line.startswith("#"):
        if line.startswith("# Nodes"):
            reading_nodes = True
            reading_cells = False
        elif line.startswith("# Cells"):
            reading_nodes = False
            reading_cells = True
        else : 
            reading_nodes = False
            reading_cells = False

    if reading_cells :
        if not line_strip or line_strip.startswith("#"):
            new_lines.append(line)
            continue

        parts = line_strip.split()
        cell_id = int(parts[0])
        first_node_id = int(parts[1])
        second_node_id = int(parts[2])

        x_first_node = nodes[first_node_id][0]  # récupérer x du premier noeud
        x_second_node = nodes[second_node_id][0] # récupérer x du second noeud
        x = 0.5 * (x_first_node + x_second_node) # position x du centre de la cellule
    
        h = h_ex(x, lx)

        if first_node_id == 1:
            u0 = q0 / h
            new_alt = z0
            c0 = 0.5*u0**2 + g*h + g*z0
        else :
            new_alt = c0/g - q0**2/(2*g*h**2) - h - (q0**2)*(n**2)*int_friction(x, lx) 

        # Remplacer la dernière colonne
        parts[-1] = f"{new_alt:.7E}"

        # Reconstituer la ligne
        new_line = "   ".join(parts)
        new_lines.append(new_line)
    else:
        new_lines.append(line)

print(c0)

with open("new_channel.geo", "w") as f:
    for line in new_lines:
        f.write(line.rstrip() + "\n")