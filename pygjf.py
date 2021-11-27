import math
import numpy as np

def cart2sph(cart):
    XsqPlusYsq = cart[0]**2 + cart[1]**2
    r = math.sqrt(XsqPlusYsq + cart[2]**2)
    elev = math.atan2(cart[2], math.sqrt(XsqPlusYsq))
    az = math.atan2(cart[1], cart[2])
    return np.array([r, elev, az])

def sph2cart(sph):
    radius = sph[0]
    elevation = sph[1]
    azimuth = sph[2]
    ax = np.cos(azimuth) * np.sin(np.pi/2 - elevation) * radius
    ay = np.sin(azimuth) * np.sin(np.pi/2 - elevation) * radius
    az = np.cos(np.pi/2 - elevation) * radius
    return np.array([ax, ay, az])

def get_hexagonal_cart(distance, layer_num, dimension=2):
    if dimension == 2:
        z = []
    elif dimension == 3:
        z = [0.0]
    else:
        raise("dimension must be 2 or 3")
    coords = []
    d = distance * layer_num
    sp = math.sin(math.pi / 3.0)
    cp = math.cos(math.pi / 3.0)
    sm = -1 * math.sin(math.pi / 3.0)
    cm = -1 * math.cos(math.pi / 3.0)
    vertexes = [
        np.array([d *  1.0, 0.0] + z), 
        np.array([d * cp, d * sp] + z), 
        np.array([d * cm, d * sp] + z), 
        np.array([d * -1.0, 0.0] + z), 
        np.array([d * cm, d * sm] + z), 
        np.array([d * cp, d * sm] + z), 
    ]
    for i, vertex in enumerate(vertexes):
        for m in range(layer_num):
            next_idx = (i + 1) % 6
            edge = vertexes[next_idx] - vertex
            coords.append(vertex + m / layer_num * edge)
    return coords


class GaussianInput():

    def __init__(self, headers, title, total_charge, multiplicity, atoms):
        self.headers = headers
        self.title = title
        self.total_charge = total_charge
        self.multiplicity = multiplicity
        self.atoms = atoms
    
    @classmethod
    def load_from_file(cls, filepath):
        with open(filepath, "r", encoding="utf-8") as f:
            string = f.read()
        return cls.load_from_string(string)

    @classmethod
    def load_from_string(cls, string):
        lines = string.split("\n")
        seps = [i for i, n in enumerate(lines) if not n]
        sep_1st = seps[0]
        sep_2nd = seps[1]
        headers = lines[0:sep_1st]
        title = lines[sep_1st + 1]
        charge_multiplicity = lines[sep_2nd + 1].split()
        total_charge = int(charge_multiplicity[0])
        multiplicity = int(charge_multiplicity[1])
        _atoms = [line.split() for line in lines[sep_2nd + 2:] if line]
        atoms = []
        for _atom in _atoms:
            atom = {}
            atom["element"] = _atom[0]
            atom["coord"] = np.array([float(coord) for coord in _atom[1:4]])
            atoms.append(atom)
        return cls(headers, title, total_charge, multiplicity, atoms)

    @classmethod
    def is_valid(cls):
        pass

    def to_string(self):
        atom_lines = []
        for atom in self.atoms:
            x, y, z = atom["coord"][0], atom["coord"][1], atom["coord"][2]
            element = atom["element"]
            atom_line = " {}                 {: .8f}   {: .8f}   {: .8f}".format(element, x, y, z)
            atom_lines.append(atom_line)
        return "\n".join(self.headers) \
            + "\n\n" + self.title \
            + "\n\n" + str(self.total_charge) + " " + str(self.multiplicity) \
            + "\n" + "\n".join(atom_lines) \
            + "\n"

    def to_file(self, filepath):
        string = self.to_string()
        with open(filepath, "w") as f:
            f.write(string)

    def rotate(self, origin, azimuth=0.0, elevation=0.0):
        for atom in self.atoms:
            cart = atom["coord"]
            rel_cart = cart - origin
            rel_sph = cart2sph(rel_cart)
            rot_sph = rel_sph + np.array([0.0, elevation, azimuth])
            rot_cart = sph2cart(rot_sph)
            atom["coord"] = rot_cart
    
    def clone(self, atom_idx, cart):
        center_atom = self.atoms[atom_idx]
        delta = cart - center_atom["coord"]
        atoms = []
        for atom in self.atoms[:]:
            atom["coord"] = atom["coord"] + delta
            atoms.append(atom)
        self.atoms += atoms

    def move(self, atom_idx, cart):
        delta = cart - self.atoms[atom_idx]["coord"]
        new_atoms = []
        for atom in self.atoms:
            new_atoms.append({"element": atom["element"], "coord": (atom["coord"] + delta)})
        self.atoms = new_atoms
    
    def clean_atoms(self):
        self.atoms = []

    def hex_repeat(self, distance, layer, origin_idx, perp_origin_idx):
        oc = self.atoms[origin_idx]["coord"]
        pc = self.atoms[perp_origin_idx]["coord"]
        origins = get_hexagonal_cart(distance, layer, dimension=3)
        rot_v = cart2sph(pc - oc)
        # rot_el = math.pi / 2 - rot_v[1]
        _atoms = self.atoms[:]
        self.clean_atoms()
        for i, origin_cart in enumerate(origins):
            _tmp_atoms = _atoms
            _tmp_mol = GaussianInput([], "", 0, 1, _tmp_atoms)
            _tmp_mol.move(origin_idx, origin_cart)
            self.atoms.extend(_tmp_mol.atoms)
        