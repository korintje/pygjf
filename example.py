import pygjf

mol = pygjf.GaussianInput.load_from_file("EtN3_ESP.gjf")
mol.hex_repeat(4.5, 1, 2, 0)
# for atom in mol.atoms:
#     print(atom["element"])
#     print(atom["coord"])
mol.to_file("test.gjf")

