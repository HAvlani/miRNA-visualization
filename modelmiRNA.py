import delegator
import pandas as pd
from Bio.PDB import PDBParser, PDBIO

name = ""
sequence = ""
fastaFileName = f"{name}_fasta.txt"
pdbFileName = f"{name}.pdb"

with open(fastaFileName, mode="w") as fastaFile:
    fastaFile.write(f">{name}\n{sequence}")

model = delegator.run(f"obabel -ifasta {fastaFileName} -O {pdbFileName}")
print(model.out)

with open(pdbFileName, mode="r") as pdbFile:
    pdb = pdbFile.read()

columns = [
    "record_type",
    "serial_number",
    "atom_name",
    "residue_name",
    "chain_id",
    "residue_number",
    "x",
    "y",
    "z",
    "occupancy",
    "temperature_factor",
    "element",
    "charge",
]
data = []

# Read each line and extract information
for line in pdb.strip().split("\n"):
    record_type = line[:6].strip()
    if record_type in ("ATOM", "HETATM"):
        data.append(
            [
                record_type,
                int(line[6:11]),
                line[12:16].strip(),
                line[17:20].strip(),
                line[21],
                int(line[22:26]),
                float(line[30:38]),
                float(line[38:46]),
                float(line[46:54]),
                float(line[54:60]),
                float(line[60:66]),
                line[76:78].strip(),
                line[78:80].strip(),
            ]
        )

# Create the DataFrame
df = pd.DataFrame(data, columns=columns)

basestring = df["atom_name"].to_list()
basestring = "".join(basestring)

templist = basestring.split("H" * 11)
baselist = []
for i in templist:
    x = i.split("H" * 10)
    if len(x[0]) > 4:
        baselist.append(x)

baselist = sum(baselist, [])
namelist = [*name]

newbaselist = []
for i in range(len(namelist)):
    repeats = len(baselist[i])
    newbaselist.append(namelist[i] * repeats)

for i in range(len(baselist)):
    if baselist[i] in basestring:
        basestring = basestring.replace(baselist[i], newbaselist[i])

finalbaselist = [*basestring]
io = PDBIO()
p = PDBParser()
structure = p.get_structure(name, pdbFileName)  # input files
i = 0
for model in structure:
    for chain in model:
        for residue in chain:
            residue.resname = finalbaselist[i]
            i += 1
io.set_structure(structure)
io.save("out.pdb", preserve_atom_numbering=True)

