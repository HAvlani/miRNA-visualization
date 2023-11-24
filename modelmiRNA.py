import delegator

name = "hsa-miR-155-5p"
sequence = "UUAAUGCUAAUCGUGAUAGGGGUU"
fastaFileName = f'{name}_fasta.txt'
outputFileName = f'{name}.pdb'

with open(fastaFileName, mode="w") as file:
    file.write(f">{name}\n{sequence}")

model = delegator.run(f'obabel -ifasta {fastaFileName} -O {outputFileName}')
print(model.out)
