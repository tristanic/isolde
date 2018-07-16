import clipper

f = clipper.MMDBfile()
f.read_file ( "test.pdb")
mmol = clipper.MiniMol ()
print mmol
f.import_minimol ( mmol )


f.write_file("silly.pdb",0)
f.write_file("silly.cif",clipper.MMDBManager.CIF)

print dir(f)
print dir(mmol)

atoms = mmol.atom_list()

print dir(atoms)
print atoms[0].coord_orth().x()
print len(atoms)

for i in range(len(atoms)):
  print atoms[i]

for at in atoms:
  c = at.coord_orth()
  print c.x(), c.y(), c.z()

mod = mmol.model()

for poly in mod:
  for mono in poly:
    for atom in mono:
      print atom.coord_orth().x(), atom.coord_orth().y(), atom.coord_orth().z()

print mmol[0][0][0]
