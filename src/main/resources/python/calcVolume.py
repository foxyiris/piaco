import sys
sys.path.insert(0, sys.argv[2])

import os
import chimeraInit
# https://gist.github.com/jaimergp/834935666066515246ca
basepath = os.path.dirname(os.path.abspath(sys.argv[0]))
# This is the magic!
# If a script needs to be executed, Chimera does not launch
# the in-house command interface! So, just passing a `pass`
# statement does the trick.
chimera_pass = os.path.join(basepath, 'pass.py')
chimeraInit.init(['@jaimergp', '--script', chimera_pass], nogui=True, eventloop=False, exitonquit=False)

from chimera import runCommand as rc
from Surfnet import show_surfnet
from MeasureVolume import surface_volume_and_area
import _surfnet
import Midas
import chimera

def interface_surfnet(receptorOSL, ligandsOSL, useMesh=True,
			interval=1.0, cutoff=10, density='Gaussian',
			color=None):
	try:
		receptorAtomList = Midas._selectedAtoms(receptorOSL)
	except (Midas.MidasError, chimera.oslParser.OSLSyntaxError), e:
		return str(e)
	if not receptorAtomList:
		return 'No receptor atoms selected'
	try:
		ligandsAtomList = Midas._selectedAtoms(ligandsOSL)
	except (Midas.MidasError, chimera.oslParser.OSLSyntaxError), e:
		return str(e)
	if not ligandsAtomList:
		return 'No ligand atoms selected'
	surfmol, mol = _surfnet.interface_surfnet(receptorAtomList,
						ligandsAtomList,
						cutoff=cutoff,
						threshold=interval)
	#print '%d spheres' % len(surfmol.atoms)
	mList = show_surfnet(surfmol, useMesh, interval, density,
						color, sameAs=mol)
	for m in mList:
		m.name = 'Surfnet for "%s" - "%s"' % (receptorOSL, ligandsOSL)
	surfmol.destroy()
	return mList

configuration = sys.argv[1]
f = open(configuration)
lines=f.readlines()
f.close

filename = lines[0]
receptor = lines[1]
ligand   = lines[2]
rc("open" + filename)

mList = interface_surfnet(receptor, ligand, interval=0.8, cutoff=5.5)

for m in mList:
    volume, area, hole_count = surface_volume_and_area(m)
    print "Volume=", volume
