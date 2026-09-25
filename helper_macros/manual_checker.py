from ij import ImagePlus, IJ
from ij.gui import NonBlockingGenericDialog
from ij.io import DirectoryChooser
from ij.plugin import LutLoader
import os

dc = DirectoryChooser("Select the folder with your images")
datadir = dc.getDirectory()
filenames = []

for obj in os.listdir(datadir):
	if os.path.isfile(datadir + "/" + obj):
		filenames.append(obj)

for f in filenames:
	imp = IJ.openImage(os.path.join(datadir, f))
	gbd = LutLoader.getLut("glasbey on dark")
	imp.getProcessor().setLut(gbd)
	imp.show()
		
	check = NonBlockingGenericDialog("Confirmation")
	check.addMessage("Check segmentation, then press okay")
	check.showDialog()
	
	imp.close()

print("Finished running")