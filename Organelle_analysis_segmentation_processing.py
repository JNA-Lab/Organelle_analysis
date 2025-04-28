from ij import IJ
from ij.gui import GenericDialog
from ij.io import DirectoryChooser
import os
import re

#get folder and file information
dc = DirectoryChooser("Select the folder with your segmented images")
dir = dc.getDirectory()
f = []
f_ext = set()
for root, dirnames, filenames in os.walk(dir):
	for filename in filenames:
		f.append(filename) #f.append(os.path.join(root, filename))
		f_ext.add(re.search(r"\.([A-Za-z0-9]+)$", filename).group(1))
#-----

#get processing options
organelles = ['nucleus', 'Golgi', 'peroxisomes', 'ER', 'mitochondria', 'lysosomes', 'bacteria']

options = GenericDialog('Options')
options.addChoice('File extension', list(f_ext), "tif")
options.addMessage("Select which organelles to analyse,and how they \nare represented in file names.");
options.addCheckbox('Cell boundaries (required)', True)
options.addToSameRow()
options.addStringField('', 'cell')
#for o in organelles:
	#options.addCheckbox(o, True)
	#options.addToSameRow()
	#options.addStringField('', o.lower())
#***TEMPORARY***
options.addCheckbox("nucleus", False)
options.addToSameRow()
options.addStringField("", "")
options.addCheckbox("Golgi", True)
options.addToSameRow()
options.addStringField("", "-Best_1dpi__2_golgi")
options.addCheckbox("peroxisomes", True)
options.addToSameRow()
options.addStringField("", "-Best_so_far_perox")
options.addCheckbox("ER", True)
options.addToSameRow()
options.addStringField("", "-Best_1dpi_ER")
options.addCheckbox("mitochondria", True)
options.addToSameRow()
options.addStringField("", "-Best_1dpi_2_mito")
options.addCheckbox("lysosomes", False)
options.addToSameRow()
options.addStringField("", "")
options.addCheckbox("bacteria", True)
options.addToSameRow()
options.addStringField("", "-Best_1dpi_Ot_LD")
#***************
options.addMessage("\n\n")
options.addCheckbox("Calculate pairwise overlaps?", True)
options.showDialog()

if options.wasOKed():
	ext = options.getNextChoice()
	cell_bool = options.getNextBoolean()
	if cell_bool != True:
		raise Exception("Cell boundaries are required for processing organelle segmentations.")
	cell_regex = options.getNextString()
	org_bool = dict()
	org_regex = dict()
	for o in organelles:
		org_bool[o] = options.getNextBoolean()
		org_regex[o] = options.getNextString()
	contacts_bool = options.getNextBoolean()
	
	#org_regex = [o for o, b in zip(org_regex, org_bool) if b]#removing match for non-selected organelles
	
else:
	print("Cancelled")#TODO actually handle something here
#-----

conditions = set()
print(len(filenames))
for f in filenames:#TODO - filter to selected extension
	f = f.replace("." + ext, "")
	print(f)
	for o in org_regex.values():
		f = f.replace(o, "")#not actually regex matching - probably best for now
	print(f)
	conditions.add(f)
print(len(conditions))


