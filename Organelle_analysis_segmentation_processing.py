from ij import IJ
from ij.gui import GenericDialog
from ij.io import DirectoryChooser
import os
import re

#get folder and file information
dc = DirectoryChooser("Select the folder with your segmented images")
datadir = dc.getDirectory()
f_ext = set()
filenames = []

for obj in os.listdir(datadir):
	print(obj)
	print(type(obj))
	if os.path.isfile(datadir + "/" + obj):
		print("found a file!")
		filenames.append(obj)
		f_ext.add(re.search(r"\.([A-Za-z0-9]+)$", obj).group(1))

	
#for root, dirnames, filenames in os.walk(dir):
	#for filename in filenames:
		#f_ext.add(re.search(r"\.([A-Za-z0-9]+)$", filename).group(1))

#create analysis folder for output
if not os.path.exists(datadir + "/analysis/"):
	os.mkdir(datadir + "/analysis/")
#-----

#get processing options
organelles = ['nucleus', 'Golgi', 'peroxisomes', 'ER', 'mitochondria', 'lysosomes', 'bacteria']

options = GenericDialog('Options')
options.addChoice('File extension', list(f_ext), "tif")
options.addMessage("Select which organelles to analyse,and how they \nare represented in file names.");
options.addCheckbox('Cell boundaries (required)', True)
options.addToSameRow()
#options.addStringField('', 'cells')
options.addStringField('', "_cp_masks")#TEMPORARY
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
	org_regex = dict()#not actually regex - just using str.replace for now
	for o in organelles:
		org_bool[o] = options.getNextBoolean()
		org_regex[o] = options.getNextString()
	contacts_bool = options.getNextBoolean()
	
	org_regex = {k:v for (k,v), b in zip(org_regex.items(), list(org_bool.values())) if b == True}

	
#-----

filenames_filtered = [f for f in filenames if re.search(ext + "$", f)]#filter to files with the correct extension

conditions = set()
for f in filenames_filtered:#TODO - filter to selected extension
	f = f.replace("." + ext, "")
	for o in org_regex.values():
		f = f.replace(o, "")#not actually regex matching - probably best for now
	f = f.replace(cell_regex, "")
	conditions.add(f)

print(len(org_regex))

for c in conditions:
	print(c)
	#open images
	c_filenames = [f for f in filenames if c in f]#find filenames matching condition
	print(len(c_filenames))
	for cf in c_filenames:
		IJ.open(os.path.join(datadir, cf))
		imp = IJ.getImage()
		
		try:
			organelle = [k for k,v in org_regex.items() if v in cf][0]
		except:
			if cell_regex in cf:
				organelle = "cells"
			else:
				imp.close()
		
		imp.setTitle(organelle)
		
