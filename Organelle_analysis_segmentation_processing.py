from ij import IJ
from ij.gui import GenericDialog, NonBlockingGenericDialog, Roi
from ij.io import DirectoryChooser
from ij.measure import ResultsTable
from ij.plugin import ImageCalculator, ImagesToStack, StackEditor
from ij.plugin.frame import RoiManager
import os
import re


#get folder and file information
dc = DirectoryChooser("Select the folder with your segmented images")
datadir = dc.getDirectory()
f_ext = set()
filenames = []

for obj in os.listdir(datadir):
	if os.path.isfile(datadir + "/" + obj):
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
organelles = ['nuclei', 'Golgi', 'peroxisomes', 'ER', 'mitochondria', 'lysosomes', 'bacteria']

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
options.addCheckbox("nuclei", False)
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

if options.wasOKed():#is this necessary?
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
	print(org_bool)
	print(org_regex)
	contacts_bool = options.getNextBoolean()
	
	org_regex = {k:v for (k,v), b in zip(org_regex.items(), list(org_bool.values())) if b == True}
	organelles_selected = []
	for o in organelles:
		if org_bool[o] == True:
			organelles_selected.append(o)
	print(organelles_selected)
	
	#ROI groups
	ROI_groups = dict(zip(organelles_selected, list(range(1, len(organelles_selected) + 1, 1))))
	print(ROI_groups)
	
	#pairwise contact groups and ROI groups
	if contacts_bool == True:
		pairwise_groups = {"ng": ("nuclei", "Golgi"),
							"np": ("nuclei", "peroxisomes"),
							"ne": ("nuclei", "ER"),
							"nm": ("nuclei", "mitochondria"),
							"nl": ("nuclei", "lysosomes"),
							"nb": ("nuclei", "bacteria"),
							"gp": ("Golgi", "peroxisomes"),
							"ge":  ("Golgi", "ER"),
							"gm": ("Golgi", "mitochondria"),
							"gl": ("Golgi", "lysosomes"),
							"gb": ("Golgi", "bacteria"),
							"pe": ("peroxisomes", "ER"),
							"pm": ("peroxisomes", "mitochondria"),
							"pl": ("peroxisomes", "lysosomes"),
							"pb": ("peroxisomes", "bacteria"),
							"em": ("ER", "mitochondria"),
							"el": ("ER", "lysosomes"),
							"eb": ("ER", "bacteria"),
							"ml": ("mitochondria", "lysosomes"),
							"mb": ("mitochondria", "bacteria"),
							"lb": ("lysosomes", "bacteria")}#TODO - ...
		pairwise_ROI_groups = dict(zip(pairwise_groups.keys(), range(max(ROI_groups.values()) + 1, max(ROI_groups.values()) + 1 + len(pairwise_groups), 1)))
		print(pairwise_ROI_groups)
	
#-----

filenames_filtered = [f for f in filenames if re.search(ext + "$", f)]#filter to files with the correct extension

#find conditions from filenames
conditions = set()
for f in filenames_filtered:#TODO - filter to selected extension
	f = f.replace("." + ext, "")
	for o in org_regex.values():
		f = f.replace(o, "")#not actually regex matching - probably best for now
	f = f.replace(cell_regex, "")
	conditions.add(f)
	



#set up ROI manager
RM = RoiManager()
rm = RM.getRoiManager()

##MAIN LOOP
for c in conditions:
	#open and rename images
	c_filenames = [f for f in filenames if c in f]#find filenames matching condition
	c_cell_image = str()
	for cf in c_filenames:
		IJ.open(os.path.join(datadir, cf))
		imp = IJ.getImage()		
		try:
			organelle = [k for k,v in org_regex.items() if v in cf][0]
		except:#index out of range - no match in organelles
			if cell_regex in cf:
				c_cell_image = cf
			imp.close()#close for now		
		imp.setTitle(organelle)
	
	images_to_stack = []
	#ORGANELLE THRESHOLDING
	for org in organelles_selected:
		print(org + " thresholding")
		IJ.selectWindow(org)
		org_img = IJ.getImage()
		IJ.setRawThreshold(org_img, 1, (2^org_img.getBitDepth())-1)#will work for 8 or 16-bit; might break on 32-bit or RGB
		IJ.run(org_img, "Convert to Mask", "")
		images_to_stack.append(org_img)
		
	#PAIRWISE OVERLAPS
	for combo in pairwise_groups.keys():
		if set(pairwise_groups[combo]).issubset(organelles_selected):
			IJ.selectWindow(pairwise_groups[combo][0])
			img1 = IJ.getImage()
			IJ.selectWindow(pairwise_groups[combo][1])
			img2 = IJ.getImage()
			img3 = ImageCalculator.run(img1, img2, "AND create")
			img3.setTitle(combo)
			images_to_stack.append(img3)
	
	#STACK
	orgstack = ImagesToStack.run(images_to_stack)
		
	#LOAD CELL IMAGE
	IJ.open(os.path.join(datadir, c_cell_image))
	cells = IJ.getImage()
	
	#Check max value
	cells_max = int(cells.getStatistics().max) #converting double to int - beware of rounding bugs	
	
	#Per step (cell) loop:
	for i in range(1, cells_max, 1):
		#duplicate cell image
		cells_copy = cells.duplicate()
		#step threshold
		IJ.setRawThreshold(cells_copy, i, i)
		#analyse particles
		IJ.run(cells_copy, "Analyze Particles...", "size=0-Infinity add composite") #should just be one cell
		#duplicate stack by ROI
		rm.select(0)#CHECK THIS
		cell_crop = orgstack.crop("stack")#CHECK if cropped or duplicated
		#re-add ROI, delete original
		rm.runCommand("Add")#how to specify auto-generated selection as ROI?
		rm.select(0)
		rm.runCommand("Delete")
		#select cell ROI, rename, clear outside
		rm.rename(0, "cell_" + str(i))#CHECK INDEXING
		rm.select(cell_crop, 0)
		IJ.run(cell_crop, "Clear Outside", "stack");
		#analyse particles per organelle/pair
		StackEditor.convertStackToImages(cell_crop)
		for org in organelles_selected:
			print(org + " thresholding")
			Roi.setDefaultGroup(ROI_groups[org])
			IJ.selectWindow(org)
			org_img = IJ.getImage()
			IJ.run(org_img, "Analyze Particles...", "size=0-Infinity add composite")
			#rename by slice/organelle/cell
			rm.selectGroup(ROI_groups[org])
			roi_start = rm.getSelectedIndex()
			roi_count = rm.selected()
			for i in range(roi_start, roi_start + roi_count, 1):
				rm.rename(i, org + "_" + str(i - roi_start + 1))
		for combo in pairwise_groups.keys():
			if set(pairwise_groups[combo]).issubset(organelles_selected):
				Roi.setDefaultGroup(pairwise_ROI_groups[combo])
				IJ.selectWindow(combo)
				combo_img = IJ.getImage()
				IJ.run(combo_img, "Analyze Particles...", "size=0-Infinity add composite")
				rm.selectGroup(pairwise_ROI_groups[combo])
				combo_roi_start = rm.getSelectedIndex()
				combo_roi_count = rm.selected()
				for i in range(combo_roi_start, combo_roi_start + combo_roi_count, 1):
					rm.rename(i, combo + "_" + str(i - combo_roi_start + 1))
		
		#save ROIs
		rm.runCommand("Select All")
		rm.save(datadir + "/analysis/" + c + "_cell" + i + "_ROIs.zip")
		#measure
		IJ.run("Set Measurements...", "area mean min centroid center shape feret's skewness kurtosis display redirect=None decimal=3")
		rm.runCommand("Select All")#should be redundant
		rm.runCommand("Measure")
		#save measurements
		results = ResultsTable.getResultsTable()
		results.save(datadir + "/analysis/" + c + "_cell" + i + "_results.csv")
		
		#clear ROIs, results, leave stack and cell image open (close duplicates)
		cell_crop.close()
		rm.reset()
		results.reset()
		
		
		
		
	
	
	
	
##******************OLD ORGANELLE PROCESSING CODE*******************************	
#	#TODO - loop this properly
#	if "nuclei" in organelles_selected:
#		print("nuclei block")
#		Roi.setDefaultGroup(ROI_groups["nuclei"])
#		IJ.selectWindow("nuclei")
#		nuc = IJ.getImage()
#		IJ.setRawThreshold(nuc, 1, 65535)#specific for 16-bit images
#		IJ.run(nuc, "Convert to Mask", "")
#		IJ.run(nuc, "Analyze Particles...", "size=0-Infinity add composite")
#		rm.selectGroup(ROI_groups["nucleus"])
#		nuc_roi_start = rm.getSelectedIndex()
#		nuc_roi_count = rm.selected()
#		for i in range(nuc_roi_start, nuc_roi_start + nuc_roi_count, 1):
#			rm.rename(i, "nucleus_" + str(i - nuc_roi_start + 1))
#	
#	if "Golgi" in organelles_selected:#TODO - lowercase?
#		print("golgi block")
#		Roi.setDefaultGroup(ROI_groups["Golgi"])
#		IJ.selectWindow("Golgi")
#		gol = IJ.getImage()
#		IJ.setRawThreshold(gol, 1, 65535)#specific for 16-bit images
#		IJ.run(gol, "Convert to Mask", "")
#		IJ.run(gol, "Analyze Particles...", "size=0-Infinity add composite")
#		rm.selectGroup(ROI_groups["Golgi"])
#		gol_roi_start = rm.getSelectedIndex()
#		gol_roi_count = rm.selected()
#		for i in range(gol_roi_start, gol_roi_start + gol_roi_count, 1):
#			rm.rename(i, "Golgi_" + str(i - gol_roi_start + 1))
#	
#	if "peroxisomes" in organelles_selected:
#		print("perox block")
#		Roi.setDefaultGroup(ROI_groups["peroxisomes"])
#		IJ.selectWindow("peroxisomes")
#		per = IJ.getImage()
#		IJ.setRawThreshold(per, 1, 65535)#specific for 16-bit images
#		IJ.run(per, "Convert to Mask", "")
#		IJ.run(per, "Analyze Particles...", "size=0-Infinity add composite")
#		rm.selectGroup(ROI_groups["peroxisomes"])
#		per_roi_start = rm.getSelectedIndex()
#		per_roi_count = rm.selected()
#		for i in range(per_roi_start, per_roi_start + per_roi_count, 1):
#			rm.rename(i, "perox_" + str(i - per_roi_start + 1))
#			
#	if "ER" in organelles_selected:#TODO - lowercase?
#		print("ER block")
#		Roi.setDefaultGroup(ROI_groups["ER"])
#		IJ.selectWindow("ER")
#		er = IJ.getImage()
#		IJ.setRawThreshold(er, 1, 65535)#specific for 16-bit images
#		IJ.run(er, "Convert to Mask", "")
#		IJ.run(er, "Analyze Particles...", "size=0-Infinity add composite")
#		rm.selectGroup(ROI_groups["ER"])
#		er_roi_start = rm.getSelectedIndex()
#		er_roi_count = rm.selected()
#		for i in range(er_roi_start, er_roi_start + er_roi_count, 1):
#			rm.rename(i, "ER_" + str(i - er_roi_start + 1))
#
#	if "mitochondria" in organelles_selected:
#		print("mito block")
#		Roi.setDefaultGroup(ROI_groups["mitochondria"])
#		IJ.selectWindow("mitochondria")
#		mit = IJ.getImage()
#		IJ.setRawThreshold(mit, 1, 65535)#specific for 16-bit images
#		IJ.run(mit, "Convert to Mask", "")
#		IJ.run(mit, "Analyze Particles...", "size=0-Infinity add composite")
#		rm.selectGroup(ROI_groups["mitochondria"])
#		mit_roi_start = rm.getSelectedIndex()
#		mit_roi_count = rm.selected()
#		for i in range(mit_roi_start, mit_roi_start + mit_roi_count, 1):
#			rm.rename(i, "mito_" + str(i - mit_roi_start + 1))
#
#	if "lysosomes" in organelles_selected:
#		print("lyso block")
#		Roi.setDefaultGroup(ROI_groups["lysosomes"])
#		IJ.selectWindow("lysosomes")
#		lys = IJ.getImage()
#		IJ.setRawThreshold(lys, 1, 65535)#specific for 16-bit images
#		IJ.run(lys, "Convert to Mask", "")
#		IJ.run(lys, "Analyze Particles...", "size=0-Infinity add composite")
#		rm.selectGroup(ROI_groups["lysosomes"])
#		lys_roi_start = rm.getSelectedIndex()
#		lys_roi_count = rm.selected()
#		for i in range(lys_roi_start, lys_roi_start + lys_roi_count, 1):
#			rm.rename(i, "lyso_" + str(i - lys_roi_start + 1))
#
#	if "bacteria" in organelles_selected:
#		print("bac block")
#		Roi.setDefaultGroup(ROI_groups["bacteria"])
#		IJ.selectWindow("bacteria")
#		bac = IJ.getImage()
#		IJ.setRawThreshold(bac, 1, 65535)#specific for 16-bit images
#		IJ.run(bac, "Convert to Mask", "")
#		IJ.run(bac, "Analyze Particles...", "size=0-Infinity add composite")
#		rm.selectGroup(ROI_groups["bacteria"])
#		bac_roi_start = rm.getSelectedIndex()
#		bac_roi_count = rm.selected()
#		for i in range(bac_roi_start, bac_roi_start + bac_roi_count, 1):
#			rm.rename(i, "bac_" + str(i - bac_roi_start + 1))
#		
#	
#	if contacts_bool == True:
#		pairwise_groups = {"ng": ("nuclei", "Golgi"),
#							"np": ("nuclei", "peroxisomes"),
#							"ne": ("nuclei", "ER"),
#							"nm": ("nuclei", "mitochondria"),
#							"nl": ("nuclei", "lysosomes"),
#							"nb": ("nuclei", "bacteria"),
#							"gp": ("Golgi", "peroxisomes"),
#							"ge":  ("Golgi", "ER"),
#							"gm": ("Golgi", "mitochondria"),
#							"gl": ("Golgi", "lysosomes"),
#							"gb": ("Golgi", "bacteria"),
#							"pe": ("peroxisomes", "ER"),
#							"pm": ("peroxisomes", "mitochondria"),
#							"pl": ("peroxisomes", "lysosomes"),
#							"pb": ("peroxisomes", "bacteria"),
#							"em": ("ER", "mitochondria"),
#							"el": ("ER", "lysosomes"),
#							"eb": ("ER", "bacteria"),
#							"ml": ("mitochondria", "lysosomes"),
#							"mb": ("mitochondria", "bacteria"),
#							"lb": ("lysosomes", "bacteria")}#TODO - ...
#		pairwise_ROI_groups = dict(zip(pairwise_groups.keys(), range(max(ROI_groups.values()) + 1, max(ROI_groups.values()) + 1 + len(pairwise_groups), 1)))
#		print(pairwise_ROI_groups)
#		
#		for combo in pairwise_groups.keys():
#			if set(pairwise_groups[combo]).issubset(organelles_selected):
#				Roi.setDefaultGroup(pairwise_ROI_groups[combo])
#				IJ.selectWindow(pairwise_groups[combo][0])
#				img1 = IJ.getImage()
#				IJ.selectWindow(pairwise_groups[combo][1])
#				img2 = IJ.getImage()
#				img3 = ImageCalculator.run(img1, img2, "AND create")
#				IJ.run(img3, "Analyze Particles...", "size=0-Infinity add composite")
#				rm.selectGroup(pairwise_ROI_groups[combo])
#				combo_roi_start = rm.getSelectedIndex()
#				combo_roi_count = rm.selected()
#				for i in range(combo_roi_start, combo_roi_start + combo_roi_count, 1):
#					rm.rename(i, combo + "_" + str(i - combo_roi_start + 1))
#		
#	#************************************************
#	
#	
#	#Save ROIs
#	rm.runCommand("Select All")
#	rm.save(datadir + "/analysis/" + c + "_ROIs.zip")
#	#Measure
#	IJ.run("Set Measurements...", "area mean min centroid center shape feret's skewness kurtosis display redirect=None decimal=3")
#	rm.runCommand("Select All")#should be redundant
#	rm.runCommand("Measure")
#	results = ResultsTable.getResultsTable()
#	results.save(datadir + "/analysis/" + c + "_results.csv")
#	#TODO - clean up table
#	
#	pause = NonBlockingGenericDialog('Pause')
#	pause.addMessage('Click OK when ready')
#	pause.showDialog()
#	IJ.run("Close All", "")
#	rm.reset()
#	results.reset()
#	
#	if "nucleus" in organelles_selected:
#		rm.setGroup(ROI_groups["nucleus"])
#		IJ.selectWindow("nucleus")
#		nuc = IJ.getImage()
#		IJ.setRawThreshold(nuc, 1, 65535)#specific for 16-bit images
#		IJ.run(nuc, "Convert to Mask", "")
#		IJ.run(nuc, "Analyze Particles...", "size=0-Infinity add composite")
#		rm.selectGroup(ROI_groups["nucleus"])
#		nuc_roi_start = rm.getIndex()
		
