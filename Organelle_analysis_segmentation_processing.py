from ij import IJ
from ij.gui import GenericDialog, NonBlockingGenericDialog, Roi
from ij.io import DirectoryChooser
from ij.measure import ResultsTable
from ij.plugin import ImageCalculator, StackEditor, ImagesToStack, SubstackMaker, ZProjector
from ij.plugin.frame import RoiManager
from ij.plugin.filter import ParticleAnalyzer
from java.awt import Color
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
#***TEMPORARY***#TODO - variables for saving in custom values
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
	contacts_bool = options.getNextBoolean()
	
	org_regex = {k:v for (k,v), b in zip(org_regex.items(), list(org_bool.values())) if b == True}
	organelles_selected = []
	for o in organelles:
		if org_bool[o] == True:
			organelles_selected.append(o)
	
	#ROI groups
	ROI_groups = dict(zip(organelles_selected, list(range(2, len(organelles_selected) + 2, 1))))
	ROI_groups["cells"] = 1
	
	#pairwise contact groups and ROI groups
	if contacts_bool == True:
		pairwise_groups = {"ng": ("nuclei", "Golgi"),
							"np": ("nuclei", "peroxisomes"),
							"ne": ("nuclei", "ER"),
							"nm": ("nuclei", "mitochondria"),
							"nl": ("nuclei", "lysosomes"),
							"nb": ("nuclei", "bacteria"),
							"gp": ("Golgi", "peroxisomes"),
							"ge": ("Golgi", "ER"),
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
SE = StackEditor()
I2S = ImagesToStack()
SM = SubstackMaker()
PA = ParticleAnalyzer()

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
			imp.setTitle(organelle)
		except:#index out of range - no match in organelles
			if cell_regex in cf:
				c_cell_image = cf
			imp.close()#close for now		
		
	
	images_to_stack = []
	#ORGANELLE THRESHOLDING
	for org in organelles_selected:
		IJ.selectWindow(org)
		org_img = IJ.getImage()
		IJ.setRawThreshold(org_img, 1, (2**org_img.getBitDepth())-1)#will work for 8 or 16-bit; might break on 32-bit or RGB
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
#	NonBlockingGenericDialog("BREAK").showDialog()#debugging
	orgstack = I2S.run(images_to_stack)
	orgstack.setTitle("Organelles")
	orgstack.show()
		
	#LOAD CELL IMAGE
	IJ.open(os.path.join(datadir, c_cell_image))
	cells = IJ.getImage()
	cells.setTitle("Cells")
	
	#Check max value
	cells_max = int(cells.getStatistics().max) #converting double to int - beware of rounding bugs
	
	#Per step (cell) loop:
	for i in range(1, cells_max + 1, 1):
		#duplicate cell image
		cells_copy = cells.duplicate()
		cells_copy.setTitle("cells_DUPLICATE")
		cells_copy.show()
		cells.hide()
		#step threshold
		IJ.setRawThreshold(cells_copy, i, i)
		#analyse particles
		Roi.setDefaultGroup(ROI_groups["cells"])#TODO - dynamic?
		IJ.run(cells_copy, "Analyze Particles...", "size=0-Infinity add composite") #should just be one cell
		#SOME CELLS RETURNING >1 - combine
		if rm.getCount() > 1:
			rm.selectGroup(ROI_groups["cells"])			
			rm.runCommand("Combine")
			rm.addRoi(cells_copy.getRoi())
			#delete non-combined ROIs
			rm.setSelectedIndexes(range(0, rm.getCount()-1, 1))#all but last
			rm.runCommand("Delete")
		#duplicate stack and clear outside ROI
		rm.select(0)
		cell_id = "cell_" + str(i)
		cell_crop = orgstack.crop([rm.getRoi(0)], "stack")[0]
		cell_crop.show()
		orgstack.hide()
		rm.addRoi(cell_crop.getRoi())#cell boundary in cropped image
		rm.select(0)
		rm.runCommand("Delete")#deleting original cell boundary - now in wrong place in cropped image
		rm.rename(0, cell_id)
		rm.select(cell_crop, 0)
		IJ.run(cell_crop, "Clear Outside", "stack")
		
		#analyse particles per organelle/pair
		SE.convertStackToImages(cell_crop)
		for org in organelles_selected:
			Roi.setDefaultGroup(ROI_groups[org])
			IJ.selectWindow(org)
			org_img = IJ.getImage()
			IJ.run(org_img, "Analyze Particles...", "size=0-Infinity add composite")
			#rename by slice/organelle/cell
			rm.deselect()
			rm.selectGroup(ROI_groups[org])
			roi_start = rm.getSelectedIndex()
			roi_count = rm.selected()
			if roi_count > 0:
				for i in range(roi_start, roi_start + roi_count, 1):
					rm.rename(i, cell_id + "_" + org + "_" + str(i - roi_start + 1))
			org_img.close()#TODO - move below?
		combo_present = []
		for combo in pairwise_groups.keys():
			if set(pairwise_groups[combo]).issubset(organelles_selected):
				combo_present.append(combo)
				Roi.setDefaultGroup(pairwise_ROI_groups[combo])
				IJ.selectWindow(combo)
				combo_img = IJ.getImage()
				IJ.run(combo_img, "Analyze Particles...", "size=0-Infinity add composite")
				rm.deselect()
				rm.selectGroup(pairwise_ROI_groups[combo])
				combo_roi_start = rm.getSelectedIndex()
				combo_roi_count = rm.selected()
				if combo_roi_count > 0:
					for i in range(combo_roi_start, combo_roi_start + combo_roi_count, 1):
						rm.rename(i, cell_id + "_" + combo + "_" + str(i - combo_roi_start + 1))
				combo_img.close()#TODO - move below?	
				
		#save cropped stack
		cell_crop.show()#for pooled ROI measurements
		IJ.saveAs(cell_crop, "Tiff", datadir + "/analysis/" + c + "_" + cell_id + "_stack.tif")
		
		#save ROIs
		rm.runCommand("Select All")
		rm.save(datadir + "/analysis/" + c + "_" + cell_id + "_ROIs.zip")
		
		#measure
		IJ.run("Set Measurements...", "area mean min centroid center shape feret's skewness kurtosis display redirect=None decimal=3")
		rm.runCommand("Select All")#should be redundant
		rm.runCommand("Measure")				
		results = ResultsTable.getResultsTable()
		#add groups for processing in R
		full_ROI_groups = dict(ROI_groups, **pairwise_ROI_groups)#combine
		groups_to_type = {v:k for k, v in full_ROI_groups.items()}
		for i in range(0, results.size(), 1):
			results.setValue("Type", i, groups_to_type[int(results.getValue("Group", i))])
		#save measurements
		results.save(datadir + "/analysis/" + c + "_" + cell_id + "_results.csv")
		
		
		#POOLED ROI ANALYSIS
		all_ROI_groups = ROI_groups#includes only selected organelles
		for combo in combo_present:#from section above - checks whether both organelles in pair are selected
			all_ROI_groups[combo] = pairwise_ROI_groups[combo]	
		for r in all_ROI_groups.keys():#includes pairwise overlap groups
			Roi.setDefaultGroup(all_ROI_groups[r])
			rm.selectGroup(all_ROI_groups[r])
			r_nselected = len(rm.getSelectedIndexes())
			if (r_nselected != rm.getCount()) and (r_nselected > 1):
				rm.runCommand("Combine")
				rm.addRoi(cell_crop.getRoi())#IMPORTANT - depends on cell_crop.show() above
				rm.selectGroup(all_ROI_groups[r])
				#delete non-combined ROIs
				r_all = rm.getSelectedIndexes()
				r_last = r_all.pop(-1)
				r_individual = r_all[:-1]
				rm.rename(r_last, cell_id + "_" + r)
				rm.setSelectedIndexes(r_individual)#all but last
				rm.runCommand("Delete")	
				rm.selectGroup(all_ROI_groups[r])
			elif (r_nselected != rm.getCount()) and (r_nselected == 1):#rename singlets from cell_r_1get
				rm.rename(rm.getSelectedIndex(), cell_id + "_" + r)
		
				#save ROIs
		rm.runCommand("Select All")
		rm.save(datadir + "/analysis/" + c + "_" + cell_id + "_POOLED_ROIs.zip")
		
		#MEASUREMENTS
		results.reset()#might need to close for next line to work properly
		IJ.run("Set Measurements...", "area centroid center shape feret's skewness kurtosis display redirect=None decimal=3")
		rm.runCommand("Select All")#should be redundant
		rm.runCommand("Measure")
		#add groups for processing in R - copied from above
		full_ROI_groups = dict(ROI_groups, **pairwise_ROI_groups)#combine
		groups_to_type = {v:k for k, v in full_ROI_groups.items()}
		for i in range(0, results.size(), 1):
			results.setValue("Type", i, groups_to_type[int(results.getValue("Group", i))])
		
		#save measurements
		results.save(datadir + "/analysis/" + c + "_" + cell_id + "_POOLED_results.csv")
		
						
		#max intensity stacks for total area covered
		results.reset()
		_, _, _, nSlices, _ = cell_crop.getDimensions()
		slicenames = dict()
		for i in range(1, nSlices + 1, 1):
			slicenames[cell_crop.getImageStack().getSliceLabel(i)] = i
		org_slicenames = dict((k, slicenames[k]) for k in organelles_selected)
		combo_slicenames = dict((k, slicenames.get(k)) for k in combo_present)
		cell_crop_mainorg = SM.makeSubstack(cell_crop, ','.join([str(i) for i in org_slicenames.values()]))
		cell_crop_pairwise = SM.makeSubstack(cell_crop, ','.join([str(i) for i in combo_slicenames.values()]))
		mainorg_max = ZProjector.run(cell_crop_mainorg, "max")
		mainorg_max.setTitle("organelles")
		pairwise_max = ZProjector.run(cell_crop_pairwise, "max")
		pairwise_max.setTitle("pairwise_contact_sites")
		PA.setSummaryTable(results)
		IJ.run(mainorg_max, "Analyze Particles...", "size=0-Infinity summarize composite")
		PA.setSummaryTable(results)
		IJ.run(pairwise_max, "Analyze Particles...", "size=0-Infinity summarize composite")
		results.save(datadir + "/analysis/" + c + "_" + cell_id + "_summary_results.csv")
		
				
		#clear ROIs, results, leave stack and cell image open (close duplicates)
		cells_copy.close()
		cell_crop.close()
		
		cells.show()
		orgstack.show()
		
		rm.reset()
		results.reset()
		
		
	
	cells.close()
	orgstack.close()
	IJ.run("Close All", "")
	rm.reset()
	results.reset()
	
rm.close()
IJ.selectWindow("Results")
IJ.run("Close")
