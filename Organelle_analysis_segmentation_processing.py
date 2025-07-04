from ij import ImagePlus, IJ
from ij.gui import GenericDialog, NonBlockingGenericDialog, Roi
from ij.io import DirectoryChooser
from ij.measure import ResultsTable
from ij.plugin import ImageCalculator, StackEditor, ImagesToStack, SubstackMaker, ZProjector, HyperStackConverter
from ij.plugin.frame import RoiManager
from ij.plugin.filter import ParticleAnalyzer
from ij.process import ImageProcessor
from java.awt import Color
import os
import re

organelles = ['nuclei', 'Golgi', 'peroxisomes', 'ER', 'mitochondria', 'lysosomes', 'bacteria']

#***USER DEFAULTS***
#change these values to set your default cell/organelle suffixes
default_cells = "_cp_masks"
default_org = dict()
default_org["nuclei"] = "_composite_seg_nucleus_mask"
default_org["Golgi"] = "-Best_1dpi__2_golgi"
default_org["peroxisomes"] = "-Best_so_far_perox"
default_org["ER"] = "-Best_1dpi_ER"
default_org["mitochondria"] = "-Best_1dpi_2_mito"
default_org["lysosomes"] = "_lysosomes"
default_org["bacteria"] = "-Best_1dpi_Ot_LD"
default_checkboxes = [False, True, True, True, True, False, True]#in same organelle order as above
#*********************


#get folder and file information
dc = DirectoryChooser("Select the folder with your segmented images")
datadir = dc.getDirectory()
f_ext = set()
filenames = []

for obj in os.listdir(datadir):
	if os.path.isfile(datadir + "/" + obj):
		filenames.append(obj)
		f_ext.add(re.search(r"\.([A-Za-z0-9]+)$", obj).group(1))

#create analysis folder for output
if not os.path.exists(datadir + "/analysis/"):
	os.mkdir(datadir + "/analysis/")


#get processing options
options = GenericDialog('Options')
options.addChoice('File extension', list(f_ext), "tif")
options.addMessage("Select which organelles to analyse,and how they \nare represented in file names.");
options.addCheckbox('Cell boundaries (required)', True)
options.addToSameRow()
options.addStringField('', default_cells)
for i in range(0, len(organelles), 1):
	o = organelles[i]
	options.addCheckbox(o, default_checkboxes[i])
	options.addToSameRow()
	options.addStringField('', default_org[o])
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
filenames_key = dict()
for f in filenames_filtered:
	g = f.replace("." + ext, "")
	for o in org_regex.values():
		g = g.replace(o, "")#not actually regex matching - probably best for now
	g = g.replace(cell_regex, "")
	conditions.add(g)
	if g in filenames_key.keys():
		filenames_key[g].append(f)
	else:
		filenames_key[g] = [f]


#set up ROI manager
RM = RoiManager(False)#don't display
rm = RM.getRoiManager()
rm.hide()
SE = StackEditor()
I2S = ImagesToStack()
SM = SubstackMaker()
PA = ParticleAnalyzer()

##MAIN LOOP
progress = 0
for c in conditions:
	print("Starting condition " + c)
	progress += 1
	#open and rename images
	c_filenames = filenames_key[c]#find filenames matching condition
	c_cell_image = str()
	imp_key = {}
	for cf in c_filenames:
		imp = IJ.openImage(os.path.join(datadir, cf))#IJ.getImage()
		try:
			organelle = [k for k,v in org_regex.items() if v in cf][0]
			imp.setTitle(organelle)
			imp_key[organelle] = imp
		except:#index out of range - no match in organelles
			if cell_regex in cf:
				c_cell_image = cf
			imp.close()#close for now		
	
	images_to_stack = []
	#ORGANELLE THRESHOLDING
	for org in organelles_selected:
		org_img = imp_key[org]
		IJ.setRawThreshold(org_img, 1, (2**org_img.getBitDepth())-1)#will work for 8 or 16-bit; might break on 32-bit or RGB
		IJ.run(org_img, "Convert to Mask", "")
		images_to_stack.append(org_img)
		
	#PAIRWISE OVERLAPS
	for combo in pairwise_groups.keys():
		if set(pairwise_groups[combo]).issubset(organelles_selected):
			img1 = imp_key[pairwise_groups[combo][0]]
			img2 = imp_key[pairwise_groups[combo][1]]
			img3 = ImageCalculator.run(img1, img2, "AND create")
			img3.setTitle(combo)
			images_to_stack.append(img3)
	
	#STACK
	orgstack = I2S.run(images_to_stack)
	orgstack.setTitle("Organelles")
		
	#LOAD CELL IMAGE
	cells = IJ.openImage(os.path.join(datadir, c_cell_image))
	cells.setTitle("Cells")
	
	#Check max value
	cells_max = int(cells.getStatistics().max) #converting double to int - beware of rounding bugs
	
	#create dict for saving ROIs (aligned to original image)
	cells_ROIs_orig = dict()
	pooled_ROIs_per_cell = dict()
	
	#Per step (cell) loop:
	for i in range(1, cells_max + 1, 1):
		print("Processing cell " + str(i) + "...")
		cell_id = "cell_" + str(i)
		#duplicate cell image
		cells_copy = cells.duplicate()
		cells_copy.setTitle("cells_DUPLICATE")
		#step threshold
		IJ.setRawThreshold(cells_copy, i, i)
		#analyse particles
		Roi.setDefaultGroup(ROI_groups["cells"])#TODO - dynamic?
		IJ.run(cells_copy, "Analyze Particles...", "size=0-Infinity add composite") #should just be one cell
		#SOME CELLS RETURNING >1 - combine
		if rm.getCount() > 1:
			rm.selectGroup(ROI_groups["cells"])
			rm.runCommand(cells_copy, "Combine")
			rm.addRoi(cells_copy.getRoi())
			#delete non-combined ROIs
			rm.setSelectedIndexes(range(0, rm.getCount()-1, 1))#all but last
			rm.runCommand("Delete")
		rm.rename(0, cell_id)#renaming for image-wide list
		cells_ROIs_orig[cell_id] = rm.getRoi(0)
		#duplicate stack and clear outside ROI
		rm.select(0)#probably redundant
		cell_crop = orgstack.crop([rm.getRoi(0)], "stack")[0]
		cell_crop.setTitle(cell_id)
		rm.addRoi(cell_crop.getRoi())#cell boundary in cropped image
		rm.select(0)
		rm.runCommand("Delete")#deleting original cell boundary - now in wrong place in cropped image
		rm.rename(0, cell_id)
		rm.select(cell_crop, 0)		
		IJ.run(cell_crop, "Clear Outside", "stack")
		cell_crop_mask = cell_crop.createRoiMask()
		cell_crop.getImageStack().addSlice('cell', cell_crop_mask)
		
		#analyse particles per organelle/pair
		#manual implementation of SE.convertStackToImages for background running:
		cell_crop_slices = cell_crop.getStackSize()
		cell_crop_stack = cell_crop.getImageStack()
		crop_key = {}
		for i in range(1, cell_crop_slices + 1, 1):
			label = cell_crop_stack.getShortSliceLabel(i)
			ip = cell_crop_stack.getProcessor(i)
			i = ImagePlus(label, ip)
			crop_key[label] = i
		#select and close cell slice image - ROI already added
		cell_img = crop_key['cell']
		cell_img.close()
		#iterate over organelles and contact types
		for org in organelles_selected:
			Roi.setDefaultGroup(ROI_groups[org])
			org_img = crop_key[org]
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
				combo_img = crop_key[combo]
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
		cell_crop_copy = cell_crop.duplicate()
		HyperStackConverter.toHyperStack(cell_crop_copy, len(images_to_stack) + 1, 1, 1)#+1 for added cell mask
		IJ.saveAs(cell_crop_copy, "Tiff", datadir + "/analysis/" + c + "_" + cell_id + "_stack.tif")
		
		#save ROIs
		rm.runCommand("Select All")
		rm.save(datadir + "/analysis/" + c + "_" + cell_id + "_ROIs.zip")
		
		#measure
		IJ.run("Set Measurements...", "area mean min centroid center shape feret's skewness kurtosis display redirect=None decimal=3")
		rm.deselect()
		rm.runCommand("Select All")#should be redundant
		rm.runCommand(cell_img, "Measure")#TODO - limit slices?
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
				rm.runCommand(cell_crop, "Combine")
				rm.addRoi(cell_crop.getRoi())
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
		rm.deselect()
		rm.save(datadir + "/analysis/" + c + "_" + cell_id + "_POOLED_ROIs.zip")
		for i in range(0, rm.getCount(), 1):
			r = rm.getRoi(i)
			cb = cells_ROIs_orig[cell_id].getBounds()
			r.translate(cb.x, cb.y)
		pooled_ROIs_per_cell[cell_id] = rm.getRoisAsArray()
			
		
		#MEASUREMENTS
		results.reset()#might need to close for next line to work properly
		IJ.run("Set Measurements...", "area centroid center shape feret's skewness kurtosis display redirect=None decimal=3")
		rm.runCommand("Select All")#should be redundant
		rm.deselect()
		rm.runCommand(cell_img, "Measure")
		#add groups for processing in R - copied from above
		full_ROI_groups = dict(ROI_groups, **pairwise_ROI_groups)#combine
		groups_to_type = {v:k for k, v in full_ROI_groups.items()}
		for i in range(0, results.size(), 1):
			results.setValue("Type", i, groups_to_type[int(results.getValue("Group", i))])
		
		#save measurements
		results.save(datadir + "/analysis/" + c + "_" + cell_id + "_POOLED_results.csv")
		
						
		#max intensity stacks for total area covered
		results.reset()#?
		_, _, _, nSlices, _ = cell_crop.getDimensions()
		slicenames = dict()
		for i in range(1, nSlices + 1, 1):
			slicenames[cell_crop.getImageStack().getShortSliceLabel(i)] = i
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

		
		rm.reset()
		results.reset()
		
		
		#save stack labels for interpreting tiff stacks in R
		with open(datadir + "analysis\\" + c + "_" + cell_id + "_slice_labels.csv", "w") as f:#different path format than ImageJ
			for i in range(1, len(slicenames) + 1, 1):
				f.write(str(i) + "," + slicenames.keys()[list(slicenames.values()).index(i)] + "\n")
		#TODO - move combo_present and this outside of loop - should be the same for all cells and conditions in a batch
		
	IJ.showProgress(progress/len(conditions))
	print("Progress - " + str(progress) + " of " + str(len(conditions)) + " conditions done")		
	
	cells.close()
	orgstack.close()
	rm.reset()
	
	for i in range(0, len(cells_ROIs_orig), 1):
		k = cells_ROIs_orig.keys()[i]
		v = cells_ROIs_orig[k]
		rm.addRoi(v)
		rm.rename(i, k)
	rm.deselect()
	rm.save(datadir + "/analysis/" + c + "_ALL_CELL_ROIs.zip")
	rm.reset()
	for i in range(0, len(pooled_ROIs_per_cell), 1):
		k = pooled_ROIs_per_cell.keys()[i]
		vl = pooled_ROIs_per_cell[k]
		for v in vl:
			rm.addRoi(v)
	rm.deselect()
	rm.save(datadir + "/analysis/" + c + "_ALL_ORGANELLE_ROIs.zip")
	rm.reset()
	
	
	results.reset()
	
rm.close()
IJ.selectWindow("Results")#can't select when not displayed
IJ.run("Close")

