from ij import IJ #ImagePlus
from ij.gui import NonBlockingGenericDialog #GenericDialog, Roi
from ij.io import DirectoryChooser
from ij.measure import ResultsTable
#from ij.plugin import Duplicator, ImageCalculator, ImagesToStack, SubstackMaker, ZProjector, HyperStackConverter
from ij.plugin.frame import RoiManager
#from ij.plugin.filter import ParticleAnalyzer
import os
import re

#get folder and file information
img_dc = DirectoryChooser("Select the folder with your fluorescence images")
imgdir = img_dc.getDirectory()

ROI_dc = DirectoryChooser("Select the folder with your analysis files")
ROIdir = ROI_dc.getDirectory()

#identify image and ROI files
print("Identifying files...")
all_ROI_files = [ f for f in os.listdir(ROIdir) if os.path.isfile(ROIdir + "/" + f) and re.match('.+_ALL_ORGANELLE_ROIs\.zip$', f) and not re.match('Thumbs.db', f) ]#includes cell ROIs
img_files = [ f for f in os.listdir(imgdir) if os.path.isfile(imgdir + "/" + f) and not re.match('^\.', f) and not re.match('Thumbs.db', f) ]

#create analysis folder for output
if not os.path.exists(ROIdir + "/fluorescence_quantification/"):
	os.mkdir(ROIdir + "/fluorescence_quantification/")


#get number of channels from first image - might break if different images have different numbers of channels
print("Processing channels...")
img1 = IJ.openImage(os.path.join(imgdir, img_files[0]))
channels = img1.getDimensions()[2]#width, height, nChannels, nSlices, nFrames
#calling nChannels field directly doesn't work - BioFormats issue; see https://forum.image.sc/t/loci-plugin-in-colorizer-when-importing-through-bioformats/27183/2

#get processing options
options = NonBlockingGenericDialog('Options')
options.addMessage("Select the channels to quantify.")
for n in range(1, channels + 1, 1):
	options.addCheckbox("Channel " + str(n), False)
	options.addToSameRow()
	options.addStringField('Name', '')
options.showDialog()

to_quantify = dict()
if options.wasOKed():
	for n in range(1, channels + 1, 1):
		include = options.getNextBoolean()
		name = options.getNextString()
		if include == True:
			to_quantify[n] = name

#matching image files to ROI files
ROI_dict = dict()
for f in img_files:
	ROI_files = [r for r in all_ROI_files if f in r]
	if len(ROI_files) > 0:
		ROI_dict[f] = ROI_files[0]#should only be one


#set up ROI manager, don't display
RM = RoiManager(False)
rm = RM.getRoiManager()
IJ.run("Set Measurements...", "mean min integrated display redirect=None decimal=3")#TODO - restore afterwards?


for f in ROI_dict.keys():#will skip any image files without matching ROIs - e.g. excluded from analysis
	img = IJ.openImage(imgdir + '/' + f)#os.path.join not working - string vs list input?
	rm.open(ROIdir + '/' + ROI_dict[f])
	for c,n in to_quantify.items():
		img.setC(c)
		rm.runCommand("Select All")
		rm.runCommand(img, "Measure")#TODO - calculate without bringing up results table
		results = ResultsTable.getResultsTable()
		results.save(ROIdir + "/fluorescence_quantification/" + f + "_" + n + "_quantification.csv")
	img.close()
	rm.reset()
	results.reset()

rm.close()
IJ.selectWindow("Results")
#can't select when not displayed
IJ.run("Close")



