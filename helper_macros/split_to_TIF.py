from ij import ImagePlus, IJ
from ij.gui import GenericDialog
from ij.io import DirectoryChooser
from ij.plugin import ChannelSplitter
import os
import re

default_labels = ['perox', 'ER', 'mito', 'lyso', 'Golgi', 'residuals']

dc = DirectoryChooser("Select the folder your images")
datadir = dc.getDirectory()
filenames = []

for obj in os.listdir(datadir):
	if os.path.isfile(datadir + "/" + obj):
		filenames.append(obj)

#create analysis folder for output
if not os.path.exists(datadir + "/split_TIFs/"):
	os.mkdir(datadir + "/split_TIFs/")

options = GenericDialog('Options')
options.addMessage("Select which channels to include,and how they \nshould be labelled.");
for i in range(0, 10, 1):
	options.addCheckbox(str(i+1), i+1 <= 6)
	options.addToSameRow()
	if i+1 <= len(default_labels):
		options.addStringField('', default_labels[i])
	else:
		options.addStringField('', '')
options.showDialog()

ch_bool = dict()
ch_label = dict()
for i in range(0, 10, 1):#match above
	ch_bool[i] = options.getNextBoolean()
	ch_label[i] = options.getNextString()
	if ch_bool[i] is False:
		del ch_label[i]

CS = ChannelSplitter()

for f in filenames:
	imp = IJ.openImage(os.path.join(datadir, f))
	chs = CS.split(imp)
	#cnames = []
	#for img in chs:
	#	cnames.append(img.getName())
	#for c in ch_labels.keys():
	#	r = re.compile("*" + c.zfill(4) + "$")
	#	cimg_name = list(filter(r.match, cnames))[0]#TODO - check if >1 match
	for i in range(0, len(chs), 1):
		IJ.saveAs(chs[i], "Tiff", datadir + "/split_TIFs/" + f + "_" + ch_label[i])
		chs[i].close()

print("Finished running")