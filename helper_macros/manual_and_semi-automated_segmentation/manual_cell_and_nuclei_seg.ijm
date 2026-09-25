setBatchMode(false);//interactive

function main(dir, nch){

	initial_name = getInfo("image.title");
	rename("orig");
	run("Duplicate...", "title=["+initial_name+"_duplicate] duplicate");
	close("orig");
	
	getDimensions(width, height, channels, slices, frames);
	
	run("Z Project...", "projection=[Max Intensity]");
	rename("max");
	autoBC();
	
	run("ROI Manager...");
	roiManager("Show None");
	
	setTool("polygon");
	run("Brightness/Contrast...");//open tool for adjustments
	
	Dialog.createNonBlocking("Cells");
	Dialog.addMessage("Trace cell boundaries, add to ROI manager, then click 'okay'.");
	Dialog.show();
	
	//make new empty image, iterate over ROIs, and set pixel values
	selectImage("max");
	//duplicating and clearing instead of creating new to preserve scale
	roiManager("deselect");
	run("Duplicate...", "title=[cells] ignore");
	run("Select All");
	setBackgroundColor(0, 0, 0);
	run("Clear", "slice");
	run("8-bit");
	//select ROIs and set values
	roi_count = roiManager("count");
	if(roi_count >= 1){
		for (i = 0; i < roi_count; i++) {
			setForegroundColor(i+1, i+1, i+1);//WRONG - SET EXACT VALUE
			roiManager("select", i);
			roiManager("Fill");
		}	
	}
	saveAs("Tiff", dir + "cell_and_nuclei_masks/" + initial_name + "_cells.tif");
	close();
	close("cells");
	close("max");
	roiManager("deselect");
	roiManager("delete");
	
	Stack.setChannel(nch);
	
	Dialog.createNonBlocking("Nuclei");
	Dialog.addMessage("Trace nuclei, add to ROI manager, then click 'okay'.");
	Dialog.show();
	
	//deselect (select all), combine, and create mask
	roiManager("deselect");
	roiManager("Combine");
	run("Create Mask", "nuclei");
	saveAs("Tiff", dir + "cell_and_nuclei_masks/" + initial_name + "_nuclei.tif");
	
	close("*");
	close("Roi Manager");
}

//---batch runner---
Dialog.create("Batch selection");
Dialog.addCheckbox("Run as batch?", true);
Dialog.addNumber("Nuclei channel", 1);
Dialog.show();
batch = Dialog.getCheckbox();
nch = Dialog.getNumber();

if(batch){
	//TODO - directory selection and iteration
	//main();//placeholder
	dir = getDirectory("Choose a Directory");
	File.makeDirectory(dir + "cell_and_nuclei_masks");
	
	files = getFileList(dir);
	//exclude directories
	temp = newArray(0);
	for(i = 0; i < files.length; i++){
		if(!endsWith(files[i], "/")&&!endsWith(toLowerCase(files[i]), ".ini")){
			temp = Array.concat(temp,files[i]);
		}
	}
	files = temp;
	print(files.length);
	print(files[0]);
		
	for(i = 0; i < files.length; i++){
		open(dir + files[i]);
		print("Processing " + getTitle());
		main(dir, nch);
		//saving and closing files in main loop
	}
}
else{
	dir = getInfo("image.directory");//directory of image - for saving results
	File.makeDirectory(dir + "cell_and_nuclei_masks");
	main(dir, nch);
}




/*-------------------------------------------------*/

function autoBC(){
//from Kota Miura and Damien Guimond 
//http://imagej.1557.x6.nabble.com/Auto-Brightness-Contrast-
//and-setMinAndMax-td4968628.html
	AUTO_THRESHOLD = 5000; 
	getRawStatistics(pixcount); 
	limit = pixcount/10; 
	threshold = pixcount/AUTO_THRESHOLD; 
	nBins = 256; 
	getHistogram(values, histA, nBins); 
	i = -1; 
	found = false; 
	do { 
	        counts = histA[++i]; 
	        if (counts > limit) counts = 0; 
	        found = counts > threshold; 
	} while ((!found) && (i < histA.length-1)) 
	hmin = values[i]; 
	i = histA.length; 
	do { 
	        counts = histA[--i]; 
	        if (counts > limit) counts = 0; 
	        found = counts > threshold; 
	} while ((!found) && (i > 0)) 
	hmax = values[i]; 
	setMinAndMax(hmin, hmax); 
	run("Apply LUT"); 
}
