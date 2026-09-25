setBatchMode(false);//interactive

function main(dir, basename){
	
	selectImage("base");
	autoBC();
	
	selectImage("seg");
	run("Add Image...", "image=[base] x=0 y=0 opacity=50 zero");
	
	Color.set("#000000")
	//can't manually select drawing tools - ensure brush or pencil is selected with width = 2
	
	Dialog.createNonBlocking("Proceed");
	Dialog.addMessage("Draw a line to separate the cells, then click 'okay'.");
	Dialog.show();
	
	selectImage("seg");
	run("Remove Overlay");
	
	saveAs("Tiff", dir + basename + "_cut.tif");
	close;
	
	selectImage("base");
	close;
	
	
}

//---batch runner---
Dialog.create("Select options");
Dialog.addCheckbox("Run as batch?", true);
Dialog.addString("Segmentation suffix", "_cellbinary");
Dialog.addChoice("File extension", newArray("tif"), "tif");
Dialog.show();
batch = Dialog.getCheckbox();
suffix = Dialog.getString();
ext = Dialog.getChoice();
ext = "." + ext

if(batch){
	dir = getDirectory("Choose a Directory");
	File.makeDirectory(dir + "separated");
		
	files = getFileList(dir);

	temp = newArray(0);
	for(i = 0; i < files.length; i++){
		if(!endsWith(files[i], "/")&&!endsWith(toLowerCase(files[i]), ".ini")&&!endsWith(files[i], suffix + ext)){
			temp = Array.concat(temp,files[i]);
		}
	}
	Array.print(temp);
	files = temp;
	
	
		
	for(i = 0; i < files.length; i++){
		print("BASE - " + files[i]);
		open(dir + files[i]);
		rename("orig");
		run("Duplicate...", "title=[base] duplicate");
		close("orig");
		print("SEG - " + files[i].replace(ext, "") + suffix + ext);
		open(dir + files[i].replace(ext, "") + suffix + ext);
		rename("orig");
		run("Duplicate...", "title=[seg] duplicate");
		close("orig");
		main(dir + "/separated/", replace(files[i], "\\" + ext, ""));
		//saving and closing files in main loop
	}
}
else{
	dir = getInfo("image.directory");//directory of image - for saving results
	File.makeDirectory(dir + "separated");
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
