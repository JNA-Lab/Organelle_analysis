setBatchMode(true);
//-----Select directory and list files-----
dir = getDirectory("Choose the directory containing segmentation masks.");
img_list =  getFileList(dir);
img_list_ext = newArray();
for (i = 0; i < img_list.length; i++) {
	if (img_list[i].indexOf("/") == -1){
		img_list_ext = Array.concat(img_list_ext, substring(img_list[i], img_list[i].lastIndexOf(".")));
	}
}
img_list_ext = ArrayUnique(img_list_ext);
//---------------

//-----Options Dialog-----
Dialog.create("Options");
//Dialog.addChoice("File extension", img_list_ext);
Dialog.addMessage("Select which organelles to analyse,and how they \nare represented in file names.");
Dialog.addCheckbox("Cell boundaries (required)", true);
Dialog.addToSameRow();
Dialog.addString("", "cell");
Dialog.addCheckbox("Nucleus", true);
Dialog.addToSameRow();
Dialog.addString("", "nucleus");
Dialog.addCheckbox("Golgi", true);
Dialog.addToSameRow();
Dialog.addString("", "golgi");
Dialog.addCheckbox("mitochondria", true);
Dialog.addToSameRow();
Dialog.addString("", "mito");
Dialog.addCheckbox("ER", true);
Dialog.addToSameRow();
Dialog.addString("", "ER");
Dialog.addCheckbox("peroxisomes", true);
Dialog.addToSameRow();
Dialog.addString("", "perox");
Dialog.addCheckbox("lysosomes", true);
Dialog.addToSameRow();
Dialog.addString("", "lyso");
Dialog.addCheckbox("bacteria", true);
Dialog.addToSameRow();
Dialog.addString("", "bac");
Dialog.addMessage("\n\n");
Dialog.addCheckbox("Calculate pairwise overlaps?", true);
Dialog.show();

//ext = Dialog.getChoice();
cell_inc = Dialog.getCheckbox();//verify true
cell_match = Dialog.getString();
nuc_inc = Dialog.getCheckbox();
nuc_match = Dialog.getString();
golgi_inc = Dialog.getCheckbox();
golgi_match = Dialog.getString();
mito_inc = Dialog.getCheckbox();
mito_match = Dialog.getString();
ER_inc = Dialog.getCheckbox();
ER_match = Dialog.getString();
perox_inc = Dialog.getCheckbox();
perox_match = Dialog.getString();
lyso_inc = Dialog.getCheckbox();
lyso_match = Dialog.getString();
bac_inc = Dialog.getCheckbox();
bac_match = Dialog.getString();
contacts = Dialog.getCheckbox();
//img_list = Array.filter(img_list, ext);
//---------------





//---------------------------------------------------------------------
//-----FUNCTIONS-------------------------------------------------------
//---------------------------------------------------------------------
//From ImageJ example macro Array_Functions.txt
function ArrayUnique(array) {
	array 	= Array.sort(array);
	array 	= Array.concat(array, 999999);
	uniqueA = newArray();
	i = 0;	
   	while (i<(array.length)-1) {
		if (array[i] == array[(i)+1]) {
			//print("found: "+array[i]);			
		} else {
			uniqueA = Array.concat(uniqueA, array[i]);
		}
   		i++;
   	}
	return uniqueA;
}
