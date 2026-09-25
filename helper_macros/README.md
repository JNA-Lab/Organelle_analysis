# ImageJ macros for image processing and segmentation

Helper macros to make image segmentation tasks easier.  None of these are required to use the organelle analysis pipeline, but they might make your life easier.

### `manual_checker.py`
For checking the quality of segmentations.  Sequentially opens images in a folder, sets the LUT to 'glasbey on dark', and waits for the user to press okay before moving to the next image.  I use this to check all my segmented data before running the analysis, and note down which images should be excluded for segmentation errors.

### `split_to_TIF.py`
For splitting multi-channel images into individual TIFs with descriptive suffixes.  I use this to separate out organelle channels for segmentation.

## Manual and semi-automated segmentation

Cell segmentation has consistently been the most challenging step of our analysis workflow.  These macros are designed to assist with two different approaches — fully manual segmentation in ImageJ (for both cell and nuclei masks), or a semi-automated approach, in which cells are initially segmented with an Ilastik pixel classifier, abutting cells are manually separated in ImageJ, and the resulting segmentation files cleaned up in Gimp or a similar image editing programme.

### Manual segmentation - `manual_cell_and_nuclei_seg.ijm`
Designed to step through a directory of multi-channel fluorescence images and guide users to trace out cell and nucleus ROIs.  Maximum intensity projections are generated for cell tracing, and the multi-channel image (opened to a user-configured channel by default) is provided for nucleus tracing.  Segmented TIF files are automatically saved before moving to the next step or image.

### Semi-automated segmentation
1. Initial segmentation (pixel classification) with Ilastik, trained and processed on maximum-intensity projections of all the collected organelle channels.
2. Manual separation in ImageJ using `manual_segmentation_cutter.ijm`, which will overlay the Ilastik segmentation with the fluorescence data so that abutting cells can be manually separated by drawing a line of black pixels between them.  Make sure to set the pencil tool width to 2px to ensure that cells are fully separated.  TIF files are automatically saved as binary masks, without the overlaid fluorescence data.
3. Post-processing in Gimp.  
	1. Flood fill each complete cell in the image with a different colour or shade of gray (it doesn't matter what colours these are, as long as each cell gets a different pixel value.  Check that the cells do not bleed into each other, and that parts of the cell haven't been cut off.
	2. Select any remaining white pixels (select by colour) and delete them.  This will remove loose pixels at the edge of the cell, partial cells around the edges of the image, etc.
	3. Carefully fill in any segmentation holes (black pixels) in the interior of each cell with that cell's colour.  I do this by selecting a very large brush (at hardness 100%) and drawing around the middle of the cell, and then switching to progressively smaller crushes as I get closer to the periphery.  The analysis pipeline will detect any black pixels within the cell as part of the cell periphery, which will interfere with MDR calculation.