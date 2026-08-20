# Organelle analysis

Software tools for extracting quantitative organelle metrics from segmented imaging data.  A companion to PRISM — a Plasmid-based Reporter for Intracellular Spectral Microscopy — but can be used for any organelle images.  These are designed to be easy to use with minimal technical expertise, but can also be readily customised to suit your data.



## Getting started

The organelle analysis pipeline runs on segmented, 2D image data of organelles within cells.  At a minimum, it requires only cell masks and at least one organelle, although we highly recommend including nucleus segmentations, as some metrics will be less accurate without this information.  As a default, this pipeline will process Golgi, ER, mitochondria, peroxisome, and lysosomes (all included in the PRISM plasmid), as well as nucleus and an optional additional organelle channel.  That said, all organelles (with some exceptions for the nucleus) are processed in the same way, so you can swap any of these out for your favourite organelle!

### Installation
No installation necessary!  Once you've downloaded the files, they can be run by dragging and dropping them directly into ImageJ (for .py and .ijm files) or RStudio (for .Rmd).

<!-- If you don't have the required R packages installed, RStudio should prompt you to install them.  If it doesn't, you can run this script:

``` ```-->




### Prerequisites
#### Segmentation
The organelle analysis pipeline runs on _segmented_ images, not directly on your fluorescence data.  You can use whatever tools you like for this, as long as the final output is segmented TIF files — cell masks must be _instance_ segmentation, i.e. a different pixel value for each cell in the image, and organelle channels can be binary masks or instance segmentaion (instance segmentations will be converted to binary masks by the processing macro).

We find that Ilastik works well for organelle segmentation, and can be trained to segment organelles even where the morphology changes significantly between experimental conditions.  The `split_to_TIF.py` macro can split your multi-channel images into single-channel TIFs suitable for Ilastik training and segmentation.

Cell segmentation can a bit more complicated, and we'd recommend trying a variety of tools to see what works best on your data.  We've included some helper scripts for manual or semi-manual segmentation.

Finally, all the segmentations should have a background value of 0 — i.e. white organelles on a black background.  If your segmentation tool outputs black organelles on a white background, tick the 'invert' checkbox when setting up the ImageJ macro.

#### File naming
File names are used to match up segmented images, so all segmented images should share a base name with a distinct suffix for each organelle and cell mask.  Furthermore, file names can be split at each underscore and different parts mapped to specific variables, so any metadata included in file names (e.g. `replicate1_control_slide1_cell3`) can be used directly in later analysis.

### Process
Once your data is segmented, the basic process is:
1) Run the `Organelle_analysis_segmentation_processing.py` macro in ImageJ.
2) (Optionally) run the `Organelle_analysis_fluorescence_quantification.py` macro in ImageJ, if you want to quantify a fluorescent marker in the context of your organelle masks.
3) Run the `Organelle_analysis.Rmd` notebook in RStudio.

User-customisable parameters (where present) are in a clearly labelled box at the top of each file.  Most should be set to sensible defaults; the one you are mostly likely to need to change is the one in `Organelle_analysis.Rmd` that controls how file names are mapped to variables in the final data table.  There is guidance within the file on setting this up.

### Outputs
The main output of this pipeline is a table called `full_metrics.filtered` with one row for each cell in your data, containing a variety of metrics related to the abundance, shape, size, and intracellular distribution of each type of organelle and organelle contact in that cell.  Subject to the options selected, this data is appropriately normalised, filtered, and ready for your downstream analyses.

The outputs of the intermediate ImageJ macros are also saved in standard formats, in case you want to use them for other applications.  These include ROI files for your segmented cells and organelles (on a per-image and per-cell basis), a TIF stack for each cell containing all included organelle masks and calculated organelle contacts, and various CSV files with information on each segmented organelle within your cells.




