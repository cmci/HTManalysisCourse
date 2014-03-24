# Analysis of High-Throughput Microscopy Image Data

BioIT-EMBL Centres Joint Training Week, Mar. 26, 2014

Kota Miura (CMCI, EMBL)

- [Course Web Page](https://bio-it.embl.de/joint-training) (only internal access)

##Contents

- [Instructors](#instructors)
- [Course Outline](#course-outline)
- [Introduction](#introduction)
- [Configuration for March 24](#configuration-for-march-24)
  - [Install Fiji on your desktop](#install-fiji-on-your-desktop)
  - [Configure access to the Fiji in the network](#configure-access-to-the-fiji-in-the-network)
  - [ Check if Fiji is working](#Check if Fiji is working)
- [Methods (Image Analysis)](#methods-image-analysis)
- [Resources](#resources)
  - [Fiji](#fiji)
  - [Data](#data)
  - [Scripts](#scripts)  
- [Workflow: GUI](#workflow-gui)
- [Workflow: Python Primer](#workflow-python-primer)
- [Workflow: Fiji Jython Scripting](#workflow:-fiji-jython-scripting)


## Instructors

- Kota Miura
  - Centre for Molecular and Cellular Imaging (CMCI), EMBL
  - <miura@embl.de>

- Christoph Schiklenk
  - Haering Lab, Cell Biology & Biophysics, EMBL
  - Christoph has been using Fiji/Jython for a year and he now is a quite literate Fiji scripting guy. He is volunteering to assist your learning. 

- Clemens Lakner
  - Bio-IT, EMBL
  - Clemens wrote the job-array generator scripts for image analysis. He will explain to you how his shell script is working.
  
- Christian Tischer
  - ALMF, EMBL
  - Tischi will be partially appearing to assist you. Imaging expert and not only an Jython-Fiji expert, but also a Cell Profiler expert.

- Aliaksandr Halavatyi
  - Pepperkok Lab, Cell Biology & Biophysics, EMBL
  - Aliaksandr will be partially appearing to assist you. He has been coding image analysis tools in Matlab, but knows ImageJ/Fiji pretty well and coding in Jython.   
  
- Acknowledgements
  - Thanks to Fatima Verissimo, Faba Neumann, Rainer Pepperkok and Jean-Karim Hériché for providing data, sharing their experiences and giving advices on these data. 
  - Thanks to Bernd Klaus, Matt Rogon and Aidan Budd (the course team) for ther effort in starting up this course and working together.   

##Course Outline

1. **Aim: Learn how to process and analyze high content image data by Fiji scripting**
2. Introduction: the background of data, setting up.
3. Processing and analyzing using Fiji GUI
   - We first try to segment nucleus image and do measurement using GUI, to know how the algorithm works.    
4. Python Crash Course
   - We want to write a script that does processing and measurements: for this we first need to learn  basics of python syntax.   
5. Fiji scripting in Jython
   - Jython is a Java-implemented Python, allows to use full functionality of Fiji by scripting. We learn how to write Jython script to do our task automatically.  
6. Submitting Jobs (executing Jython script) to the cluster.
   - As there are many images, we use cluster to parallelize the computation, to get results in a short time. 

## Introduction

Image based screening (high content screening) is one of the modern methods for screening genes involved in a biological process. The main idea is to sytematically acquire a huge amount of data with various treatments, and try to isolate treatments that affects the target biological process. 

Dataset we use in this course is from a scretory pathway screening. Proteins synthesized in endoplasmic reticulum are transported first to the Golgi, where proteins are added with modifications, and then to the plasma membrane or other desitned intracellular structures. 

VSVG is a model protein used for studying this transport system and is a heat-sensitive protein.  This protein could beaccumulated in / released from ER by temperature shifting, to experimentally analyze the process of intracellular transport from ER to the plasma membrane. 

To study the transport process, it is desirable to take time lapse sequence to follow the changes in VSVG localization within cell. However, it is sufficient to assess the transport process by taking a snapshot at a defined time point after releasing the protein from ER and check where the protein is located. In control condition, proteins should be localized in the plasma membrane after successful transports. If the transport out from the ER is interferred due to suppression of the component of transport, proteins should be stably located in the ER. Likewise, finding proteins in the Golgi suggests that transport machinary from Golgi to the plasma membrane is blocked. 

With this snapshot strategy, it is possible to systematically prepare a set of siRNAs and drug treatments to test their effect on vesicle transport system: the high content image screening. 

In this course, we will try to analyze all images from different types of treatments automatically, evaluate and compare the difference in the effeciency of transport. This COULD be done manually, but we have huge number of various treatments: approximately 700,000 of images to analyze. To not to end the rest of our life clicking images, we will completely automate the image processing and  analysis to extract transport parametes and output results as a huge list of parameter table. 

This table then is studied using statistical techniques to distinguish treatments that are largely affecting the transport system. 

Followings are the literatures published with this dataset. 

- Liebel, U., Starkuviene, V., Erfle, H., Simpson, J.C., Poustka, A., Wiemann, S., Pepperkok, R., 2003. A microscope-based screening platform for large-scale functional protein analysis in intact cells. FEBS Lett. 554, 394–8. <http://www.ncbi.nlm.nih.gov/pubmed/17275941>

- Simpson, J.C., Cetin, C., Erfle, H., Joggerst, B., Liebel, U., Ellenberg, J., Pepperkok, R., 2007. An RNAi screening platform to identify secretion machinery in mammalian cells. J. Biotechnol. 129, 352–65. <http://www.ncbi.nlm.nih.gov/pubmed/17275941>

- Simpson, J.C., Joggerst, B., Laketa, V., Verissimo, F., Cetin, C., Erfle, H., Bexiga, M.G., Singan, V.R., Hériché, J.-K., Neumann, B., Mateos, A., Blake, J., Bechtel, S., Benes, V., Wiemann, S., 
Ellenberg, J., Pepperkok, R., 2012. Genome-wide RNAi screening identifies human proteins with a regulatory function in the early secretory pathway. Nat. Cell Biol. 14, 764–74. <http://www.ncbi.nlm.nih.gov/pubmed/22660414>

## Configuration for March 24


### Install Fiji on your desktop

Access the following page

```
http://fiji.sc/Downloads
```

then install "Life-Line Version", Linux 64-bit. Place the Fiji folder on your desktop. 

###  Configure access to the Fiji in the network 

For using Fiji in the network from command line, do the following. 

#### Add an alias in your local bash profile.

Open ~/.bash_profile

```
open ~/.bash_profile
```

Add the following line in your bash profile (.bash_profile). 

```
alias cluster='ssh -t <YOUR UNIX USERNAME>@submaster '\''bsub -M 5000 -Is bash'\'''
```

By the way, this will allow you to execute the command below as an alias'cluster'

```
ssh -t <YOUR UNIX USERNAME>@submaster \bsub -M 5000 -Is bash
```
If you have time, try the above command directly and see that it works. 

Reload the .bash_profle by 

```
source .bash_profile
```

to validate your changes. 

#### Add several paths to your .bash_profile in the network.

If your are ready with the above new alias, try

```
cluster
```

You will then be asked for your password, and then will be connected to one of the cluster nodes in an interactive mode. 

Open your .bash_profile:

```
vim ~/.bash_profile
```

If you are not use to vi editor, you could use others such as


```
nano ~/.bash_profile

```

Otherwise, it might be able to edit the .bash_profile directly from your course linux machine. 


Add the follwoing lines in your .bash_profile

```
export PATH=/g/software/bin:$PATH
export PATH=/g/almf/software/bin2:$PATH
export PATH=/g/almf/software/bin/java/jdk1.6.0_45/bin/:$PATH
export JAVA_HOME=/g/almf/software/bin/java/jdk1.6.0_45
```

This is just to run Fiji from your command line (the one that I am maintaining) and to use the right Java version. You could delete these lines after the course. 

Reload the profile

```
source .bash_profile
```

### Check if Fiji is working 

If you are still connected to one of the nodes, do

```
java -version
```

The output should look like

```
[miura@compute083]~ $
java -version
java version "1.6.0_45"
Java(TM) SE Runtime Environment (build 1.6.0_45-b06)
Java HotSpot(TM) 64-Bit Server VM (build 20.45-b01, mixed mode)
```

If your java version is OK (same as above), then 

```
fiji --help
```
 The output should be
 
 ```
 [miura@compute083]~ $
fiji --help
Usage: fiji [<Java options>.. --] [<ImageJ options>..] [<files>..]

Java options are passed to the Java Runtime, ImageJ
options to ImageJ (or Jython, JRuby, ...).

In addition, the following options are supported by ImageJ:
General options:
--help, -h
	show this help
--dry-run
	show the command line, but do not run anything
--debug
	verbose output
--system
	do not try to run bundled Java
--java-home <path>
	specify JAVA_HOME explicitly
--print-java-home
	print ImageJ's idea of JAVA_HOME
--print-ij-dir
	print where ImageJ thinks it is located
--headless
	run in text mode
--ij-dir <path>
	set the ImageJ directory to <path> (used to find
	 jars/, plugins/ and macros/)
--heap, --mem, --memory <amount>
	set Java's heap size to <amount> (e.g. 512M)
--class-path, --classpath, -classpath, --cp, -cp <path>
	append <path> to the class path
--jar-path, --jarpath, -jarpath <path>
	append .jar files in <path> to the class path
--pass-classpath
	pass -classpath <classpath> to the main() method
--full-classpath
	call the main class with the full ImageJ class path
--dont-patch-ij1
	do not try to runtime-patch ImageJ1
--ext <path>
	set Java's extension directory to <path>
--default-gc
	do not use advanced garbage collector settings by default
	(-Xincgc -XX:PermSize=128m)
--gc-g1
	use the G1 garbage collector
--debug-gc
	show debug info about the garbage collector on stderr
--debugger=<port>
	start the JVM in a mode so Eclipse's debugger can attach to it
--no-splash
	suppress showing a splash screen upon startup

Options for ImageJ:
--ij2
	start ImageJ2 instead of ImageJ1
--ij1
	start ImageJ1
--allow-multiple
	do not reuse existing ImageJ instance
--plugins <dir>
	use <dir> to discover plugins
--run <plugin> [<arg>]
	run <plugin> in ImageJ, optionally with arguments
--compile-and-run <path-to-.java-file>
	compile and run <plugin> in ImageJ
--edit [<file>...]
	edit the given file in the script editor

Options to run programs other than ImageJ:
--update
	start the command-line version of the ImageJ updater
--javah
	start javah instead of ImageJ
--javap
	start javap instead of ImageJ
--javadoc
	start javadoc instead of ImageJ
--main-class <class name> (this is the
	default when called with a file ending in .class)
	start the given class instead of ImageJ
--retrotranslator
	use Retrotranslator to support Java < 1.6
 ```
 
## Methods (Image Analysis)

As we are interested in the attenuation or the enhancement of the protein transport machinary from ER to the plasma membrane via the Golgi apparatus, we compare the protein density in ER, the Golgi and the plasma membrane. For this, we use three different images. 

The minimal detaset consists of three images:

1. Dapi signal: nucleus image
2. cy3 signal: VSVG signal
3. Immunostaining signal: VSVG signal exposed to the extracellular space. 

Dapi signal is used for locating individual cells and marking them. VSVG signal tells us where the VSVG protein is located within the cell. When they are at the Golgi, they are concentrated close to the nucleus. When they are in ER and the plasma membrane, they are evenly distributed within the cell. To check the VSVG protein specifically localized to the plasma membrane, the third image from immunostaining of the VSVG protein is used: the signal only appears when the protein is exposed to the extracellular space.

###Prescreening

There are bad images just because of simple failures during acquisition, such as

- out of focus
- over exposure

We exclude these images by prescreening images. Standard deviation of pixel intensity in each image is measured, and those with extreme values will be excluded from the analysis. For this, a black list of file names which will be used in the main analysis is created.

We probably will not explain this part in detail, just show you how it was done.  

###Nucleus Segmentation

To measure these values, we first segment each nucleus using atutomatic local intensity thresholding (Bernsen algorithm). Then the rim of the nucleus is dilated outwards to define an approximate area to define the region of **Cell Area**. Within this region, we call the dilated region (Cell Area - Nucelus Area) as the **Juxta-Nuclear** cytoplasmic region. 

(Figure required here)

By measuring the intensity ratio of those at the plasma membrane and those at the Golgi (the transport ratio **T**), we assess the effeciency of transport.

### Transport Ratio
The transport ratio (T) calculated in Simpson (2007) is as follows:

Mean intensity of PM signal / Mean Intensity of VSVG signal. 

If **T** is smaller than control cells, transport is blocked. If it is larger than control, then the transport is accelerated. 

## Resources

###Fiji

Fiji is located at

	/g/almf/software/bin2/fiji

Please put the following lines in your .bashrc or .bash_profile

	export PATH=/g/almf/software/bin2:$PATH
	export PATH=/g/almf/software/bin/java/jdk1.6.0_45/bin:$PATH

Then 

	source .bashrc

or

	source .bash_profile
	
###Data

All data are located under

	/g/data/bio-it_centres_course/data/VSVG/

Plates are numbered such that

	plateID-Number--year-month-date

For exmaple:

	0001-03--2005-08-01

This means that this is the plate ID 1, third experiment (same plate configuration is repeated at least three times each). We analyze plateID from 1 to 159. There are 703 folders so the cluster job array will be composed of 703 jobs.  

Within each of these folders, images are stored under 
	
	data/

###Scripts

The script for processing images for a single plate is

	/g/almf/software/scripts2/measTransportBatch3.py
	
To process a plate, for example the first plate '0001-03--2005-08-01', the plate folder name is given as an argument to a command like below:

	fiji --headless --mem=1000m /g/almf/software/scripts2/measTransportBatch3.py '0001-03--2005-08-01'
   
Prescreening script is 

	/g/almf/software/scripts2/listFocusedImagesV2.py
	
A script for generating job array written by Clemens Lankner is
	
	
## Workflow: GUI

The protocol below is described in Simpson (2007).

###Background subtraction 

Background subtraction will only be done with PM and VSVG images. In the final script, the funciton name is **backgroundSubtraction**

Get minimum and mean intensity of the image.

	[Analyze > Set measurements…] : min and max, mean
	[Analyze> Measure]
		
Results will be shown in the Results Table. We use min and mean values. We call them min1, mean1. 

Then calculate mean of the pixels those with values greater than min1 and less than the mean1. 

	[Analyze > Set measurements…]
		:select min and max, mean, limit to threshold
	[Image > Adjust > threshold…]
		click "Set" and put 
			Lower Threshold Level: min1
			Upper Threshold level: mean1
	[Analyze> Measure]
		results: mean2
	Click "Reset" in Threshold GUI.      

Subtract mean2 value from the image. 
     
	[Process > Math > Subtract]: use "mean2". 

###nucleus segmentation

Nucleus segmentration. We only do this with Dapi image. Corresponding function in the script is **nucleusSegmentation(imp2)**.
 
	[Image > Type > 8bit]
	[Process > Filters > Gaussian Blur...]
		: radius=2.0
	[Image . Adjust > Auto Local Threshold]
		: algorithm, Bernsen, other parameters stay same
   
Pre-filtering by size and circularity

	[Analyze > Set Measurements]
          :select  area, standard deviation, shape descriptors, integrated density
	[Analyze > Analyze Particles…]
		:set area 20 - 1000
		:set circularity 0.8 - 1.0
		:show: Masks
		:check Display Results, Clear results, Exclude on Edges, Include Holes
Rename the resulted image as impNuc

	[Image > Rename]

Duplicate the resulted image and then invert LUT

	[Image > Duplicate…]
	[Image > Lookup Tables > Invert LUT]
          
Dilate segmented nucleus by 15 pixels (to check if nucleus will merge)
	
	[Process > Binary > Options] set iterations to 15, Check "Black Background"
	[Process > Binary > Dilate]

(called imp2)


Particle analysis 1 (with imp2) 

	[Analyze > Analyze Particles…]
		:set area 20 - 10000
		:set circularity 0.8 - 1.0
		:show: None
		:check Display Results, Clear results, Exclude on Edges, Include Holes

q1, q3, offset calculation: Copy and paste <https://gist.github.com/cmci/9641262>

     Example results: 2926.0 3854.0 1392.0
     min2 = q1 - offset (= 2462)
     max2 = q3 + offset (= 5246)


Particle analysis 2 (remove outlier)

	with imp2, 
	[Analyze > Analyze Particles…]
		:set area min2 - max2
		:set circularity 0.8 - 1.0
		:show: Mask
		:check Display Results, Clear results, Exclude on Edges, Include Holes
     
	output: the mask. rename it as impDilatedNuc

Filter original nulcues

     [Process > Image Calculator …]
          AND operation for impNuc and impDilatedNuc


###Check Segmentation Results. 

     Visual inspection. We terminate analysis if there is no nucleus segmented. 


###Generate nuclus ROIs from segmented nucleus. Store them as a list. 

     [Image > Selection > Create Selection]
     [Image > Selection > Make Inverse]
     [Image > Selection > To ROImanager]

###Roi modifications

     In ROI manager
     a. Generate All Cell ROIs (enlarge)
         Activate the ROI above, then 
         [Image > Selection > Enlarge…] : 15 pixels
          click "Add". Second ROI appears in the ROI manager list. 
          this will be the "All cell"
     b. Generate JuxtaNuclear ROIs (enlarged minus nucleus)
          select both ROIs, and then [More >> XOR]. 
          click "Add". Third ROI appears in the listing. 
          This ail be the zone surrounding nucleus. 


### Measurements
Do measurements using modified ROIs (measureROIs).

     select corresponding ROI and 
     
     [Analyze > Measure]
     
     a. enlarged: all cell plasma membrane
     b. ring rois: juxtanuclear signal for VSVG
     c. all cell VSVG intensity


### Transport Effeciency

Calculate transport efficiency and add the results to results table. All results from single plate are merged and placed in a file, but we do not do this manually. 

## Workflow: Python Primer

Jython is Java-implemented Python, allowing us to access Java classes in Python syntax. To write Jython scripts in Fiji, we use Script Editor. 

	[File > New > Script]

This command launches the editor and a new script file. 

Script Editor has ist own menu, and one of them is 

	[Language > Python]
	
To set the language to Python (Jython).

 The editor has two panels. Upper one is the input field, and the lower one is the output field. Output field has two types, standard output and error output. You could toggle this by clicking "error" or "output".
 
We learn the basics of Python commands interactively. Below is a list of Python funtions we go over.

````
print
print with comma separated

variables

functions
     return
     multiple returns
     
list
     initialization
     range
     sorted(list)
     slicing
     l.append()
     len

for

     allcellA = [roiEnlarger(r) for r in nucroiA]

     enumerate
     extend()

if, else, break

String
     concatenation
     f.endswith()
     
dictionary
     initialization
     d['a'] = 3

Accessing files
     path / file managements
          os.path.listdir
          os.path.isfile
          os.path.join

Python regex
          re.compile
          re.search
          re.group
          
I/O 
f = open(path, 'rb')
sys.argv

csv.reader

````

Codes used in this interactive session is in 

<https://gist.github.com/cmci/9687511>

For more information, see

<http://www.jython.org/docs/library/indexprogress.html><http://www.jython.org/docs/library/functions.html>


## Workflow: Fiji Jython Scripting

### IJ.log

Let's first try accessing Java classes. Try the following. 

````
IJ.log("Hello World!")
````

This will print the argument "Hello World!" in the Log window (in ImageJ macro, similar job can be done by print(string) command).

We could add another line

````
IJ.log("Hello World!")
IJ.log("\\Clear")
````
Run this code, and compare with the output of 

````
IJ.log("\\Clear")
IJ.log("Hello World!")
````

The double back-slash

	\\ 

is called escape sequence, and is a command send to the interpretor. Insted of printing the following string in the Log window, it is interpreted as a command. In case of "\\Clear", texts in the Log window will be cleared. 

### Javadoc

Now, let's see a bit more of details about "IJ.log". Click the following link and see the page

[Javadoc - IJ](http://rsbweb.nih.gov/ij/developer/api/ij/IJ.html)

This is a documentation called "Javadoc" and is a reference for all the functionality that an Java applicaiton has. It's **VERY IMPORTANT** to understand how to use this documentation, so I wll explain in details. 

ImageJ is a collection of many classes, and you could find all those classes in this javadoc. The web page has side bar in the left side, and it is separated to the upper and the lower panels. The upper panel has its title "Packages" and the bottle panel is "Classes". As I have already said, ImageJ is a colleciton of many classes - and they are all listed in the bottom panel. 

Since there are many classes, they are organized in hierarchical way by "packages". You could imagine each packages as a folder in file system. For example, you could find a package named "ij". It means something like a folder named "ij", and for a package "ij.plugin.filter", its a folder "filter" under a folder "plguin" within a folder "ij" (in fact, class files are located in file system in this way). 

In the bottom panel, all classes are listed as long as you have not filtered them by clicking one of the packages in the upper panel. Try to look for the class **IJ** by scrolling down the lower panel, and click the link "IJ". 

You will then see a new page in the right side: There is a big text saying "Class IJ". IJ from IJ.log is a **class**, which contains many **methods**. To see those methods, scroll down the page and you will find a table titled "Method Summary". This is an alphabetical list of all methods that IJ contains. Among them, you could find 

	log(java.lang.String s) 

To use a method from IJ class, you just need a dot after the class name and then add the name of the method. In case of the method "log", you need an argument and it should be a String type. In the example we tried above, it was "Hello World!". 

How about other methods? Try the method beep(). 

	IJ.beep()
	
When this code is executed, your machine will make a beep sound. In such a way all the methods are usable. 

###Loading and showing Image

	imp = IJ.openImage('/Users/miura/image.tif')
	
	imp.show()
	
###Class ImagePlus


imp.getTitle()
imp.deleteRoi()
imp.setRoi()
imp.duplicate()
imp.getProcessor()

###Class ImageProcessor


impstats = ip.getStatistics()
ip.getMax()
ip.setThreshold(min, max, lut)
ip.resetThreshold()

ip.subtract(val)
ip.invertLut()
r = ip.restRoi()
ip.setRoi()

###Class ImageStatistics

impstats.mean


###Constructing an Object


	imp = IJ.getImage()
	pf = ProfilePlot(imp)
	profile = pf.getProfile()
	for val in profile:
		print val
		
We combine this with the python package csv to write the output to a comma separated value text file. 

	import csv
	
	imp = IJ.getImage()
	pf = ProfilePlot(imp)
	profile = pf.getProfile()
	for val in profile:
		print val
	f = open('/Users/miura/Desktop/prof.csv', 'wb')
	writer = csv.writer(f)
	for index, val in enumerate(profile):
		writer.writerow([index, val])
	f.close()


## Workflow: Implementing Image Analysis

### Prescreening of Images.

<https://gist.github.com/cmci/9684453>

Prescreen results are saved under

	/g/data/bio-it_centres_course/data/VSVG/prescreen

### Nucleus segmentation

A script for processing single image dataset is:

<https://gist.github.com/cmci/9644827>

This script segments nucleus. 

### ROI generator

### Measurements

#### Transport Ratio

#### Haralick Features

## Final Code

