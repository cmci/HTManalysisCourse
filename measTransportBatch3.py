"""
Measure with  Simpson (2007) primary screening method, 
additionally with 
1. Haralick features
2. (readial gradient) ... not implemented yet

example command:
  fiji --mem=2000m --headless measTransportbatch.py 0137-17--2006-05-06 out4

20140200 miura@embl.de
http://cmci.embl.de
"""
from ij import IJ, Prefs
from ij.process import ImageProcessor, ImageStatistics, ImageConverter
from ij.measure import ResultsTable
from ij.gui import Roi, ShapeRoi
from ij.plugin import ImageCalculator, Duplicator, RoiEnlarger
from ij.plugin.filter import GaussianBlur
from ij.plugin.filter import ParticleAnalyzer as PA
from ij.plugin.filter import RankFilters
from ij.plugin.filter import Analyzer
from ij.plugin.filter import ThresholdToSelection
from ij.plugin.frame import RoiManager
from fiji.threshold import Auto_Local_Threshold as ALT
import jarray
from org.apache.commons.math3.stat.descriptive import DescriptiveStatistics as DSS
from org.jfree.data.statistics import BoxAndWhiskerCalculator
from java.util import ArrayList, Arrays
from emblcmci.glcm import GLCMtexture
import os, csv, re, sys

# size of juxtanuclear region. In pixels.
RIMSIZE = 15
# image background is expected to be black. 
Prefs.blackBackground = True
# verbose output
VERBOSE = False
# visual outputs, such as Images, Results and ROI lists. Only possible with graphic card. 
GUIMODE = False
# test mode limits execution only to first 10 image sets.
TESTMODE = False
# Root folder of the data.
#rootfolder = '/Volumes/data/bio-it_centres_course/data/VSVG'
rootfolder = '/g/data/bio-it_centres_course/data/VSVG'


class InstBWC(BoxAndWhiskerCalculator):
    def __init__(self):
        pass

def backgroundSubtraction(imp):
  """ subtract background, Cihan's method.
    see Simpson(2007)
  """
  impstats = imp.getProcessor().getStatistics()
#    backlevel = impstats.min + (impstats.mean - impstats.min)/2
  imp.getProcessor().setThreshold(impstats.min, impstats.mean, ImageProcessor.RED_LUT)
  measOpt = ImageStatistics.MEAN + ImageStatistics.LIMIT
  impstats = ImageStatistics.getStatistics(imp.getProcessor(), measOpt, None)
  backlevel = impstats.mean
  imp.getProcessor().resetThreshold()
  imp.getProcessor().subtract(backlevel)
  print imp.getTitle(), " : background intensity - ", backlevel
  return backlevel

def roiRingGenerator(r1):
  """ Create a band of ROI outside the argument ROI.
  See Liebel (2003) Fig. 1
  """
  #r1 = imp.getRoi()
  r2 = RoiEnlarger.enlarge(r1, RIMSIZE)
  sr1 = ShapeRoi(r1)
  sr2 = ShapeRoi(r2)
  return sr2.not(sr1)

def roiEnlarger(r1):
  """ Enlarges ROI by a defined iterations.
  """
  return ShapeRoi(RoiEnlarger.enlarge(r1, RIMSIZE))

def getOutlierBound(rt):
  """ Analyzes the results of the 1st partcile analysis.
  Since the dilation of nuclear perimeter often causes 
  overlap of neighboring neculeus 'terrirories', such nucleus 
  are discarded from the measurements. 

  Small nucelei are already removed, but since rejection of nuclei depends on 
  standard outlier detection method, outliers in both smaller and larger sizes
  are discarded. 
  """
  area = rt.getColumn(rt.getColumnIndex('Area'))
  circ = rt.getColumn(rt.getColumnIndex("Circ."))
  arealist = ArrayList(Arrays.asList(area.tolist()))
  circlist = ArrayList(Arrays.asList(circ.tolist()))
  bwc = InstBWC()
  ans = bwc.calculateBoxAndWhiskerStatistics(arealist)
  #anscirc = bwc.calculateBoxAndWhiskerStatistics(circlist)
  if (VERBOSE):
    print ans.toString()
    print ans.getOutliers()
  q1 = ans.getQ1()
  q3 = ans.getQ3()
  intrange = q3 - q1 
  outlier_offset = intrange * 1.5
  # circularity better be fixed. 
  #circq1 = anscirc.getQ1()
  #circq3 = anscirc.getQ3()
  #circintrange = circq3 - circq1 
  #circoutlier_offset = circintrange * 1.5
  return q1, q3, outlier_offset

def nucleusSegmentation(imp2):
    """ Segmentation of nucleus image. 
    Nucleus are selected that:
    1. No overlapping with dilated regions
    2. close to circular shape. Deformed nuclei are rejected.
    Outputs a binary image.
    """
#Convert to 8bit
    ImageConverter(imp2).convertToGray8()
#blur slightly using Gaussian Blur 
    radius = 2.0
    accuracy = 0.01
    GaussianBlur().blurGaussian( imp2.getProcessor(), radius, radius, accuracy)
# Auto Local Thresholding
    imps = ALT().exec(imp2, "Bernsen", 15, 0, 0, True)
    imp2 = imps[0]


#ParticleAnalysis 0: prefiltering by size and circularity
    rt = ResultsTable()
    paOpt = PA.CLEAR_WORKSHEET +\
                    PA.SHOW_MASKS +\
                    PA.EXCLUDE_EDGE_PARTICLES +\
                    PA.INCLUDE_HOLES #+ \
#		PA.SHOW_RESULTS 
    measOpt = PA.AREA + PA.STD_DEV + PA.SHAPE_DESCRIPTORS + PA.INTEGRATED_DENSITY
    MINSIZE = 20
    MAXSIZE = 10000
    pa0 = PA(paOpt, measOpt, rt, MINSIZE, MAXSIZE, 0.8, 1.0)
    pa0.setHideOutputImage(True)
    pa0.analyze(imp2)
    imp2 = pa0.getOutputImage() # Overwrite 
    imp2.getProcessor().invertLut()
#impNuc = imp2.duplicate()	## for the ring. 
    impNuc = Duplicator().run(imp2)

#Dilate the Nucleus Area
## this should be 40 pixels in Cihan's method, but should be smaller. 
#for i in range(20):
#	IJ.run(imp2, "Dilate", "")
    rf = RankFilters()
    rf.rank(imp2.getProcessor(), RIMSIZE, RankFilters.MAX)

#Particle Analysis 1: get distribution of sizes. 

    paOpt = PA.CLEAR_WORKSHEET +\
                    PA.SHOW_NONE +\
                    PA.EXCLUDE_EDGE_PARTICLES +\
                    PA.INCLUDE_HOLES #+ \
#		PA.SHOW_RESULTS 
    measOpt = PA.AREA + PA.STD_DEV + PA.SHAPE_DESCRIPTORS + PA.INTEGRATED_DENSITY
    rt1 = ResultsTable()
    MINSIZE = 20
    MAXSIZE = 10000
    pa = PA(paOpt, measOpt, rt1, MINSIZE, MAXSIZE)
    pa.analyze(imp2)
    #rt.show('after PA 1')
#particle Analysis 2: filter nucleus by size and circularity. 
    #print rt1.getHeadings()
    if (rt1.getColumnIndex('Area') > -1):
      q1, q3, outlier_offset = getOutlierBound(rt1)
    else:
      q1 = MINSIZE
      q3 = MAXSIZE
      outlier_offset = 0
      print imp2.getTitle(), ": no Nucleus segmented,probably too many overlaps"

    paOpt = PA.CLEAR_WORKSHEET +\
                    PA.SHOW_MASKS +\
                    PA.EXCLUDE_EDGE_PARTICLES +\
                    PA.INCLUDE_HOLES #+ \
#		PA.SHOW_RESULTS 
    rt2 = ResultsTable()
    #pa = PA(paOpt, measOpt, rt, q1-outlier_offset, q3+outlier_offset, circq1-circoutlier_offset, circq3+circoutlier_offset)
    pa = PA(paOpt, measOpt, rt2, q1-outlier_offset, q3+outlier_offset, 0.8, 1.0)
    pa.setHideOutputImage(True)
    pa.analyze(imp2)
    impDilatedNuc = pa.getOutputImage() 

#filter original nucleus

    filteredNuc = ImageCalculator().run("AND create", impDilatedNuc, impNuc)
    return filteredNuc    

def measureTexture(imp, thisrt, roiA):
  """ texture measurement
  """
  ip8 = imp.duplicate().getProcessor().convertToByte(True)
  glmc = GLCMtexture(1, 45, True, False)
  for index, r in enumerate(roiA):
    ip8.resetRoi()
    #br = Roi(r.getBounds())
    if r.getBounds().getWidth() != 0 and r.getBounds().getHeight() != 0:
      ip8.setRoi(r)
      glmc.calcGLCM(ip8)
      resmap = glmc.getResultsArray()
      pa = glmc.paramA
      #print 'cell', index, r
      for p in pa:
        thisrt.setValue(p, index, resmap.get(p))
        #print p, resmap.get(p)
def measureROIs(imp, measOpt, thisrt, roiA, backint, doGLCM):    
  """ Cell-wise measurment using ROI array. 
  """
  analObj = Analyzer(imp, measOpt, thisrt)
  for index, r in enumerate(roiA):
    imp.deleteRoi()
    imp.setRoi(r)
    analObj.measure()
    maxint = thisrt.getValue('Max', thisrt.getCounter()-1)
    saturation = 0
    if ( maxint + backint) >= 4095:
      saturation = 1
      if (VERBOSE):
        print 'cell index ', index, 'maxint=', maxint
    thisrt.setValue('CellIndex', thisrt.getCounter()-1, index)
    thisrt.setValue('Saturation', thisrt.getCounter()-1, saturation)
  if (doGLCM):
    imp.deleteRoi()
    measureTexture(imp, thisrt, roiA)

def checkPrescreenResult(folderpath, aplate, wnumber):
  """ Accesses the prescreening result csv file, check if it was decided to be
  abnormal (failed in capturing, out of focus, none even illuminations...)
  """
  filename = aplate + '.csv'
  csvfilepath = os.path.join(folderpath, 'prescreen', filename )
  f = open(csvfilepath, 'rb')
  data = csv.reader(f, delimiter=',')
  healthy = True
  for row in data:
    f = row[0]
    if (f.startswith("'--W" + wnumber)):
      #print row
      if row[1] == '1':
        healthy = False
      break
  return healthy

def procOneImage(pathpre, wnumber, endings):
  """ Analyzes a single image set (Dapi, VSVG, PM images)
  pathpre: fullpath prefix, down till "endings". 
  endings: a dictionary with signiture for three different channels. 
  wnumber: a number in string, indicating the spot ID.
  Returns three results tables. 
  """
  imp = IJ.openImage(pathpre + endings['dapi'] + '.tif')
  impVSVG = IJ.openImage(pathpre + endings['vsvg'] + '.tif')
  impPM = IJ.openImage(pathpre + endings['pm'] + '.tif')
  imp2 = imp.duplicate()

  rtallcellPM = ResultsTable()
  rtjnucVSVG = ResultsTable()
  rtallcellVSVG = ResultsTable()

  backVSVG = backgroundSubtraction(impVSVG)
  backPM = backgroundSubtraction(impPM)
  impfilteredNuc = nucleusSegmentation(imp2)

  intmax = impfilteredNuc.getProcessor().getMax()
  if intmax == 0:
    return rtallcellPM, rtjnucVSVG, rtallcellVSVG

  impfilteredNuc.getProcessor().setThreshold(1, intmax, ImageProcessor.NO_LUT_UPDATE)
  nucroi = ThresholdToSelection().convert(impfilteredNuc.getProcessor())
  nucroiA = ShapeRoi(nucroi).getRois()
#print nucroiA
  allcellA = [roiEnlarger(r) for r in nucroiA]
  jnucroiA = [roiRingGenerator(r) for r in nucroiA]
#print allcellA
  print 'Detected Cells: ', len(jnucroiA)  
  if len(jnucroiA) <2:
      print "measurement omitted, as there is only on nucleus detected"
      return  rtallcellPM, rtjnucVSVG, rtallcellVSVG
  if (GUIMODE):
    rm = RoiManager()
    for r in jnucroiA:
      rm.addRoi(r)
    rm.show()
    impfilteredNuc.show()
  

  measOpt = PA.AREA + PA.MEAN + PA.CENTROID + PA.STD_DEV + PA.SHAPE_DESCRIPTORS + PA.INTEGRATED_DENSITY + PA.MIN_MAX +\
    PA.SKEWNESS + PA.KURTOSIS + PA.MEDIAN + PA.MODE

## All Cell Plasma Membrane intensity
  measureROIs(impPM, measOpt, rtallcellPM, allcellA, backPM, True)
  meanInt_Cell = rtallcellPM.getColumn(rtallcellPM.getColumnIndex('Mean'))
  print "Results Table rownumber:", len(meanInt_Cell)
# JuxtaNuclear VSVG intensity 
  measureROIs(impVSVG, measOpt, rtjnucVSVG, jnucroiA, backVSVG, False)    
  meanInt_jnuc = rtjnucVSVG.getColumn(rtjnucVSVG.getColumnIndex('Mean'))

# AllCell VSVG intensity 
  measureROIs(impVSVG, measOpt, rtallcellVSVG, allcellA, backVSVG, True)    
  meanInt_vsvgall = rtallcellVSVG.getColumn(rtallcellVSVG.getColumnIndex('Mean'))
  
#Calculation of Transport Ratio JuxtaNuclear VSVG intensity / All Cell Plasma Membrane intensity results will be appended to PM results table.
  for i in range(len(meanInt_Cell)):
    if meanInt_Cell[i] != 0.0:
      transportR = meanInt_jnuc[i] / meanInt_Cell[i]
      transportRall = meanInt_vsvgall[i] / meanInt_Cell[i]
    else:
      transportR = float('inf')
      transportRall = float('inf')
    rtjnucVSVG.setValue('TransportRatio', i, transportR)
    rtallcellVSVG.setValue('TransportRatio', i, transportRall)
    rtjnucVSVG.setValue('WellNumber', i, int(wnumber)) 
    rtallcellVSVG.setValue('WellNumber', i, int(wnumber))
    rtallcellPM.setValue('WellNumber', i, int(wnumber)) 
  return rtallcellPM, rtjnucVSVG, rtallcellVSVG


def concatResultsTable(rt, rtAll):
  """ Appending measurement results to the previous measurements. 
  Intended for collecting all measurement results from a plate to a single csv file. 
  """
  heads = rt.getHeadings()
  # there should be 32 columns. This number is better be generalized. 
  #if len(heads) >= 32:
  if len(rtAll) == 0:
    for head in heads:
        col = rt.getColumnAsDoubles(rt.getColumnIndex(head))
        rtAll.append(col)
  else:
    for i, head in enumerate(heads):
      col = rt.getColumnAsDoubles(rt.getColumnIndex(head))
      rtAll[i].extend(col)
  #else: 
    #print "--- a mismatch in results table header length! ---"
    #print " omitted from merging the results from this spot"
    #print "   Parent headers:", len(rtAll)
    #print "   current headers:", len(heads), heads
  ##return heads

def outputResultsTable(outpath, rt, headings):
  """ CSV writer
  """
  trt = map(list, zip(*rt))
  f = open(outpath, 'wb')
  writer = csv.writer(f)
  writer.writerow(headings)
  writer.writerows(trt)
  f.close()

def measurePlate(rootfolder, aplate, outfolder):
  """ List files in a plate directory and run the measurement for each images. 
  Image set is detected by image names. 

  """
  platepath = os.path.join(rootfolder, aplate, 'data')
  print platepath
  pattern = re.compile('(.*)--W(.*)--P(.*)--Z(.*)--T(.*)--(.*)\.(.*)')
  pat2 = re.compile('(.*--W.*--P.*--Z.*--T.*--)(.*)\.tif')
  files = []
  # dictionary to keep channel signature strings
  endings = {'dapi':'', 'vsvg':'', 'pm':''}
  for f in os.listdir(platepath):
    if os.path.isfile(os.path.join(platepath, f)):
      res = re.search(pat2, f)
      if f.endswith('dapi.tif'):
        files.append(f)
        if endings['dapi'] == '': endings['dapi'] = res.group(2)
      elif f.endswith('cfp.tif'):
        if endings['vsvg'] == '': endings['vsvg'] = res.group(2)
      elif f.endswith('647.tif'):
        if endings['pm'] == '': endings['pm'] = res.group(2)
  files = sorted(files)

  allcellPMA = []
  jnucVSVGA = []
  allcellVSVGA = []

  if (TESTMODE):
    listoffiles = files[0:13]
  else:
    listoffiles = files

  for f in listoffiles: 
    res = re.search(pattern, f)
    wnumber = res.group(2)
    res2 = re.search(pat2, f)
    filepre = res2.group(1)
    pathpre = os.path.join(rootfolder, aplate, 'data', filepre)
    if checkPrescreenResult(rootfolder, aplate, wnumber):
      print '=======', wnumber
      rtallcellPM, rtjnucVSVG, rtallcellVSVG= procOneImage(pathpre, wnumber, endings)
      if rtallcellPM.getColumnIndex('Area') > -1:
        concatResultsTable(rtallcellPM, allcellPMA)
        concatResultsTable(rtjnucVSVG, jnucVSVGA)
        concatResultsTable(rtallcellVSVG, allcellVSVGA)    
      
      if (GUIMODE):
        rtallcellPM.show('PM')
        rtjnucVSVG.show('VSVG_RING')
        rtallcellVSVG.show('VSVG_ALL')  
    else:
      print "Rejected in Prescreen:", f

  outPMallPath = os.path.join(rootfolder, outfolder, aplate + '--PMall.csv')
  outputResultsTable(outPMallPath, allcellPMA, rtallcellPM.getHeadings())

  outVSVGallPath = os.path.join(rootfolder, outfolder, aplate + '--VSVGall.csv')
  outputResultsTable(outVSVGallPath, allcellVSVGA, rtallcellVSVG.getHeadings())

  outVSVGjnucPath = os.path.join(rootfolder, outfolder, aplate + '--VSVGjnuc.csv')
  outputResultsTable(outVSVGjnucPath, jnucVSVGA, rtjnucVSVG.getHeadings())

if len(sys.argv) > 1:
  aplate = sys.argv[1]
else:
  aplate = '0001-03--2005-08-01'
if len(sys.argv) > 2:
  outfolder = sys.argv[2]
else:
  outfolder = 'out4'
measurePlate(rootfolder, aplate, outfolder)

