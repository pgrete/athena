import sys
import glob
import os
import os.path
import shutil
import commands
import visit
import time

def viewit():
   tl = visit.CreateAnnotationObject("Text2D")
   tl.position = (.08, .93)
   tl.height = .026
   tl.fontFamily = 0
   tl.fontBold = 1
   tl.fontShadow = 1

   visit.DeleteAllPlots()

   k = 30000
   kMax = 40001

   while k < kMax:
      filename = "scalar" + str(k) + ".Point3D"
      print "opening database file:  " + filename
      status = visit.OpenDatabase(filename, 0, "Point3D")

      if status != 1:
	print "Could not open " + filename
	return

      status = visit.AddPlot("Pseudocolor", "density")
      visit.DrawPlots()
      pca = visit.GetPlotOptions()
      pca.colorTableName = "density"
      pca.invertColorTable = 1
      pca.minFlag = 1
      pca.maxFlag = 1
      pca.min = -.21
      pca.max =  .21
#     visit.SetPlotOptions(pca)

      names = visit.GetAnnotationObjectNames()
      print names
      lastName = names[-1]
      print "lastName = " + lastName

      legend = visit.GetAnnotationObject(lastName)
      legend.numTicks = 3
      legend.managePosition = 0
      legend.position = (.1, .88)
      legend.orientation = 1
      legend.drawMinMax = 0
      legend.drawTitle = 0
      legend.drawLabels = "Labels"
      legend.suppliedLabels = (" -.15", " 0.", " .15")
      legend.fontHeight = .024
      legend.numberFormat = "%#-9.2g"
      legend.fontFamily = 2
      legend.fontBold = 1

      tl.text = "Density Fluctuation  -  Step " + str(k)
      tl.visible = 1
      visit.SaveWindow()
      k = k + 1000

      if k < kMax:
      	tl.visible = 0
   	visit.DeleteAllPlots()
        visit.CloseDatabase(filename)

#  tl.visible = 0

# Here is the main
windowA = visit.SaveWindowAttributes()
windowA.format = windowA.JPEG
windowA.progressive = 1
windowA.width = 807
windowA.height = 908
windowA.fileName = "fdensity"
#windowA.outputDirectory = "/global/u2/e/efeibush/gtc/flow/flow3/compare/densityimages"
windowA.outputToCurrentDirectory = 1
#windowA.screenCapture = 1
visit.SetSaveWindowAttributes(windowA)

#a = visit.ListPlots()
b = visit.SetActivePlots(0)

anot = visit.GetAnnotationAttributes()
anot.databaseInfoFlag = 0
visit.SetAnnotationAttributes(anot)

viewit()
