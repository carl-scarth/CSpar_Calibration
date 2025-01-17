# state file generated using paraview version 5.12.0-RC3
import os
import paraview

paraview.compatibility.major = 5
paraview.compatibility.minor = 12

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'

paraview.simple._DisableFirstRenderCameraReset()


#-----------------------------------------------------------------

# Set up parameters used to get the input files
filestr = "LHSDesign100x8"
vtk_path = "C:\\Users\\cs2361\\Documents\\CSpar_Calibration\\outputs\\gp_predictions_nonlinear_" + filestr + ""
filenames = os.listdir(vtk_path)
#filenos = [int(''.join(filter(str.isdigit, file))) for file in filename]  # Get the numeric part

# Sort out the string into natural counting order using base python
file0, ext = filenames[0].split(".")
prefix = file0.strip("0123456789")
sorted_suffix = sorted([int(''.join(filter(str.isdigit,file))) for file in filenames])
filenames = [prefix+str(suffix)+"."+ext for suffix in sorted_suffix]
filenames = [os.path.join(vtk_path, file) for file in filenames]

output_list = [{"name":"eta_sam_mu_u", "range":[-3.8, 0], "title":"Mean $u_1$ (mm)", "Labels": [0.0, -1.0, -2.0, -3.0, -3.8], "LabelFormat":"%-#6.1f"},
               {"name":"eta_sam_mu_w", "range":[-1.1, 0], "title":"Mean $u_3$ (mm)", "Labels": [0.0, -0.25, -0.5, -0.75, -1.1], "LabelFormat":"%-#6.2f"},
               {"name":"eta_sam_sigma_u", "range":[0, 0.016], "title":"$u_1$ standard \n deviation (mm)", "Labels": [0.00, 0.004, 0.008, 0.012, 0.016], "LabelFormat":"%-#6.1e"},
               {"name":"eta_sam_sigma_w", "range":[0, 0.006], "title":"$u_3$ standard \n deviation (mm)", "Labels": [0.00, 0.002, 0.004, 0.006], "LabelFormat":"%-#6.1e"}
          ]


# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1571, 794]
renderView1.AxesGrid = 'Grid Axes 3D Actor'
renderView1.OrientationAxesInteractivity = 1
renderView1.CenterOfRotation = [30.5, 75.0, 210.0]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [-440.3718938869471, 309.7178349580358, -215.2347379949773]
renderView1.CameraFocalPoint = [80.13445297892076, -42.35593989511255, 397.36779614683184]
renderView1.CameraViewUp = [0.272463763275882, 0.9158702835120002, 0.29486458159838863]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 227.13707315187452
renderView1.LegendGrid = 'Legend Grid Actor'
renderView1.BackEnd = 'OSPRay raycaster'
renderView1.OSPRayMaterialLibrary = materialLibrary1

SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.AssignView(0, renderView1)
layout1.SetSize(1028, 463)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'Legacy VTK Reader'
frame_0vtk = LegacyVTKReader(registrationName='frame_0.vtk*', FileNames=filenames)

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from frame_0vtk
frame_0vtkDisplay = Show(frame_0vtk, renderView1, 'UnstructuredGridRepresentation')


for output in output_list:
    # get 2D transfer function
    TF2D = GetTransferFunction2D(output["name"])
    TF2D.ScalarRangeInitialized = 1
    TF2D.Range = output["range"] + [0.0, 1.0]

    # get color transfer function/color map
    LUT = GetColorTransferFunction(output["name"])
    LUT.TransferFunction2D = TF2D
    LUT.RGBPoints = [output["range"][0], 0.231373, 0.298039, 0.752941, (output["range"][0]+output["range"][1])/2, 0.865003, 0.865003, 0.865003, output["range"][1], 0.705882, 0.0156863, 0.14902]
    LUT.ScalarRangeInitialized = 1.0

    # get opacity transfer function/opacity map
    PWF = GetOpacityTransferFunction(output["name"])
    PWF.Points = [output["range"][0], 0.0, 0.5, 0.0, output["range"][1], 1.0, 0.5, 0.0]
    PWF.ScalarRangeInitialized = 1

    # setup the color legend parameters for each legend in this view
    LUTColorBar = GetScalarBar(LUT, renderView1)
    LUTColorBar.WindowLocation = 'Any Location'
    LUTColorBar.Position = [0.7, 0.127]
    LUTColorBar.Title = output["title"]
    LUTColorBar.ComponentTitle = ''
    LUTColorBar.HorizontalTitle = 1
    LUTColorBar.TitleFontFamily = 'Times'
    LUTColorBar.TitleFontSize = 40
    LUTColorBar.LabelFontFamily = 'Times'
    LUTColorBar.LabelFontSize = 30
    LUTColorBar.ScalarBarThickness = 54
    LUTColorBar.ScalarBarLength = 0.675
    LUTColorBar.AutomaticLabelFormat = 0
    LUTColorBar.LabelFormat = output["LabelFormat"]
    LUTColorBar.UseCustomLabels = 1
    LUTColorBar.CustomLabels = output["Labels"]
    LUTColorBar.AddRangeLabels = 0
    LUTColorBar.RangeLabelFormat = output["LabelFormat"]
    LUTColorBar.Visibility = 0


# trace defaults for the display properties.
frame_0vtkDisplay.Representation = 'Surface'
frame_0vtkDisplay.ColorArrayName = ['POINTS', output_list[-1]["name"]]
frame_0vtkDisplay.LookupTable = LUT
frame_0vtkDisplay.SelectTCoordArray = 'None'
frame_0vtkDisplay.SelectNormalArray = 'None'
frame_0vtkDisplay.SelectTangentArray = 'None'
frame_0vtkDisplay.OSPRayScaleArray = output_list[-1]["name"]
frame_0vtkDisplay.OSPRayScaleFunction = 'Piecewise Function'
frame_0vtkDisplay.Assembly = ''
frame_0vtkDisplay.SelectOrientationVectors = 'None'
frame_0vtkDisplay.ScaleFactor = 42.0
frame_0vtkDisplay.SelectScaleArray = 'None'
frame_0vtkDisplay.GlyphType = 'Arrow'
frame_0vtkDisplay.GlyphTableIndexArray = 'None'
frame_0vtkDisplay.GaussianRadius = 2.1
frame_0vtkDisplay.SetScaleArray = ['POINTS', output_list[-1]["name"]]
frame_0vtkDisplay.ScaleTransferFunction = 'Piecewise Function'
frame_0vtkDisplay.OpacityArray = ['POINTS', output_list[-1]["name"]]
frame_0vtkDisplay.OpacityTransferFunction = 'Piecewise Function'
frame_0vtkDisplay.DataAxesGrid = 'Grid Axes Representation'
frame_0vtkDisplay.PolarAxes = 'Polar Axes Representation'
frame_0vtkDisplay.ScalarOpacityFunction = PWF
frame_0vtkDisplay.ScalarOpacityUnitDistance = 15.10230564779437
frame_0vtkDisplay.OpacityArrayName = ['POINTS', output_list[-1]["name"]]
frame_0vtkDisplay.SelectInputVectors = [None, '']
frame_0vtkDisplay.WriteLog = ''

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
frame_0vtkDisplay.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 0.0, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
frame_0vtkDisplay.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 0.0, 1.0, 0.5, 0.0]

# show color legend
frame_0vtkDisplay.SetScalarBarVisibility(renderView1, True)
# set color bar visibility
LUTColorBar.Visibility = 1

# ----------------------------------------------------------------
# setup animation scene, tracks and keyframes
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get time animation track
timeAnimationCue1 = GetTimeTrack()

# initialize the animation scene

# get the time-keeper
timeKeeper1 = GetTimeKeeper()

# initialize the timekeeper

# initialize the animation track

# get animation scene
animationScene1 = GetAnimationScene()

# initialize the animation scene
animationScene1.ViewModules = renderView1
animationScene1.Cues = timeAnimationCue1
animationScene1.AnimationTime = 15.0
animationScene1.EndTime = 15.0
animationScene1.PlayMode = 'Snap To TimeSteps'

# ----------------------------------------------------------------
# restore active source
SetActiveSource(frame_0vtk)
# ----------------------------------------------------------------


##--------------------------------------------
## You may need to add some code at the end of this python script depending on your usage, eg:
#
## Render all views to see them appears
# RenderAllViews()
#
## Interact with the view, usefull when running from pvpython
# Interact()
#
## Save a screenshot of the active view
# SaveScreenshot("path/to/screenshot.png")
#
## Save a screenshot of a layout (multiple splitted view)
# SaveScreenshot("path/to/screenshot.png", GetLayout())
#
## Save all "Extractors" from the pipeline browser
# SaveExtracts()
#
## Save a animation of the current active view
# SaveAnimation()
#
## Please refer to the documentation of paraview.simple
## https://kitware.github.io/paraview-docs/latest/python/paraview.simple.html
##--------------------------------------------