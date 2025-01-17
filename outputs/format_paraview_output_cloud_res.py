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
# csv_path = "E:\\Calibration_outputs_for_paper\\prior_residuals_LC\\"
csv_path = "E:\\Calibration_outputs_for_paper\\cloud_output_dataset_LC_100x8\\"

filenames = os.listdir(csv_path)

# Sort out the string into natural counting order using base python
file0, ext = filenames[0].split(".")
prefix = file0.strip("0123456789")
sorted_suffix = sorted([int(''.join(filter(str.isdigit,file))) for file in filenames])
filenames = [prefix+str(suffix)+"."+ext for suffix in sorted_suffix]
filenames = [os.path.join(csv_path, file) for file in filenames]


""" output_list = [{"name":"w_res", "range":[-0.25, 0.25], "title":"Residual $u_3$ (mm)", "Labels": [-0.25, -0.2, -0.1, 0, 0.1, 0.2, 0.25], "LabelFormat":"%-#6.2f"},
               {"name":"w_absres", "range":[0, 0.25], "title":"Abs residual $u_3$ (mm)", "Labels": [0.0, 0.05, 0.1, 0.15, 0.2, 0.25], "LabelFormat":"%-#6.2f"},
               {"name":"u_res", "range":[-0.7, 0.7], "title":"Residual $u_1$ (mm)", "Labels": [-0.7, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.7], "LabelFormat":"%-#6.1f"},
               {"name":"u_absres", "range":[0, 0.7], "title":"Abs residual $u_1$ (mm)", "Labels": [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7], "LabelFormat":"%-#6.1f"}
          ]
 """

output_list = [{"name":"eta_sam_y_mu_w_res", "range":[-0.25, 0.25], "title":"Residual $u_3$ (mm)", "Labels": [-0.25, -0.2, -0.1, 0, 0.1, 0.2, 0.25], "LabelFormat":"%-#6.2f"},
               {"name":"eta_sam_y_mu_w_absres", "range":[0, 0.25], "title":"Abs residual $u_3$ (mm)", "Labels": [0.0, 0.05, 0.1, 0.15, 0.2, 0.25], "LabelFormat":"%-#6.2f"},
               {"name":"eta_sam_y_mu_u_res", "range":[-0.7, 0.7], "title":"Residual $u_1$ (mm)", "Labels": [-0.7, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.7], "LabelFormat":"%-#6.1f"},
               {"name":"eta_sam_y_mu_u_absres", "range":[0, 0.7], "title":"Abs residual $u_1$ (mm)", "Labels": [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7], "LabelFormat":"%-#6.1f"}
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

# create a new 'CSV Reader'
prior_res = CSVReader(registrationName='prior_res.csv*', FileName=filenames)
# create a new 'Table To Points'
tableToPoints1 = TableToPoints(registrationName='TableToPoints1', Input=prior_res)
tableToPoints1.XColumn = 'x_proj'
tableToPoints1.YColumn = 'y_proj'
tableToPoints1.ZColumn = 'z_proj'

#adsadasd
# create a new 'CSV Reader'
#image_Inc_0csv_1 = CSVReader(registrationName='Image_Inc_0.csv*', FileName=filenames2)

# create a new 'Table To Points'
#tableToPoints2 = TableToPoints(registrationName='TableToPoints2', Input=image_Inc_0csv_1)
#tableToPoints2.XColumn = 'x_proj'
#tableToPoints2.YColumn = 'y_proj'
#tableToPoints2.ZColumn = 'z_proj'

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

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
    LUTColorBar.Position = [0.75, 0.127]
    LUTColorBar.Title = output["title"]
    LUTColorBar.ComponentTitle = ''
    LUTColorBar.HorizontalTitle = 1
    LUTColorBar.TitleFontFamily = 'Times'
    LUTColorBar.TitleFontSize = 40
    LUTColorBar.LabelFontFamily = 'Times'
    LUTColorBar.LabelFontSize = 30
    LUTColorBar.ScalarBarThickness = 36
    LUTColorBar.ScalarBarLength = 0.675
    LUTColorBar.AutomaticLabelFormat = 0
    LUTColorBar.LabelFormat = output["LabelFormat"]
    LUTColorBar.UseCustomLabels = 1
    LUTColorBar.CustomLabels = output["Labels"]
    LUTColorBar.AddRangeLabels = 0
    LUTColorBar.RangeLabelFormat = output["LabelFormat"]
    LUTColorBar.Visibility = 0


# show data from tableToPoints1
tableToPoints1Display = Show(tableToPoints1, renderView1, 'GeometryRepresentation')
#tableToPoints2Display = Show(tableToPoints2, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
tableToPoints1Display.Representation = 'Surface'
tableToPoints1Display.ColorArrayName = ['POINTS', output_list[-1]["name"]]
tableToPoints1Display.LookupTable = LUT
tableToPoints1Display.PointSize = 8.0
# tableToPoints1Display.RenderPointsAsSpheres = 0 # BRING BACK IN MAYBE IF NEEDED
tableToPoints1Display.SelectTCoordArray = 'None'
tableToPoints1Display.SelectNormalArray = 'None'
tableToPoints1Display.SelectTangentArray = 'None'
tableToPoints1Display.OSPRayScaleArray = 'Compressive Force'
tableToPoints1Display.OSPRayScaleFunction = 'Piecewise Function'
tableToPoints1Display.Assembly = ''
tableToPoints1Display.SelectOrientationVectors = 'None'
tableToPoints1Display.ScaleFactor = 41.07055345227981
tableToPoints1Display.SelectScaleArray = 'Compressive Force'
tableToPoints1Display.GlyphType = 'Arrow'
tableToPoints1Display.GlyphTableIndexArray = 'Compressive Force'
tableToPoints1Display.GaussianRadius = 2.0535276726139906
tableToPoints1Display.SetScaleArray = ['POINTS', 'Compressive Force']
tableToPoints1Display.ScaleTransferFunction = 'Piecewise Function'
tableToPoints1Display.OpacityArray = ['POINTS', 'Compressive Force']
tableToPoints1Display.OpacityTransferFunction = 'Piecewise Function'
tableToPoints1Display.DataAxesGrid = 'Grid Axes Representation'
tableToPoints1Display.PolarAxes = 'Polar Axes Representation'
tableToPoints1Display.SelectInputVectors = [None, '']
tableToPoints1Display.WriteLog = ''

#tableToPoints2Display.Representation = 'Surface'
#tableToPoints2Display.ColorArrayName = ['POINTS', output_list[-1]["name"]]
#tableToPoints2Display.LookupTable = LUT
#tableToPoints2Display.PointSize = 8.0
#tableToPoints2Display.RenderPointsAsSpheres = 1
#tableToPoints2Display.SelectTCoordArray = 'None'
#tableToPoints2Display.SelectNormalArray = 'None'
#tableToPoints2Display.SelectTangentArray = 'None'
#tableToPoints2Display.OSPRayScaleArray = 'Compressive Force'
#tableToPoints2Display.OSPRayScaleFunction = 'Piecewise Function'
#tableToPoints2Display.Assembly = ''
#tableToPoints2Display.SelectOrientationVectors = 'None'
#tableToPoints2Display.ScaleFactor = 41.07055345227981
#tableToPoints2Display.SelectScaleArray = 'Compressive Force'
#tableToPoints2Display.GlyphType = 'Arrow'
#tableToPoints2Display.GlyphTableIndexArray = 'Compressive Force'
#tableToPoints2Display.GaussianRadius = 2.0535276726139906
#tableToPoints2Display.SetScaleArray = ['POINTS', 'Compressive Force']
#tableToPoints2Display.ScaleTransferFunction = 'Piecewise Function'
#tableToPoints2Display.OpacityArray = ['POINTS', 'Compressive Force']
#tableToPoints2Display.OpacityTransferFunction = 'Piecewise Function'
#tableToPoints2Display.DataAxesGrid = 'Grid Axes Representation'
#tableToPoints2Display.PolarAxes = 'Polar Axes Representation'
#tableToPoints2Display.SelectInputVectors = [None, '']
#tableToPoints2Display.WriteLog = ''

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
tableToPoints1Display.ScaleTransferFunction.Points = [150.0, 0.0, 0.5, 0.0, 150.03125, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
tableToPoints1Display.OpacityTransferFunction.Points = [150.0, 0.0, 0.5, 0.0, 150.03125, 1.0, 0.5, 0.0]

# show data from tableToPoints2
#tableToPoints2Display = Show(tableToPoints2, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
#tableToPoints2Display.Representation = 'Surface'
#tableToPoints2Display.ColorArrayName = ['POINTS', output_list[-1]["name"]]
#tableToPoints2Display.LookupTable = LUT
#tableToPoints2Display.SelectTCoordArray = 'None'
#tableToPoints2Display.SelectNormalArray = 'None'
#tableToPoints2Display.SelectTangentArray = 'None'
#tableToPoints2Display.OSPRayScaleArray = 'Compressive Force'
#tableToPoints2Display.OSPRayScaleFunction = 'Piecewise Function'
#tableToPoints2Display.Assembly = ''
#tableToPoints2Display.SelectOrientationVectors = 'None'
#tableToPoints2Display.ScaleFactor = 40.812430235347854
#tableToPoints2Display.SelectScaleArray = 'Compressive Force'
#tableToPoints2Display.GlyphType = 'Arrow'
#tableToPoints2Display.GlyphTableIndexArray = 'Compressive Force'
#tableToPoints2Display.GaussianRadius = 2.0406215117673927
#tableToPoints2Display.SetScaleArray = ['POINTS', 'Compressive Force']
#tableToPoints2Display.ScaleTransferFunction = 'Piecewise Function'
#tableToPoints2Display.OpacityArray = ['POINTS', 'Compressive Force']
#tableToPoints2Display.OpacityTransferFunction = 'Piecewise Function'
#tableToPoints2Display.DataAxesGrid = 'Grid Axes Representation'
#tableToPoints2Display.PolarAxes = 'Polar Axes Representation'
#tableToPoints2Display.SelectInputVectors = [None, '']
#tableToPoints2Display.WriteLog = ''

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
#tableToPoints2Display.ScaleTransferFunction.Points = [150.0, 0.0, 0.5, 0.0, 150.03125, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
#tableToPoints2Display.OpacityTransferFunction.Points = [150.0, 0.0, 0.5, 0.0, 150.03125, 1.0, 0.5, 0.0]

# setup the color legend parameters for each legend in this view

# set color bar visibility
LUTColorBar.Visibility = 1

# show color legend
tableToPoints1Display.SetScalarBarVisibility(renderView1, True)
# show color legend
#tableToPoints2Display.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup animation scene, tracks and keyframes
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get the time-keeper
timeKeeper1 = GetTimeKeeper()

# initialize the timekeeper

# get time animation track
timeAnimationCue1 = GetTimeTrack()

# initialize the animation track

# get animation scene
animationScene1 = GetAnimationScene()

# initialize the animation scene
animationScene1.ViewModules = renderView1
animationScene1.Cues = timeAnimationCue1
animationScene1.AnimationTime = 15.0
animationScene1.EndTime = 15.0
animationScene1.PlayMode = 'Snap To TimeSteps'

# initialize the animation scene


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
