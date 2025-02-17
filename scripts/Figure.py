from Compute import *
import cmaps

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

SmoothingIterations = 1

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

def figure(file, method, size=500):

    # Create a new 'Render View'
    renderView = CreateView('RenderView')
    renderView.ViewSize = [size, size]
    renderView.InteractionMode = '2D'
    renderView.AxesGrid = 'GridAxes3DActor'
    renderView.OrientationAxesVisibility = 0
    SetActiveView(renderView)


    # ----------------------------------------------------------------
    # read input and output
    # ----------------------------------------------------------------

    input = CSVReader(registrationName='input', FileName=[inputpath + file + ".csv"])
    input.HaveHeaders = header(file)

    output = CSVReader(registrationName='output', FileName=[savedatapath + "{}_{}.csv".format(file, method)])


    # -----------------------------------------------------
    # setup point cloud
    # ----------------------------------------------------------------

    # create a new 'Table To Points'
    tableToPoints = TableToPoints(registrationName='TableToPoints1', Input=output)
    tableToPoints.XColumn = 'Component_0'
    tableToPoints.YColumn = 'Component_1'
    tableToPoints.a2DPoints = 1

    # show data from tableToPoints1
    tableToPointsDisplay = Show(tableToPoints, renderView, 'GeometryRepresentation')

    # trace defaults for the display properties.
    tableToPointsDisplay.Representation = 'Surface'
    tableToPointsDisplay.PointSize = 6.0
    tableToPointsDisplay.RenderPointsAsSpheres = 1
    tableToPointsDisplay.Ambient = 1.0
    tableToPointsDisplay.Diffuse = 0.0

    if file == "3Clusters":
        clusterIdTF2D = GetTransferFunction2D('ClusterId')
        clusterIdLUT = GetColorTransferFunction('ClusterId')
        clusterIdLUT.TransferFunction2D = clusterIdTF2D
        clusterIdLUT.RGBPoints = cmaps.cmap_3blobs
        clusterIdLUT.ScalarRangeInitialized = 1.0
        tableToPointsDisplay.ColorArrayName = ['POINTS', 'ClusterId']
        tableToPointsDisplay.LookupTable = clusterIdLUT
        tableToPointsDisplay.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.0, 1.0, 0.5, 0.0]
        tableToPointsDisplay.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.0, 1.0, 0.5, 0.0]
    else:
        tableToPointsDisplay.AmbientColor = [0.0, 0.0, 0.0]
        tableToPointsDisplay.ColorArrayName = [None, '']


    # ----------------------------------------------------------------
    # setup generators
    # ----------------------------------------------------------------

    if file not in ["3Clusters"]:
        # create a new 'TTK RipsPersistenceGenerators'
        tTKRipsPersistenceGenerators = TTKRipsPersistenceGenerators(registrationName='TTKRipsPersistenceGenerators1',
                                                                    Inputsource=input,
                                                                    Input3Drepresentation=tableToPoints)
        setColumns(file, tTKRipsPersistenceGenerators)
        tTKRipsPersistenceGenerators.Simplexmaximumdiameter = 1000.0

        if file == "K5":
            plot_K5(tTKRipsPersistenceGenerators, renderView)
        else:
            # create a new 'Threshold'
            threshold = Threshold(registrationName='Threshold1', Input=tTKRipsPersistenceGenerators)
            threshold.Scalars = ['CELLS', 'ClassPersistence']
            threshold.LowerThreshold = threshold_v(file)
            threshold.UpperThreshold = 1000.0

            # create a new 'TTK GeometrySmoother'
            tTKGeometrySmoother = TTKGeometrySmoother(registrationName='TTKGeometrySmoother1', Input=threshold)
            tTKGeometrySmoother.IterationNumber = SmoothingIterations

            # show data from threshold1
            thresholdDisplay = Show(tTKGeometrySmoother, renderView, 'UnstructuredGridRepresentation')

            # get 2D transfer function for 'ClassIdentifier'
            classIdentifierTF2D = GetTransferFunction2D('ClassIdentifier')

            # get color transfer function/color map for 'ClassIdentifier'
            classIdentifierLUT = GetColorTransferFunction('ClassIdentifier')
            classIdentifierLUT.TransferFunction2D = classIdentifierTF2D
            if file == "K4":
                classIdentifierLUT.RGBPoints = cmaps.cmap_K4
            elif file == "Twist":
                classIdentifierLUT.RGBPoints = cmaps.roll_cmap_cyclic(45)
            else:
                classIdentifierLUT.RGBPoints = cmaps.cmap_cyclic
            classIdentifierLUT.ColorSpace = 'HSV'
            classIdentifierLUT.ScalarRangeInitialized = 1

            # get opacity transfer function/opacity map for 'ClassIdentifier'
            classIdentifierPWF = GetOpacityTransferFunction('ClassIdentifier')
            if file == "K4":
                classIdentifierPWF.Points = [28.0, 0.0, 0.5, 0.0, 32.0, 1.0, 0.5, 0.0]
            if file == "K5":
                classIdentifierPWF.Points = [44.0, 0.0, 0.5, 0.0, 56.0, 1.0, 0.5, 0.0]
            else:
                classIdentifierPWF.Points = [0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]
            classIdentifierPWF.ScalarRangeInitialized = 1

            # trace defaults for the display properties.
            thresholdDisplay.Representation = 'Surface'
            thresholdDisplay.Opacity = 0.4
            thresholdDisplay.LineWidth = 10.0
            if file == "K4":
                thresholdDisplay.ColorArrayName = ['CELLS', 'ClassIdentifier']
            else:
                thresholdDisplay.ColorArrayName = ['CELLS', 'GeneratorParametrization']
            thresholdDisplay.LookupTable = classIdentifierLUT

    # ----------------------------------------------------------------
    # print image
    # ----------------------------------------------------------------

    GetActiveView().ResetCamera()
    SaveScreenshot(savefigurepath + "{}_{}.png".format(file, method), TransparentBackground=True)



def plot_K5(tTKRipsPersistenceGenerators, renderView):
    # create a new 'Threshold'
    threshold4 = Threshold(registrationName='Threshold4', Input=tTKRipsPersistenceGenerators)
    threshold4.Scalars = ['CELLS', 'ClassIdentifier']
    threshold4.LowerThreshold = 47.0
    threshold4.UpperThreshold = 47.0

    # create a new 'Connectivity'
    connectivity2 = Connectivity(registrationName='Connectivity2', Input=threshold4)

    # create a new 'Threshold'
    threshold3 = Threshold(registrationName='Threshold3', Input=connectivity2)
    threshold3.Scalars = ['POINTS', 'RegionId']

    # create a new 'TTK GeometrySmoother'
    tTKGeometrySmoother1 = TTKGeometrySmoother(registrationName='TTKGeometrySmoother1', Input=threshold3)
    tTKGeometrySmoother1.IterationNumber = 20
    tTKGeometrySmoother1.InputMaskField = ['POINTS', 'RegionId']

    # create a new 'Threshold'
    threshold2 = Threshold(registrationName='Threshold2', Input=tTKRipsPersistenceGenerators)
    threshold2.Scalars = ['CELLS', 'ClassIdentifier']
    threshold2.LowerThreshold = 44.0
    threshold2.UpperThreshold = 44.0

    # create a new 'Threshold'
    threshold9 = Threshold(registrationName='Threshold9', Input=tTKRipsPersistenceGenerators)
    threshold9.Scalars = ['CELLS', 'ClassIdentifier']
    threshold9.LowerThreshold = 54.0
    threshold9.UpperThreshold = 54.0

    # create a new 'Connectivity'
    connectivity5 = Connectivity(registrationName='Connectivity5', Input=threshold9)

    # create a new 'Threshold'
    threshold10 = Threshold(registrationName='Threshold10', Input=connectivity5)
    threshold10.Scalars = ['POINTS', 'RegionId']

    # create a new 'Transform'
    transform4 = Transform(registrationName='Transform4', Input=threshold10)
    transform4.Transform = 'Transform'

    # init the 'Transform' selected for 'Transform'
    transform4.Transform.Translate = [0.03, 0.03, 0.0]

    # create a new 'TTK GeometrySmoother'
    tTKGeometrySmoother5 = TTKGeometrySmoother(registrationName='TTKGeometrySmoother5', Input=transform4)
    tTKGeometrySmoother5.IterationNumber = 20
    tTKGeometrySmoother5.InputMaskField = ['POINTS', 'RegionId']

    # create a new 'Threshold'
    threshold5 = Threshold(registrationName='Threshold5', Input=tTKRipsPersistenceGenerators)
    threshold5.Scalars = ['CELLS', 'ClassIdentifier']
    threshold5.LowerThreshold = 48.0
    threshold5.UpperThreshold = 48.0

    # create a new 'Threshold'
    threshold11 = Threshold(registrationName='Threshold11', Input=tTKRipsPersistenceGenerators)
    threshold11.Scalars = ['CELLS', 'ClassIdentifier']
    threshold11.LowerThreshold = 56.0
    threshold11.UpperThreshold = 56.0

    # create a new 'Connectivity'
    connectivity6 = Connectivity(registrationName='Connectivity6', Input=threshold11)

    # create a new 'Threshold'
    threshold12 = Threshold(registrationName='Threshold12', Input=connectivity6)
    threshold12.Scalars = ['POINTS', 'RegionId']

    # create a new 'TTK GeometrySmoother'
    tTKGeometrySmoother6 = TTKGeometrySmoother(registrationName='TTKGeometrySmoother6', Input=threshold12)
    tTKGeometrySmoother6.IterationNumber = 20
    tTKGeometrySmoother6.InputMaskField = ['POINTS', 'RegionId']

    # create a new 'Threshold'
    threshold7 = Threshold(registrationName='Threshold7', Input=tTKRipsPersistenceGenerators)
    threshold7.Scalars = ['CELLS', 'ClassIdentifier']
    threshold7.LowerThreshold = 53.0
    threshold7.UpperThreshold = 53.0

    # create a new 'Connectivity'
    connectivity4 = Connectivity(registrationName='Connectivity4', Input=threshold7)

    # create a new 'Threshold'
    threshold8 = Threshold(registrationName='Threshold8', Input=connectivity4)
    threshold8.Scalars = ['POINTS', 'RegionId']

    # create a new 'Transform'
    transform2 = Transform(registrationName='Transform2', Input=threshold8)
    transform2.Transform = 'Transform'

    # init the 'Transform' selected for 'Transform'
    transform2.Transform.Translate = [-0.03, -0.03, 0.0]

    # create a new 'TTK GeometrySmoother'
    tTKGeometrySmoother4 = TTKGeometrySmoother(registrationName='TTKGeometrySmoother4', Input=transform2)
    tTKGeometrySmoother4.IterationNumber = 20
    tTKGeometrySmoother4.InputMaskField = ['POINTS', 'RegionId']

    # create a new 'Connectivity'
    connectivity3 = Connectivity(registrationName='Connectivity3', Input=threshold5)

    # create a new 'Threshold'
    threshold6 = Threshold(registrationName='Threshold6', Input=connectivity3)
    threshold6.Scalars = ['POINTS', 'RegionId']

    # create a new 'Transform'
    transform3 = Transform(registrationName='Transform3', Input=threshold6)
    transform3.Transform = 'Transform'

    # init the 'Transform' selected for 'Transform'
    transform3.Transform.Translate = [0.03, -0.03, 0.]

    # create a new 'TTK GeometrySmoother'
    tTKGeometrySmoother3 = TTKGeometrySmoother(registrationName='TTKGeometrySmoother3', Input=transform3)
    tTKGeometrySmoother3.IterationNumber = 20
    tTKGeometrySmoother3.InputMaskField = ['POINTS', 'RegionId']

    # create a new 'Connectivity'
    connectivity1 = Connectivity(registrationName='Connectivity1', Input=threshold2)

    # create a new 'Threshold'
    threshold1 = Threshold(registrationName='Threshold1', Input=connectivity1)
    threshold1.Scalars = ['POINTS', 'RegionId']
    threshold1.LowerThreshold = 1.0
    threshold1.UpperThreshold = 1.0

    # create a new 'Transform'
    transform1 = Transform(registrationName='Transform1', Input=threshold1)
    transform1.Transform = 'Transform'

    # init the 'Transform' selected for 'Transform'
    transform1.Transform.Translate = [-0.03, 0.03, 0.]

    # create a new 'TTK GeometrySmoother'
    tTKGeometrySmoother2 = TTKGeometrySmoother(registrationName='TTKGeometrySmoother2', Input=transform1)
    tTKGeometrySmoother2.IterationNumber = 20
    tTKGeometrySmoother2.InputMaskField = ['CELLS', 'ClassBirth']

    # show data from tTKGeometrySmoother2
    tTKGeometrySmoother2Display = Show(tTKGeometrySmoother2, renderView, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    tTKGeometrySmoother2Display.Representation = 'Surface'
    tTKGeometrySmoother2Display.AmbientColor = [0.9764705882352941, 0.5176470588235295, 0.9607843137254902]
    tTKGeometrySmoother2Display.ColorArrayName = ['POINTS', '']
    tTKGeometrySmoother2Display.DiffuseColor = [0.9764705882352941, 0.5176470588235295, 0.9607843137254902]
    tTKGeometrySmoother2Display.Opacity = 0.8
    tTKGeometrySmoother2Display.LineWidth = 5.0
    tTKGeometrySmoother2Display.RenderLinesAsTubes = 0
    tTKGeometrySmoother2Display.Ambient = 0.2

    # show data from tTKGeometrySmoother1
    tTKGeometrySmoother1Display = Show(tTKGeometrySmoother1, renderView, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    tTKGeometrySmoother1Display.Representation = 'Surface'
    tTKGeometrySmoother1Display.AmbientColor = [0.8, 0.6666666666666666, 0.1803921568627451]
    tTKGeometrySmoother1Display.ColorArrayName = ['POINTS', '']
    tTKGeometrySmoother1Display.DiffuseColor = [0.8, 0.6666666666666666, 0.1803921568627451]
    tTKGeometrySmoother1Display.Opacity = 0.8
    tTKGeometrySmoother1Display.LineWidth = 5.0
    tTKGeometrySmoother1Display.RenderLinesAsTubes = 0
    tTKGeometrySmoother1Display.Ambient = 0.2

    # show data from tTKGeometrySmoother3
    tTKGeometrySmoother3Display = Show(tTKGeometrySmoother3, renderView, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    tTKGeometrySmoother3Display.Representation = 'Surface'
    tTKGeometrySmoother3Display.AmbientColor = [0.7137254901960784, 0.10588235294117647, 0.08627450980392157]
    tTKGeometrySmoother3Display.ColorArrayName = ['POINTS', '']
    tTKGeometrySmoother3Display.DiffuseColor = [0.7137254901960784, 0.10588235294117647, 0.08627450980392157]
    tTKGeometrySmoother3Display.Opacity = 0.8
    tTKGeometrySmoother3Display.LineWidth = 5.0
    tTKGeometrySmoother3Display.RenderLinesAsTubes = 0
    tTKGeometrySmoother3Display.Ambient = 0.2

    # show data from tTKGeometrySmoother4
    tTKGeometrySmoother4Display = Show(tTKGeometrySmoother4, renderView, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    tTKGeometrySmoother4Display.Representation = 'Surface'
    tTKGeometrySmoother4Display.AmbientColor = [0.45098039215686275, 0.5686274509803921, 0.5215686274509804]
    tTKGeometrySmoother4Display.ColorArrayName = ['POINTS', '']
    tTKGeometrySmoother4Display.DiffuseColor = [0.45098039215686275, 0.5686274509803921, 0.5215686274509804]
    tTKGeometrySmoother4Display.Opacity = 0.8
    tTKGeometrySmoother4Display.LineWidth = 5.0
    tTKGeometrySmoother4Display.RenderLinesAsTubes = 0
    tTKGeometrySmoother4Display.Ambient = 0.2

    # show data from tTKGeometrySmoother5
    tTKGeometrySmoother5Display = Show(tTKGeometrySmoother5, renderView, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    tTKGeometrySmoother5Display.Representation = 'Surface'
    tTKGeometrySmoother5Display.AmbientColor = [0.13333333333333333, 0.44313725490196076, 0.7333333333333333]
    tTKGeometrySmoother5Display.ColorArrayName = ['POINTS', '']
    tTKGeometrySmoother5Display.DiffuseColor = [0.13333333333333333, 0.44313725490196076, 0.7333333333333333]
    tTKGeometrySmoother5Display.Opacity = 0.8
    tTKGeometrySmoother5Display.LineWidth = 5.0
    tTKGeometrySmoother5Display.RenderLinesAsTubes = 0
    tTKGeometrySmoother5Display.Ambient = 0.2

    # show data from tTKGeometrySmoother6
    tTKGeometrySmoother6Display = Show(tTKGeometrySmoother6, renderView, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    tTKGeometrySmoother6Display.Representation = 'Surface'
    tTKGeometrySmoother6Display.AmbientColor = [0.5294117647058824, 0.39215686274509803, 0.9803921568627451]
    tTKGeometrySmoother6Display.ColorArrayName = ['POINTS', '']
    tTKGeometrySmoother6Display.DiffuseColor = [0.5294117647058824, 0.39215686274509803, 0.9803921568627451]
    tTKGeometrySmoother6Display.Opacity = 0.8
    tTKGeometrySmoother6Display.LineWidth = 5.0
    tTKGeometrySmoother6Display.RenderLinesAsTubes = 0
    tTKGeometrySmoother6Display.Ambient = 0.2