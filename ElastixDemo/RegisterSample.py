#%% #!/usr/bin/env python3

Description = """

"""

__author__ = ['Mathieu Simon']
__date_created__ = '03-10-2025'
__license__ = 'GPL'
__version__ = '1.0'

#%% Imports

import vtk
import struct
import numpy as np
import pyvista as pv
import SimpleITK as sitk
from pathlib import Path
from vtk.util.numpy_support import vtk_to_numpy # type: ignore

#%% Define functions

def Get_AIM_Ints(File):

    """
    Function by Glen L. Niebur, University of Notre Dame (2010)
    reads the integer data of an AIM file to find its header length
    """

    nheaderints = 32
    File.seek(0)
    binints = File.read(nheaderints * 4)
    header_int = struct.unpack("=32i", binints)

    return header_int

def ReadAIM(File):

    """
    Reads an AIM file and provides
    the corresponding itk image additional
    data (i.e. spacing, calibration data, 
    and header)

    from Denis hFE pipeline
    """

    # read header
    with open(File, 'rb') as f:
        AIM_Ints = Get_AIM_Ints(f)
        # check AIM version
        if int(AIM_Ints[5]) == 16:
            if int(AIM_Ints[10]) == 131074:
                Format = "short"
            elif int(AIM_Ints[10]) == 65537:
                Format = "char"
            elif int(AIM_Ints[10]) == 1376257:
                Format = "bin compressed"
                print("     -> format " + Format + " not supported! Exiting!")
                exit(1)
            else:
                Format = "unknown"
                exit(1)
            Header = f.read(AIM_Ints[2])
            Header_Length = len(Header) + 160
            Extents = (0, AIM_Ints[14] - 1, 0, AIM_Ints[15] - 1, 0, AIM_Ints[16] - 1)
        else:
            print("     -> version 030")
            if int(AIM_Ints[17]) == 131074:
                Format = "short"
                print("     -> format " + Format)
            elif int(AIM_Ints[17]) == 65537:
                Format = "char"
                print("     -> format " + Format)
            elif int(AIM_Ints[17]) == 1376257:
                Format = "bin compressed"
                print("     -> format " + Format + " not supported! Exiting!")
                exit(1)
            else:
                Format = "unknown"
                print("     -> format " + Format + "! Exiting!")
                exit(1)
            Header = f.read(AIM_Ints[8])
            Header_Length = len(Header) + 280
            Extents = (0, AIM_Ints[24] - 1, 0, AIM_Ints[26] - 1, 0, AIM_Ints[28] - 1)

    # collect data from header if existing
    # header = re.sub('(?i) +', ' ', header)
    Header = Header.split('\n'.encode())
    Header.pop(0)
    Header.pop(0)
    Header.pop(0)
    Header.pop(0)
    Scaling = None
    Slope = None
    Intercept = None
    IPLPostScanScaling = 1
    for Line in Header:
        if Line.find(b"Orig-ISQ-Dim-p") > -1:
            origdimp = ([int(s) for s in Line.split(b" ") if s.isdigit()])

        if Line.find("Orig-ISQ-Dim-um".encode()) > -1:
            origdimum = ([int(s) for s in Line.split(b" ") if s.isdigit()])

        if Line.find("Orig-GOBJ-Dim-p".encode()) > -1:
            origdimp = ([int(s) for s in Line.split(b" ") if s.isdigit()])

        if Line.find("Orig-GOBJ-Dim-um".encode()) > -1:
            origdimum = ([int(s) for s in Line.split(b" ") if s.isdigit()])

        if Line.find("Scaled by factor".encode()) > -1:
            Scaling = float(Line.split(" ".encode())[-1])
        if Line.find("Density: intercept".encode()) > -1:
            Intercept = float(Line.split(" ".encode())[-1])
        if Line.find("Density: slope".encode()) > -1:
            Slope = float(Line.split(" ".encode())[-1])
        # if el_size scale was applied, the above still takes the original voxel size. This function works
        # only if an isotropic scaling was applied!!!
        if Line.find('downscaled'.encode()) > -1:
            pass
        elif Line.find("scale".encode()) > -1:
            IPLPostScanScaling = float(Line.split(" ".encode())[-1])
    # Spacing is calculated from Original Dimensions. This is wrong, when the images were coarsened and
    # the voxel size is not anymore corresponding to the original scanning resolution!

    try:
        Spacing = IPLPostScanScaling * (
            np.around(np.asarray(origdimum) / np.asarray(origdimp) / 1000, 5)
        )
    except:
        pass

    # read AIM with vtk
    Reader = vtk.vtkImageReader2()
    Reader.SetFileName(File)
    Reader.SetDataByteOrderToLittleEndian()
    Reader.SetFileDimensionality(3)
    Reader.SetDataExtent(Extents)
    Reader.SetHeaderSize(Header_Length)
    if Format == "short":
        Reader.SetDataScalarTypeToShort()
    elif Format == "char":
        Reader.SetDataScalarTypeToChar()
    Reader.SetDataSpacing(Spacing)
    Reader.Update()
    VTK_Image = Reader.GetOutput()


    # Convert VTK to numpy
    Data = VTK_Image.GetPointData().GetScalars()
    Dimension = VTK_Image.GetDimensions()
    Numpy_Image = vtk_to_numpy(Data)
    Numpy_Image = Numpy_Image.reshape(Dimension[2], Dimension[1], Dimension[0])

    # Y symmetry (thanks Michi for notifying this!)
    Numpy_Image = Numpy_Image[:,::-1,:]

    AdditionalData = {'Spacing':Spacing,
                      'Scaling':Scaling,
                      'Slope':Slope,
                      'Intercept':Intercept,
                      'Header':Header}

    return Numpy_Image.T[::-1,::-1,::-1], AdditionalData

def Resample(Image, Factor=None, Size=[None], Spacing=[None], Order=1):

    """
    Resample a SimpleITK image by either a given factor, a new size, or
    a new voxel spacing. Order stands for interpolation order e.g.
    Order = 1: Linear interpolation 
    """

    if type(Image) == np.array:
        Image = sitk.GetImageFromArray(Image)

    Dimension = Image.GetDimension()
    OriginalSpacing = np.array(Image.GetSpacing())
    OriginalSize = np.array(Image.GetSize())
    PhysicalSize = OriginalSize * OriginalSpacing

    Origin = Image.GetOrigin()
    Direction = Image.GetDirection()

    if Factor:
        NewSize = [round(Size/Factor) for Size in Image.GetSize()] 
        NewSpacing = [PSize/(Size-1) for Size,PSize in zip(NewSize, PhysicalSize)]
    
    elif Size[0]:
        NewSize = Size
        NewSpacing = [PSize/Size for Size, PSize in zip(NewSize, PhysicalSize)]
    
    elif Spacing[0]:
        NewSpacing = Spacing
        NewSize = [np.floor(Size/Spacing).astype('int') + 1 for Size,Spacing in zip(PhysicalSize, NewSpacing)]

    NewArray = np.zeros(NewSize[::-1],'int')
    NewImage = sitk.GetImageFromArray(NewArray)
    NewImage.SetOrigin(Origin - OriginalSpacing/2)
    NewImage.SetDirection(Direction)
    NewImage.SetSpacing(NewSpacing)
  
    Transform = sitk.TranslationTransform(Dimension)
    Resampled = sitk.Resample(Image, NewImage, Transform, Order)
    
    return Resampled

def PlotImage(Image, Spacing):
    """Visualize a 3D volume using PyVista.

    Parameters
    ----------
    Image : np.ndarray or SimpleITK.Image
    Spacing : sequence(float)
    Sampling : int
        Subsampling factor to speed rendering.
    """
    if type(Image) == sitk.SimpleITK.Image:
        Array = sitk.GetArrayFromImage(Image).T
    else:
        Array = Image

    Shape = np.round(np.array(Array.shape) * Spacing,1)

    Args = dict(font_family='times', 
                font_size=30,
                location='outer',
                show_xlabels=False,
                show_ylabels=False,
                show_zlabels=False,
                all_edges=True,
                fmt='%i',
                xtitle=f'{Shape[0]} mm',
                ytitle=f'{Shape[1]} mm',
                ztitle=f'{Shape[2]} mm',
                use_3d_text=False
                )

    # Plot using pyvista
    Plot = pv.Plotter(off_screen=True)
    Actors = Plot.add_volume(Array,cmap='bone',show_scalar_bar=False, opacity='sigmoid_7')
    Actors.prop.interpolation_type = 'linear'
    Plot.camera_position = 'xz'
    Plot.camera.roll = 0
    Plot.camera.elevation = 30
    Plot.camera.azimuth = -60
    Plot.camera.zoom(1)
    Plot.show_bounds(**Args)
    Plot.add_axes()
    Plot.show()


    return

def PlotDice(Image, FixedMask, MovingMask, Spacing):
    """Visualize a 3D volume using PyVista.

    Parameters
    ----------
    Image : np.ndarray or SimpleITK.Image
    Spacing : sequence(float)
    Sampling : int
        Subsampling factor to speed rendering.
    """

    if type(Image) == sitk.SimpleITK.Image:
        Image = sitk.GetArrayFromImage(Image).T
    else:
        Image = Image
    
    if type(FixedMask) == sitk.SimpleITK.Image:
        Array1 = sitk.GetArrayFromImage(FixedMask).T
    else:
        Array1 = FixedMask

    if type(MovingMask) == sitk.SimpleITK.Image:
        Array2 = sitk.GetArrayFromImage(MovingMask).T
    else:
        Array2 = MovingMask

    Shape = np.round(np.array(Array1.shape) * Spacing,1)
    Array1 = Array1 > 0
    Array2 = Array2 > 0
    Array = (Array1 & Array2) * 1
    Array1[Array == 1] = 0
    Array2[Array == 1] = 0
    Image[Array != 1] = 0

    Args = dict(font_family='times', 
                font_size=30,
                location='outer',
                show_xlabels=False,
                show_ylabels=False,
                show_zlabels=False,
                all_edges=True,
                fmt='%i',
                xtitle=f'{Shape[0]} mm',
                ytitle=f'{Shape[1]} mm',
                ztitle=f'{Shape[2]} mm',
                use_3d_text=False
                )

    # Plot using pyvista
    Plot = pv.Plotter(off_screen=True)
    Actors = Plot.add_volume(Array1*1.,cmap='kr',show_scalar_bar=False, opacity=[0.0,1.0])
    Actors = Plot.add_volume(Array2*1.,cmap='kb',show_scalar_bar=False, opacity=[0.0,1.0])
    Actors = Plot.add_volume(Array,cmap='bone',show_scalar_bar=False, opacity=[0.0,1.0])
    Actors.prop.interpolation_type = 'linear'
    Plot.camera_position = 'xz'
    Plot.camera.roll = 0
    Plot.camera.elevation = 30
    Plot.camera.azimuth = -60
    Plot.camera.zoom(1)
    Plot.show_bounds(**Args)
    Plot.add_axes()
    Plot.show()


    return

def PlotVectors(Vectors, Spacing):
    """Visualize a 3D vector field as arrows using PyVista.

    Parameters
    ----------
    Vectors : np.ndarray or SimpleITK.Image
        Expected shape (X, Y, Z, 3) or transposed variants handled by the function.
    Spacing : sequence(float)
    Sampling : int
    """
    if type(Vectors) == sitk.SimpleITK.Image:
        Array = sitk.GetArrayFromImage(Vectors).T
    else:
        Array = Vectors

    Shape = np.round(np.array(Array.shape[:-1]) * Spacing,1)

    nX, nY, nZ = Array.shape[:-1]
    X, Y, Z = np.meshgrid(np.arange(nX), np.arange(nY), np.arange(nZ), indexing="ij")

    Points = np.column_stack((X.ravel(), Y.ravel(), Z.ravel()))
    Array = Array.reshape(-1, 3)

    PointCloud = pv.PolyData(Points, force_float=False)
    PointCloud['Vectors'] = Array
    PointCloud['Norm'] = np.linalg.norm(Array, axis=-1)
    Arrows = PointCloud.glyph(orient='Vectors', scale=True, factor=1.0)

    Args = dict(font_family='times', 
                font_size=30,
                location='outer',
                show_xlabels=False,
                show_ylabels=False,
                show_zlabels=False,
                all_edges=True,
                fmt='%i',
                xtitle=f'{Shape[0]} mm',
                ytitle=f'{Shape[1]} mm',
                ztitle=f'{Shape[2]} mm',
                use_3d_text=False
                )

    sArgs = dict(font_family='times', 
                 width=0.05,
                 height=0.75,
                 vertical=True,
                 position_x=0.88,
                 position_y=0.125,
                 title_font_size=30,
                 label_font_size=20,
                 title='u (mm)',
                 fmt="%.3f")
    
    # Plot using pyvista
    Plot = pv.Plotter(off_screen=True)
    Actors = Plot.add_mesh(Arrows, scalars='Norm', cmap='jet', show_scalar_bar=True, scalar_bar_args=sArgs)
    Plot.camera_position = 'xz'
    Plot.camera.roll = 0
    Plot.camera.elevation = 30
    Plot.camera.azimuth = -60
    Plot.camera.zoom(0.8)
    Plot.show_bounds(**Args)
    Plot.add_axes()
    Plot.show()

    return

def PlotJ(Image, J, Spacing):
    """Visualize a 3D volume using PyVista.

    Parameters
    ----------
    Image : np.ndarray or SimpleITK.Image
    Spacing : sequence(float)
    Sampling : int
        Subsampling factor to speed rendering.
    """
    if type(Image) == sitk.SimpleITK.Image:
        Array = sitk.GetArrayFromImage(Image).T
    else:
        Array = Image
    if type(J) == sitk.SimpleITK.Image:
        J = sitk.GetArrayFromImage(J).T
    else:
        J = J


    Shape = np.round(np.array(Array.shape) * Spacing,1)
    J[Array == 0] = 1.0

    Args = dict(font_family='times', 
                font_size=30,
                location='outer',
                show_xlabels=False,
                show_ylabels=False,
                show_zlabels=False,
                all_edges=True,
                fmt='%i',
                xtitle=f'{Shape[0]} mm',
                ytitle=f'{Shape[1]} mm',
                ztitle=f'{Shape[2]} mm',
                use_3d_text=False
                )
    
    sArgs = dict(font_family='times', 
                 width=0.05,
                 height=0.75,
                 vertical=True,
                 position_x=0.9,
                 position_y=0.125,
                 title_font_size=30,
                 label_font_size=20,
                 title='J (-)')
    

    # Plot using pyvista
    Plot = pv.Plotter(off_screen=True)
    Actors = Plot.add_volume(Array,cmap='bone', show_scalar_bar=False, opacity='sigmoid_7')
    Actors = Plot.add_volume(J,cmap='jet', show_scalar_bar=True, scalar_bar_args=sArgs, opacity=[1.,0.,0.,0.,1.], clim=[0.75,1.25])
    Actors.prop.interpolation_type = 'linear'
    Plot.camera_position = 'xz'
    Plot.camera.roll = 0
    Plot.camera.elevation = 30
    Plot.camera.azimuth = -60
    Plot.camera.zoom(0.8)
    Plot.show_bounds(**Args)
    Plot.add_axes()
    Plot.show()

    return

def PlotFTilde(Image, FTilde, Spacing):
    """Visualize a 3D volume using PyVista.

    Parameters
    ----------
    Image : np.ndarray or SimpleITK.Image
    Spacing : sequence(float)
    Sampling : int
        Subsampling factor to speed rendering.
    """
    if type(Image) == sitk.SimpleITK.Image:
        Array = sitk.GetArrayFromImage(Image).T
    else:
        Array = Image
    if type(FTilde) == sitk.SimpleITK.Image:
        FTilde = sitk.GetArrayFromImage(FTilde).T
    else:
        FTilde = FTilde

    Shape = np.round(np.array(Array.shape) * Spacing,1)

    Args = dict(font_family='times', 
                font_size=30,
                location='outer',
                show_xlabels=False,
                show_ylabels=False,
                show_zlabels=False,
                all_edges=True,
                fmt='%i',
                xtitle=f'{Shape[0]} mm',
                ytitle=f'{Shape[1]} mm',
                ztitle=f'{Shape[2]} mm',
                use_3d_text=False
                )
    
    sArgs = dict(font_family='times', 
                 width=0.05,
                 height=0.75,
                 vertical=True,
                 position_x=0.9,
                 position_y=0.125,
                 title_font_size=30,
                 label_font_size=20,
                 title='||F|| (-)')
    

    # Plot using pyvista
    Min, Max = FTilde.min(), FTilde.max()
    # Max = Max - (Max-Min)*0.2
    FTilde[Array == 0] = Min
    Plot = pv.Plotter(off_screen=True)
    Actors = Plot.add_volume(Array,cmap='bone', show_scalar_bar=False, opacity='sigmoid_7')
    Actors = Plot.add_volume(FTilde, cmap='seismic', show_scalar_bar=True, scalar_bar_args=sArgs, opacity=[0.,0.,1.], clim=[Min,Max])
    Actors.prop.interpolation_type = 'linear'
    Plot.camera_position = 'xz'
    Plot.camera.roll = 0
    Plot.camera.elevation = 30
    Plot.camera.azimuth = -60
    Plot.camera.zoom(0.8)
    Plot.show_bounds(**Args)
    Plot.add_axes()
    Plot.show()

    return

#%% Define image files

# List files
FilePath = Path(__file__).parent / 'Sample'
Files = sorted([F for F in FilePath.iterdir() if F.name.endswith('.AIM')])

# Defines images
ReferenceFile = Files[0]
DeformedFile = Files[1]

#%% Read uCT reference sample
Reference, AdditionalData = ReadAIM(ReferenceFile)
Spacing = AdditionalData['Spacing']
PlotImage(Reference[300:], Spacing)

#%% Read uCT deformed sample
Deformed, AdditionalData = ReadAIM(DeformedFile)
Spacing = AdditionalData['Spacing']
PlotImage(Deformed[300:], Spacing)

#%% Rigid registration

# Cast our numpy arrays to SimpleITK images
Fixed = sitk.GetImageFromArray(Reference.T)
Moving = sitk.GetImageFromArray(Deformed.T)
Fixed.SetSpacing(Spacing)
Moving.SetSpacing(Spacing)

# Define parameter map
ParameterMap = sitk.GetDefaultParameterMap('rigid')
ParameterMap['SP_alpha'] = ['0.6']
ParameterMap['SP_A'] = ['100']
ParameterMap['MaximumNumberOfIterations'] = ['256']
ParameterMap['NewSamplesEveryIteration'] = ['true']

# Perform rigid registration
Elastix = sitk.ElastixImageFilter()
Elastix.SetFixedImage(Fixed)
Elastix.SetMovingImage(Moving)
Elastix.SetParameterMap(ParameterMap)
Elastix.SetOutputDirectory('Sample')
Rigid = Elastix.Execute()

#%% Measure registration quality

# Segment images
Otsu = sitk.OtsuMultipleThresholdsImageFilter()
Otsu.SetNumberOfThresholds(2)
ReferenceMask = Otsu.Execute(Fixed)
RigidMask = Otsu.Execute(Rigid)

# Measure goodness of registration
Measure = sitk.LabelOverlapMeasuresImageFilter()
Measure.Execute(ReferenceMask, RigidMask)
print(f'Dice coefficient {Measure.GetDiceCoefficient():.3f}')
Factor = 5
rReferenceMask = Resample(ReferenceMask, Factor=Factor)
rRigidMask = Resample(RigidMask, Factor=Factor)
rFixed = Resample(Fixed, Factor=Factor)
PlotDice(rFixed, rReferenceMask, rRigidMask, rRigidMask.GetSpacing())

#%% B-Spline registration

# Register reference image with displaced image
ParameterMap = sitk.GetDefaultParameterMap('bspline')
ParameterMap['SP_alpha'] = ['0.6']
ParameterMap['SP_A'] = ['100']
ParameterMap['MaximumNumberOfIterations'] = ['2000']
ParameterMap['NewSamplesEveryIteration'] = ['true']

Schedule = np.repeat([32, 16, 8, 4, 2], Fixed.GetDimension())
ParameterMap['FixedImagePyramidSchedule'] = [str(S) for S in Schedule]
ParameterMap['MovingImagePyramidSchedule'] = [str(S) for S in Schedule]

# Perform bspline registration
Elastix = sitk.ElastixImageFilter()
Elastix.SetFixedImage(Fixed)
Elastix.SetMovingImage(Rigid)
Elastix.SetParameterMap(ParameterMap)
Elastix.SetOutputDirectory('Sample')
BSpline = Elastix.Execute()

#%% Measure registration quality

BSplineMask = Otsu.Execute(BSpline)
Measure.Execute(ReferenceMask, BSplineMask)
print(f'Dice coefficient {Measure.GetDiceCoefficient():.3f}')
rReferenceMask = Resample(ReferenceMask, Factor=Factor)
rBSplineMask = Resample(BSplineMask, Factor=Factor)
PlotDice(rFixed, rReferenceMask, rBSplineMask, rBSplineMask.GetSpacing())

#%% Compute deformation field and gradient

ResultParameterMap = Elastix.GetTransformParameterMap()[0]
Transformix = sitk.TransformixImageFilter()
Transformix.SetTransformParameterMap(ResultParameterMap)
Transformix.ComputeDeterminantOfSpatialJacobianOff()
Transformix.ComputeDeformationFieldOn()
Transformix.ComputeSpatialJacobianOn()
Transformix.SetMovingImage(Rigid)
Transformix.SetOutputDirectory('Sample')
Transformed = Transformix.Execute()

#%% Load Deformation field
File = Path('Sample/deformationField.nii')
DefField = sitk.ReadImage(File)
rDefField = Resample(DefField, Factor=33)
rSpacing = rDefField.GetSpacing()
rDefField = sitk.GetArrayFromImage(rDefField)
rDefField = rDefField.transpose((2,1,0,3))
PlotVectors(rDefField, rSpacing)

#%% Load Jacobian
File = Path('Sample/fullSpatialJacobian.nii')
Jacobian = sitk.ReadImage(File)
Jacobian = sitk.GetArrayFromImage(Jacobian)
Jacobian = Jacobian.transpose((2,1,0,3))

#%% Build deformation gradient field
F_Field = np.zeros(Jacobian.shape[:-1] + (3,3))
for i in range(3):
    for j in range(3):
        F_Field[...,i,j] = Jacobian[...,-(i*3+j)-1]

#%% Subsample array
rF_Field = F_Field[::Factor,::Factor,::Factor]
rSpacing = rReferenceMask.GetSpacing()
rReferenceMask = sitk.GetArrayFromImage(ReferenceMask).T
rReferenceMask = rReferenceMask[::Factor,::Factor,::Factor]

#%% Decompose F
J = np.linalg.det(rF_Field)
F_tilde = np.reshape(J, J.shape + (1,1)) ** (-1 / 3) * rF_Field
NormF_tilde = np.linalg.norm(F_tilde, axis=(-1, -2))

# Resample in higher resolution
rJ = sitk.GetImageFromArray(J.T)
rJ.SetSpacing(rSpacing)
rNormF_tilde = sitk.GetImageFromArray(NormF_tilde.T)
rNormF_tilde.SetSpacing(rSpacing)
rReference = sitk.GetImageFromArray(Reference.T)
rReference.SetSpacing(rSpacing)

Size = [round(s/2) for s in rReference.GetSize()]
rJ = Resample(rJ, Size=Size)
rNormF_tilde = Resample(rNormF_tilde, Size=Size)
rReference = Resample(rReference, Size=Size)

#%% Plot change of volume
PlotJ(rReference, rJ, rJ.GetSpacing())

#%% Plot shear and rotations
PlotFTilde(rReference, rNormF_tilde, rNormF_tilde.GetSpacing())


# %%
