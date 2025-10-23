# Elastix Demo for Mechanical-Testing Scan Registration

Repository by MathieuSm

## 🎯 Project Overview

This repository demonstrates how to use the open-source toolbox elastix for image registration in the context of mechanical testing. Specifically, we show how to take a scanned image of a specimen before mechanical testing and a scan after failure, and align them via registration. This helps visualise deformation, damage or failure-induced changes in material or structure.

While elastix is widely used in medical imaging, this demo repurposes it for materials/mechanical-testing workflows. It’s intended as a friendly, practical tutorial rather than an exhaustive research pipeline.

## 📂 What’s in this repository

Here’s a breakdown of the key files/folders you’ll find:

**Elastix**/
├── **ElastixDemo**/
│ ├── **Sample**/ – Folder containing your scans and results of RegisterSample.py
│ │ └── *TransformParameters.0.txt* – Resulting parameter map
│ ├── *Elastix_Demo.ipynb* – General demo of Elastix
│ ├── *RegisterSample.py* – Elastix registration applied to experimental case
│ ├── *TransformParameters.0.txt* – First parameter map resulting from Elastix_Demo.ipynb
│ └── *TransformParameters.1.txt* – Second parameter map resulting from Elastix_Demo.ipynb
├── **STBio_Presentation**/ – Presentation done during the STBio meeting
├── *.gitignore* – Ignoring large files
├── *LICENSE* – MIT license file
└── *README.md* – this file

## 🛠 Getting started / Installation

By running the *Elastix_Demo.ipynb* file, the first cell will install the required packages

    # Install required packages
    !pip install SimpleITK-SimpleElastix pyvista xmltodict


## 🧭 What you’ll see / Typical results

The *RegisterSample.py* will show the following:
- 3D images of “before” and “after” scans, showing how the specimen failed
- A simple visualisation of the difference (e.g., thresholded scans) and the Dice coefficient
- A transformation / deformation field (especially from the Bspline step) that quantifies how each voxel moved
- A change of volume field that quantifies local expansion or compression
- A kind of amount of shear + rotations field

## 💡 Tips and notes

Pay careful attention to voxel spacing, orientation, and image origin in your scans: mismatches can cause registration to fail or produce poor results.

Always start with the simpler rigid/affine step before jumping into non-linear warping.

Tweak the parameter files (especially grid spacing, number of iterations) for your material domain – what works for medical imaging may need adaptation for mechanical testing and the opposite as well.

Use visual overlays (e.g., colour blend, checkerboards) to qualitatively assess alignment before trusting quantitative outputs.

If you have very large volumes / high resolution, consider downsampling or using fewer resolution levels in the parameter file to reduce computation time.

## 🧾 License & Attribution

This repository is provided under the MIT license.
The registration engine elastix itself is open source (Apache-2.0) and should be cited if you publish results obtained using it. 

If you use this demo in your work, you’re welcome (but not required) to credit the author

## 📚 Reference

An example of elastix registration and comparison with homogenized finite element analysis can be found here:

M. Simon, et al.
“Homogenized finite element analysis of distal tibia sections: Achievements and limitations.”
Bone Reports, 21, 101752, 2024.
https://doi.org/10.1016/j.bonr.2024.101752