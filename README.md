# Vtrak: Image Analysis Software for Vesicle Aspiration Experiments

![gif1](ex_raw.gif "title-1") ![gif2](ex_process.gif "title-2") <img src="rot2_area.png" width=360px>

---

# Process

The bulk of the mage analysis is handled by the sci-kit image libray
(skimage). The process consists of canny edge detection filters followed by
Hough transforms to isolate and extract the circular components of the
vesicle and protrusion. The protrustion length is calculated from these
circles and the relative area change produced following the approximations
used by Rawicz et. al.
```
W Rawicz, KC Olbrich, T McIntosh, D Needham, and E Evans.
Effect of chain length and unsaturation on elasticity of lipid bilayers.
Biophysical Journal. 79(1). 328-339.
2000.
```

# Usage

Input files must be multimage TIFF stacks (foo.tif). Parameters must be
specified in a valid YAML file (foo.yml) in the same directory. All data
logs and output graphs will be saved in the  current directory with the
root name `foo`.

To use the program, simply call `python vtrak.py foo.tif`. If multiple files
are to be analyzed, the local buildlisting dialog (WXpython for Windows
and Linux Dialog for Linux/OSX) will be called to return a datalist from
the current directory by supplying the option `--fs`. Multiple files will
be analyzed and graphed together. The user will be prompted to supply a
name for the multi-experiment graph, and it will be saved in the current
directory. If no name is supplied, the graph will not be saved. 

# Dependencies:

* numpy
* skimage
* WXpython (Windows)
* pydialog (Linux/OSX)


