# Progressive 3D Reconstruction All the Way

**Author:** [Alex Locher](http://www.vision.ee.ethz.ch/~alocher)

This is an implementation of the method presented in the paper 'Progressive 
3D Reconstruction All the Way'. The code bridges the gap between progressive
sparse 3D modelling and progressive dense 3D modelling by adapting an existing
dense 3D model to an updated sparse scene. This allows to reuse the parts of the
dense model which did not changed while still incorporate and account for all the 
changes in internal and external camera calibration. 

#### Input and Output of the algorithm
<img src="http://www.vision.ee.ethz.ch/~alocher/pdf/3dv16/overview.png" 
alt="PATW - Input / Output" width="50%" border="10" />


## Related Publication
[1] Alex Locher, Michal Havlena and Luc Van Gool. Progressive 3D Modeling All the Way. *3DV 2016*. [pdf](http://people.ee.ethz.ch/~alocher/pdf/locher_3dv16_progressive_all_the_way.pdf)

[2] Alex Locher, Michal Perdoch and Luc Van Gool. Progressive prioritized multi-view stereo. *CVPR 2016*. [pdf](http://people.ee.ethz.ch/~alocher/pdf/locher_cvpr16_progressive_prioritized_mvs.pdf)


# Licence
PATW is realeased under the [GPLv3 Licence](https://www.gnu.org/licenses/gpl-3.0.txt). 


If you use the algortihm in your academic work, please cite:

    @article{locher163dv,
      title={Progressive 3D Modeling All the Way},
      author={Locher, Alex and Havlena, Michal and Van Gool, Luc},
      journal={3DV},
      year={2016}
     }

# Disclaimer
The software is research code and should not be taken as an example for good coding style ;-)


# Build and Install
The software can be build with the cmake framework and has some dependencies 

 - [TheiaSfM](http://www.theia-sfm.org/)
 - [Pointcloud Library](http://pointclouds.org/)

If you don't have the dependencies installed yet, you can try to execute the 
`install_deps.sh` scritp included in the repository. 


```
git clone https://github.com/alexlocher/patw
cd patw
git submodule update
# optinally install dependencies
./install_deps.sh
mkdir build && cd build
cmake ..
make -j4
```

# Usage
The binary *patw* takes two nvm files [nvm-file](http://ccwu.me/vsfm/doc.html#nvm) for the sparse models a 
ply file for the dense model before as input and the output can be specified:

```
./patw --nvm1=<nvm-file-before> --nvm2=<nvm-file-after> --patches=<dense-patches-before> --outdir=/tmp/patw
```

