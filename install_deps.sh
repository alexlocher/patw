#!/bin/bash
# #####################################

SCRIPTPATH=$( cd $(dirname $0) ; pwd -P )
INSTALL_DIR=${SCRIPTPATH}/thirdparty
DOWNLOAD_DIR=${SCRIPTPATH}/thirdparty_src
mkdir $INSTALL_DIR
mkdir $DOWNLOAD_DIR

# install the dependencies:
# ###########################################

# Gflags
# ------
cd $DOWNLOAD_DIR
git clone https://github.com/gflags/gflags.git
cd gflags
git checkout release
mkdir build && cd build
cmake \
	-DBUILD_SHARED_LIBS=ON \
	-DBUILD_STATIC_LIBS=ON \
	-DGFLAGS_NAMESPACE=google \
	-DCMAKE_INSTALL_PREFIX=$INSTALL_DIR ..
make -j8 install

# GLOG
# ----
cd $DOWNLOAD_DIR
git clone https://github.com/google/glog.git
cd glog
git checkout tags/v0.3.3
./configure  --prefix=$INSTALL_DIR --with-gflags=$INSTALL_DIR
make -j8 install


# Ceres
# -----
cd $DOWNLOAD_DIR
git clone https://github.com/ceres-solver/ceres-solver.git
cd ceres-solver
git checkout tags/1.11.0
mkdir build && cd build
cmake \
	-DCMAKE_INSTALL_PREFIX=$INSTALL_DIR \
	-DGFLAGS_INCLUDE_DIR_HINTS=${INSTALL_DIR}/include \
	-DGFLAGS_LIBRARY_DIR_HINTS=${INSTALL_DIR}/lib \
	-DGLOG_INCLUDE_DIR_HINTS=${INSTALL_DIR}/include \
	-DGLOG_LIBRARY_DIR_HINTS=${INSTALL_DIR}/lib \
	..
make -j8 install

# OpenImageIO
# -----------
cd $DOWNLOAD_DIR
git clone https://github.com/OpenImageIO/oiio.git
cd oiio
mkdir build && cd build
cmake \
	-DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} \
	-DEMBEDPLUGINS=OFF \
	-DOIIO_BUILD_TOOLS=OFF \
	-DOIIO_BUILD_TESTS=OFF \
	-DUSE_OPENGL=OFF \
	-DUSE_QT=OFF \
	-DUSE_PYTHON=OFF \
	-DUSE_FFMPEG=OFF \
	..
make -j4 install



# Theia
# -----

cd $DOWNLOAD_DIR
git clone https://github.com/sweeneychris/TheiaSfM.git
cd TheiaSfM
mkdir build
cd build
cmake \
	-DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} \
	-DGFLAGS_INCLUDE_DIR_HINTS=${INSTALL_DIR}/include \
	-DGFLAGS_LIBRARY_DIR_HINTS=${INSTALL_DIR}/lib \
	-DGLOG_INCLUDE_DIR_HINTS=${INSTALL_DIR}/include \
	-DGLOG_LIBRARY_DIR_HINTS=${INSTALL_DIR}/lib \
	-DCeres_DIR=${INSTALL_DIR}/share/Ceres \
	-DBUILD_UI=OFF \
	-DBUILD_TESTING=OFF \
	..
make -j4 install


# install pointcloud library
# ##########################################
pcl_name=pcl-1.8.0
cd $DOWNLOAD_DIR
wget https://github.com/PointCloudLibrary/pcl/archive/${pcl_name}.tar.gz
tar -a -xf ${pcl_name}.tar.gz
rm ${pcl_name}.tar.gz
cd pcl-${pcl_name} && mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} -DCMAKE_BUILD_TYPE=Release
make -j8 install


