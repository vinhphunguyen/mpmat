MATLAB="/Applications/MATLAB_R2015a.app"
Arch=maci64
ENTRYPOINT=mexFunction
MAPFILE=$ENTRYPOINT'.map'
PREFDIR="/Users/vinhphunguyen/.matlab/R2015a"
OPTSFILE_NAME="./setEnv.sh"
. $OPTSFILE_NAME
COMPILER=$CC
. $OPTSFILE_NAME
echo "# Make settings for mlsLinearBasis2D" > mlsLinearBasis2D_mex.mki
echo "CC=$CC" >> mlsLinearBasis2D_mex.mki
echo "CFLAGS=$CFLAGS" >> mlsLinearBasis2D_mex.mki
echo "CLIBS=$CLIBS" >> mlsLinearBasis2D_mex.mki
echo "COPTIMFLAGS=$COPTIMFLAGS" >> mlsLinearBasis2D_mex.mki
echo "CDEBUGFLAGS=$CDEBUGFLAGS" >> mlsLinearBasis2D_mex.mki
echo "CXX=$CXX" >> mlsLinearBasis2D_mex.mki
echo "CXXFLAGS=$CXXFLAGS" >> mlsLinearBasis2D_mex.mki
echo "CXXLIBS=$CXXLIBS" >> mlsLinearBasis2D_mex.mki
echo "CXXOPTIMFLAGS=$CXXOPTIMFLAGS" >> mlsLinearBasis2D_mex.mki
echo "CXXDEBUGFLAGS=$CXXDEBUGFLAGS" >> mlsLinearBasis2D_mex.mki
echo "LD=$LD" >> mlsLinearBasis2D_mex.mki
echo "LDFLAGS=$LDFLAGS" >> mlsLinearBasis2D_mex.mki
echo "LDOPTIMFLAGS=$LDOPTIMFLAGS" >> mlsLinearBasis2D_mex.mki
echo "LDDEBUGFLAGS=$LDDEBUGFLAGS" >> mlsLinearBasis2D_mex.mki
echo "Arch=$Arch" >> mlsLinearBasis2D_mex.mki
echo OMPFLAGS= >> mlsLinearBasis2D_mex.mki
echo OMPLINKFLAGS= >> mlsLinearBasis2D_mex.mki
echo "EMC_COMPILER=Xcode with Clang" >> mlsLinearBasis2D_mex.mki
echo "EMC_CONFIG=optim" >> mlsLinearBasis2D_mex.mki
"/Applications/MATLAB_R2015a.app/bin/maci64/gmake" -B -f mlsLinearBasis2D_mex.mk
