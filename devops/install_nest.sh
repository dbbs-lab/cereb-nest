echo ""
echo "Empty cache, starting NEST v$NEST_VERSION installation."
echo ""
export MY_BEFORE_DIR=$PWD
cd /home/travis
git clone --depth 1 https://github.com/dbbs-lab/nest-simulator nest-simulator-$NEST_VERSION
mkdir nest-simulator-$NEST_VERSION-build
mkdir nest-install-$NEST_VERSION
cd nest-simulator-$NEST_VERSION-build
cmake -Dwith-python=3 \
  -DPYTHON_EXECUTABLE=/home/travis/ \
  -DPYTHON_LIBRARY=/usr/lib/x86_64-linux-gnu/libpython3.6m.so \
  -DPYTHON_INCLUDE_DIR=/usr/include/python3.6m/ \
  -DCMAKE_INSTALL_PREFIX:PATH=/home/travis/nest-$NEST_VERSION \
  /home/travis/nest-simulator-$NEST_VERSION
make
make install
cd $MY_BEFORE_DIR
