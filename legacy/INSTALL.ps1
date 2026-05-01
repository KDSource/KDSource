# PowerShell script for installing KDSource on Windows 10

# Requirements: Git 2.14+, MinGW-GCC 11+, CMake 3+, Pip 22+ (Python 3.8+).

$SOURCEDIR = "$pwd"
$INSTALLDIR = "$pwd\..\KDinstall"

# Install dependencies
mkdir dependencies
cd dependencies
Invoke-Webrequest http://xmlsoft.org/sources/win32/64bit/libxml2-2.9.3-win32-x86_64.7z -Outfile libxml2-2.9.3-win32-x86_64.7z
Invoke-Webrequest http://xmlsoft.org/sources/win32/64bit/iconv-1.14-win32-x86_64.7z -Outfile iconv-1.14-win32-x86_64.7z
Invoke-Webrequest http://xmlsoft.org/sources/win32/64bit/zlib-1.2.8-win32-x86_64.7z -Outfile zlib-1.2.8-win32-x86_64.7z
7z x libxml2-2.9.3-win32-x86_64.7z -o"libxml2-2.9.3-win32-x86_64" -y
7z x iconv-1.14-win32-x86_64.7z -o"iconv-1.14-win32-x86_64" -y
7z x zlib-1.2.8-win32-x86_64.7z -o"zlib-1.2.8-win32-x86_64" -y
$env:PATH = "$SOURCEDIR\dependencies\libxml2-2.9.3-win32-x86_64;" + $env:PATH
$env:PATH = "$SOURCEDIR\dependencies\iconv-1.14-win32-x86_64;" + $env:PATH
$env:PATH = "$SOURCEDIR\dependencies\zlib-1.2.8-win32-x86_64;" + $env:PATH
$env:C_INCLUDE_PATH="$SOURCEDIR\dependencies\iconv-1.14-win32-x86_64\include"
cd ..

# Build and test
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX="$INSTALLDIR" -G "MinGW Makefiles"
mingw32-make
mingw32-make test
mingw32-make install
cd ..
$env:PATH += ";$INSTALLDIR\bin;$INSTALLDIR\lib"
kdtool-resample -h
