echo "Unpacking data"

unzip ./data.zip

echo "Configuring and building RCD and RCDL"

cd src
mkdir bin
cd bin
cmake .. 
make -j8

cd ..
