#!/usr/local/bin/bash

for MESH in 2 3 4 5 6 
do
    # Run make geometry
    cd ../../build
    cmake .. -DMESHDEPTH=$MESH
    make
    ./makeGeometry

    # Rename file
    cd ../outputs
    mv laplace_out.txt laplace-$MESH.txt

    # Go back to starting dir
    cd solver_analysis
done

# Move to appropriate folder
cd ..
mv laplace* solver_analysis/laplace_mesh_refine/circle_2rad2