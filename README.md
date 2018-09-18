# Memory-Aware-LBM 

Current work on 2D grid, using D2Q9 model.
The paper is accepted in SBAC-PAD 2018. https://graal.ens-lyon.fr/sbac-pad/index.php/cfp/accepted-papers
Our future work will extend and apply the memory-aware idea to 3D and distributed memory HPC systems.

## Modules
```
module load cmake
module load papi
```

## Compile
```
mkdir build
cd build/
cmake -DCMAKE_BUILD_TYPE=Bridges ..
cmake -DCMAKE_BUILD_TYPE=Stampede ..
```
