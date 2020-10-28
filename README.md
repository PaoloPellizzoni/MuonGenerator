# MuonGenerator
Monte carlo simulations and data analysis algorithms for a high-energy physics experiments on positron-induced muon generation.

### Compilation
Compiling this software requires g++9 or newer, or a equivalent c++ compiler.
To compile the binaries under Linux environment use

    c++ -o testXXX testXXX.cpp utils.cpp reader.cpp
and under Windows use

    c++ -o testXXX.exe testXXX.cpp utils.cpp reader.cpp
    
### Project overview
The production of a high brillance muon beam is one of the most important challenge for the future of Particle Physics. A particularly interesting idea consists of shooting high energy positrons on a target, aiming at the production of muons by means of the process $e^+ + e^- \rightarrow \mu^+ + \mu^-$. To mimize the divergence of the resulting "muon beam", the positrons energy is chosen so that the reaction occurs close to threshold (assuming the electrons in the target to be at rest). We study the cross section oh the scattering and properties of the outgoing muon beam. We present an algorithm for reconstructing the trajectory of the muons emerging from a magnetic field region.
For a detailed description of the methods employed in this software, please check the file `report.pdf`.
   
    
