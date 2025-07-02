# Fission-Yeast-Rec8

This code base is using the [Julia Language](https://julialang.org/) and
[DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project for the paper  

**"Meiotic cohesin Rec8 imposes fitness costs on fission yeast gametes to drive evolution of parental bias in gene expression"**  

authored by:  

Celso Martins (1,5) , Harry Booth (2,5), Clàudia Salat-Canela (1), Zena Hadjivasiliou (2,3,4), Aleksandar Vještica (1)  

1 Center for Integrative Genomics, Génopode, University of Lausanne, Lausanne, Switzerland  
2 Mathematical and Physical Biology Laboratory, The Francis Crick Institute, London, United Kingdom  
3 Department of Physics and Astronomy, University College London, London, United Kingdom  
4 Institute for the Physics of Living Systems, University College London, London, United Kingdom  

To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.

You may notice that most scripts start with the commands:
```julia
using DrWatson
@quickactivate "FissionYeast"
```
which auto-activate the project and enable local path handling from DrWatson.
