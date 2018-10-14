# Lorenz 96 with Ginelli et al Algorithm for computing Covariant Lyapunov Vectors (CLVs)
by Sebastian Schubert (https://www.mi.uni-hamburg.de/arbeitsgruppen/theoretische-meteorologie/personen/schubert-sebastian.html)

E-Mail: sebastian-schubert(at)posteo.de

# General References #
# On the Lorenz 96 

Lorenz, E. N. (1996). Predictability - A problem partially solved. In Proc. Seminar on predictability.

# On CLVs

Trevisan, A., & Pancotti, F. (1998). Periodic Orbits, Lyapunov Vectors, and Singular Vectors in the Lorenz System. Journal of the Atmospheric Sciences, 55(3), 390–398. http://doi.org/10.1175/1520-0469(1998)055<0390:POLVAS>2.0.CO;2

Kuptsov, P. V., & Parlitz, U. (2012). Theory and Computation of Covariant Lyapunov Vectors. Journal of Nonlinear Science, 22(5), 727–762. http://doi.org/10.1007/s00332-012-9126-5

# Ginelli algorithm

Ginelli, F., Poggi, P., Turchi, A., Chaté, H., Livi, R., & Politi, A. (2007). Characterizing dynamics with covariant Lyapunov vectors. Physical Review Letters, 99(13), 130601. http://doi.org/10.1103/PhysRevLett.99.130601

# CLVs in a quasi-geostrophic two-layer model

Schubert, S., & Lucarini, V. (2015). Covariant Lyapunov vectors of a quasi-geostrophic baroclinic model: Analysis of instabilities and feedbacks. Quarterly Journal of the Royal Meteorological Society, 141(693), 3040–3055. http://doi.org/10.1002/qj.2588

Schubert, S., & Lucarini, V. (2016). Dynamical analysis of blocking events: spatial and temporal fluctuations of covariant Lyapunov vectors. Quarterly Journal of the Royal Meteorological Society, 142(698), 2143–2158. http://doi.org/10.1002/qj.2808

# On the Model Code

# Compilation

export MKL_NUM_THREADS=8 # choose number of mkl threads

ifort -r8 -O3 -assume byterecl -mkl lorenz96v5.90 -o l96.x

# Namelist parameters meaning

< to be included>

# Execute Programm

./l96.x < namelist
