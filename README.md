**Packages to install:**
Gurobi, CSV, DataFrames, DelimitedFiles, JuMP, LibCURL, LinearAlgebra

1)Type 'include("Stoich_Matrix.jl")' to see how the stoichiometric matrix for the reactions were generated

2)Type 'include("Ecoli_Optimizer_v3.jl")' to see the data returned from the
optimizer. 

**For 2:**
Call alc_data to receive the tuple of arrays specified by the function optimizer. 

Including the alcohol of interest:
alc_data[1] yields the array of delGs rxns calculated from the optimizer for each rxn for all 3 concentrations. 
alc_data[2] yields the array of the final moles of each compound for all 3 concentrations
alc_data[3] yields the array containing the error from the equilibrium for all 3 concentrations 

no_alc_data is essentially the same except the data does not include the alcohols of interest within the reaction network. 

See documentation for other functions.
