using JuMP
using DelimitedFiles
using DataFrames
using CSV
using Gurobi
using LinearAlgebra
 

# Get S_matrix and delta Gs precomputed and processed
s_matrix = readdlm("s_matrix.csv") 
no_alc_s_matrix = readdlm("s_matrix_NoAlc.csv")
delG_data = CSV.read("Delta_G_Data.csv")
no_alc_delG_data = CSV.read("Delta_G_Data_NoAlc.csv")
  
function zero_cols_matrix(s_matrix, delG_data)
	""" 
	Zero out the columns with compounds with missing delG data 
			
	Parameters:
	-----------
	s_matrix: n_compound by n_rxn matrix 
		The stoich matrix for the reaction network
	delG_data: DataFrame 
		The delG data for the kegg compound. Missing data indicated by delG value of 
		1.01e11 kJ/mol
	
	Returns:
	The s_matrix with rxn columns zeroed out if missing chemical compound data.

	""" 
	n_chemicals,n_reactions = size(s_matrix)
	for i in 1:n_chemicals 
		if delG_data[i,2] > 10^6 
			for j in 1:n_reactions 
				if s_matrix[i,j] != 0
						s_matrix[:,j] .= 0
				end
			end
		end
	end
	return s_matrix
end
			 


function optimizer(v_cell,s_matrix,delG_data)
	"""
	Solve the optimization problem of the project 

	Parameters:
	-----------
	v_cell: float
		The volume of the E. coli. cells

	s_matrix: matrix
		Stoichiometric matrix for the reaction

	delG_data:
		The delG values for the chemical reactions 

	Returns:
	-----------
	Three arrays containing the delG of each rxn, final moles of each compound, 
	and error from the equilibrium condition. 
	
	"""
	all_DG, all_ne, all_error = [],[],[]
	# Concentrations in Molarity for the cell solution 
	for conc in [1*10^-3,5*10^-3,10*10^-3]
		# Remove all rxns with missing data 
		s_matrix = zero_cols_matrix(s_matrix,delG_data)	

		# Getting parameters 
		n_chemicals,n_reactions = size(s_matrix)	

		# Number of Moles initially for optimizer arbritrary
		n_mol_init= conc*v_cell*ones(n_chemicals,1)	
		s_matrix = zero_cols_matrix(s_matrix, delG_data) 

		# Calculating deltaGrxn of all rxns in units of J/mol
		conversion = 1000
		delG_formation = [delG_data[i,2] for i in 1:n_chemicals] 
		delG_rxn = [conversion*dot(Transpose(delG_formation),s_matrix[:,i]) for i in 1:n_reactions]

		# Optimization 
		EF_Model = Model(with_optimizer(Gurobi.Optimizer))	
		@variable(EF_Model, extent[i=1:n_reactions])	

		# STP conditions 1 atm and 298.15K
		R = 8.314 	# Units of J/(K*Mol)
		TK = 298.15  # Units of K
		summation = sum(extent[i]*delG_rxn[i] for i=1:n_reactions)
		to_min_term = summation/(R*TK)
		@objective(EF_Model, Min, to_min_term^2)	

		# Number of moles can not go negative
		for i in 1:n_chemicals
			@constraint(EF_Model, n_mol_init[i] + sum(s_matrix[i,j]*extent[j] for j=1:n_reactions) >= 0) 
		end 	
		optimize!(EF_Model)
		error = objective_value(EF_Model)	

		# Building solution vector and matrix
		extent = JuMP.value.(extent[:])	

		# Getting the n_e vector		
		ne = abs.([n_mol_init[i] + sum(s_matrix[i,j]*extent[j] for j in 1:n_reactions)
						for i=1:n_chemicals])

		# Equilibrium constants for each reaction
		K = [prod((ne[i]/v_cell)^s_matrix[i,r] for i in 1:n_chemicals) 
				 for r in 1:n_reactions]
		
		# DG value from equilibrium constant 
		DG=-TK*R*log.(K)./conversion	

		# Error from equilibrium condition
		error=sum(extent[i]*delG_rxn[i] for i=1:n_reactions)/(R*TK)	

		all_DG = append!(all_DG, [DG])
		all_ne = append!(all_ne, [ne])
		all_error = append!(all_error, [error])
	end
	return all_DG,all_ne,all_error
end

# Data including alcohol of interest
alc_data = optimizer(1e-15, s_matrix, delG_data)

# Data not including alcohol of interest
no_alc_data = optimizer(1e-15, no_alc_s_matrix, no_alc_delG_data)
