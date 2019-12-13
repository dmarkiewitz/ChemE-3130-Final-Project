using JuMP
using DelimitedFiles
using DataFrames
using CSV
using LinearAlgebra
using Gurobi
 

# Get S_matrix and delta Gs precomputed and processed
s_matrix = readdlm("s_matrix.csv") 
no_alc_s_matrix = readdlm("s_matrix_NoAlc.csv")
delG_data = CSV.read("Delta_G_Data.csv")
no_alc_delG_data = CSV.read("Delta_G_Data_NoAlc.csv")
  
# Storing col whose reactions contain compounds missing delG data
function zero_cols_matrix(s_matrix, delG_data)
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
			 


function optimizer(conc, v_cell,s_matrix,delG_data) 
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

	# number of moles can not go below 0! Possibly mole fraction criteria too?
	for i in 1:n_chemicals
		@constraint(EF_Model, n_mol_init[i] + sum(s_matrix[i,j]*extent[j] for j=1:n_reactions) >= 0) 
	end 

	optimize!(EF_Model)
	error = objective_value(EF_Model)

	# Building solution vector and matrix
	extent = JuMP.value.(extent[:])

	#getting ne vector		
	ne = abs.([n_mol_init[i] + sum(s_matrix[i,j]*extent[j] for j in 1:n_reactions)
					for i=1:n_chemicals])
						
		 
	K = [prod((ne[i]/v_cell)^s_matrix[i,r] for i in 1:n_chemicals) 
			 for r in 1:n_reactions]
	
	DG=-TK*R*log.(K)./conversion

	error=sum(extent[i]*delG_rxn[i] for i=1:n_reactions)/(R*TK)

	return DG, ne, error
end
#.+0.0084
#n_mol_init[30]=10
#nMol[288]=111
#nMol[28]=0
#nMol[287]=0
#nMol[6]=10
#nMol[16]=10 
