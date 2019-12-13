using JuMP
using DelimitedFiles
using DataFrames
using CSV
using LinearAlgebra
#using Gurobi
 

# Get S_matrix and delta Gs precomputed and processed
s_matrix = readdlm("S_matrix.csv") 
delG_data = CSV.read("Delta_G_Data.csv")

# S_matrix is the stoichiometric matrix with rows representing the chemical 
# and columns representing the reactions 
n_chemicals,n_reactions = size(s_matrix)

# Storing col whose reactions contain compounds missing delG data
function zero_cols_matrix(s_matrix)
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
			 
# Update the s_matrix 
s_matrix = zero_cols_matrix(s_matrix)

# Define volume of cell 
Vcell = 1*50*10^-6

# Number of Moles
nMol= Vcell*ones(n_chemicals,1)
  
#.+0.0084
#nMol[30]=10
#nMol[288]=111
#nMol[28]=0
#nMol[287]=0
#nMol[6]=10
#nMol[16]=10


# Calculating deltaGrxn of all rxns in units of J/mol
conversion = 1000
delG_formation = [delG_data[i,2] for i in 1:n_chemicals] 
delG_rxn = [conversion*dot(Transpose(delG_formation),s_matrix[:,i]) 
												 for i in 1:n_reactions]
					 					
# Optimization 
EF_Model = Model(with_optimizer(Gurobi.Optimizer))

@variable(EF_Model, extent[i=1:n_reactions])


# STP conditions 1 atm and 298.15K
R = 8.314 	# Units of J/(K*Mol)
T = 298.15  # Units of K
summmation = sum(extent[i]*delG_rxn[i] for i=1:n_reactions)
to_min_term = summation/(R*T)
@objective(EF_Model, Min, to_min_term^2)

for i in 1:n_chemicals
	@constraint(EF_Model, nMol[i] + sum(S_matrix[i,j]*extent[j] for j=1:n_reactions) >= 0) 
end 

optimize!(EF_Model)
fudge=objective_value(EF_Model)

#prebuilding vectors
extent=zeros(n_reactions,1)
# Cheaty way of building ne 
ne=nMol
# Building solution vector and matrix
extent[:]=JuMP.value.(extent[:])
#getting ne vector
for i in 1:n_chemicals		
	ne[i]=nMol[i] + sum(S_matrix[i,j]*extent[j] for j=1:n_reactions)
end
#delta_G=dot(Transpose(delta_G_nMolt_compound),nMol)-dot(Transpose(delta_G_nMolt_compound),ne)


#Prebuilding and calculating the K's
K=zeros(n_reactions,1)

K=[prod((ne[i]/Vcell)^(S_matrix[i,r]) for i in 1:n_chemicals) for r in 1:n_reactions]		
  =#