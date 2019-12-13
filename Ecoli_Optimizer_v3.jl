using JuMP
using GLPKMathProgInterface
using GLPK
using DelimitedFiles
using DataFrames
using CSV
using PyPlot
using LinearAlgebra
using Gurobi
 

# Get S_matrix and delta Gs precomputed and processed
S_matrix = CSV.read("S_matrix.csv") 
deltaG_list=CSV.read("Thermo_Delta_G_Data.csv")

# S_matrix is the stoichiometric matrix with rows representing the chemical 
# and columns representing the reactions 
n_chemicals,n_reactions = size(S_matrix)

# Storing indexes of chemicals whose reactions are missing delG data
bad_rxn= [j for i in 1:n_chemicals 
            if deltaG_list[i,2] >= 10^6 for j in 1:n_reactions 
              if abs(S_matrix[i,j]) > 0.5]
	
#= # 	# Storing indexes of bad reactions
  bad_rxn=[k for j in bad_comp for k in 1:n_reactions if abs(S_matrix[j,k]) > 0.5]
 # =#

  a,b=size(S_matrix)
  S_matrix = hcat([S[:,i] for i in 1:b if !(i in bad_rxn)])

  updatedchem,misc=size(deltaG_list)
  
	n_chemicals, n_reactions=size(S_matrix)

  # Define volume of cell 
  Vcell = 1*50*10^-6

  # Number of Moles
  no= Vcell*ones(n_chemicals,1)
  
	#.+0.0084
	#no[30]=10
	#no[288]=111
	#no[28]=0
	#no[287]=0
	#no[6]=10
	#no[16]=10

	# Selecting out the ordered standard delta G formation of chemical compounds
	delta_G_not_compound=deltaG_list[:,2]

	#prebuilding the delta_G_not_reaction vector
	delta_G_not_reaction=zeros(n_reactions,1)
	
	#Calculating and building delta G not of reaction array

  delta_G_not_reaction= [dot(Transpose(delta_G_not_compound),S_matrix[:,i]) 
                           for i in 1:n_reactions]

	#Model Construction

	EF_Model=Model(with_optimizer(Gurobi.Optimizer))

	#Variables

	@variable(EF_Model, excess[i=1:n_reactions])

	#Objective

	@objective(EF_Model, Min, (sum(excess[i]*delta_G_not_reaction[i] for i=1:n_reactions)/(8.314*298.15))^2)

	#Constraints
	for i in 1:n_chemicals
		@constraint(EF_Model, no[i] + sum(S_matrix[i,j]*excess[j] for j=1:n_reactions) >= 0) #If ni finals are known then can put those in as a constraint
	end	
		
	#Output
	optimize!(EF_Model)

	fudge=objective_value(EF_Model)
	
	#prebuilding vectors
  Excess=zeros(n_reactions,1)
  # Cheaty way of building ne 
	ne=no

	# Building solution vector and matrix
	Excess[:]=JuMP.value.(excess[:])

	#getting ne vector
	for i in 1:n_chemicals		
		ne[i]=no[i] + sum(S_matrix[i,j]*Excess[j] for j=1:n_reactions)
	end

	#delta_G=dot(Transpose(delta_G_not_compound),no)-dot(Transpose(delta_G_not_compound),ne)
	
	
	#Prebuilding and calculating the K's
	K=zeros(n_reactions,1)
	
  K=[prod((ne[i]/Vcell)^(S_matrix[i,r]) for i in 1:n_chemicals) for r in 1:n_reactions]		
