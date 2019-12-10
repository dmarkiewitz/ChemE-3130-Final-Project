using JuMP
using GLPKMathProgInterface
using GLPK
using DelimitedFiles
using DataFrames
using CSV
using PyPlot
using LinearAlgebra
	
	include("Stoichiometric_array_and_reference_vector_constructor.jl")
	
	nochemicals,noreactions=size(S_matrix)

	#selecting out the ordered delta G not of formation of chemical compounds
	delta_G_not_compound="comes from the spread sheet"

	#prebuilding the delta_G_not_reaction vector
	delta_G_not_reaction=zeros(length(chemicals_raw),1)
	
	#Calculating and building delta G not of reaction array
	for i in 1:noreactions
		delta_G_not_reaction[i]=dot(Transpose(delta_G_not_compound),S_matrix[:,i])
	end

		#Model Construction

		EF_Model=Model(with_optimizer(GLPK.Optimizer))

		#Variables

		@variable(EF_Model, excess[i=1:noreactions])

		#Objective

		@objective(EF_Model, Min, sum(excess[i]*delta_G_not_reaction[i] for i=1:noreactions)/(8.314*298.15))

		#Constraints
		for i in 1:nochemicals
			@constraint(EF_Model, no[i] + sum(S_matrix[i,j]*excess[j] for j=1:noreactions) >= 0) #If ni finals are known then can put those in as a constraint
		end
		#@constraint(EF_Model, S_matrix*V_vector[:,t] .== FBA_constraint)

	end	

		#Output
		optimize!(EF_Model)



