using JuMP
using GLPKMathProgInterface
using GLPK
using DelimitedFiles
using DataFrames
using CSV
using PyPlot
using LinearAlgebra
	
	include("Stoichiometric_array_and_reference_vector_constructor.jl")

	deltaG_list=CSV.read("Thermo_Delta_G_Data.csv")

	nochemicals,noreactions=size(S_matrix)

	#storing indexes of chemicals whose reactions 
	bad_comp=[]

	for i in 1:nochemicals
	global bad_comp
		if deltaG_list[i,2]>= 10^6
			bad_comp=push!(bad_comp,i)
		end

	end
	
	#Storing indexes of bad reactions
	bad_rxn=[]

	for j in bad_comp
	global bad_rxn
		for k in 1:noreactions
			if S_matrix[j,k]>0.5
				bad_rxn=push!(bad_rxn,k)
			end		
		end
	end

	bad_rxn=sort(bad_rxn,rev=true)

	for l in bad_rxn
	global S_matrix
		a,b=size(S_matrix)
		S_matrix=hcat(S_matrix[:,1:(l-1)],S_matrix[:,(l+1):b])		
	end

	bad_rxn_v2=[]
	for k in [34,147,40,53,60,64,67,73,102,111,130,141,144,147,161,168,208,209,220,225]
	global bad_rxn_v2
		bad_rxn_v2=push!(bad_rxn_v2,k)
	end
	
	bad_rxn_v2=sort(bad_rxn_v2,rev=true)

	for l in bad_rxn_v2
	global S_matrix
		a,b=size(S_matrix)
		S_matrix=hcat(S_matrix[:,1:(l-1)],S_matrix[:,(l+1):b])		
	end

	nochemicals,noreactions=size(S_matrix)

	no=0.002*ones(length(chemicals_raw))

	#selecting out the ordered delta G not of formation of chemical compounds
	delta_G_not_compound=deltaG_list[:,2]

	#prebuilding the delta_G_not_reaction vector
	delta_G_not_reaction=zeros(length(chemicals_raw),1)
	
	#Calculating and building delta G not of reaction array
	for i in 1:noreactions
		#if i!=34 && i!=147 && i!=40 && i!=53 && i!=60 && i!=64 && i!=67 && i!=73 && i!=102 && i!=111 && i!=130 && i!=141 && i!=144 && i!=161 && i!=168 && i!=208 && i!=209 && i!=220 && i!=225 #53 is 2x error 64 is 3x error 73 is a factor of 100x error 130 is a factor of 3x error
			delta_G_not_reaction[i]=dot(Transpose(delta_G_not_compound),S_matrix[:,i])
		#end	
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

		#Output
		optimize!(EF_Model)

		fudge=objective_value(EF_Model)

