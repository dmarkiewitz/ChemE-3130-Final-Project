using JuMP
using GLPKMathProgInterface
using GLPK
using DelimitedFiles
using DataFrames
using CSV
using PyPlot
using LinearAlgebra
using Gurobi
	
	#include("Stoichiometric_array_and_reference_vector_constructor.jl")

	deltaG_list=CSV.read("Thermo_Delta_G_Data.csv")
	aaa,bbb=size(deltaG_list)
	deltaG_list=deltaG_list[1:aaa-2,:]

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
	bad_rxns=Set{Int}()

	for j in bad_comp
	global bad_rxns
		for k in 1:noreactions
			if abs(S_matrix[j,k])>0.5
				bad_rxns=union(bad_rxns,k)
			end		
		end
	end

	bad_rxn=collect(bad_rxns)
	
	bad_rxn=sort(bad_rxn,rev=true)

	for l in bad_rxn
	global S_matrix
		a,b=size(S_matrix)
		S_matrix=hcat(hcat(S_matrix[:,1:(l-1)],zeros(a,1)),S_matrix[:,(l+1):b])		
	end

	Volume=50*10^-3

	nochemicals,noreactions=size(S_matrix)

	updatedchem,misc=size(deltaG_list)

	#adding the zero rows first
	#zrow=zeros(convert(Int,updatedchem-nochemicals),noreactions)

	#S_matrix=vcat(S_matrix,zrow)

	#The is	adding the reactions required for the diol
	#radd=zeros(updatedchem,1)
	#radd[17]=1
	#radd[28]=-1
	#radd[30]=1
	#radd[287]=-1
	#radd[288]=1	
	#S_matrix=hcat(S_matrix,radd)
	
	#radd2=zeros(updatedchem,1)
	#radd2[6]=1
	#radd2[277]=-1
	#radd2[288]=1
	#S_matrix=hcat(S_matrix,radd2)

	nochemicals,noreactions=size(S_matrix)

	no=(10*50*10^-6)*ones(nochemicals,1)
	#.+0.0084
	#no[30]=10
	#no[288]=111
	#no[28]=0
	#no[287]=0
	#no[6]=10
	#no[16]=10

	#selecting out the ordered delta G not of formation of chemical compounds
	delta_G_not_compound=deltaG_list[:,2]

	#prebuilding the delta_G_not_reaction vector
	delta_G_not_reaction=zeros(noreactions,1)
	
	#Calculating and building delta G not of reaction array
	for i in 1:noreactions
			delta_G_not_reaction[i]=dot(Transpose(delta_G_not_compound),S_matrix[:,i])	
	end

		#Model Construction

		EF_Model=Model(with_optimizer(Gurobi.Optimizer))

		#Variables

		@variable(EF_Model, excess[i=1:noreactions])

		#Objective

		@objective(EF_Model, Min, (sum(excess[i]*delta_G_not_reaction[i] for i=1:noreactions)/(8.314*298.15))^2)

		#Constraints
		for i in 1:nochemicals
			@constraint(EF_Model, no[i] + sum(S_matrix[i,j]*excess[j] for j=1:noreactions) >= 0) #If ni finals are known then can put those in as a constraint
		end	
		
		#Output
		optimize!(EF_Model)

		fudge=objective_value(EF_Model)
	
		#prebuilding vectors
		Excess=zeros(noreactions,1)
		ne=no

		#Building solution vector and matrix
		Excess[:]=JuMP.value.(excess[:])

		#getting ne vector
		for i in 1:nochemicals		
			ne[i]=no[i] + sum(S_matrix[i,j]*Excess[j] for j=1:noreactions)
		end

		#delta_G=dot(Transpose(delta_G_not_compound),no)-dot(Transpose(delta_G_not_compound),ne)
		
		
		#Prebuilding and calculating the K's
		K=ones(noreactions,1)
		
		for r in 1:noreactions
			for i in 1:nochemicals
				if S_matrix[i,r]!=0			
					k=(ne[i]/(Volume))^(S_matrix[i,r])
					K[r]=K[r]*k
					println(" $i ")
				end	
			end
		end		
