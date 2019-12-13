using DelimitedFiles

function Therm_Integration()
	#obtaining delta g information in kj/mol
	eqinfo=readdlm("kegg_compounds_Alberty_ph6.5.csv")
	chemicals_raw = readdlm("raw_chemicals.csv", '\n', header = false)

	#intializing storgage arrays
	DGMF=[]

	#gathering the delta g's of formation and the indexes of the reactions and the metabolites involved
	for p in 1:length(chemicals_raw)
		DGMF=push!(DGMF,101010101010) #This will be the number that represents to Delta G of formation was not found
			for k in 1:length(eqinfo[:,1])
				if occursin(chemicals_raw[p],eqinfo[k,1])
					DGMF[p]=eqinfo[k,3]
				end						
			end					

	end

	#Building a coresponding value array
	DGPCN=hcat(chemicals_raw,DGMF)
	DGPCNP=vcat(["Kegg Name" "DeltaG of formation (KJ/mol)"],DGPCN)
	
	#Saved
	writedlm("JustinLee.csv",DGPCNP)	

	#Returning wanted information
	return [DGPCNP]

end




