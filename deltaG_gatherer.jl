using DelimitedFiles

#= Collect the delG values for reactions at pH = 6.5 and equilibrium. =#
function Therm_Integration()
	""" 
	Get the delG values from the website for all the chemicals of interest

	Parameters:
	-----------
	None
	
	Returns:
	-----------
	Array of containing the delG values for each molecule. The molecule has a 
	delG value of 101010101010 if the value is missing from the website. 

	"""
	# Read chemical info.
	eqinfo = readdlm("kegg_compounds_Alberty_ph6.5.csv")
	chemicals_raw = readdlm("raw_chemicals.csv", '\n', header = false)

	DGMF=[]

	# Gathering the delG values.
	for p in 1:length(chemicals_raw)
		DGMF=push!(DGMF,101010101010) 
			for k in 1:length(eqinfo[:,1])
				if occursin(chemicals_raw[p],eqinfo[k,1])
					DGMF[p]=eqinfo[k,3]
				end						
			end					

	end
	
	DGPCN=hcat(chemicals_raw,DGMF)
	DGPCNP=vcat(["Kegg Name" "DeltaG of formation (KJ/mol)"],DGPCN)
	
	# Arbritrary save file 
	writedlm("JustinLee.csv",DGPCNP)	
	return [DGPCNP]
end




