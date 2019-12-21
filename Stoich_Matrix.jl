using LibCURL
using DelimitedFiles

#= Develop the stoich matrix for our reaction network =#

function parse_from_web() 
    """ 
    Writes the reaction numbers categorized within the e.c. to a csv. 
    
    Parameters
    -----------
    None 
    
    Returns
    ------------
    None

    """
    #= 
    User input (need list of e.c's) put as vector only the form "x.x.x.x" 
    x's can be more than one digit but the x's represent a number.
    Results will not contrain exchange reactions those must be manually inputed
    Example: ecs = ["3.5.3.1";"2.1.3.3";"4.3.2.1";"6.3.4.5";"1.14.13.39"] 
    =#
    ecs = readdlm("ec.csv")
    rxns=[] 
    for ec in ecs 
        # Get molecular data from website 
        data =read(`curl -X GET http://rest.kegg.jp/link/rn/$ec/`, String) 

        # Parse data to get reaction number of the form of R$$$$$ where $ are ints
        rem_newline = split.(data,"\n")
        rem_tab = split.(rem_newline,"\t")
        rem_empty = filter(x -> x != [""],rem_tab)  # Filter empty web results
        rem_rn = map(x -> SubString(x[2],4),rem_empty) # Get the reaction
        rxns = union(rxns,rem_rn)                   # Could just be append 
    end 
    # Contains duplicate reactions that have different reaction numbers 
    writedlm("Reaction.csv", rxns, "\n")
end

function parse_rxn_exp()  
    """ 
    Writes the reaction expressions categorized by their reaction numbers 
    to a csv. 
    
    Parameters
    -----------
    None 
    
    Returns
    ------------
    None
    
    """
    exprxns=[] 
    # Removed rxns classified as duplicates (manual task because only a few)
    rxns = readdlm("Reaction_Rem_Dupe.csv",'\n',header=false)
    for rxn in rxns
        # Get the reaction expression from Kegg
        data=read(`curl -X GET http://rest.kegg.jp/get/reaction:$rxn/`, String)

        # Parse out the reaction expression
        split_eq=split(data,"\nEQUATION    ")[2]
        rem_rxn=split(split_eq,"\nRCLASS")[1]
        rem_comment=split(rem_rxn,"\nCOMMENT")[1]
        exprxns = append!(exprxns,[rem_comment])
    end
    # Contains messages about duplicate reaction numbers
    writedlm("ReactionExp.csv", exprxns, "\n")
end

function parse_rxns(rxn_file) 
    
    """ 
    Converts the reaction expressions to the respective Kegg chemical values and
    coefficients.

    Parameters
    -----------
    None 

    Returns
    ------------
    all_reactants: Array-like
        Reactant kegg values with their coefficients
    all_products: Array-like
        Product kegg values with their coefficients

    """
    # Remove dupe messages manually
    all_reactants = []
    all_products = []
    exprxns = readdlm(rxn_file, '\n', header=false)
    for exprxn in exprxns 
        # Split reactants and products to two arrays
        reactants, products = split(exprxn, r" -> | <=> ")   
        reactants, products = split(reactants, " + "), split(products, " + ")

        # Add 1 coefficients for easier parsing and negate reactant coeff.
        reactants = map(x -> if startswith(x,"C") "-1 "*x 
                             else "-"*x end, reactants)
        products = map(x -> if startswith(x,"C") "1 "*x else x end, products)
        all_reactants = append!(all_reactants,[split.(reactants, " ")])
        all_products = append!(all_products,[split.(products, " ")])
    end 
    return all_reactants, all_products
end

function raw_chemicals(new_file,rxn_file) 
    """
    Writes the unique list of kegg identifiers for the chemicals within all the 
    reactions we analyzed to a csv.

    Parameters
    -----------
    None

    Returns
    -----------
    None

    """
    # Implement unique sort
    reactants, products = parse_rxns(rxn_file)
    all_comp = append!(reactants, products)
    flat_comp = collect(Iterators.flatten(all_comp))
    all_chem = [flat_comp[i][2] for i in 1:length(flat_comp)]
    unq_sort = sort(unique(all_chem))
    writedlm(new_file, unq_sort, '\n',)
end


function reaction_matrix(chem_file,rxn_file,s_matrix_file)
    """
    Constructs the reaction matrix for all the reactions we are analyzing
    and writes it to a CSV
    
    Parameters
    -----------
    None

    Returns
    -----------
    None

    """
    # Important rxn data 
    chemicals = readdlm(chem_file,'\n')
    reactants,products = parse_rxns(rxn_file) 
    n_chemicals = length(chemicals)
    n_rxns = length(reactants)

    # Instantiate matrix = n_chemicals * n_rxns
    S = zeros(n_chemicals,n_rxns)
    
    # S[i][j]
    for j in 1:n_rxns    # n_cols
        for i in 1:n_chemicals  # n_rows
           # Find the stoich coeff of the ith chemical within the matrix  
            stoich_coeff = 0 
            num_prod = length(products[j])
            num_reac = length(reactants[j]) 
            chemical = chemicals[i]
            for k in 1:num_prod                      # Products of jth rxn 
                product = products[j][k]             # kth product in jth rxn
                product_co = product[1]              # Product coefficient

                # Update stoich coeff for the chemical in reaction j
                if chemical in product
                   stoich_coeff = stoich_coeff + parse(Int, product_co)
                end
            end  
            for k in 1:num_reac                      # Reactants of jth rxn 
                reactant = reactants[j][k]
                reactant_co = reactant[1]

                # Update stoich coeff for the chemical in reaction j
                if chemical in reactant
                    stoich_coeff = stoich_coeff + parse(Int, reactant_co)
                end
            end
            # Finalize the stoich coeff in the S matrix
            S[i,j] = stoich_coeff 
        end
    end
    writedlm(s_matrix_file,S,header=false)
end

# What we called to develop the two S matrices #

# No alcohol compounds of interest 
reaction_matrix("raw_chemicals_NoAlc.csv",
                "ReactionExp_Rem_Dupe_NoAlc.csv",
                "s_matrix_NoAlc.csv")

# Contains alcohol compounds of interest
reaction_matrix("raw_chemicals.csv",
                "ReactionExp_Rem_Dupe.csv",
                "s_matrix.csv")