using LibCURL
using DelimitedFiles

#= 
User input (need list of e.c's) put as column vector only the form "x.x.x.x" 
x's can be more than one digit but the x's represent a number
results will not contrain exchange reactions those must be manuelly inputed
Example ec=["3.5.3.1";"2.1.3.3";"4.3.2.1";"6.3.4.5";"1.14.13.39"] 
=#
ecs=["3.1.3.10";"5.4.2.2";"2.7.1.199";"3.1.3.9";"2.7.1.1";"2.7.1.2";"2.7.1.63";
    "2.7.1.147";"5.1.3.3";"5.1.3.15";"5.3.1.9";"2.7.1.-";"3.2.1.86";"3.1.3.11"; 
    "2.7.1.11";"2.7.1.146";"2.7.1.90";"4.1.2.13";"5.3.1.1";"1.2.1.12";
    "1.2.1.59";"1.2.1.9";"1.2.7.6";"1.2.1.90";"5.4.2.4";"2.7.2.3";"5.4.2.11";
    "5.4.2.12";"3.1.3.80";"4.2.1.11";"2.7.1.40";"4.1.1.32";"4.1.1.49";"1.2.7.1";
    "1.2.7.11";"1.1.1.27";"1.2.4.1";"2.3.1.12";"1.1.1.1";"1.1.2.7";"1.1.1.2";
    "1.1.5.5";"1.1.2.8";"6.2.1.1";"6.2.1.13";"1.8.1.4";"1.2.1.3";"1.2.1.5";
    "1.2.1.-";"6.4.1.1";"4.2.1.3";"1.1.1.42";"1.1.1.41";"1.1.1.286";"1.2.4.2";
    "2.3.1.61";"2.3.3.3";"2.3.3.8";"2.3.3.1";"1.2.7.3";"1.1.1.37";"1.1.5.4";
    "4.2.1.2";"1.3.5.4";"1.3.5.1";"6.2.1.4";"6.2.1.5";"2.8.3.18";"1.1.1.47";
    "1.1.3.4";"1.1.3.5";"3.1.1.17";"1.1.1.360";"1.1.1.359";"1.1.5.9";"1.1.5.2";
    "1.1.99.3";"1.1.1.215";"2.7.1.13";"1.1.1.43";"2.7.1.12";"1.1.1.49";
    "1.1.1.363";"1.1.1.1388";"5.3.1.27";"2.2.1.1";"2.2.1.2";"4.1.2.9";
    "4.1.2.4";"2.7.1.15";"5.4.2.7";"4.2.1.39";"4.2.1.40";"4.3.1.9";"4.1.2.55";
    "4.1.2.51";"2.7.1.45";"2.7.1.178";"3.1.1.31";"4.2.1.12";"1.1.1.44";
    "1.1.1.343";"2.7.1.203";"4.3.1.29";"4.1.2.43";"5.1.3.1";"5.3.1.6";
    "2.7.6.1";"2.7.4.23";"4.1.2.14";"1.2.99.8";"1.2.1.89";"1.2.7.5";"2.7.1.165"]

""" 
    Writes the reaction numbers categorized within the e.c. to a csv. 
    
    Parameters
    -----------
    None 
    
    Returns
    ------------
    None
"""
function parse_from_web() 
  
    rxns=[] 
    for ec in ecs 
        # Get molecular data from website 
        data =read(`curl -X GET http://rest.kegg.jp/link/rn/$ec/`, String) 
        rem_newline = split.(data,"\n")
        rem_tab = split.(rem_newline,"\t")
        rem_empty = filter(x -> x != [""],rem_tab)
        rem_rn = map(x -> SubString(x[2],4),rem_empty)
        # Remove possible duplicate numeraically classified reactions
        rxns = union(rxns,rem_rn)
    end 
    # Contains duplicates
    writedlm("Reaction.csv", rxns, "\n")
end

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
function parse_rxn_exp()  
    exprxns=[] 
    # Removed Duplicate rxns (manual task because aliases)
    rxns = readdlm("Reaction_Rem_Dupe.csv",'\n',header=false)
    for rxn in rxns
        data=read(`curl -X GET http://rest.kegg.jp/get/reaction:$rxn/`, String)
        split_eq=split(data,"\nEQUATION    ")[2]
        rem_rxn=split(split_eq,"\nRCLASS")[1]
        rem_comment=split(rem_rxn,"\nCOMMENT")[1]
        exprxns = append!(exprxns,[rem_comment])
    end
    # Contains duped messages
    writedlm("ReactionExp.csv", exprxns, "\n")
end

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
function parse_rxns() 
    # Remove dupe messages manually
    all_reactants = []
    all_products = []
    exprxns = readdlm("ReactionExp_Rem_Dupe.csv", '\n', header=false)
    for exprxn in exprxns 
        reactants, products = split(exprxn, r" -> | <=> ")   
        reactants, products = split(reactants, " + "), split(products, " + ")

        # Add 1 coefficients for easier parsing and negate reactant coeff.
        reactants = map(x -> if startswith(x,"C") "-1 "*x 
                             else "-"*x end, reactants)
        products = map(x -> if startswith(x,"C") "1 "*x else x end, products)
        all_reactants = append!(all_reactants,[split.(reactants, " ")])
        all_products = append!(all_products,[split.(products, " ")])
    end 
    return all_reactants,all_products
end

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
function raw_chemicals() 
    reactants, products = parse_rxns()
    all_comp = append!(reactants, products)
    flat_comp = collect(Iterators.flatten(all_comp))
    all_chem = [flat_comp[i][2] for i in 1:length(flat_comp)]
    writedlm("raw_chemicals.csv", sort(unique!(all_chem)), '\n',)
end

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
function reaction_matrix()
    # Important rxn data 
    chemicals = readdlm("raw_chemicals.csv", '\n', header=false)
    reactants,products = parse_rxns() 
    n_chemicals = length(chemicals)
    n_rxns = length(reactants)

    # Instantiate matrix = n_chemicals * n_rxns
    S = zeros(n_chemicals,n_rxns)

    # S[i][j]
    for j in 1:n_rxns   # n_cols
        for i in 1:n_chemicals  # n_rows
            stoich_coeff = 0 
            num_prod = length(products[j])
            num_reac = length(reactants[j]) 
            chemical = chemicals[i]
            for k in 1:num_prod                      # Products of jth rxn 
                product = products[j][k]              
                product_co = product[1]

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
    writedlm("s_matrix.csv",S,header=false)
end