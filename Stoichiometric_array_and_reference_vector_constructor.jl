using LibCURL
using DelimitedFiles

#User input (need list of e.c's) put as column vector only the form "x.x.x.x" x's can
#be more than one digit but the x's represent a number
#results will not contrain exchange reactions those must be manuelly inputed
#Example ec=["3.5.3.1";"2.1.3.3";"4.3.2.1";"6.3.4.5";"1.14.13.39"]
ec=["3.1.3.10";"5.4.2.2";"2.7.1.199";"3.1.3.9";"2.7.1.1";"2.7.1.2";"2.7.1.63"
;"2.7.1.147";"5.1.3.3";"5.1.3.15";"5.3.1.9";"2.7.1.-";"3.2.1.86";"3.1.3.11";"2.7.1.11";"2.7.1.146"
;"2.7.1.90";"4.1.2.13";
"5.3.1.1";"1.2.1.12";"1.2.1.59";"1.2.1.9";"1.2.7.6";"1.2.1.90";"5.4.2.4";"2.7.2.3";"5.4.2.11"
;"5.4.2.12";"3.1.3.80";"4.2.1.11";"2.7.1.40";"4.1.1.32";"4.1.1.49"
;"1.2.7.1";"1.2.7.11";"1.1.1.27";"1.2.4.1";"2.3.1.12";"1.1.1.1"
;"1.1.2.7";"1.1.1.2";"1.1.5.5";"1.1.2.8";"6.2.1.1";"6.2.1.13";"1.8.1.4";"1.2.1.3";"1.2.1.5";"1.2.1.-";
"6.4.1.1";"4.2.1.3";"1.1.1.42";"1.1.1.41";"1.1.1.286"
;"1.2.4.2";"2.3.1.61";"2.3.3.3";"2.3.3.8";"2.3.3.1";"1.2.7.3";"1.1.1.37";"1.1.5.4";"4.2.1.2";
"1.3.5.4";"1.3.5.1";"6.2.1.4";"6.2.1.5";"2.8.3.18";"1.1.1.47"
;"1.1.3.4";"1.1.3.5";"3.1.1.17";"1.1.1.360";"1.1.1.359";"1.1.5.9";"1.1.5.2";"1.1.99.3";"1.1.1.215";"2.7.1.13";
"1.1.1.43";"2.7.1.12";"1.1.1.49";"1.1.1.363";"1.1.1.1388"
;"5.3.1.27";"2.2.1.1";"2.2.1.2";"4.1.2.9";"4.1.2.4";"2.7.1.15";"5.4.2.7";"4.2.1.39";"4.2.1.40";"4.3.1.9";"4.1.2.55";
"4.1.2.51";"2.7.1.45";"2.7.1.178";"3.1.1.31"
;"4.2.1.12";"1.1.1.44";"1.1.1.343";"2.7.1.203";"4.3.1.29";"4.1.2.43";"5.1.3.1";"5.3.1.6";"2.7.6.1";"2.7.4.23";"4.1.2.14";"1.2.99.8";"1.2.1.89";"1.2.7.5";"2.7.1.165"]

#rxns will contain a list of the chemical reactions
#creactants contains a list of chemical reactants in terms of
#kegg nomincalature from which matrix EB short for elemental balence is
#constructed and with the rows correctioning to the elements (can be expanded to entire periodic table if wanted)
#but here it is [C;H;N;O;P;S], as for the columns they are ordered by the metabolites with row equal to column number

reb=["C";"H";"N";"O";"P";"S"]#rows of elemental balance matrix

rxns=Set{String}() # used to prevent duplicates
creactants=Set{String}() #no used here but could be used for larger reaction networks
exprxns=Set{String}() #used to preven duplicates

#h's are temporary vectors
for i in 1:length(ec)
    global rxns
    global stp
    h=ec[i]
    h1=read(`curl -X GET http://rest.kegg.jp/link/rn/$h/`, String)

        s=split(h1,"\t")
        st=split.(s,"\n")
        stp=st[2:length(st)]
        for i in 1:length(stp)
            global stp
            global rxns
            rxns=union(rxns,[split.(stp[i][1],"rn:")[2]])
        end
end

#change set to vector for easy access
rxnv=collect(rxns)

for j in 1:length(rxnv)
    global rxnv
    global exprxns

    h=rxnv[j]
    v=read(`curl -X GET http://rest.kegg.jp/get/reaction:$h/`, String)
    t=split(v,"\nEQUATION    ")[2]
    d=split(t,"\nRCLASS")[1]
    c=split(d,"\nCOMMENT")[1]
    exprxns=union(exprxns,[c])

end

#vectorizing
exprxnv=collect(exprxns)


#adding on ex rxn givens
#exprxnv=push!(exprxnv," -> C00169")



#adding on ex rxn to balance elements, Added after running the first tome
#needed an ATP exchange rxn and an ADP exchange rxn to have a full elemental balance
#assume this is an exchange required and a H20 export rxn. the none exchange rxns were already balenced
# example exprxnv=push!(exprxnv,"C00001 -> ")


chemicals_raw=[]
for i in 1:length(exprxnv)#length(rxns) #column i of matrix
    global exprxnv
    global chemicals_raw
    bigsplit=split(exprxnv[i],r" -> | <=> ")#splits rxn into left and right half
    left_rxn=bigsplit[1]
    left_chem_raw=split(left_rxn," + ")#subdivides the left half into chemicals and their coeffiecent
    left_chem_and_co=split.(left_chem_raw," ")#splits up chemicals and their coeffiecent
    for j in 1:length(left_chem_and_co)#collects all unique chemicals on left side
        if length(left_chem_and_co[j])==2 && (!(in(left_chem_and_co[j][2],chemicals_raw)) && left_chem_and_co[j][2]!="")#if chemical has coeffiecent + uniqueness constraints + protection against "" inculsion
            chemicals_raw=push!(chemicals_raw,left_chem_and_co[j][2])
        elseif length(left_chem_and_co[j])==1 && (!(in(left_chem_and_co[j][1],chemicals_raw)) && left_chem_and_co[j][1]!="")#if chemical has no coeffiecent + uniqueness constraints + protection against "" inculsion
            chemicals_raw=push!(chemicals_raw,left_chem_and_co[j][1])
        end
    end
    right_rxn=bigsplit[2]
    right_chem_raw=split(right_rxn," + ")#subdivides the right half into chemicals and their coeffiecent
    right_chem_and_co=split.(right_chem_raw)#splits up chemicals and their coeffiecent
    for j in 1:length(right_chem_and_co)#collects all unique chemicals on right side
        if length(right_chem_and_co[j])==2 && !(in(right_chem_and_co[j][2],chemicals_raw)) && !(right_chem_and_co[j][2]=="")#if chemical has coeffiecent + uniqueness constraints + protection against "" inculsion
            chemicals_raw=push!(chemicals_raw,right_chem_and_co[j][2])
        elseif length(right_chem_and_co[j])==1 && !(in(right_chem_and_co[j][1],chemicals_raw)) && !(right_chem_and_co[j][1]=="")#if chemical has no coeffiecent + uniqueness constraints + protection against "" inculsion
            chemicals_raw=push!(chemicals_raw,right_chem_and_co[j][1])
        end
    end
end
chemicals_raw

#building S_matrix

S_matrix=zeros(length(chemicals_raw),length(exprxnv))

#processing of rxns into stoichometric matrix
for j in 1:length(exprxnv)#length(rxns) #column j of matrix
    global S_matrix
    global exprxnv
    global chemicals_raw
    bigsplit=split(exprxnv[j],r" -> | <=> ")#splits rxn into left and right half
    left_rxn=bigsplit[1]
    left_chem_raw=split(left_rxn," + ")#subdivides the left half into chemicals and their coeffiecent
    left_chem_and_co=split.(left_chem_raw," ")#splits up chemicals and their coeffiecent
    for k in 1:length(left_chem_and_co)#Places coeffiecents in s_matrix and since its the left side the entries recieve an negative sign
        if length(left_chem_and_co[k])==2 && left_chem_and_co[k][2]!=""#grabs coeffiecent of chemicals if present
            for i in 1:length(chemicals_raw) #possible row number
                if occursin(left_chem_and_co[k][2],chemicals_raw[i])#finds correct row number
                    S_matrix[i,j]=-1*parse(Int,left_chem_and_co[k][1])
                end
            end
        elseif length(left_chem_and_co[k])==1 && left_chem_and_co[k][1]!=""#grabs present chemicals and assigns -1 to their S_matrix possition
            for i in 1:length(chemicals_raw) #possible row number
                if occursin(left_chem_and_co[k][1],chemicals_raw[i])#finds correct row number
                    S_matrix[i,j]=-1
                end
            end
        end
    end
    right_rxn=bigsplit[2]
    right_chem_raw=split(right_rxn," + ")#subdivides the right half into chemicals and their coeffiecent
    right_chem_and_co=split.(right_chem_raw)#splits up chemicals and their coeffiecent
    for k in 1:length(right_chem_and_co)
        if length(right_chem_and_co[k])==2 && !(right_chem_and_co[k][2]=="")#grabs coeffiecent of chemicals if present
            for i in 1:length(chemicals_raw) #possible row number
                if occursin(right_chem_and_co[k][2],chemicals_raw[i])#finds correct row number
                    S_matrix[i,j]=parse(Int,right_chem_and_co[k][1])
                end
            end
        elseif length(right_chem_and_co[k])==1 && right_chem_and_co[k][1]!=""#grabs present chemicals and assigns 1 to their S_matrix position
            for i in 1:length(chemicals_raw) #possible row number
                if occursin(right_chem_and_co[k][1],chemicals_raw[i])#finds correct row number
                    S_matrix[i,j]=1
                end
            end
        end
    end
end

S_matrix

chemicals_raw
