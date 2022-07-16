


# setdiff peut être intéressant
# example : compB=setdiff(cI[i],B)

# A[1,:] different A[[1],:]

struct IBvu{U,V}      #Structure for I with, the B_I, the v_B^i and u^B_i
    I::Vector{Int64}        #Set of indices
    boolI::BitArray{1}      #Boolean does i in I
    B::Vector{Int64}        #Extracted bases B_I
    boolB::BitArray{1}      #Boolean does i in B
    IlessB::Vector{Int64}   #I \ B_I
    compI::Vector{Int64}    #[nb_const]\backslash I
    v::Matrix{U}            #Matrix of (v^B_i)
    u::Vector{U}            #Vector of scalar (u^B_i)
    invWBTB::Matrix{U}      #Matrix W_{B_I}^{-1}T_{B_I}
    invWBhB::Vector{U}      #Vector W_{B_I}^{-1}h_{B_I}
    espq::Vector{V}         #Vector EspCostInd
    α::Vector{V}            #Vector T_{B_I}^\top W_{B_I}^{-1\top}EspCostInd
    β::V                    #Scalar h_{B_I}^\top W_{B_I}^{-1\top}EspCostInd
    
    #The following are added not to recompute the subdivisions in the vertex to edge phase
    IlessBnotempty::Bool
    subdiv::Vector{Vector{IBvu}}    #Vector of the same size of I
                                    #The cIBvu corresponding to I\{i} for all i in I
    subdivα::Vector{Vector{U}}      #sum(cJvunew[i].α for i in 1:lengcJnew)-cIvu[J].α
    subdivβ::Vector{V}              #sum(cJvunew[i].β for i in 1:lengcJnew)-cIvu[J].β
    subdiv_computed::BitArray{1}    #List of Booleanto know if the subdivisions are computed
end


function IBvu(cons::CoupCons,I::Vector{Int64},cost::CostDistribution,typecons::DataType=Float64,typecost::DataType=Float64) #Constructor of IBvu
    boolI=falses(cons.nb_const)
    for i in I
        boolI[i]=true
    end
    B=BasisExtraction(cons.W,I)
    boolB=falses(cons.nb_const)
    for i in B
        boolB[i]=true
    end
    U=typeof(cons.h[1])
    u=typecons.(zeros(cons.nb_const))
    v=typecons.(zeros(cons.nb_const,cons.dimx))
    IlessB=Vector{Int64}(undef,0)
    compI=Vector{Int64}(undef,0)
    invWB=inv(cons.W[B,:])   #Store the inverse not to compute it many times
    invWBTB=typecons.(invWB*cons.T[B,:])
    invWBhB=typecons.(invWB*cons.h[B])
    for j in 1:cons.nb_const
        if !boolB[j]
        	# auxWBj=cons.W[[j],:]*invWB   #W_i W_B^{-1} 
            # v[[j],:]=cons.T[[j],:]-auxWBj*cons.T[B,:] #Compute the definition of v^B_j
            # u[j]=cons.h[j]-dot(auxWBj,cons.h[B]) #Compute the definition of u^B_j
            v[[j],:]=cons.T[[j],:]-cons.W[[j],:]*invWBTB #Compute the definition of v^B_j
            u[j]=cons.h[j]- dot(cons.W[j,:],invWBhB)   #Compute the definition of u^B_j
            if boolI[j]      #If j is in I
                push!(IlessB,j) #Add j to I less B since i is in I and  not in B
            else
                push!(compI,j)  #Add j to the complementary of I
            end
        end
    end
    espq=EspCostInd(cons,cost,I)
    #Call the good oracle with cost I
    α=-transpose(invWBTB)*espq 
    #Implement -bespcost_I* W_B_I^{-1}*T_B_I etc
    β=transpose(invWBhB)*espq
    #Implement bespcost_I* W_B_I^{-1}*h_B_I etc   
    IlessBnotempty=(IlessB!=[])
    if IlessBnotempty
        lI=length(I)
    else
        lI=0
    end
    subdiv=Vector{Vector{IBvu}}(undef,lI)  
    subdivα=Vector{Vector{Float64}}(undef,lI)
    subdivβ=Vector{Float64}(undef,lI)
    subdiv_computed=falses(lI)
    (espq,α,β)=compute_espq_αβ(cons,cost,I,invWBTB,invWBhB)
    return IBvu(I,boolI,B,boolB,IlessB,compI,v,u,invWBTB,invWBhB,espq,α,β,IlessBnotempty,subdiv,subdivα,subdivβ,subdiv_computed)
end

function IBvu(cons::CoupCons,boolI::Vector{Bool},cost::CostDistribution,typecons::DataType=Float64,typecost::DataType=Float64) #Constructor of IBvu
    I=Vector{Int64}(undef,0)
    for i in 1:cons.nb_const
        if boolI[i]
            push!(I,i)
        end
    end
    return IBvu(cons,I,cost,typecons,typecost)
end



function IBvu(cons::CoupCons,I::Vector{Int64},typecons::DataType=Float64,typecost::DataType=Float64) #Constructor of IBvu
    return IBvu(cons,I,ZeroCost(cons.dimy),typecons,typecost)
end





# Subfunction to compute IBvu parameters espq α and β
function compute_espq_αβ(cons::CoupCons,cost::CostDistribution,I::Vector{Int64},invWBTB::Matrix{<:Real},invWBhB::Vector{<:Real})
    espq=EspCostInd(cons,cost,I)
    #Call the good oracle with cost I
    α=-transpose(invWBTB)*espq 
    #Implement -bespcost_I* W_B_I^{-1}*T_B_I etc
    β=transpose(invWBhB)*espq
    #Implement bespcost_I* W_B_I^{-1}*h_B_I etc
    return (espq,α,β)
end






#Function to plot easily a vector of IBvu
function vecI(cI::Vector{IBvu})
	return [cI[i].I for i in 1:length(cI)]
end
function veccompI(cI::Vector{IBvu})
    return [cI[i].compI for i in 1:length(cI)]
end
function vecIlessB(cI::Vector{IBvu})
    return [cI[i].IlessB for i in 1:length(cI)]
end
function vecα(cI::Vector{IBvu})
    return [cI[i].α for i in 1:length(cI)]
end
function vecβ(cI::Vector{IBvu})
    return [cI[i].β for i in 1:length(cI)]
end



#Function to evaluate the value of x directly without simplex
function value_at_x(cons::CoupCons,cost::CostDistribution,x::Vector{<:Real},c::Vector{<:Real})
    cIx=cIBvu(cons,x,cost)
    l=length(cIx)
    α=sum(cIx[i].α for i in 1:l)
    β=sum(cIx[i].β for i in 1:l)
    println(α)
    println(β)
    return dot(α+c,x)+β
end










#Compute ValidityDomain i.e. projection of face of P_I thanks to the v_B_I and u_B_I
# Return a polyhedron from Polyhedral.jl type
function ValidityDomain(basevu::IBvu)
    VD=polyhedron(hrep(basevu.v,basevu.u,BitSet(basevu.I)), CDDLib.Library(:exact))
    # removehredundancy!(VD)
    return VD
end




#Functions without CostDistribution
function cIBvu(cons::CoupCons,cI::Vector{Vector{Int64}},typecons::DataType=Float64,typecost::DataType=Float64)
    cIvu=Vector{IBvu}(undef,0)
    for i in 1:length(cI)
        push!(cIvu,IBvu(cons,cI[i],typecons,typecost))
    end
    return cIvu
end

function cIBvu(cons::CoupCons,x::Vector{<:Real},typecons::DataType=Float64,typecost::DataType=Float64)
    return cIBvu(cons,ActiveConsColl(cons,x),typecons,typecost)
end


#Function with CostDistribution
function cIBvu(cons::CoupCons,x::Vector{<:Real},cost::CostDistribution,typecons::DataType=Float64,typecost::DataType=Float64)
    return cIBvu(cons,ActiveConsColl(cons,x),cost,typecons,typecost)
end

function cIBvu(cons::CoupCons,cI::Vector{Vector{Int64}},cost::CostDistribution,typecons::DataType=Float64,typecost::DataType=Float64)
    cIvu=Vector{IBvu}(undef,0)
    for i in 1:length(cI)
        myIBvu=IBvu(cons,cI[i],cost,typecons,typecost)
        push!(cIvu,myIBvu)
    end
    return cIvu
end


function cIBvu(cons::CoupCons,cboolI::Vector{Vector{Bool}},cost::CostDistribution,typecons::DataType=Float64,typecost::DataType=Float64)
    cIvu=Vector{IBvu}(undef,0)
    for i in 1:length(cboolI)
        myIBvu=IBvu(cons,cboolI[i],cost,typecons,typecost)
        push!(cIvu,myIBvu)
    end
    return cIvu
end










function cIBvu_to_polymakefan(cons::CoupCons,cIBvu::Vector{IBvu})
    return fan.PolyhedralFan(INPUT_RAYS=vcat(zeros(1,cons.dimy),cons.W),INPUT_CONES=vecI(cIBvu))
end

function IBvu_to_polymakecone(cons::CoupCons,myIBvu::IBvu)
    return polytope.Cone(INPUT_RAYS=cons.W[myIBvu.I,:],INPUT_RAY_LABELS=string.(myIBvu.I))
end

function I_to_polymakecone(cons::CoupCons,I::Vector{Int64})
    return polytope.Cone(INPUT_RAYS=cons.W[I,:],INPUT_RAY_LABELS=string.(I))
end







# OLD Function to copy
#  Not necessary now that we have a non mutable Structure
# function copyIBvu(myIBvu::IBvu,withcost::Bool=true)
#     I=copy(myIBvu.I)
#     boolI=copy(myIBvu.boolI)
#     B=copy(myIBvu.B)
#     boolB=copy(myIBvu.boolB)
#     IlessB=copy(myIBvu.IlessB)
#     compI=copy(myIBvu.compI)
#     v=copy(myIBvu.v)
#     u=copy(myIBvu.u)
#     invWBTB=copy(myIBvu.invWBTB)
#     invWBhB=copy(myIBvu.invWBhB)
#     espq=copy(myIBvu.espq)
#     α=copy(myIBvu.α)
#     β=copy(myIBvu.β)
#     myIBvu=IBvu(I,boolI,B,boolB,IlessB,compI,v,u,invWBTB,invWBhB,espq,α,β)
#     return myIBvu
# end

# function copycIBvu(mycIBvu::Vector{IBvu})
#     l=length(mycIBvu)
#     cIvu=Vector{IBvu}(undef,l)
#     for i in 1:l
#         cIvu[i]=copyIBvu(mycIBvu[i])
#     end
#     return cIvu
# end


