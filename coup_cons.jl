
using Polyhedra, Polymake, LinearAlgebra, CDDLib

using JuMP

using Distributions
using GLPK
using InvertedIndices
# using Exceptions

# ../../../../Logiciels_packages/julia-1.5.1/bin/julia


# Library for GLPK and Float64
# solver = GLPK.Optimizer
# lib = DefaultLibrary{Float64}(solver)



function egal(a::Float64,b::Float64)
    return abs(a-b)<=(abs(a)+abs(b)+10^(-10))*10^(-10)
end

function egal(a::Float64,b::Int64)
    return abs(a-b)<=(abs(a)+abs(b)+10^(-10))*10^(-10)
end


function egal(a::BigFloat,b::Float64)
    return abs(a-b)<=(abs(a)+abs(b)+10^(-10))*10^(-10)
end

function egal(a::Rational,b::Rational)
    return  abs(a-b)<=(abs(a)+abs(b)+10^(-10))*10^(-10)
    # a==b
end

function egal(a::BigFloat,b::Rational)
    return abs(a-b)<=(abs(a)+abs(b)+10^(-10))*10^(-10)
end
function egal(a::Float64,b::Rational)
    return abs(a-b)<=(abs(a)+abs(b)+10^(-10))*10^(-10)
end



mutable struct CoupCons{U} #Structure for the coupling linear constraint in inequational form Tx+Wy\leq h or (x,y) \in poly
    dimx::Int64             #Dimension of first state
    dimy::Int64             #Dimension of second state
    nb_const::Int64         #Number of second stage constraint
    T::Matrix{U}           #Constraint matrix of x
    W::Matrix{U}           #Constraint matrix of y
    h::Vector{U}            #Right term in the inequality constraint
    A::Matrix{U}            #First stage constrain matrix
    b::Vector{U}            #First stage right term
    nb_const_1st::Int64     #Number of first stage constraints
end






mutable struct CoupConsPoly{U} #Structure with Polyhedron for the coupling linear constraint in inequational form Tx+Wy\leq h or (x,y) \in poly 
    poly
    dimx::Int64             #Dimension of first state
    dimy::Int64             #Dimension of second state
    nb_const::Int64         #Number of constraint
    T::Array{U,2}           #Constraint matrix of x
    W::Array{U,2}           #Constraint matrix of y
    h::Vector{U}            #Right term in the inequality constraint
    prjx                    #Projection of the polyhedron on x
end






function CoupCons(T::Array{<:Real},W::Array{<:Real},h::Vector{<:Real},A::Matrix{<:Real}=Matrix{Float64}(undef,0,0),b=Vector{Float64}(undef,0)::Vector{<:Real},typecons::DataType=Float64)    #Create with the matrix of constraints T, W and h
    T=typecons.(T)
    W=typecons.(W)
    h=typecons.(h)
    M=hcat(T,W)
    nb_const=length(h)
    nb_const_1st=length(b)
    if nb_const!=length(T[:,1])
        error("T,W and h must have the same height")
    end
    
    dimx=length(T[1,:])
    dimy=length(W[1,:])
    if size(A)==(0,0)
        A=Matrix{Float64}(undef,0,dimx)
    end
    # poly=polyhedron(hrep(M,h), lib)
    # prjx=eliminate(poly,[i for i in dimx+1:1:dimx+dimy])
    # removehredundancy!(prjx)
    # return CoupCons(poly,dimx,dimy,nb_const,T,W,h,prjx)
    return CoupCons(dimx,dimy,nb_const,T,W,h,A,b,nb_const_1st)
end

function CoupCons(poly,dimx::Int64,dimy::Int64,A::Matrix{<:Real}=Matrix{Float64}(undef,0,0),b=Vector{Float64}(undef,0)::Vector{<:Real}) #thanks to a polyhedron object
    if dimx+dimy !=fulldim(poly)
        error("The sum of dimensions must be equal to the dimension of the ambient space of the polyhedron")
    end
    removehredundancy!(poly)
    hr=MixedMatHRep(hrep(poly))  #Taking the H-representation
    T=Rational.(hr.A[:,1:1:dimx])           #T is the first dimx columns
    W=Rational.(hr.A[:,dimx+1:1:dimx+dimy])      #W is the second dimy columns
    h=Rational.(hr.b)
    # prjx=eliminate(poly,[i for i in dimx+1:1:dimx+dimy])
    # removehredundancy!(prjx)
    # return CoupCons(poly,dimx,dimy,length(h),T,W,h,prjx)
    return CoupCons(T,W,h,A,b)
end


function CoupCons(mod::Model,dimx::Int64,dimy::Int64,A::Matrix{<:Real}=Matrix{Float64}(undef,0,0),b=Vector{Float64}(undef,0)::Vector{<:Real}) #Create thanks to a JuMP model
    poly = polyhedron(mod,CDDLib.Library())
    return CoupCons(poly,dimx,dimy,A,b)
end


function perturbate_coupcons_gaussian(cons::CoupCons,stand_dev::Real,typecons::DataType=Float64)
    W=cons.W + stand_dev*randn(cons.nb_const,cons.dimy)
    T=cons.T + stand_dev*randn(cons.nb_const,cons.dimx)
    h=cons.h + stand_dev*randn(cons.nb_const)
    b=cons.b + stand_dev*randn(cons.nb_const_1st)
    A=cons.A + stand_dev*randn(cons.nb_const_1st,cons.dimx)
    return CoupCons(T,W,h,A,b)
end

function perturbate_coupcons_det(cons::CoupCons,pert::Real,typecons::DataType=Float64)
    k=0
    T=Matrix{typecons}(undef,cons.nb_const,cons.dimx)
    W=Matrix{typecons}(undef,cons.nb_const,cons.dimy)
    h=Vector{typecons}(undef,cons.nb_const)
    for i in cons.nb_const
        for j in cons.dimx
            T[i,j]=cons.T[i,j]+k*pert
            k=k+1
        end
    end
    for i in cons.nb_const
        for j in cons.dimy
            W[i,j]=cons.W[i,j]+k*pert
            k=k+1
        end
    end
    for i in cons.nb_const
        h[i]=cons.h[i]+k*pert
        k=k+1
    end
    return CoupCons(T,W,h,typecons)
end




function Fiber(cons::CoupCons,x::Vector)
    # if !in(x,cons.prjx)
    #     error("Empty fiber, the input vector is not in the projection")
    # end
    if cons.dimx!=length(x)
        error("The length of the input vector is different from the first state dimension of the constraint")
    end
    F=polyhedron(hrep(cons.W,cons.h-cons.T*x), CDDLib.Library(:exact))
    # removehredundancy!(F)
    return F
end


function FiberPolymake(cons::CoupCons,x::Vector)
    if cons.dimx!=length(x)
        error("The length of the input vector is different from the first state dimension of the constraint")
    end
    M=hcat(cons.h-cons.T*x,-cons.W)
    return polytope.Polytope(INEQUALITIES=M)
end



#Permit not to recompute the fiber if it is already computed
function ActiveConsColl(cons::CoupCons,x::Vector,F) #Compute the collection of active constraint
    # removevredundancy!(F)
    rep=vrep(F)
    cI=Vector{Vector{Int64}}(undef,0)                               #Collection of maximal active constraints set initialized
    for pi in eachindex(points(rep))    #Iteration on the vertices of the Fiber
        v=get(rep, pi)                  #Take v the current vertex
        I=Vector{Int64}(undef,0)                       #New constraint set initialized
        for i in 1:1:cons.nb_const
            if egal(dot(cons.T[i,:],x)+dot(cons.W[i,:],v),cons.h[i]) #If the constraint is active on this set
                push!(I,i)                           #Add the index to the current constraint set
            end
        end
        if !in(I,cI)
            push!(cI,I)         #Add the new constraint set to the collection
        end
    end
    return cI
end

function PrintVerticesCoupCons(cons::CoupCons)
    removevredundancy!(cons.poly)
    rep=vrep(cons.poly)
    cI=Vector{Vector{Int64}}(undef,0)                               #Collection of maximal active constraints set initialized
    for pi in eachindex(points(rep))    #Iteration on the vertices of the Fiber
        v=get(rep, pi)                  #Take v the current vertex
        I=Vector{Int64}(undef,0)                       #New constraint set initialized
        println(Float64.(v))
    end
end



insert_and_dedup!(v::Vector, x) = (splice!(v, searchsorted(v,x), [x]); v)

function remove_from_sorted!(v::Vector,x::Int64)
    u=length(v)
    l=1
    while l<u-1
        i=(l+u)÷2
        if v[i]<x
            l=i
        elseif v[i]>x
            u=i
        else
            l=i
            u=i
            deleteat!(v,i)
        end
    end
    if v[l]==x
        deleteat!(v,l)
    elseif v[u]==x
        deleteat!(v,u)
    end
         # println("Impossible to delete, value not in the vector")
end

#Recompute the fiber
function ActiveConsColl(cons::CoupCons,x::Vector)
    F=Fiber(cons,x)
    return ActiveConsColl(cons,x,F)
end

function BasisExtraction(A::Array)  #Extract a square inversible submatrix from a full rank rectangular matrix
    d=length(A[1,:])
    if rank(A)!=d
        error("Cannot extract a basis from a non full rank matrix")
    end
    currentAB=zeros(1,d)
    B=Vector{Int64}(undef,0)
    r=0
    i=1
    while r<d
        M=vcat(currentAB,transpose(A[i,:]))
        if rank(M)==r+1
            currentAB=M
            push!(B,i)
            r=r+1
        end
        i=i+1
    end
    return B
end


function BasisExtraction(A::Array,I::Vector) #Extract a square inversible submatrix from a full rank rectangular submatrix A[I,:] of A defined thanks to a set of indices I 
    d=length(A[1,:])
    if Rank(A[I,:])!=d
        error("Cannot extract a basis from a non full rank matrix")
    end
    currentAB=zeros(1,d)
    B=Vector{Int64}(undef,0)
    r=0
    i=1
    while r<d
        M=vcat(currentAB,transpose(A[I[i],:]))
        if Rank(M)==r+1
            currentAB=M
            push!(B,I[i]) #Add the new index of the great matrix
            r=r+1
        end
        i=i+1
    end
    if cond(A[B,:])>10^10
        println("Warning : basis matrix is ill-conditionned")
    end
    return B
end


function Rank(A::Array) #Rank for catching if we have a BigInt Rational matrix
    try 
        return rank(A)
    catch e
        return rank(Float64.(A))
    end
end


function Chamber(cons::CoupCons,x::Vector,cI::Array) #Permit not to recompute the ActiveConsColl if already computed
    q=[i for i in 1:1:cons.nb_const]
    A=transpose(zeros(cons.dimx))
    b=[0]
    for k in 1:1:length(cI)
        I=cI[k]
        B=BasisExtraction(cons.W,I)
        for j in q[Not(I)] #Constraint for the indices j not in the active constraint set I
            v= transpose(cons.W[j,:])*inv(cons.W[B,:])
            A=vcat(A,transpose(cons.T[j,:])-v*cons.T[B,:])
            b=vcat(b,cons.h[j]-v*cons.h[B])
        end
        for i in setdiff(I,B) #Constraint for the indices i in I (if i in B then the cons) #Constraint for the indices i in I (if i in B then the constraint is trivial 0=0)
            v= transpose(cons.W[i,:])*inv(cons.W[B,:])
            a=transpose(cons.T[i,:])-v*cons.T[B,:]
            r=cons.h[i]-v*cons.h[B]
            A=vcat(A,a,-a)
            b=vcat(b,r,-r)
        end
    end
    P = polyhedron(hrep(A,b),CDDLib.Library(:exact))
    # removehredundancy!(P)
    return P
end

function Chamber(cons::CoupCons,x::Vector) #Recompute the ActiveConsColl
    cI=ActiveConsColl(cons,x)
    return Chamber(cons,x,cI)
end



function αβCoef(cons::CoupCons,cost::CostDistribution,x::Vector)
    #For the moment only works for generic x i.e. in the relint of a maximal chamber
    α=Rational.(zeros(cons.dimx))
    β=0
    cI=ActiveConsColl(cons,x)
    for k in 1:1:length(cI)
        I=cI[k]
        espcost=EspCostInd(cons,cost,I)
        minusλI = transpose(inv(cons.W[I,:]))*espcost
        α = α - transpose(cons.T[I,:])*minusλI
        β =β + dot(cons.h[I],minusλI)
    end
    return (α,β)
end



#Iteration of αβCoef on representants
function αβCoefIt(cons::CoupCons,cost::CostDistribution,rep::Array)
    k=length(rep)
    A=Rational.(zeros(k,cons.dimx))
    println(A)
    b=Rational.(zeros(k))
    i=1
    for x in rep
        println(x)
        (α,β)=αβCoef(cons,cost,x)
        println(A[i,:])
        A[i,:]=α
        b[i]=β
        i=i+1
    end
    return (A,b)
end





mutable struct oneSLP
	cost::Distribution
	cons::CoupCons
end





mutable struct MSLP
	horizon::Int64
    multi::Vector{oneSLP}
end




