using JuMP

using Distributions


mutable struct MSLP1dim
	horizon::Int64                      #One for each time step
    a::Vector{Vector{Float64}}          #The contraints are of the form x_t a_t +x_{t-1} b_t \leq h_t
    b::Vector{Vector{Float64}}
    h::Vector{Vector{Float64}}
    c::Vector{UnivariateDistribution}   #cost given as a distribution
    m::Vector{Int64}                    #Number of constraints, i.e. size of a_t,b_t and h_t
end





function MSLP1dimConst(hzn::Int64,a::Vector{Float64},b::Vector{Float64},h::Vector{Float64},c::UnivariateDistribution)
    return MSLP1dim(hzn,[a for i in 1:hzn],[b for i in 1:hzn],[h for i in 1:hzn],[c for i in 1:hzn],length(h))
end


mutable struct MSLP1dimUnif
    horizon::Int64                      #One for each time step
    a::Vector{Vector{Float64}}          #The contraints are of the form x_t a_t +x_{t-1} b_t \leq h_t
    b::Vector{Vector{Float64}}
    h::Vector{Vector{Float64}}
    cup::Vector{Float64}                #cost given as a uniform distribution with upper bound
    clow::Vector{Float64}               #and lower bound
    m::Vector{Int64}                    #Number of constraints, i.e. size of a_t,b_t and h_t
end





function MSLP1dimUnifConst(hzn::Int64,a::Vector{Float64},b::Vector{Float64},h::Vector{Float64},cup::Float64,clow::Float64)
    return MSLP1dimUnif(hzn,[a for i in 1:hzn],[b for i in 1:hzn],[h for i in 1:hzn],[cup for i in 1:hzn],[clow for i in 1:hzn],[length(h) for i in 1:hzn])
end