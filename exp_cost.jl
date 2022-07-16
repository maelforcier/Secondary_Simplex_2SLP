abstract type CostDistribution end

mutable struct ExpCost<:CostDistribution
	cone
	vertex::Vector{<:Real}
	rays::Array{<:Real}
	fulldim::Int64
	dim::Int64
	param::Vector{<:Real}
	norm_factor::Real
end

function ExpCost(cone,param::Vector)
	dim=Polyhedra.dim(cone)
	fulldim=Polyhedra.fulldim(cone)
	#TODO
end

function ExpCost(vertex::Vector,rays::Array)
	fulldim=length(vertex)
	if length(rays[1,:])!=fulldim
		error("The dimension of rays must be equal to the dimension of the vertex")
	end
	#TODO
end