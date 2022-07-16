
abstract type CostDistribution end

struct DeterministicCost<:CostDistribution
	value::Vector{<:Real}
	fulldim::Int64
end

#Constructor with probabilities
function DeterministicCost(value::Vector{<:Real})
	return DiscreteCost(value,length(value))
end



#Method of IntegrateCostPoly for deterministic cost, returns ∫_relint(poly) c d\PP(c)
function IntegrateCostPoly(cost::DeterministicCost,poly) 
	removehredundancy!(poly)
	if inrelativeinterior(cost.scenarios[i,:],poly)
		return(value)
	else
		return zeros(cost.fulldim)
	end
end


#Method of EspCostInd for deterministic cost
function EspCostInd(cons::CoupCons,cost::DeterministicCost,indexset::Array) #Input a coupling constraint object, the polyhedron where c has an uniform distribution on, the basis I where we want to compute E(c 1 c \in W_I^\top)
	if cost.fulldim!=cons.dimy
		error("The ambient space of the cost and the recourse must be of same dimension")
	end
	coneI=polyhedron(vrep([0 0],-cons.W[indexset,:]),CDDLib.Library(:exact)) #We put a minus before W because we minimize
	return IntegrateCostPoly(cost,coneI)
end



struct ZeroCost<:CostDistribution
	fulldim::Int64
end





#Method of IntegrateCostPoly for deterministic cost, returns ∫_relint(poly) c d\PP(c)
function IntegrateCostPoly(cost::ZeroCost,poly) 
	return zeros(cost.fulldim)
end


#Method of EspCostInd for deterministic cost
function EspCostInd(cons::CoupCons,cost::ZeroCost,indexset::Array) #Input a coupling constraint object, the polyhedron where c has an uniform distribution on, the basis I where we want to compute E(c 1 c \in W_I^\top)
	return zeros(cost.fulldim)
end


