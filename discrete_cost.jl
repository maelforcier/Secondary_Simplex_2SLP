
abstract type CostDistribution end

struct DiscreteCost<:CostDistribution
	scenarios::Array{<:Real}
	prob::Vector{<:Real}
	fulldim::Int64
	nb_scenarios::Int64
	esp::Vector{<:Real}
end

#Constructor with probabilities
function DiscreteCost(scen::Array,prob::Vector)
	d=length(scen[1,:])
	n=length(scen[:,1])
	if length(prob)!=n
		error("Wrong dimensions : The length of the vector of probabilities must be equal to the number of scenarios")
	end
	for i in 1:1:n
		if prob[i]<0
			error("Probilities must be non-negative")
		end
	end
	s = sum(prob[i] for i in 1:1:n)
	if s!=1
		println("WARNING : The probabilities do not sum to one, we rescaled the probas")
		prob=prob/n
	end
	esp=sum(cost.prob[i]*cost.scenarios[i,:] for i in 1:cost.nb_scenarios)/cost.nb_scenarios
	return DiscreteCost(scen,prob,d,n,esp)
end

#Constructor withouth probabilities -> uniform probabilities by default
function DiscreteCost(scen::Array)
	d=length(scen[1,:])
	n=length(scen[:,1])
	prob=ones(n)/n
	esp=sum(scen[i,:] for i in 1:n)/n
	return DiscreteCost(scen,prob,d,n,esp)
end




#Method of IntegrateCostPoly for discrete cost, returns ∫_relint(poly) c d\PP(c)
function IntegrateCostPoly(cost::DiscreteCost,poly) 
	# removehredundancy!(poly) #Takes too long!
	eci=zeros(cost.fulldim)
	for i in 1:1:cost.nb_scenarios
		#if in(cost.scenarios[i,:],poly) #Check the i-th cost scenario is in poly #use ininterior or inrelativeinterior if we want to generalize 
		if inrelativeinterior(cost.scenarios[i,:],poly)
			eci=eci + cost.prob[i]*cost.scenarios[i,:]
		end
	end
	return eci
end




#Method of EspCostInd for uniform cost
function EspCostInd(cons::CoupCons,cost::DiscreteCost,indexset::Array) #Input a coupling constraint object, the polyhedron where c has an uniform distribution on, the basis I where we want to compute E(c 1 c \in W_I^\top)
	if cost.fulldim!=cons.dimy
		error("The ambient space of the cost and the recourse must be of same dimension")
	end
	coneI=polyhedron(vrep(zeros(1,cons.dimy),-cons.W[indexset,:]),CDDLib.Library(:exact)) #We put a minus before W because we minimize
	return IntegrateCostPoly(cost,coneI)
end


scens=[0 1;1 0.1;0 -1; -1 0]

disccost=DiscreteCost(scens)