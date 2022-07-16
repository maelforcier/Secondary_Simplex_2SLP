
abstract type CostDistribution end

mutable struct UniformCost<:CostDistribution
	supp
	vol::Real
	dim::Int64
	fulldim::Int64
end

function UniformCost(supp)
	vol=volume(supp)
	dim=Polyhedra.dim(supp)
	fulldim=Polyhedra.fulldim(supp)
	return UniformCost(supp,vol,dim,fulldim)
end




function UnscaledCentroid(poly) #Returns∫x dx integrated on the polyhedron P with dx the Lebesgue measure
	d=Polyhedra.fulldim(poly)
	Δ=polyhedron.(triangulation(poly),CDDLib.Library(:exact)) #We do the ponderated sum of centroid of the triangulation
	centr=Rational.(zeros(d))
	for i in 1:1:length(Δ)					#For all simplices in the triangulation
		smplxcentr=Rational.(zeros(d))
		for p in points(Δ[i])				#Compute the centroid of the simplex
			smplxcentr=smplxcentr+p
		end
		smplxcentr=smplxcentr/(d+1)
		volsmplx=volume_simplex(Δ[i])
		centr = centr + volsmplx*smplxcentr
	end
	return centr   		#For the real (or scaled) centroid, we should divide by the volume of the input polytope
						#but it is unecessary because it simplify in the computation of Ecind
end

#Method of IntegrateCostPoly for uniform cost, returns ∫_poly c d\PP(c)
function IntegrateCostPoly(cost::UniformCost,poly) 
	if cost.dim==cost.fulldim
		polyintersupp=intersect(poly,cost.supp)
		removevredundancy!(polyintersupp)
		removehredundancy!(polyintersupp)
		return UnscaledCentroid(polyintersupp)//cost.vol
	else
		error("Non full-dimensional support not coded yet")
	end
end




#Method of EspCostInd for uniform cost
function EspCostInd(cons::CoupCons,cost::UniformCost,indexset::Array) #Input a coupling constraint object, the polyhedron where c has an uniform distribution on, the basis I where we want to compute E(c 1 c \in W_I^\top)
	if cost.fulldim!=cons.dimy
		error("The ambient space of the cost and the recourse must be of same dimension")
	end
	coneI=polyhedron(vrep([0 0],-cons.W[indexset,:]),CDDLib.Library(:exact)) #We put a minus before W because we minimize
	return IntegrateCostPoly(cost,coneI)
end



rayon=1

ballinfty=polyhedron(vrep([rayon rayon;rayon -rayon;-rayon rayon;-rayon -rayon]),CDDLib.Library(:exact))
costballinfty=UniformCost(ballinfty)
ballone=polyhedron(vrep([rayon 0;0 rayon;-rayon 0;0 -rayon]),CDDLib.Library(:exact))
costballone=UniformCost(ballone)
