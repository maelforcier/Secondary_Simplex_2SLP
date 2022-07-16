using Plots
using MeshCat

function PlotScatterPoly(poly)
    plot!(poly,aspect_ratio=:1)
    scatter!([x[1] for x in points(poly)], [x[2] for x in points(poly)])
end

function PlotPoly3d(poly)
	if fulldim(poly)!=3
		error("The Polyhedra is not of dimension 3")
	end
	vis=Visualizer()
	m=Polyhedra.Mesh(F)
	setobject!(vis,m)
	open(vis)
end


#Before using this function, run plotChamber=plot() first 
function PlotChamberFiber(cons::CoupCons,x::Array)
	F=Fiber(cons,x)
 	C=Chamber(cons,x)
 	println(ActiveConsColl(cons,x))
	plotFiber=plot(F,aspect_ratio=:1);
 	scatter!([x[1] for x in points(F)], [x[2] for x in points(F)]); #add to the internal variable plotFiber

 	plot(plotChamber,plotFiber,layout=(1,2),legend=false)
	plot!(C,aspect_ratio=:1)					#add only to the external variable plotChamber
 	#Since when there is a layout the adding of plots only apply to the first subplot
 	scatter!([x[1] for x in points(C)], [x[2] for x in points(C)])
end

#Iterated on representants
function PlotChamberFiberIt(cons::CoupCons,rep::Array)
	for x in rep
		PlotChamberFiber(cons,x)
	end
end


function PlotChamber(cons::CoupCons,x::Vector)
	PlotScatterPoly(Chamber(cons,x))
end

function PlotChamberIt(cons::CoupCons,rep::Vector)
	for x in rep
		PlotChamber(cons,x)
	end
end


function PlotValidityDomain(Ivu::IBvu)
	PlotScatterPoly(ValidityDomain(Ivu))
end

function PlotValidityDomain(cons::CoupCons,I::Vector{Int64})
	PlotValidityDomain(IBvu(cons,I))
end

