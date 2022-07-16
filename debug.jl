

function IBvu_coherent(myIBvu::IBvu)
	if !(myIBvu.boolB<=myIBvu.boolI)
		error("BoolB is not <= to BoolI")
	end
	for i in myIBvu.I
		if !myIBvu.boolI[i]
			error("I and boolI disagree on index "*string(i))
		end
	end
	for i in myIBvu.compI
		if myIBvu.boolI[i]
			error("compI and boolI disagree on index "*string(i))
		end
	end
	for i in myIBvu.B
		if !myIBvu.boolB[i]
			error("B and boolB disagree on index "*string(i))
		end
	end
	for i in myIBvu.IlessB
		if myIBvu.boolB[i]
			error("IlessB and boolB disagree on index "*string(i))
		end
	end
end

function cIBvu_coherent(cIvu::Vector{IBvu})
	for I in 1:length(cIvu)
		println("Checking term at index "*string(I))
		IBvu_well_defined(cIvu[I])
	end
end


# function αβ_coherent(ab::αβ)
# 	if length(ab.listα)!=ab.n_terms
# 		error("Length of listα different from n_terms")
# 	end
# 	if length(ab.listβ)!=ab.n_terms
# 		error("Length of listβ different from n_terms")
# 	end
# 	if length(ab.listβ)!=length(ab.listα)
# 		error("Length of listα different from length of listβ")
# 	end
# 	println(sum(ab.listα))
# 	println(ab.α)
# 	println(sum(ab.listβ))
# 	println(ab.β)
# end



function check_edgetovertex(cons::CoupCons,x::Vector{<:Real},d::Vector{<:Real},cost::CostDistribution,mytype::DataType=Float64)
	# Compare \cI(W,h-T(x+λd)) computed directly versus updated from the computation of x
	# zerocost=ZeroCost(cons.dimy)
	xe=x+10^-4*d
	cIedge=cIBvu(cons,xe,cost)
	αedge=sum(cIedge[i].α for i in 1:length(cIedge))
	βedge=sum(cIedge[i].β for i in 1:length(cIedge))
	(y,cJupdate,α,β)=edge_to_vertex(cons,cost,cIedge,x,d,αedge,βedge)
	println(" ")
	println("New point")
	println(y)
	cJdirect=cIBvu(cons,y,cost)
	println(" ")
	println("cI of edge")
	println(vecI(cIedge))
	println(" ")
	println("cI of vertex")
	println(" by direct computation")
	println(vecI(cJdirect))
	println(" by update")
	println(vecI(cJupdate))
	println(" ")
	# println("Length of cI vertex")
	# println(" by direct computation")
	# println(length(cJdirect))
	# println(" by update")
	# println(length(cJupdate))
	# println(" ")
	println("α")
	println(" by direct computation")
	αdirect=sum(cJdirect[i].α for i in 1:length(cJdirect))
	println(αdirect)
	println(" by update")
	αupdate=sum(cJupdate[i].α for i in 1:length(cJupdate))
	println(αupdate)
	println("β")
	println(" by direct computation")
	βdirect=sum(cJdirect[i].β for i in 1:length(cJdirect))
	println(βdirect)
	println(" by update")
	βupdate=sum(cJupdate[i].β for i in 1:length(cJupdate))
	println(βupdate)
	println(" ")
	println("α list")
	println("edge")
	println(vecα(cIedge))
	println(" by direct computation")
	println(vecα(cJdirect))
	println(" by update")
	println(vecα(cJupdate))
	println(" ")
	println("β list")
	# println("edge")
	# println(vecβ(cIedge))
	println(" by direct computation")
	println(vecβ(cJdirect))
	println(" by update")
	println(vecβ(cJupdate))
	println(" ")
	println("Value one point of vertex")
	println(" by direct computation")
	println(dot(αdirect,y)+βdirect)
	println(" by update")
	println(dot(αupdate,y)+βupdate)
	return (y,cJdirect,cJupdate,cIedge);
end






#TODO REWRITE WITH THE NEW IBvu including α β and espq
function check_vertextoedge(cons::CoupCons,x::Vector{<:Real},J::Int64,j::Int64,cost::CostDistribution,mytype::DataType=Float64)

	# Compare α and β at the edge computed directly or via the update of vertex to edge
	# J is the index of active constraint we want to modify
	# j is the index in J we want to "delete/split" into J\{j} U {j}Uvis_face
	
	# Compute the αβ for the "old" vertex
	cIvuvrtx=cIBvu(cons,x,cost)
	αβvrtx=compute_αβ(cons,cost,cIvuvrtx,mytype)

	#Compute the update for the edge 
	# Piece of code we cant to debug for the real vertex to edge function
	(d,cJvunew)=find_edge_J_less_j(cons,cIvuvrtx,J,j,true)
	Jlessjvu=remove_one_index_from_IBvu(cons,cIvuvrtx[J],j)
	push!(cJvunew,Jlessjvu)
	αβupdate=update_αβ_vertextoedge(C22,αβvrtx,cost,cIvuvrtx,[J],cJvunew)
	cJvu_update=copycIBvu(cIvuvrtx)
	deleteat!(cJvu_update,J)
	append!(cJvu_update,cJvunew)

	# Compute directly for the edge
	cJvu_direct=cIBvu(cons,x+10^-4*d)
	αβdir=compute_αβ(cons,cost,cJvu_direct,mytype)

	println(" ")
	println("Direction")
	println(d)
	println(" ")
	println("cI of edge")
	println(" by direct computation")
	println(vecI(cJvu_direct))
	println(" by update")
	println(vecI(cJvu_update))
	println(" ")
	println("α by direct computation")
	println(αβdir.α)
	println("α by update")
	println(αβupdate.α)
	println(" ")
	println("β by direct computation")
	println(αβdir.β)
	println("β by update")
	println(αβupdate.β)
	println(" ")
	println("Length of cI edge")
	println(" by direct computation")
	println(length(cJvu_direct))
	println(" by update")
	println(length(cJvu_update))
	println(" ")
	println("Length term α by direct computation")
	println(length(αβdir.listα))
	println("Length term α by update")
	println(length(αβupdate.listα))
	println(" ")
	println("Length term β by direct computation")
	println(length(αβdir.listβ))
	println("Length term β by update")
	println(length(αβupdate.listβ))
	println(" ")
	println("List of term α by direct computation")
	println(αβdir.listα)
	println("List of term α by update")
	println(αβupdate.listα)
	println(" ")
	println("List of term β by direct computation")
	println(αβdir.listβ)
	println("List of term β by update")
	println(αβupdate.listβ)
end


