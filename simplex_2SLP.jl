using LinearAlgebra
using JuMP
using Cbc

abstract type CostDistribution end


include("coup_cons.jl")
include("cIBvu.jl")
include("edgetovertex.jl")
include("vertextoedge.jl")
include("deterministic_cost.jl")


#Runner avec julia 1.5 
#../../../../Logiciels_packages/julia-1.5.1/bin/julia




function simplex_2SLP_generic(cons::CoupCons,cost::CostDistribution,xstart::Vector{<:Real},c::Vector{<:Real})
	cIvu=cIBvu(cons,xstart,cost)
	(K,I)=ActiveCons(cons,xstart,zeros(cons.dimy))
	return simplex_2SLP_generic(cons,cost,K,cIvu,xstart,c)
end



function simplex_2SLP_generic(cons::CoupCons,cost::CostDistribution,Kstart::Vector{Int64},cIvu::Vector{IBvu},xstart::Vector{<:Real},c::Vector{<:Real})
	#cons is the constraint structure
	#cost is the structure defining the recours cost distribution
	#cIvu is the vector of IBvu at start
	#xstart is the starting vertex of the simplex
	#c is the first stage cost
	#Kstart is the first stage active constraints set at start
	x=xstart
	K=Kstart
	l=length(cIvu)
	α=sum(cIvu[i].α for i in 1:l)
	β=sum(cIvu[i].β for i in 1:l)
	i=1
	stop=false
	value=dot(c+α,x)+β
	time_vertex_to_edge=0
	time_edge_to_vertex=0
	while !stop
		println(i)
		# println(" ")
		# println("Vertex")
		# println(x)
		# println("cI of vertex and first stage constraints")
		# println(vecI(cIvu))
		# println(K)
		# println("IlessB of vertex")
		# println(vecIlessB(cIvu))
		# println("Slacks Ax-b")
		# println(max.(cons.A*x-cons.b,zeros(cons.nb_const_1st)))
		# println("α")
		# println(α)
		# println("β")
		# println(β)
		println("Current value")
		new_value=dot(c+α,x)+β
		println(new_value)
		if value<new_value-10^-9
			error("The value went higher")
		end
		value=new_value
		# push!(value_it,dot(c+α,x)+β)
		# it_time=@elapsed (d,K,cJvu,αedge,βedge,stop)=vertex_to_edge_greedy(cons,x,K,cIvu,α,β,cost,c)
		it_time=@elapsed (d,K,cJvu,αedge,βedge,stop)=vertex_to_edge_steepest(cons,x,K,cIvu,α,β,cost,c)
		time_vertex_to_edge=time_vertex_to_edge+it_time
		println("Total time vertex to edge")
		println(time_vertex_to_edge)
		if !stop
			# println(" ")
			# println("Direction")
			# println(d)
			# println("cI of edge")
			# println(vecI(cJvu))
			# println("cIlessB of edge")
			# println(vecIlessB(cJvu))
			# println(K)
			# println("α")
			# println(αedge)
			# println("β")
			# println(βedge)
			it_time= @elapsed (x,K,cIvu,α,β,stop)=edge_to_vertex(cons,cost,K,cJvu,x,d,αedge,βedge)
			time_edge_to_vertex=time_edge_to_vertex+it_time
			println("Total time edge to vertex")
			println(time_edge_to_vertex)
			i=i+1
			println(" ")
		end
	end
end




function simplex_2SLP_generic(cons::CoupCons,cost::CostDistribution,c::Vector{<:Real})
	# q=sum(cost.prob[i]*cost.scenarios[i,:] for i in 1:cost.nb_scenarios)
	println("Initialization:")
	if length(c)!=cons.dimx
		error("The length of the cost and the first state dimension must be the same size.")
	end
	if cost.fulldim!=cons.dimy
		error("The recourse cost and the second state must have the same dimension.")
	end
	println("Time spent on the first deterministic linear program")
	time_first_LP = @elapsed (xfl,yfl)=find_coupling_vertex_jump(cons,cost,c)
	println(time_first_LP)
	(K,I)=ActiveCons(cons,xfl,yfl)
	M=vcat(hcat(cons.T,cons.W),hcat(cons.A,zeros(cons.nb_const_1st,cons.dimy)))
	right=vcat(cons.h,cons.b)
	B=vcat(I,cons.nb_const.+K)
	# B=BasisExtraction(A,I) #Not necessary in generic form
	xystart=M[B,:]\right[B]
	xstart=xystart[1:cons.dimx]
 	
 	time_first_cIBvu= @elapsed cIvustart=cIBvu(cons,xstart,cost)
 	 println("Time spent on Initialization of first cIBvu")
 	println(time_first_cIBvu)
 	if length(cIvustart)==0
 		cIvustart=cIBvu(cons,[I],cost)
 	end
 	cIlessB=vecIlessB(cIvustart)
 	println(sum(length(cIlessB[i]) for i in 1:length(cIvustart)))
 	println(" ")
    return simplex_2SLP_generic(cons,cost,K,cIvustart,xstart,c)
end


function find_coupling_vertex_jump(cons::CoupCons,cost::CostDistribution,c::Vector{<:Real})
	q=cost.esp
	model= Model(GLPK.Optimizer)
	@variable(model,x[1:cons.dimx])
	@variable(model,y[1:cons.dimy])
	@constraint(model,cons.T*x+cons.W*y.<=cons.h)
	@constraint(model,cons.A*x.<=cons.b)
	@objective(model,Min,dot(c,x)+dot(q,y))
	# @objective(model,Min,dot(c,x))
	optimize!(model)
	println(primal_status(model))
	xs=value.(x)
	ys=value.(y)
	return (value.(x),value.(y))
end

function ActiveCons(cons::CoupCons,x::Vector{<:Real},y::Vector{<:Real})
	I=Vector{Int64}(undef,0)                   #New constraint set initialized
	K=Vector{Int64}(undef,0)
        for i in 1:cons.nb_const
            if egal(dot(cons.T[i,:],x)+dot(cons.W[i,:],y),cons.h[i]) #If the constraint is active on this set
                push!(I,i)                         #Add the index to the current constraint set
            end
        end
        for k in 1:cons.nb_const_1st
        	if egal(dot(cons.A[k,:],x),cons.b[k]) #If the constraint is active on this set
                push!(K,k)                         #Add the index to the current constraint set
            end
        end
    return (K,I)
end


function ActiveCons(x::Vector{<:Real},A::Matrix{<:Real},b::Vector{<:Real})
	I=Vector{Int64}(undef,0)                   #New constraint set initialized

        for i in 1:length(b)
            if egal(dot(A[i,:],x),b[i]) #If the constraint is active on this set
                push!(I,i)                         #Add the index to the current constraint set
            end
        end
    return I
end




# OLD CODE FOR USING POLYMAKE IN Initialization of simplex

# hxs=cons.h-cons.T*xstart
 	# F=Fiber(cons,xstart)
 	# # println(F.FACETS_THRU_VERTICES)
 	# lcI=length(F.VERTICES[:,2])
 	# vrtx=Vector{Float64}(undef,cons.dimy)
 	# cI=Vector{Vector{Int64}}(undef,lcI)
 	# for i in 1:lcI
 	# 	vrtx[:]=Float64.(F.VERTICES[i,2:cons.dimy+1])
 	# 	I=ActiveCons(vrtx,cons.W,hxs)
 	# 	cI[i]=I
 	# end
 	# llcI=sum(length(cI[i]) for i in 1:lcI)
 	# println(cI)
 	# println(llcI)
 	# println(lcI)
 	# println(llcI-lcI*cons.dimy)
 	# time_first_cIBvu= @elapsed cIvustart=cIBvu(cons,cI,cost)


 	# cboolI=Vector{Vector{Bool}}(undef,lcI)
 	# for i in 1:lcI
 	# 	cboolI[i]=F.FACETS_THRU_VERTICES[i,:]
 	# end
 	# time_first_cIBvu= @elapsed cIvustart=cIBvu(cons,cboolI,cost)