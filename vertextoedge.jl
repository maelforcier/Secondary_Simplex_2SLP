

#TODO : ne pas rechecker inutilement l'arête d'où l'on vient

function vertex_to_edge_greedy(cons::CoupCons,x::Vector{<:Real},K::Vector{Int64},cIvu::Vector{IBvu},α::Vector{<:Real},β::Real,cost::CostDistribution,c::Vector{<:Real},steepest_edge::Bool=true)
	#Inputs:
	# cons : Coupling constraints structure
	# cIvu : structure corresponding to the vertex of the chamber complex
	# grad : actual gradient computed
	# term_grad : list of terms used to compute the gradient, coherent with the indices of cIvu
	# cost : structure which gives a cone valuation oracle to compute expectations of the cost
	#K : first stage active constraints

	M=cons.A[K,:] #We define the matrix of the system of equation, that will give the direction
	# println(size(M))
	for I in 1:length(cIvu)
		# Construct the matrix v_i^B_I for i in I\B_I for I in \cI
		M=vcat(M,cIvu[I].v[cIvu[I].IlessB,:])
	end
	# println(size(M))
	l=1 #Index of row of M to test
	right=Float64.(zeros(cons.dimx))
	for k in K 		#First try to remove the first (easier to do)
		right[l]=-1.
		d=M\right
		if dot(c+α,d)<0
			remove_from_sorted!(K,k)
			return (d,K,cIvu,α,β,false)
		end
		right[l]=0.
		l=l+1 #Test the next index
	end
	for J in 1:length(cIvu)
		if cIvu[J].IlessB!=[]     #Verify that we can remove indices
			lJ=length(cIvu[J].I)
			for indj in 1:lJ
				j=cIvu[J].I[indj]
				if cIvu[J].subdiv_computed[indj] #If subdivision already computed take the values
					cJvunew=cIvu[J].subdiv[indj]
					αedge=α+cIvu[J].subdivα[indj]
					βedge=β+cIvu[J].subdivβ[indj]
					lcJvunew=length(cJvunew)
					Jlessjvu=cJvunew[lcJvunew]
				else    					#Otherwise compute all the subdivision
					Jlessjvu=remove_one_index_from_IBvu(cons,cIvu[J],j,cost) #Compute the IBvu term Jlessj
					vis_indices=visible_facets_input_rays_labels(cons,Jlessjvu,cons.W[j,:])			#Complete the subdivision with visible indices
					cJvunew=cIBvu_edge_new(cons,cost,j,vis_indices)
					#Compute the new Vector{IBvu} corresponding to the new subdivision
					push!(cJvunew,Jlessjvu)
					lengcJnew=length(cJvunew)
					cIvu[J].subdiv[indj]=cJvunew
					αJlessj=sum(cJvunew[i].α for i in 1:lengcJnew)-cIvu[J].α
					cIvu[J].subdivα[indj]=αJlessj
					αedge=α+αJlessj
					βJlessj=sum(cJvunew[i].β for i in 1:lengcJnew)-cIvu[J].β
					cIvu[J].subdivβ[indj]=βJlessj
					βedge=β+βJlessj
					cIvu[J].subdiv_computed[indj]=true
				end
				d=find_direction_J_less_j(cons,cIvu,J,Jlessjvu,j,K)
				if dot(c+αedge,d)<0
					deleteat!(cIvu,J)
					append!(cIvu,cJvunew)
					return (d,K,cIvu,αedge,βedge,false)
				end
				l=l+1 #Test the next index
			end
		end
	end
	println(" ")
	println("Minimum is reached")
	println(" ")
	println("Optimal solution")
	println(x) 
	println(" ")
	# print("Active constraint sets")
	# println(" ")
	# println(vecI(cIvu))
	print("Value")
	println(dot(c+α,x)+β)
	return (zeros(cons.dimx),K,cIvu,α,β,true)
end





function vertex_to_edge_steepest(cons::CoupCons,x::Vector{<:Real},K::Vector{Int64},cIvu::Vector{IBvu},α::Vector{<:Real},β::Real,cost::CostDistribution,c::Vector{<:Real})
	#Inputs:
	# cons : Coupling constraints structure
	# cIvu : structure corresponding to the vertex of the chamber complex
	# grad : actual gradient computed
	# term_grad : list of terms used to compute the gradient, coherent with the indices of cIvu
	# cost : structure which gives a cone valuation oracle to compute expectations of the cost
	#K : first stage active constraints

	M=cons.A[K,:] #We define the matrix of the system of equation, that will give the direction
	# println(size(M))
	for I in 1:length(cIvu)
		# Construct the matrix v_i^B_I for i in I\B_I for I in \cI
		M=vcat(M,cIvu[I].v[cIvu[I].IlessB,:])
	end
	# println(size(M))
	l=1 #Index of row of M to test
	right=Float64.(zeros(cons.dimx))
	min_descent=0.
	lK=length(K)	
	indkmin=-1	
	dirmin=Vector{Float64}(undef,cons.dimx)
	for indk in 1:lK 		#Compute all of them
		right[indk]=-1.
		d=M\right
		dir=d/norm(d)
		descent_factor=dot(c+α,dir)
		if descent_factor<min_descent
			min_descent=descent_factor
			indkmin=indk
			dirmin=dir
		end
			right[indk]=0.
	end
	update_cI=false
	Jmin=-1
	cJvunewmin=[]
	αedgemin=Vector{Float64}(undef,cons.dimx)
	βedgemin=0.
	for J in 1:length(cIvu)
		if cIvu[J].IlessB!=[]     #Verify that we can remove indices
			lJ=length(cIvu[J].I)
			for indj in 1:lJ
				j=cIvu[J].I[indj]
				if cIvu[J].subdiv_computed[indj] #If subdivision already computed take the values
					cJvunew=cIvu[J].subdiv[indj]
					αedge=α+cIvu[J].subdivα[indj]
					βedge=β+cIvu[J].subdivβ[indj]
					lcJvunew=length(cJvunew)
					Jlessjvu=cJvunew[lcJvunew]
				else    					#Otherwise compute all the subdivision
					Jlessjvu=remove_one_index_from_IBvu(cons,cIvu[J],j,cost) #Compute the IBvu term Jlessj
					vis_indices=visible_facets_input_rays_labels(cons,Jlessjvu,cons.W[j,:])			#Complete the subdivision with visible indices
					cJvunew=cIBvu_edge_new(cons,cost,j,vis_indices)
					#Compute the new Vector{IBvu} corresponding to the new subdivision
					push!(cJvunew,Jlessjvu)
					lengcJnew=length(cJvunew)
					cIvu[J].subdiv[indj]=cJvunew
					αJlessj=sum(cJvunew[i].α for i in 1:lengcJnew)-cIvu[J].α
					cIvu[J].subdivα[indj]=αJlessj
					αedge=α+αJlessj
					βJlessj=sum(cJvunew[i].β for i in 1:lengcJnew)-cIvu[J].β
					cIvu[J].subdivβ[indj]=βJlessj
					βedge=β+βJlessj
					cIvu[J].subdiv_computed[indj]=true
					
	
				end
				d=find_direction_J_less_j(cons,cIvu,J,Jlessjvu,j,K)
				dir=d/norm(d)
				descent_factor=dot(c+αedge,dir)
				if descent_factor<min_descent
					update_cI=true
					min_descent=descent_factor
					Jmin=J
					cJvunewmin=cJvunew
					dirmin=dir
					αedgemin=αedge
					βedgemin=βedge
				end
				l=l+1 #Test the next index
			end
		end
	end
	if update_cI				#Check if we choose an edge by updating the second stafe collection cI
		deleteat!(cIvu,Jmin)
		append!(cIvu,cJvunewmin)
		return (dirmin,K,cIvu,αedgemin,βedgemin,false)
	elseif indkmin!=-1   #Otherwise check if we adapt the first stag constraints K
		remove_from_sorted!(K,K[indkmin])
		return (dirmin,K,cIvu,α,β,false)
	else
		#There is no descent direction, we stop the simplex
		println(" ")
		println("Minimum is reached")
		println(" ")
		println("Optimal solution")
		println(x) 
		println(" ")
		print("Active constraint sets")
		println(" ")
		println(vecI(cIvu))
		print("Value")
		println(dot(c+α,x)+β)
		return (zeros(cons.dimx),K,cIvu,α,β,true)
	end
end













function vertex_to_edge(cons::CoupCons,x::Vector{<:Real},cost::CostDistribution,c::Vector{<:Real})
	# Function for debug that recomputes everything from saturated
	cIvu=cIBvu(cons,x)
	abvrtx=compute_gradient_full_partial(cons,c,cost,cIvu)
	return vertex_to_edge(cons,x,cIvu,abvrtx,cost,c)
end


function find_direction_J_less_j(cons::CoupCons,cIvu::Vector{IBvu},J::Int64,Jlessjvu::IBvu,j::Int64,K::Vector{Int64})
	M=Jlessjvu.v[[j],:]
	for I in 1:length(cIvu)
		# Construct the matrix v_i^B_I  for I in \cI
		# and v_j'^B_{J\j} for j' in J\B_{J\j}\{j} 
		if I!=J
			M=vcat(M,cIvu[I].v[cIvu[I].IlessB,:])
		else
			M=vcat(M,Jlessjvu.v[Jlessjvu.IlessB,:])
		end
	end
	M=vcat(M,cons.A[K,:])
		b=vcat(-1,zeros(cons.dimx-1))
		# b=vcat(-1,zeros(length(M[:,1])-1))
		#Solve v_i^B_i d =0 for I in \cI and i in I\B
		#  and v_j'^B_J\j d=0 for j' in J\B_{J\j} \{j'}
		d=M\b #Solve the system to find the direction of the edge
		# println(d)
		# println(cons.W[j,:])
	return d
end

function find_edge_J_less_j(cons::CoupCons,cIvu::Vector{IBvu},cost::CostDistribution,J::Int64,j::Int64,new::Bool,K::Vector{Int64}=Vector{Int64}(undef,0))
	Jlessjvu=remove_one_index_from_IBvu(cons,cIvu[J],j,cost)
	M=Jlessjvu.v[[j],:]
	for I in 1:length(cIvu)
		# Construct the matrix v_i^B_I  for I in \cI
		# and v_j'^B_{J\j} for j' in J\B_{J\j}\{j} 
		if I!=J
			M=vcat(M,cIvu[I].v[cIvu[I].IlessB,:])
		else
			M=vcat(M,Jlessjvu.v[Jlessjvu.IlessB,:])
		end
	end
	M=vcat(M,cons.A[K,:])
		b=vcat(-1,zeros(cons.dimx-1))
		# b=vcat(-1,zeros(length(M[:,1])-1))
		#Solve v_i^B_i d =0 for I in \cI and i in I\B
		#  and v_j'^B_J\j d=0 for j' in J\B_{J\j} \{j'}
		d=M\b #Solve the system to find the direction of the edge
		# println(d)
		# println(cons.W[j,:])
		vis_indices=visible_facets_input_rays_labels(cons,Jlessjvu,cons.W[j,:])
		# println(vis_indices)
		if new 
			cIBvu_edge=cIBvu_edge_new(cons,cost,j,vis_indices)
			push!(cIBvu_edge,Jlessjvu)
		else	
			cIBvu_edge=cIBvu_of_edge(cons,cIvu,J,j,Jlessjvu,vis_indices)
		end
		return(d,cIBvu_edge)
end




function remove_one_index_from_IBvu(cons::CoupCons,myIBvu::IBvu,j::Int64,cost::CostDistribution)
    # This function make a copy a IBvu item by suppressing 
    Ilessj=copy(myIBvu.I)
    remove_from_sorted!(Ilessj,j) # We remove j from the I set and the corresponding boolean items
    boolIlessj=copy(myIBvu.boolI)
    boolIlessj[j]=false 

    if !myIBvu.boolB[j] #If we can keep the same basis
        # We only copy the IBvu item
        IlessjlessB=copy(myIBvu.IlessB)
        remove_from_sorted!(IlessjlessB,j)
        compIlessj=copy(myIBvu.compI)
        insert_and_dedup!(compIlessj,j)
        B=myIBvu.B
        boolB=myIBvu.boolB
        v=myIBvu.v
        u=myIBvu.u
        invWBTB=myIBvu.invWBTB
        invWBhB=myIBvu.invWBhB       
        (espq,α,β)=compute_espq_αβ(cons,cost,Ilessj,invWBTB,invWBhB)
        IlessjlessBnotempty=(IlessjlessB!=[])
    	if IlessjlessBnotempty
        	lIlessj=length(Ilessj)
    	else
        	lIlessj=0
    	end
    	subdiv=Vector{Vector{IBvu}}(undef,lIlessj)  
    	subdivα=Vector{Vector{Float64}}(undef,lIlessj)
    	subdivβ=Vector{Float64}(undef,lIlessj)
    	subdiv_computed=falses(lIlessj)
        #TODO check if j is in the conic hul of Jlessj and then do not recompute espq and thus α and β since they are the same
        Ilessjvu=IBvu(Ilessj,boolIlessj,B,boolB,IlessjlessB,compIlessj,v,u,invWBTB,invWBhB,espq,α,β,IlessjlessBnotempty,subdiv,subdivα,subdivβ,subdiv_computed)
    else
        #  Otherwise, we create a new IBvu item
        #If the list of indices Jlessj is not generating the whole space, 
        # there could be an error from BasisExtraction
        # need to cacth this error in the general case
        Ilessjvu=IBvu(cons,Ilessj,cost) #TODO check if j is in the conic hul of Jlessj and then do not recompute espq since it is the same
    end
    return Ilessjvu
end



function cIBvu_edge_new(cons::CoupCons,cost::CostDistribution,j::Int64,vis_indices::Vector{Vector{Int64}})
	# Compute the new active constraints sets with the visible faces
	l_vis_indices=length(vis_indices)
	cJBvunew=Vector{IBvu}(undef,l_vis_indices)
	for k in 1:l_vis_indices
		v=vis_indices[k]
		insert_and_dedup!(v,j)
		cJBvunew[k]=IBvu(cons,v,cost)
	end
	return cJBvunew
end

function cIBvu_of_edge(cons::CoupCons,cIBvu::Vector{IBvu},J::Int64,j::Int64,Jlessjvu::IBvu,vis_indices::Vector{Vector{Int64}})
	# Compute the subdivision associated to the edge by adding the new and modifying the splitted one
	cJBvu=copy(cIBvu)
	cJBvu[J]=Jlessjvu
	cJBvunew=cIBvu_edge_new(cons,j,vis_indices)
	append!(cJBvu,cJBvunew)
	return cJBvu
end





#Functions to compute the visible faces, with their indices

# Maybe we can improve these functions by keeping in memory the H-representations of cones outside of polymake objects
function visible_facets_input_rays_labels(mycone::Polymake.BigObjectAllocated,x::Vector{<:Real})
	x=Rational.(x) #Polymake only accept rational
	list_facets=polytope.visible_facet_indices(mycone,x)
	# Polymake already have a visible_facet function but only returns the indices of faces
	# We want to have the ste of indices of incident rays
	rays_labels_list=Vector{Vector{Int64}}(undef,0)
	# It is a list of list, because faces are represented by sets of indices
	# Thus a collection of faces is a collection of set of inidices
	# i.e. a list of list
	for f in list_facets
		rays_labels=Vector{Int64}(undef,0)
		for r in 1:mycone.N_INPUT_RAYS
			if mycone.INPUT_RAYS_IN_FACETS[f+1,r]
				insert_and_dedup!(rays_labels,parse(Int64,mycone.INPUT_RAY_LABELS[r]))
				# We use parse because input ray labels must be string
				# parse convert string to int
			end
		end
		append!(rays_labels_list,[rays_labels])
	end
	return rays_labels_list
end

# Convert input of the function to use it in different formats
function visible_facets_input_rays_labels(cons::CoupCons,myIBvu::IBvu,x::Vector{<:Real})
	return visible_facets_input_rays_labels(IBvu_to_polymakecone(cons,myIBvu),x)
end

function visible_facets_input_rays_labels(cons::CoupCons,I::Vector{Int64},x::Vector{<:Real})
	return visible_facets_input_rays_labels(I_to_polymakecone(cons,I),x)
end







function find_all_edges(cons::CoupCons,cost::CostDistribution,cIvu::Vector{IBvu})
	#For a vertex of the chamber complex represented by its collection of sets of active constraints
	#Returns all edges with directions and collections of sets of active constraints
	#For the moment only for generic case
 	sol=[]
	for J in 1:length(cIvu)
		if cIvu[J].IlessB!=[]     #Verify that we can remove indices
			for j in cIvu[J].I
				(d,newcIBvu)=find_edge_J_less_j(cons,cIvu,cost,J,j,false)
				append!(sol,[(d,newcIBvu)])
			end
		end
	end
	return sol
end

function find_all_edges(cons::CoupCons,cost::CostDistribution,x::Vector{<:Real})
	return find_all_edges(cons,cIBvu(cons,x))
end
