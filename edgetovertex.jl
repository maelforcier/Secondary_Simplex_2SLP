



function edge_to_vertex(cons::CoupCons,cost::CostDistribution,K::Vector{Int64},cIvu::Vector{IBvu},x::Vector{<:Real},d::Vector{<:Real},α::Vector{<:Real},β::Real)
	#Given a collection of active constraints sets, a point x and a direction d
	#Returns the lambda max and the new saturated constraints sets
	lambda=-1 #Represents - infty
	sat=Vector{Pair}(undef,0)
	lcI=length(cIvu)
	ksat=0
	for k in 1:cons.nb_const_1st
		if !in(k,K)
			Akd=dot(cons.A[k,:],d)
			if Akd>0
				# println("slack")
				slack=cons.b[k]-dot(cons.A[k,:],x)
				if !egal(slack,0.)
					temp=slack/Akd
					if temp<lambda||(lambda==-1)
						# println("temp")
						# println(temp)
						lambda=temp
						ksat=k
					end
				end
			end
		end
	end
	for I in 1:length(cIvu)
		jsat=Vector{Int64}(undef,0)
		for j in cIvu[I].compI
			vjd=dot(cIvu[I].v[j,:],d)
			if vjd>0
				temp=(cIvu[I].u[j]-dot(cIvu[I].v[j,:],x))/vjd
				if egal(temp,lambda) #Alternative check for float temp==lambda
					# Important to check the equality first because it is an approximate equality
					push!(jsat,j)
				elseif temp<lambda||(lambda==-1)
					# Here it is a strict inequality and not approximate, but one it is approximately equal it is never checked.
					sat=Vector{Pair}(undef,0)
					jsat=[j]
					lambda=temp
				end
			end
			# jsat is the list of j in compI such that we have saturation at constraints I,j
		end
		if length(jsat)!=0
			push!(sat,Pair(I,jsat))
		end
	end
	if lambda==-1
		println(" ")
		println("Infinite LP")
		println(" ")
		println("Direction of ray")
		println(d)
		println(" ")
		# println("Gradient on its ray")
		# println(c+α)
		return (x,K,cIvu,α,β,true)
	end
	# println("λ")
	# println(lambda)
	if length(sat)!=0 
		(cJBvu,αvrtx,βvrtx)=gather_saturated_constraints(cons,cIvu,cost,sat,α,β)
		# println("SAT")
		# println(sat)
		return (x+lambda*d,K,cJBvu,αvrtx,βvrtx,false)
	else
		# println(ksat)
		insert_and_dedup!(K,ksat)
		return (x+lambda*d,K,cIvu,α,β,false)
	end
end



function edge_to_vertex(cons::CoupCons,cI::Vector{Vector{Int64}},x::Vector{<:Real},d::Vector{<:Real},α::Vector{:Real},β::Real)
	# Function to debug:
	# in the simplex do not call directly cIBvu
	# but update locally the cIBvu
	return edgetovertex(cIBvu(cons,cI),x,d,α,β)
end


function edge_to_vertex(cons::CoupCons,x::Vector{<:Real},d::Vector{<:Real},α::Vector{:Real},β::Real)
	# Function to debug:
	# in the simplex do not call directly ActiveConsColl
	# but update locally the cI
	cI=ActiveConsColl(cons,x+10^-4*d)
	return edge_to_vertex(cIBvu(cons,cI),x,d,α,β)
end





# TODO Improve the way we compute expectation
# TODO? rethink the gathering cleverly
function gather_saturated_constraints(cons::CoupCons,cIvu::Vector{IBvu},cost::CostDistribution,sat::Vector{Pair},α::Vector{<:Real},β::Real)
    lcI=length(cIvu)
    cJBvu=copy(cIvu)  
	cboolJ=Vector{BitVector}(undef,0)
	 #Get rid of old saturated constraint i.e. I_old
	 #Then update the news in the IBvu object
	cJ=Vector{Vector{Int64}}(undef,0)
	lengthcJ=0
	for k in 1:length(sat)
		I=sat[k].first
		J=copy(cIvu[I].I)
		boolJ=copy(cIvu[I].boolI)
		JlessB=copy(cIvu[I].IlessB)
		compJ=copy(cIvu[I].compI)
		for j in sat[k].second
			# I update all the constraints items
			insert_and_dedup!(J,j)
			insert_and_dedup!(JlessB,j)
			remove_from_sorted!(compJ,j)
			#But not B, boolB, u and v since the same base works
		end
		# JnotincludedincJ=true
		# l=1
		# while l<=lengthcJ #Check if J was not added !in(J,cJ) and if J is not included in  another element of cJ
		# 	JnotincludedincJ=JnotincludedincJ&&!(boolJ<=cboolJ[l]) #J must not be included in J' for all J' in cJ
		# 	# We check if J \notin J'=cJ[l] i.e. the boolean is not inf or equal
		# 	if boolJ>cboolJ[l]  #If J strictly contains another set
		# 		for j in sat[k].second
		# 			insert_and_dedup!(J,j)
		# 			boolJ[j]=true
		# 		end
		# 		cJBvu[lcI+l].I=J  #We add the forgotten index
		# 		cboolJ[l]=boolJ
		# 		JnotincludedincJ=false
		# 	end
		# 	l=l+1
		# end
		# if JnotincludedincJ #Check if J was not added before 
		if !in(J,cJ) #Check if J was not added before 
			BJ=cIvu[I].B
			boolBJ=cIvu[I].boolB
			vJ=cIvu[I].v
			uJ=cIvu[I].u
			invWBTB=cIvu[I].invWBTB
			invWBhB=cIvu[I].invWBhB
			(espq,αJ,βJ)=compute_espq_αβ(cons,cost,J,invWBTB,invWBhB)
			α=α+αJ
			β=β+βJ
			# push!(cboolJ,boolJ)
			push!(cJ,J)
			# lengthcJ=lengthcJ+1
			#OLD WITHOUT REMEMBERING OF SUBDIVISION push!(cJBvu,IBvu(J,boolJ,BJ,boolBJ,JlessB,compJ,vJ,uJ,invWBTB,invWBhB,espq,αJ,βJ))
			#Remembering subdivision but naive
			#Does not add the subdivision we came from
			#TODO ADD THE SUBDIVISION WE CAME FROM
			lJ=length(J)
			subdiv=Vector{Vector{IBvu}}(undef,lJ)  
    		subdivα=Vector{Vector{Float64}}(undef,lJ)
    		subdivβ=Vector{Float64}(undef,lJ)
    		subdiv_computed=falses(lJ)
			push!(cJBvu,IBvu(J,boolJ,BJ,boolBJ,JlessB,compJ,vJ,uJ,invWBTB,invWBhB,espq,αJ,βJ,true,subdiv,subdivα,subdivβ,subdiv_computed)) 
		end
	end
	for k in 1:length(sat)
		α=α-cJBvu[sat[k].first].α
		β=β-cJBvu[sat[k].first].β
	end
	deleteat!(cJBvu,[sat[k].first for k in 1:length(sat)])
	return (cJBvu,α,β)
end