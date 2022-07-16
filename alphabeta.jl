
# OLD ALPHA BETA STRUCTURE NOW INTEGRATED IN EACH IBvu

struct αβ{U}
	α::Vector{U}
	β::U
	listα::Vector{Vector{U}}
	listβ::Vector{U}
	n_terms::Int64
end






function compute_term_αβ(cons::CoupCons,cost::CostDistribution,myIBvu::IBvu)
	bespcostI=EspCostInd(cons,cost,myIBvu.I)
	#Call the good oracle with cost I
	term_α=-transpose(myIBvu.invWBTB)*bespcostI 
	#Implement -bespcost_I* W_B_I^{-1}*T_B_I etc
	term_β=transpose(myIBvu.invWBhB)*bespcostI
	#Implement bespcost_I* W_B_I^{-1}*h_B_I etc
	return (term_α,term_β)
end




function compute_αβ(cons::CoupCons,cost::CostDistribution,cIvu::Vector{IBvu},mytype::DataType=Float64)
	#This compute a gradient and the list of terms by I in cI from scratch
	#c represents the first stage cost
	l=length(cIvu)
	list_term_α=Vector{Vector{mytype}}(undef,l)
	list_term_β=Vector{mytype}(undef,l)
	α=mytype.(zeros(cons.dimx))
	β=mytype(0)
	for I in 1:l
		(term_α,term_β)=compute_term_αβ(cons,cost,cIvu[I])		
		list_term_α[I]=term_α
		list_term_β[I]=term_β
		α=α+term_α
		β=β+term_β
	end
	return αβ(α,β,list_term_α,list_term_β,l)
end

function compute_αβ(cons::CoupCons,cost::CostDistribution,x::Vector{<:Real},mytype::DataType=Float64)
	return compute_αβ(cons,cost,cIBvu(cons,x),mytype)
end


function update_αβ_vertextoedge(cons::CoupCons,αβold::αβ,cost::CostDistribution,cIvuold::Vector{IBvu},list_old::Vector{Int64},cJvunew::Vector{IBvu})
	#Compute the gradient by changing locally some active constraint sets I
	#Allow not to recompute every expectation from scratch
	#cons :: Constraint structure
	#αβold :: structure containing the information on gradient and ordonnee a l'origine
	#list_old :: Indices of active constrain sets that will be removed
	#cJvunew :: New IBvu to be added
	α=copy(αβold.α)
	β=copy(αβold.β)
	for I in list_old #Indices of active constrain sets that will be removed
		α=α-αβold.listα[I] #Delete the old contributions from total gradient
		β=β-αβold.listβ[I]
	end
	listα_new=copy(αβold.listα)
	listβ_new=copy(αβold.listβ)
	deleteat!(listα_new,list_old) #and from the list
	deleteat!(listβ_new,list_old)
	n_terms_new=length(cJvunew)
	for J in 1:n_terms_new
		#Compute the new contribution
		(term_α,term_β)=compute_term_αβ(cons,cost,cJvunew[J])
		α=α+term_α #add to the new contribution
		β=β+term_β 
		push!(listα_new,term_α) #and the list α
		push!(listβ_new,term_β) #and the list β
	end
	n_terms=αβold.n_terms+n_terms_new-length(list_old) #Compute the number of terms of the new αβ item
	return αβ(α,β,listα_new,listβ_new,n_terms)
end




# OLD FUNCTIONS WHEN IBvu had not α and β inside

function compute_espq_αβ(cons::CoupCons,cost::CostDistribution,myIBvu::IBvu)
    espq=EspCostInd(cons,cost,myIBvu.I)
    #Call the good oracle with cost I
    α=-transpose(myIBvu.invWBTB)*espq 
    #Implement -bespcost_I* W_B_I^{-1}*T_B_I etc
    β=transpose(myIBvu.invWBhB)*espq
    #Implement bespcost_I* W_B_I^{-1}*h_B_I etc
    return (espq,α,β)
end


function IBvu(cons::CoupCons,I::Vector{Int64},cost::CostDistribution,typecons::DataType=Float64,typecost::DataType=Float64) #Constructor of IBvu with cost distribution
    myIBvu=IBvu(cons,I,typecons,typecost)
    (espq,α,β)=compute_espq_αβ(cons,cost,myIBvu)
    myIBvu.espq=espq
    myIBvu.α=α
    myIBvu.β=β
    return myIBvu
end






# OLD FUNCTIONS WITH ONLY THE TERM GRADS PREFER ALPHA BETA NOW

function compute_term_grad(cons::CoupCons,cost::CostDistribution,myIBvu::IBvu)
	bespcostI=EspCostInd(cons,cost,myIBvu.I)
	#Call the good oracle with cost I
	term_grad=transpose(myIBvu.invWBTB)*bespcostI 
	#Implement bespcost_I* W_B_I^{-1}*T_B_I etc
	return term_grad
end

function compute_gradient_full_partial(cons::CoupCons,c::Vector{<:Real},cost::CostDistribution,cIvu::Vector{IBvu})
	#This compute a gradient and the list of terms by I in cI from scratch
	#c represents the first stage cost
	l=length(cIvu)
	list_term_grad=Vector{Vector{<:Real}}(undef,l)
	grad=c
	for I in 1:l
		term_grad=compute_term_grad(cons,cost,cIvu[I])		
		list_term_grad[I]=term_grad
		grad=grad-term_grad
	end
	return (grad,list_term_grad)
end

function compute_gradient_full_partial(cons::CoupCons,c::Vector{<:Real},cost::CostDistribution,x::Vector{<:Real})
	return compute_gradient_full_partial(cons,c,cost,cIBvu(cons,x))
end


function update_gradient(cons::CoupCons,grad::Vector{<:Real},cost::CostDistribution,cIvuold::Vector{IBvu},list_old::Vector{Int64},term_grad_old::Vector{<:Vector{<:Real}},cJvunew::Vector{IBvu})
	#Compute the gradient by changing locally some active constraint sets I
	#Allow not to recompute every expectation from scratch
	for I in list_old #Indices of active constrain sets that will be removed
		grad=grad+term_grad_old[I] #Delete the old contributions from total gradient
	end
	term_grad_new=copy(term_grad_old)
	deleteat!(term_grad_new,list_old) #and from the list
	for J in 1:length(cJvunew)
		#Compute the new contribution
		term_grad=compute_term_grad(cons,cost,cJvunew[J])
		grad=grad-term_grad #add to the total gradient 
		push!(term_grad_new,term_grad) #and the list
	end
	return(grad,term_grad_new)
end