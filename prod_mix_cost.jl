
include("discrete_cost.jl")

h_prod_mix=[6000., 4000.]
T_prod_mix=[4.0  9.0  7.0  10.0; 1.0  1.0  3.0  40.0]

c_prod_mix=[-12., -20., -18., -40.]
q_prod_mix=[5.,10.]

W_prod_mix=[-1. 0.; 0. -1.]

A_prod_mix=Float64.(-Matrix(I,4,4))
b_prod_mix=Float64.(zeros(4))

hbig_prod_mix=vcat(h_prod_mix,zeros(2))
Tbig_prod_mix=vcat(vcat(T_prod_mix,zeros(2,4)))
Wbig_prod_mix=vcat(vcat(W_prod_mix,-Matrix(I,2,2)))
println(size(Tbig_prod_mix))
println(size(Wbig_prod_mix))
println(size(hbig_prod_mix))
println(hbig_prod_mix)



cost_prod_mix=DiscreteCost(transpose(q_prod_mix).+1000*randn(1000,2))


Cprod_mix=CoupCons(Tbig_prod_mix,Wbig_prod_mix,hbig_prod_mix,A_prod_mix,b_prod_mix)
prod_mix_pert=perturbate_coupcons_gaussian(Cprod_mix,0.01)