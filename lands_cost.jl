
include("discrete_cost.jl")

c_lands=[10., 7., 16., 6.]
f=[40 24 4;45 27 4.5; 32 19.2 3.2; 55 33 5.5]
fvec=[40.,24.,4.,45.,27.,4.5,32.,19.2,3.2,55.,33.,5.5]
m=12.
b=120.
d1_lands=3*0.3+5*0.4+7*0.3
d_lands=[5.,3.,2.]


A_lands=vcat([-1. -1. -1. -1.;10. 7. 16. 6.],-Matrix(I,4,4))
b_lands=vcat([-m,b],zeros(4))

T_lands=vcat(-Matrix(I,4,4),zeros(15,4))

W_lands=vcat([1. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0.;
			  0. 0. 0. 1. 1. 1. 0. 0. 0. 0. 0. 0.;
			  0. 0. 0. 0. 0. 0. 1. 1. 1. 0. 0. 0.;
			  0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 1.]
	,vcat([-1. 0. 0. -1. 0. 0. -1. 0. 0. -1. 0. 0.;
			0. -1. 0. 0. -1. 0. 0. -1. 0. 0. -1. 0.;
			0. 0. -1. 0. 0. -1. 0. 0. -1. 0. 0. -1.]
			,-Matrix(I,12,12)))
h_lands=Float64.(vcat(zeros(4),vcat(-d_lands,zeros(12))))




C_lands=CoupCons(T_lands,W_lands,h_lands,A_lands,b_lands)

C_lands_pert=perturbate_coupcons_gaussian(C_lands,0.01)

cost_lands=DiscreteCost(ones(100)*transpose(fvec)+100*randn(100,12))

