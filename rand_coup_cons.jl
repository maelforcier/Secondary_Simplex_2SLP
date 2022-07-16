include("discrete_cost.jl")


function RandPolyNormalVrep(d::Int64,q::Int64)
    return polyhedron(vrep(round.(20*rand(Normal(), q,d)+1*ones(q,d))),lib)
end

function RandPolyNormalHrep(d::Int64,q::Int64)
    return polyhedron(vrep(round.(20*rand(Normal(), q,d)+1*ones(q,d))),lib)
end

function RandCoupConsNormalVrep(dimx::Int64,dimy::Int64,q::Int64)
	return CoupCons(RandPolyNormal(dimx+dimy,q),dimx,dimy)
end


function RandPolyNormalHrep(d::Int64,q::Int64)
    return polyhedron(hrep(20*rand(Normal(), q,d),ones(q)),lib)
end

function RandCoupConsNormalHrep(dimx::Int64,dimy::Int64,nb_const::Int64,stand_dev::Float64)
	Trand=stand_dev*rand(Normal(), nb_const,dimx)
	Wrand=stand_dev*rand(Normal(), nb_const,dimy)
	hones=10*stand_dev*ones(nb_const)
	return CoupCons(Trand,Wrand,hones)
end

function RandCoupConsNormalHrep(dimx::Int64,dimy::Int64,nb_const::Int64,nb_const_1st::Int64,stand_dev::Float64)
	Trand=stand_dev*rand(Normal(),nb_const ,dimx)
	Wrand=stand_dev*rand(Normal(),nb_const ,dimy)
	hones=10*stand_dev*ones(nb_const)
	Arand=stand_dev*rand(Normal(),nb_const_1st ,dimx)
	bones=10*stand_dev*ones(nb_const_1st)
	return CoupCons(Trand,Wrand,hones,Arand,bones)
end


function saa_square_cost(dimy::Int64,nb_scen::Int64)
	return DiscreteCost(2*rand(nb_scen,dimy).-1)
end

function saa_gaussian_cost(dimy::Int64,nb_scen::Int64,stand_dev::Float64)
	return DiscreteCost(stand_dev*randn(nb_scen,dimy))
end

squarecost2=DiscreteCost(2*rand(100,2).-1)
squarecost3=DiscreteCost(2*rand(100,3).-1)
squarecost4=DiscreteCost(2*rand(100,4).-1)
squarecost5=DiscreteCost(2*rand(100,5).-1)
squarecost6=DiscreteCost(2*rand(100,6).-1)
squarecost7=DiscreteCost(2*rand(100,7).-1)

cost2=saa_gaussian_cost(2,100,0.1)
cost3=saa_gaussian_cost(3,100,0.1)
cost4=saa_gaussian_cost(4,100,0.1)
cost5=saa_gaussian_cost(5,100,0.1)
cost6=saa_gaussian_cost(6,100,0.1)
cost7=saa_gaussian_cost(7,100,0.1)
cost8=saa_gaussian_cost(8,100,0.1)
cost9=saa_gaussian_cost(9,100,0.1)
cost10=saa_gaussian_cost(10,100,0.1)
cost11=saa_gaussian_cost(11,100,0.1)
cost12=saa_gaussian_cost(12,100,0.1)
cost13=saa_gaussian_cost(13,100,0.1)

c5=0.01*ones(5)
c10=0.01*ones(10)
c20=0.01*ones(20)
c50=0.01*ones(50)
c100=0.01*ones(100)

C104=RandCoupConsNormalHrep(10,4,10,30,10.)
C504=RandCoupConsNormalHrep(50,4,10,100,10.)
C1004=RandCoupConsNormalHrep(100,4,10,200,10.)

C105=RandCoupConsNormalHrep(10,5,20,30,10.)
C505=RandCoupConsNormalHrep(50,5,20,100,10.)
C1005=RandCoupConsNormalHrep(100,5,10,200,10.)

C106=RandCoupConsNormalHrep(10,6,20,30,10.)
C506=RandCoupConsNormalHrep(50,6,20,100,10.)
C1006=RandCoupConsNormalHrep(100,6,10,200,10.)
C10010=RandCoupConsNormalHrep(100,10,20,500,10.)