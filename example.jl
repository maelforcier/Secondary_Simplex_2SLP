
lib = DefaultLibrary{Float64}(GLPK.Optimizer)
# lib= CDDLib.Library(:exact)

include("rand_coup_cons.jl")

include("plot.jl")
include("discrete_cost.jl")

dimx=1 #dimension of the initial state
dimy=2 #dimension of the terminal state

plotmodel= Model(GLPK.Optimizer)
@variable(plotmodel,x[1:dimx])
@variable(plotmodel,y[1:dimy])
@constraint(plotmodel,0.1<=x[1]<=5)
@constraint(plotmodel,0.1<=y[2]<=5)
@constraint(plotmodel,0.1<=y[1]<=5)
@constraint(plotmodel,x[1]+y[1]+y[2]<=10)
@constraint(plotmodel,x[1]+y[1]+3y[2]<=15)
@constraint(plotmodel,2x[1]+y[1]<=13)
@constraint(plotmodel,x[1]-2y[1]<=3)
@constraint(plotmodel,x[1]+y[1]+3y[2]<=15)
@constraint(plotmodel,0.5x[1]+2y[1]+3y[2]<=18);


#Example in FGL20
model1d= Model(GLPK.Optimizer)
@variable(model1d,x)
@variable(model1d,y[1:2])
@constraint(model1d,y[1]+y[2]<=1)
@constraint(model1d,y[1]-y[2]<=1)
@constraint(model1d,-y[1]-y[2]<=1)
@constraint(model1d,-y[1]+y[2]<=1)
@constraint(model1d,y[1]<=x)
@constraint(model1d,y[2]<=x)
@constraint(model1d,x<=3);

cross12= CoupCons(model1d,1,2)


#Generalisation in dimension 2+2 of FGL20
model2d= Model(GLPK.Optimizer)
@variable(model2d,x[1:2])
@variable(model2d,y[1:2])
@constraint(model2d,y[1]+y[2]<=1)
@constraint(model2d,y[1]-y[2]<=1)
@constraint(model2d,-y[1]-y[2]<=1)
@constraint(model2d,-y[1]+y[2]<=1)
@constraint(model2d,y[1]<=x[1])
@constraint(model2d,y[2]<=x[2])
# @constraint(model2d,x[1]<=3)
# @constraint(model2d,x[2]<=3);
A22=Matrix{Float64}(undef,2,2)
b22=Vector{Float64}(undef,2)
A22[:,:]=I(2)
b22[:]=3*ones(2)

cross22= CoupCons(model2d,2,2,A22,b22)
cross22pert=perturbate_coupcons_gaussian(cross22,0.01)


#An interior point in each relint of chamber of cross22
chrepcross22=[[0.25,0.25],[-0.25,0.25],[-0.25,-0.25],[0.25,-0.25],[1.5,1.5],[1.5,0.5],[1.5,-0.5],[0.75,0.75],[0.5,1.5],[-0.5,1]]


# for x in chrep
# 	println(x)
# 	println(ActiveConsColl(cross22,x))
# end



#Example with generic normal vectors
quadri=polyhedron(vrep([[-3,0],[-1,4],[0,-3],[5,-1]]))
hrep(quadri)
model2d2= Model(GLPK.Optimizer)
@variable(model2d2,x[1:2])
# @variable(model2d2,y[1:2])
y = @variable(model2d2, [1:2])
@constraint(model2d2,5y[1]+6y[2]<=19)
@constraint(model2d2,2y[1]-5y[2]<=15)
@constraint(model2d2,-y[1]-y[2]<=3)
@constraint(model2d2,-2y[1]+y[2]<=6)
# @constraint(model2d2,y in quadri)
@constraint(model2d2,y[1]<=x[1])
@constraint(model2d2,y[2]<=x[2])
@constraint(model2d2,x[1]<=8)
@constraint(model2d2,x[2]<=8);
C22=CoupCons(model2d2,2,2)

#An interior point in each relint of chamber of C22
chrepC22=[[-1,-1],[-1,1],[-2,3],[-0.5,3.8],[-0.5,5],[1,5],[1,2],[2,-1],[4,2.5],[4.5,-0.1],[4,-2.5],[6,-0.5],[6,1],[6,6]]

#Interesting vertices
 # [19/5,0.] 
 # [5,4]
 # [5,8.]
 # [5,-1.] 
 # [0,-3]
 # [-3.,0]
 # [-1,4.]
 # [0.,19/6]

#Example with generic normal vectors
quadri=polyhedron(vrep([[-3,0],[-1,4],[0,-3],[5,-1]]))
hrep(quadri)
model2d2= Model(GLPK.Optimizer)
@variable(model2d2,x[1:2])
# @variable(model2d2,y[1:2])
y = @variable(model2d2, [1:2])
@constraint(model2d2,5y[1]+6y[2]<=19)
@constraint(model2d2,2y[1]-5y[2]<=15)
@constraint(model2d2,-y[1]-y[2]<=3)
@constraint(model2d2,-2y[1]+y[2]<=6)
# @constraint(model2d2,y in quadri)
@constraint(model2d2,y[1]<=x[1])
@constraint(model2d2,y[2]<=x[2])
C22unbound=CoupCons(model2d2,2,2)
A22=Matrix{Float64}(undef,2,2)
b22=Vector{Float64}(undef,2)
A22[:,:]=I(2)
b22[:]=8*ones(2)
C22new=CoupCons(model2d2,2,2,A22,b22)


#Example with generic normal vectors in dim 2+3
model2d3= Model(GLPK.Optimizer)
@variable(model2d3,x[1:2])
y = @variable(model2d3, [1:3])
@constraint(model2d3,5y[1]+6y[2]<=19)
@constraint(model2d3,2y[1]-5y[2]<=15)
@constraint(model2d3,-y[1]-y[2]<=3)
@constraint(model2d3,-2y[1]+y[2]<=6)
@constraint(model2d3,y[1]+y[3]<=1)
@constraint(model2d3,y[3]<=2)
@constraint(model2d3,y[3]>=-2)
# @constraint(model2d3,y in quadri)
@constraint(model2d3,y[1]<=x[1])
@constraint(model2d3,y[2]<=x[2])
@constraint(model2d3,y[3]<=x[2]+x[1])
@constraint(model2d3,x[1]<=8)
@constraint(model2d3,x[2]<=8);
C23= CoupCons(model2d3,2,3)

chrepC23=
[[0.2,0.1],[0.6,0.6],[0.5,0.3],[1.5,0.3],[1.5,3],[1.5,1.5],[-0.5,1.5],[-0.5,2.8],[-0.5,-0.5],[2,-1],[2,-1.3],[0.5,-1.3],[2,-0.3],[4.1,0.3],[4.1,-0.3],[4.1,-2],[4.1,-2.3],[4.1,-1.9],[2.5,-1.9],[2.5,-2.9],[2.5,-2.4],[5.5,2.4],[-2.4,3.1],[2.5,4.5]]
exsommetofC23=[3.,2/3]

C23pert=perturbate_coupcons_gaussian(C23,0.01)
# C23pertrat=perturbate_coupcons_det(C23,1//1000,Rational{Int64})

model2= Model(GLPK.Optimizer)
@variable(model2,x[1:2])
@variable(model2,y[1:2])
@constraint(model2,y[1]+y[2]<=1)
@constraint(model2,y[1]-y[2]<=1)
@constraint(model2,-y[1]-y[2]<=1)
@constraint(model2,-y[1]+y[2]<=1)
@constraint(model2,y[1]<=x[1])
@constraint(model2,x[2]<=y[1])
@constraint(model2,x[1]<=3)
@constraint(model2,-3<=x[2]<=3);


#Example in FGL20 but x and y inversed
model21d= Model(GLPK.Optimizer)
@variable(model21d,x[1:2])
@variable(model21d,y)
@constraint(model21d,x[1]+x[2]<=1)
@constraint(model21d,x[1]-x[2]<=1)
@constraint(model21d,-x[1]-x[2]<=1)
@constraint(model21d,-x[1]+x[2]<=1)
@constraint(model21d,x[1]<=y)
@constraint(model21d,x[2]<=y)
@constraint(model21d,y<=3);

cross21= CoupCons(model21d,2,1)




# T=[0 0;0 0;0 0;0 0;-1 0;0 -1] # constraint matrix of x
# W=[1 1;1 -1;-1 -1;-1 1;1 0;0 1] #constraint matrix of y
# h=[1;1;1;1;0;0];   #right term in the inequality constraint

T=[0 0;0 0;0 0;0 0;-1 0;0 -1;1 0;0 1] # constraint matrix of x
W=[1 1;1 -1;-1 -1;-1 1;1 0;0 1;0 0;0 0] #constraint matrix of y
h=[1;1;1;1;0;0;3;3];   #right term in the inequality constraint

# T=[0. 0.;0. 0.;0. 0.;0. 0.;-1. 0.;0. -1.] # constraint matrix of x
# W=[1. 1.;1. -1.;-1. -1.;-1. 1.;1. 0.;0. 1.] #constraint matrix of y
# h=[1.;1.;1.;1.;0.;0.];   #right term in the inequality constraint

T1d=Int64.(zeros(7,1))
T1d[5,1]=-1
T1d[6,1]=-1
T1d[7,1]=1 # constraint matrix of x
W1d=[1 1;1 -1;-1 -1;-1 1;1 0;0 1;0 0] #constraint matrix of y
h1d=[1;1;1;1;0;0;3];   #right term in the inequality constraint


#Other constraints
T2=[0 0;0 0;0 0;0 0;-1 0;0 -1;1 0;0 -1;0 1] # constraint matrix of x
W2=[1 1;1 -1;-1 -1;-1 1;1 0;-1 0;0 0;0 0;0 0] #constraint matrix of y
h2=[1;1;1;1;0;0;3;3;3]   #right term in the inequality constraint





# plotChamber=plot()
# PlotChamberIt(cross22pert,chrepcross22)
# plotChamber

plotChamber=plot()
PlotChamberIt(C22,chrepC22)
plotChamber

# plotChamber=plot()
# PlotChamberIt(C23,chrepC23)
# plotChamber