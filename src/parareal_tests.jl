using Plots
using Distributed

include("parareal.jl")

#First test case
a         = 0.0
b         = 12.56
nC        = 10
nF        = 100
K         = Int(nF / nC)
y0        = 0.0
f         = ((x,y) -> sin(x))
coarseInt = rungeKutta
fineInt   = rungeKutta
showPrev  = true

function go()
   simulate(a,b,nC,nF,K,y0,f,coarseInt,fineInt,showPrev);
end

go();

#Second test case
a         = 0.0
b         = 12.56
nC        = 10
nF        = 100
K         = Int(nF / nC)
y0        = 1.57
f         = ((x,y) -> sin(y^3))
coarseInt = rungeKutta
fineInt   = rungeKutta
showPrev  = true

go();
