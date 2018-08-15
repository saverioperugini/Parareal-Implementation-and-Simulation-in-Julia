# parareal.jl: Parareal Algorithm Implementation and Simulation in Julia
# Copyright (C) 2018  Tyler M. Mastay and Saverio Perugini

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


@everywhere function parareal(a,b,nC,nF,K,y0,f,coarseIntegrator,fineIntegrator)
#initialize coarse information
xC = linspace(a,b,nC+1);
yC = zeros(size(xC,1),K);
deltaC = (b-a) / (nC + 1);
yC[1,:] = y0;

#"coarse integrator partially evaluated"
ciPEvaled = ((x1,y1) -> coarseIntegrator(deltaC,x1,y1,f));

#get initial coarse integration solution
for i=2:(nC+1)
   yC[i,1] = ciPEvaled(xC[i-1],yC[i-1,1]);
end
correctC = copy(yC);

#initialize fine information
xF = zeros(nC,nF+1);
for i=1:nC
   xF[i,:] = linspace(xC[i],xC[i+1],nF+1);
end
sub = zeros(nC,nF+1,K);
deltaF = xF[1,2] - xF[1,1];

#"fine integrator partially evaluated"
fiPEvaled = ((x1,y1) -> fineIntegrator(deltaF,x1,y1,f));

for k=2:K
   #run fine integration on each subdomain
   tic();
   @sync for i=1:nC
      sub[i,1,k] = correctC[i,k-1];
      @async for j=2:(nF+1)
         sub[i,j,k] = fiPEvaled(xF[i,j-1],sub[i,j-1,k]);
      end
   end
   toc();
   
   #predict and correct
   for i=1:nC
      yC[i+1,k] = ciPEvaled(xC[i],correctC[i,k]);
      correctC[i+1,k] = yC[i+1,k] - yC[i+1,k-1] + sub[i,nF+1,k];
   end
end

yF = zeros(nC*(nF+1),K-1);
for k=2:K
   yF[:,k-1] = reshape(sub[:,:,k]',nC*(nF+1));
end

return reshape(xF',nC*(nF+1)),reshape(sub[:,:,K]',nC*(nF+1)),yF,sub,xC,correctC,yC;
end

@everywhere function fullMethod(n,a,b,y0,f,integrator)
   #setup domain and range space
    x = linspace(a,b,n+1);
   deltaX = x[2] - x[1];
    y = ones(n+1,1);
   
   #initialize left endpoint
    y[1] = y0;
   
   #integrate each point
    for i=1:n
        y[i+1] = integrator(deltaX,x[i],y[i],f);
    end
   return x,y;
end

function simulate(a,b,N,M,K,y0,f,coarseInt,fineInt,showPrev)
   x1,y1 = fullMethod(N*(M+1),a,b,y0,f,fineInt);
   x,y,yF,sub,xC,yC,iC = parareal(a,b,N,M,K,y0,f,coarseInt,fineInt);
   xF = (reshape(x,M+1,N))';
   fine = M+1;
   for k=2:K
      display(plot(x1,y1));
      if(showPrev && k > 2 )
         display(scatter!(xC,yC[:,k-2],color="red",legend=false));
      end
      display(scatter!(xC,yC[:,k-1],color="green",legend=false));
      done = zeros(Int64,N,1);
      workingSubdomains = 1:N;
      while(done != (M+1) * ones(N,1) )
         index = Int64(ceil(size(workingSubdomains,1)*rand()));
         currThread = workingSubdomains[index];
         while( done[currThread] == M+1 )
            currThread = Int64(ceil(N * rand()));
         end
         currThreadPlot = Int64(ceil(fine*rand()));
         totalAdvance = done[currThread] + currThreadPlot;
         if(totalAdvance > fine) totalAdvance = fine; end
         newP = (done[currThread]+1):totalAdvance;
         display(plot!(xF[currThread,newP],sub[currThread,newP,k],color="black"));
         done[currThread] = totalAdvance;
         workingSubdomains = find( ((x)->x != M+1), done );
         print(join(["Working on subdomain #", currThread, "...",
            "Pending Subdomains: ", workingSubdomains', "\n"]));
      end
      display(plot!(x,yF[:,k-1],color="orange"));
      sleep(5);
   end
end

# Implementation schemes.
function euler(delta,x0,y0,f)
   return y0 + delta * f(x0,y0);
end

function rungeKutta(delta,x0,y0,f)
   k1 = f(x0,y0);
   k2 = f(x0+delta/2,y0 + (delta/2)*k1);
   k3 = f(x0+delta/2,y0 + (delta/2)*k2);
   k4 = f(x0+delta,y0+delta*k3);
   return y0 + (delta/6)*(k1+2*k2+2*k3+k4);
end
