using Distributed


#    INPUT PARAMETERS:
#        -- a -- left endpoint of domain
#        -- b -- right endpoint of domain
#        -- nC -- number of coarse grid dofs
#        -- nF -- number of find grid dofs
#        -- K  -- number of subdomains
#        -- y0 -- IVP value on left endpoint
#        -- f  -- RHS forcing function
#        -- coarseIntegrator -- integration scheme for coarse integration
#        -- fineIntegrator   -- integration scheme for fine integration
#            ** integrators take 4 arguments
#                && h && step size 
#                && x1 && current point
#                && y1 && value at x1
#                && f  && forcing function
#
#    OUTPUT PARAMETERS:
#        -- x -- x values
#        -- y -- y values
#        -- yF -- fine integration values
#        -- sub -- subdomain data
#        -- xC  -- coarse grid values
#        -- correctC -- history of corrected coarse-grid values
#        -- yC       -- history of coarse-grid y-values?
    
@everywhere function parareal(a,b,nC,nF,K,y0,f,coarseIntegrator,fineIntegrator)
#initialize coarse information
xC = LinRange(a,b,nC+1);
yC = zeros(size(xC,1),K);
deltaC = (b-a) / (nC + 1);
yC[1,:] = y0 * ones(size(yC[1,:]));

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
   xF[i,:] = LinRange(xC[i],xC[i+1],nF+1);
end
sub = zeros(nC,nF+1,K);
deltaF = xF[1,2] - xF[1,1];

#"fine integrator partially evaluated"
fiPEvaled = ((x1,y1) -> fineIntegrator(deltaF,x1,y1,f));

for k=2:K
   #run fine integration on each subdomain
   start = time();
   @sync for i=1:nC
      sub[i,1,k] = correctC[i,k-1];
      @async for j=2:(nF+1)
         sub[i,j,k] = fiPEvaled(xF[i,j-1],sub[i,j-1,k]);
      end
   end
   elapsed = time() - start;
   display(elapsed);
 
   #predict and correct
   for i=1:nC
      yC[i+1,k] = ciPEvaled(xC[i],correctC[i,k]);
      correctC[i+1,k] = yC[i+1,k] - yC[i+1,k-1] + sub[i,nF+1,k];
   end
end

#data reformatting step
yF = zeros(nC*(nF+1),K-1);
for k=2:K
   yF[:,k-1] = reshape(sub[:,:,k]',nC*(nF+1));
end

return reshape(xF',nC*(nF+1)),reshape(sub[:,:,K]',nC*(nF+1)),yF,sub,xC,correctC,yC;
end

@everywhere function fullMethod(n,a,b,y0,f,integrator)
   #setup domain and range space
   x = LinRange(a,b,n+1);
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

# CREATE DOCUMENTATION FOR SIMULATE

function simulate(a,b,N,M,K,y0,f,coarseInt,fineInt,showPrev)
   #simulation setup parameters
   sim_sleep_time      = 3

   #grab the ``correct'' solution up to the scheme used
   x1,y1               = fullMethod(N*(M+1),a,b,y0,f,fineInt);
   x,y,yF,sub,xC,yC,iC = parareal(a,b,N,M,K,y0,f,coarseInt,fineInt);

   #setup main simulation
   xF = (reshape(x,M+1,N))';
   fine = M+1;

   #main simulation loop
   for k=2:K
      #plot ``correct'' solution obtained from direct fine integration
      display(plot(x1,y1));

      #display the coarse integration points from previous step
      if(showPrev && k > 2 )
         display(scatter!(xC,yC[:,k-2],color="red",legend=false));
      end

      #display current guesses
      display(scatter!(xC,yC[:,k-1],color="green",legend=false));
      
      #done -- stores if a subdomain has completed its execution or not
      done = zeros(Int64,N,1);

      #workingSubdomains -- an array of which subdomians are still working
      workingSubdomains = 1:N;

      #as long as not everyone is done, continue the simulation
      while(done != (M+1) * ones(N,1) )
         #grab random subdomain to advance on
         index = Int64(ceil(size(workingSubdomains,1)*rand()));
         currThread = workingSubdomains[index];

         #continue grabbing thread as long as it is still on working
         #   subdomain
         while( done[currThread] == M+1 )
            currThread = Int64(ceil(N * rand()));
         end

         #advance appropriate subdomain through newP
         currThreadPlot = Int64(ceil(fine*rand()));
         totalAdvance = done[currThread] + currThreadPlot;
         if(totalAdvance > fine) totalAdvance = fine; end
         newP = (done[currThread]+1):totalAdvance;

         #display updated values on fine mesh level 
         display(plot!(xF[currThread,newP],
            sub[currThread,newP,k],color="black"));

         #update how much work has been completed by currThread
         done[currThread] = totalAdvance;

         #update to see if any subdomains have completed their work
         #   if so, remove them from the list 
         workingSubdomains = map( ((y) -> y[1]), 
            findall( ((x)->x != M+1), done ) );

         #print out info of interest
         print(join(["Working on subdomain #", currThread, "...",
            "Pending Subdomains: ", workingSubdomains', "\n"]));
      end
      display(plot!(x,yF[:,k-1],color="orange"));
      sleep(sim_sleep_time);
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
