
function load_pars(file_name::String)
   input = open(file_name)
   s = read(input, String)
   # Aggiungo anche il nome del file alla lista di parametri.
   p =JSON.parse(s)
   return p
end

function defineSystem(;sys_type::String="S=1/2",sys_istate::String="Up",chain_size::Int64=250,local_dim::Int64=6)
   
   sys = siteinds(sys_type,1);
   env = siteinds("Boson",dim = local_dim, chain_size);
   sysenv = vcat(sys,env);

   #isTrans = sys_istate == "i"
   
   #Temporary fix for initial state +1 of σ_y
   #If sys_istate == "i", then initialize to "Up"...we will then apply 1/√2 (1+2 Sy) to 
   #get the initial state.

   #stateSys = [(isTrans ? "Up" : sys_istate)];
   stateSys = [sys_istate]

   #Standard approach: chain always in the vacuum state
   stateEnv = ["0" for n=1:chain_size];

   stateSE = vcat(stateSys,stateEnv);

   psi0 = productMPS(ComplexF64,sysenv,stateSE);

   #In case Initial State is "i"
   # if isTrans
   #    sy = 2*op("Sy",sysenv[1])
   #    si = op("Id",sysenv[1])
   #    rot = 1/sqrt(2)*(si+sy)
   #    psi0 = apply(rot,psi0)
   # end
   
   return (sysenv,psi0);
end


#function setInital()
function createMPO(sysenv, eps::Float64, delta::Float64, intHsysSide::String, freqfile::String, coupfile::String)::MPO
   
   
   coups = readdlm(coupfile);
   freqs = readdlm(freqfile);

   NN::Int64  = size(sysenv)[1];
   NChain::Int64 = NN-1;
   

   #Check if the interaction operator on system side is one of the
   #spin operators.
   #ATTENTION:    
   #!Sx = 0.5 σx
   twoFact = any(map(x->intHsysSide==x,["Sx","Sy","Sz"]))

   #in this case, whenever we need to use σ_i we need to multiply
   #by 2.

   #to test
   if intHsysSide == "S+"
      intHsysSideDag = "S-"
   elseif intHsysSide == "S-"
      intHsysSideDag = "S+"
   else
      intHsysSideDag = intHsysSide
   end

   if twoFact
      println("spin-boson interaction operator: ", intHsysSide," ⊗ (A+Adag)")
   else
      println("spin-boson interaction operator: ", intHsysSide," ⊗ A +",intHsysSideDag, "⊗ Adag")
   end
   
   thempo = OpSum();
   #system Hamiltonian
   #pay attention to constant S_x/y/z = 0.5 σ_x/y/z
   thempo += 2*eps,"Sz",1;
   thempo += 2*delta,"Sx",1
   #system-env interaction
   
   thempo += (twoFact ? 2 : 1)*coups[1],intHsysSide,1,"Adag",2;
   thempo += (twoFact ? 2 : 1)*coups[1],intHsysSideDag,1,"A",2;
   #chain local Hamiltonians
   for j=2:NChain
      thempo += freqs[j-1],"N",j;
   end
   
   for j=2:NChain-1
      thempo += coups[j],"A",j,"Adag",j+1;
      thempo += coups[j],"Adag",j,"A",j+1;
   end
   return MPO(thempo,sysenv);
end


#For MC
#ATTENTION: we need to implement the 2 factor used in the  No-MC createMPO
function createMPO(sysenv, eps::Float64, delta::Float64,intHsysSide::String, freqfile::String, coupfile::String,
   MC_alphafile::String,MC_betafile::String,MC_coupfile::String,omega::Float64; kwargs...)::MPO

   #perm:  permutation of the closure oscillators
   perm = get(kwargs,:perm,nothing)

   #Chain parameters
   coups = readdlm(coupfile);
   freqs = readdlm(freqfile);
   
   #Closure parameters loaded from file
   alphas_MC = readdlm(MC_alphafile);
   betas_MC = readdlm(MC_betafile);
   coups_MC = readdlm(MC_coupfile);
   gammas = omega * alphas_MC[:,1];
   eff_freqs = [omega + 1im * g for g in gammas];
   eff_gs = omega * betas_MC[:,2];
   eff_coups = omega/2* (coups_MC[:,1]+ 1im *coups_MC[:,2]);

   NN = length(sysenv);
   MC_N = length(gammas);
   NP_Chain = NN-MC_N;


   #Check if the interaction operator on system side is one of the
   #spin operators.
   #ATTENTION:    
   #!Sx = 0.5 σx
   twoFact = any(map(x->intHsysSide==x,["Sx","Sy","Sz"]))

   #in this case, whenever we need to use σ_i we need to multiply
   #by 2.

   #to test
   if intHsysSide == "S+"
      intHsysSideDag = "S-"
   elseif intHsysSide == "S-"
      intHsysSideDag = "S+"
   else
      intHsysSideDag = intHsysSide
   end

   if twoFact
      println("spin-boson interaction operator: ", intHsysSide," ⊗ (A+Adag)")
   else
      println("spin-boson interaction operator: ", intHsysSide," ⊗ A +",intHsysSideDag, "⊗ Adag")
   end
   
   if(perm != nothing)
      
      if(length(perm)!= MC_N)
         println("The provided permutation is not correct")
      end

      pmtx=Permutation(perm)
      @show pmtx
   else
      #Identity permutation
      pmtx=Permutation(collect(1:MC_N))
      @show pmtx
   end

   #Lavoriamo qui
   thempo = OpSum();
   #system Hamiltonian
   #pay attention to constant S_x/y/z = 0.5 σ_x/y/z
   thempo += 2*eps,"Sz",1;
   thempo += 2*delta,"Sx",1;
   #system-env interaction
   #!Sx = 0.5 σx
   thempo += (twoFact ? 2 : 1)*coups[1],intHsysSide,1,"Adag",2
   thempo += (twoFact ? 2 : 1)*coups[1],intHsysSideDag,1,"A",2

   #Primary chain local Hamiltonians
   for j=2:NP_Chain
      thempo += freqs[j-1],"N",j
   end

   for j=2:NP_Chain-1
      thempo += coups[j],"A",j,"Adag",j+1
      thempo += coups[j],"Adag",j,"A",j+1
   end
   #################################

   #Markovian closure Hamiltonian
   for j=1:MC_N
      thempo += eff_freqs[j],"N",NP_Chain+pmtx(j)
   end

   for j=1:MC_N-1
      thempo += eff_gs[j],"A",NP_Chain+pmtx(j),"Adag",NP_Chain+pmtx(j+1)
      thempo += eff_gs[j],"Adag",NP_Chain+pmtx(j),"A",NP_Chain+pmtx(j+1)
   end

   #################################

   #Primary chain - MC interaction
   for j=1:MC_N
      thempo += eff_coups[j],"A",NP_Chain,"Adag",NP_Chain+pmtx(j)
      thempo += conj(eff_coups[j]),"Adag",NP_Chain,"A",NP_Chain+pmtx(j)
   end
   #################################

   return MPO(thempo,sysenv);
end
function createMPO2MC(sysenv, eps::Float64, delta::Float64, intHsysSide::String, freqfile::String, coupfile::String,
   MC_alphafile::String,MC_betafile::String,MC_coupfile::String,omega::Float64; kwargs...)::MPO

   #perm:  permutation of the closure oscillators
   perm = get(kwargs,:perm,nothing)

   #Chain parameters
   coups = readdlm(coupfile);
   freqs = readdlm(freqfile);
   
   #Closure parameters loaded from file
   alphas_MC = readdlm(MC_alphafile);
   betas_MC = readdlm(MC_betafile);
   coups_MC = readdlm(MC_coupfile);
   gammas = omega * alphas_MC[:,1];
   eff_freqs = [omega + 1im * g for g in gammas];
   eff_gs = omega * betas_MC[:,2];
   #Reduce interaction with each closure
   eff_coups = 1/sqrt(2)*omega/2* (coups_MC[:,1]+ 1im *coups_MC[:,2]);

   NN = length(sysenv);
   MC_N = length(gammas);
   NP_Chain = NN - 2 * MC_N;
   
   #Check if the interaction operator on system side is one of the
   #spin operators.
   #ATTENTION:    
   #!Sx = 0.5 σx
   twoFact = any(map(x->intHsysSide==x,["Sx","Sy","Sz"]))

   #in this case, whenever we need to use σ_i we need to multiply
   #by 2.

   #to test
   if intHsysSide == "S+"
      intHsysSideDag = "S-"
   elseif intHsysSide == "S-"
      intHsysSideDag = "S+"
   else
      intHsysSideDag = intHsysSide
   end

   if twoFact
      println("spin-boson interaction operator: ", intHsysSide," ⊗ (A+Adag)")
   else
      println("spin-boson interaction operator: ", intHsysSide," ⊗ A +",intHsysSideDag, "⊗ Adag")
   end
   

   if(perm != nothing)
      
      if(length(perm)!= MC_N)
         println("The provided permutation is not correct")
      end

      pmtx=Permutation(perm)
      @show pmtx
   else
      #Identity permutation
      pmtx=Permutation(collect(1:MC_N))
      @show pmtx
   end

   #Lavoriamo qui
   thempo = OpSum();
   #system Hamiltonian
   #pay attention to constant S_x/y/z = 0.5 σ_x/y/z
   thempo += 2*eps,"Sz",1;
   thempo += 2*delta,"Sx",1;
   #system-env interaction
   #!Sx = 0.5 σx
   thempo += (twoFact ? 2 : 1)*coups[1],intHsysSide,1,"Adag",2
   thempo += (twoFact ? 2 : 1)*coups[1],intHsysSideDag,1,"A",2

   #Primary chain local Hamiltonians
   for j=2:NP_Chain
      thempo += freqs[j-1],"N",j
   end

   for j=2:NP_Chain-1
      thempo += coups[j],"A",j,"Adag",j+1
      thempo += coups[j],"Adag",j,"A",j+1
   end
   #################################

   #First+Second Markovian closure Hamiltonian
   for j=1:MC_N
      thempo += eff_freqs[j],"N",NP_Chain+pmtx(j)
      #Second closure
      thempo += eff_freqs[j],"N",NP_Chain+MC_N+pmtx(j)
   end

   for j=1:MC_N-1
      thempo += eff_gs[j],"A",NP_Chain+pmtx(j),"Adag",NP_Chain+pmtx(j+1)
      thempo += eff_gs[j],"Adag",NP_Chain+pmtx(j),"A",NP_Chain+pmtx(j+1)
      #Second clousure
      thempo += eff_gs[j],"A",NP_Chain+ MC_N + pmtx(j),"Adag",NP_Chain+MC_N+pmtx(j+1)
      thempo += eff_gs[j],"Adag",NP_Chain+MC_N + pmtx(j),"A",NP_Chain+MC_N+pmtx(j+1)
   end

   #################################

   #Primary chain - MC interaction
   for j=1:MC_N
      thempo += eff_coups[j],"A",NP_Chain,"Adag",NP_Chain+pmtx(j)
      thempo += conj(eff_coups[j]),"Adag",NP_Chain,"A",NP_Chain+pmtx(j)
      #Second closure
      thempo += eff_coups[j],"A",NP_Chain,"Adag",NP_Chain + MC_N + pmtx(j)
      thempo += conj(eff_coups[j]),"Adag",NP_Chain,"A",NP_Chain + MC_N + pmtx(j)
   end
   #################################

   return MPO(thempo,sysenv);
end
function createObs(lookat)
   
   vobs = [];
   for a in lookat
      push!(vobs,opPos(a[1],a[2]));
   end

   return vobs;

end

function stretchBondDim(state::MPS,extDim::Int64)
   psiExt = copy(state);
   NN = length(psiExt)
   for n in 1:NN-1
      a = commonind(psiExt[n],psiExt[n+1])
      tagsa = tags(a)
      add_indx = Index(extDim, tags= tagsa)
      psiExt[n]=psiExt[n]*delta(a,add_indx)
      psiExt[n+1]=psiExt[n+1]*delta(a,add_indx)
   end
   #println("Overlap <original|extended> states: ", dot(state,psiExt));
   return psiExt, dot(state,psiExt)
end





function test_exsistence()
   return "MPSSpinBoson package does exist!"
end

