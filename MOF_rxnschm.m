function dNpdt = MOF_rxnschm(t,Np,params,n,nk,k1,k2,k3,ksf,Ksi,kssf,Kssi,a,b,all_i,Ns0,alpha,beta,gamma)

   dNpdt        = zeros(5*n+3,1);
   Np(Np<1e-30) = 0;  
   v            = params.bvf*params.dvl+(1-params.bvf)*params.dvl*(params.va*exp (-params.vb*t)+ ...
                  params.vc*exp(-params.vd*t))/(params.va+params.vc)+params.dvs;
   N1_Sol_tot   = sum(all_i(1:n).*Np(1:n));

   k1n          = k1*(1+alpha*N1_Sol_tot/params.dv+beta*Np(end-2)/Ns0+gamma*(Ns0-Np(end-2))/Ns0);
   bn           = params.b_new + (b-params.b_new)./(1+exp(0.15*(t-params.tlim)));
   
   if (v <= params.dvs) || (Np(end-1)<=0) || (Np(end)<=0)  
        dNpdt = zeros(5*n+3,1);                                    
   else
        dNpdt = rxn_solbulk(t,Np,k1n,a,b,n,nk,k2,v,bn,params.L0,dNpdt);           
        dNpdt = rxn_physisorp1(t,Np,n,ksf,Ksi,v,dNpdt);    
        dNpdt = rxn_physisorp2(t,Np,n,kssf,Kssi,v,dNpdt);          
        dNpdt = rxn_physisorp3(t,Np,n,kssf,Kssi,v,dNpdt);          
        dNpdt = rxn_growth(t,Np,n,nk,k3,v,dNpdt);                  
   end
end