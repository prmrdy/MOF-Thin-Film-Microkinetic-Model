function dNpdt = rxn_physisorp1(t,Np,n,ksf,Ksi,v,dNpdt)
    for i = 1:n
            ksNp     =   ksf*(  Np(i)*Np(end-2)/v  -  Ksi*Np(n+i)  );
            if ksNp>0
                dNpdt(i)     = dNpdt(i)    -   ksNp;   
                dNpdt(n+i)   = dNpdt(n+i)  +   ksNp;
                dNpdt(end-2) = dNpdt(end-2)-   i*ksNp;
            end
    end
end