function dNpdt = rxn_physisorp3(t,Np,n,kssf,Kssi,v,dNpdt)
    for i = 2:n 
        for j = i:n
                ksNp         =   (kssf)*(  Np(i)*Np(3*n+j)/v  -  Kssi*Np(3*n+i) );
                if ksNp >0

                dNpdt(i)       = dNpdt(i)       -  j*ksNp;
                dNpdt(3*n+j)   = dNpdt(3*n+j)   -  i*ksNp;
		        dNpdt(4*n+j)   = dNpdt(4*n+j)   +  i*ksNp;
		        dNpdt(3*n+i)   = dNpdt(3*n+i)   +  j*ksNp;
                end
        end
    end
end