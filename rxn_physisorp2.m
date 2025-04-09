function dNpdt = rxn_physisorp2(t,Np,n,kssf,Kssi,v,dNpdt)
    for i = 1:n 
        for j = i:n
                ksNp         =   (kssf)*(  Np(i)*Np(n+j)/v  -  Kssi*Np(3*n+i)  );
                if ksNp >0
                    dNpdt(i)     = dNpdt(i)     -  j*ksNp;
                    dNpdt(n+j)   = dNpdt(n+j)   -  i*ksNp;
		            dNpdt(2*n+j) = dNpdt(2*n+j) +  i*ksNp;
		            dNpdt(3*n+i) = dNpdt(3*n+i) +  j*ksNp;
                end
        end
    end
end