function dNpdt = rxn_solbulk(t,Np,k1,a,b,n,nk,k2,v,bn,L0,dNpdt)

    for i = 1:nk 
            if i==1 
                r1 = (k1/L0^(bn-b))*(v^(1-a-bn))*(Np(end-1)^a)*(Np(end)^bn); 
            else 
                r1 = 0; 
            end
            dNpdt(end)   = dNpdt(end)  -6*r1;            
            dNpdt(end-1) = dNpdt(end-1)-6*r1;
            dNpdt(i)     = dNpdt(i) + r1 - (k2/v)*Np(i)*sum(Np(i:n-i));
    
            for j = i:n 
                    if j<2*i
                        dNpdt(j) = dNpdt(j) - (k2/v)*Np(i)*Np(j);
                    end
                    if j>=2*i && j+i<=n
                        dNpdt(j) = dNpdt(j) - (k2/v)*Np(i)*(Np(j)-Np(j-i));
                    end
                    if j>=2*i && j+i>n
                        dNpdt(j) = dNpdt(j) + (k2/v)*Np(i)*Np(j-i);
                    end
            end
    end

end