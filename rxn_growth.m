function dNpdt = rxn_growth(t,Np,n,nk,k3,v,dNpdt)
% fprintf('%0.2f\t',t)

    for i = 1:round(nk/2)
            % Bulk Pi reaction
            dNpdt(i) = dNpdt(i)-k3*Np(i)*sum(Np(3*n+1:4*n-i))/v;
            % Adsorbed Pi reactions
            for j = 1:n 
                    if j<=i
                        dNpdt(3*n+j) = dNpdt(3*n+j)  - k3*Np(i)*Np(3*n+j)/v;
                    end
            
                    if j>i && j+i<=n
                        dNpdt(3*n+j) = dNpdt(3*n+j)  + k3*Np(i)*( Np(3*n+j-i)-Np(3*n+j)  )/v;
                    end
            
                    if j>i && j+i>n 
                        dNpdt(3*n+j) = dNpdt(3*n+j)  + k3*Np(i)*Np(3*n+j-i)/v;
                    end
            end
    end
end
