function [ coilMap] = coil_map_study_2d_Inati( data, ks, power )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

RO = size(data,1);
E1 = size(data,2);
CHA = size(data,3);

N = prod(size(data)) / (RO*E1*CHA);

% str_msg=sprintf('%d %d %d %d\n', N , RO , E1 ,CHA ); disp(str_msg);

%         if (ks % 2 != 1)
%         {
%             ks++;
%         }

kss = ks*ks;

halfKs = ks / 2;

% str_msg=sprintf('%d %d %d \n', kss , kss , halfKs); disp(str_msg);

D=zeros(ks*ks, CHA);
DC = zeros(ks*ks, CHA);
DH_D=zeros(CHA, CHA);
U1 = zeros(ks*ks, 1);
V1 = zeros(CHA,1);
V = zeros(CHA,1);

v1Norm=1;
u1Norm=1;

 halfKs=2;

% str_msg=sprintf('%d %d %d %d \n', halfKs+1 , E1 - halfKs ,halfKs+1 , RO - halfKs); disp(str_msg);

for e1 = 1:1:E1
    
    for ro = 1:1:RO
        
        %  fill the data matrix D
        %  if (e1 >= halfKs && e1<E1 - halfKs && ro >= halfKs && ro<RO - halfKs)
        if (e1 >= halfKs+1 && e1<E1 - halfKs && ro >= halfKs+1 && ro<RO - halfKs)
            
            for cha=1:1:CHA
                
                ind = 0;
                
                for ke1 = -halfKs:1:halfKs
                    
                    de1 = e1 + ke1;
                    
                    for kro = -halfKs:1:halfKs
                        
                        ind=ind+1;
                        
                        % str_msg=sprintf('%d %d %d %d %d\n', ind , cha*kss, ro+kro, de1 ,cha   ); disp(str_msg);
                        D(ind, cha)=data(ro+kro, de1, cha);
                        
                    end
                end
            end
            
        else
            
            for cha=1:1:CHA
                
                ind = 0;
                
                for ke1 = -halfKs:1:halfKs
                    de1 = e1 + ke1;
                    
                    if (de1 < 1)
                        de1 = de1+ E1;
                    end
                    if (de1 > E1)
                        de1 = de1- E1;
                    end
                    
                    for kro = -halfKs:1:halfKs
                        dro = ro + kro;
                        %                           str_msg=sprintf('%d = %d + %d \n',dro ,  ro , kro  ); disp(str_msg);
                        if (dro < 1)
                            %                               str_msg=sprintf('dro < 1   %d  \n',dro  ); disp(str_msg);
                            dro = dro+RO;
                        end
                        
                        if (dro > RO)
                            %  str_msg=sprintf('dro >= RO  %d\n',dro   ); disp(str_msg);
                            dro =dro-RO;
                            
                        end
                        
                        % pD[ind + cha*kss] = pDataCurr[de1*RO + dro];
                        ind=ind+1;
                        %                         str_msg=sprintf('%d %d \n', ro , e1  ); disp(str_msg);
                        %                         str_msg=sprintf('%d cha %d:  dro %d de1 %d  \n', ind , cha, dro, de1    ); disp(str_msg);
                        D(ind,cha)=data(dro,de1,cha);
                        
                    end
                end
            end
        end
        
        for cha=1:1:CHA
            
            V1(cha,1) =sum(D(:,cha));
            
        end
        
        somme=0;
        
        for cha=1:1:CHA
            
            c = V1(cha,1);
            re = real(c);
            im = imag(c);
            somme = somme +((re*re) + (im * im));
        end
        
        v1Norm = sqrt(somme);        
        v1NormInv = 1/ v1Norm;
        
        for cha=1:1:CHA
            
            V1(cha,1) = V1(cha,1)*v1NormInv;
        end
        
        DC=D;        
        % gemm(DH_D, DC, true, D, false);
        DH_D=DC'*D;
        
        % void gemm(hoNDArray<T>& C, const hoNDArray<T>& A, bool transA,  const hoNDArray<T>& B, bool transB);
        % if transA==true, C = A'*B
        % if transB==true, C=A*B'
        % if both are true, C=A'*B'
        
        for po = 1:1:power
            
            % gemm(V, DH_D, false, V1, false);
            V=DH_D*V1;
            
            somme = 0;
            
            V1=V;
            
            for cha=1:1:CHA
                c = V1(cha,1);
                re = real(c);
                im = imag(c);
                somme = somme +((re*re) + (im * im));
            end
            
            v1Norm = sqrt(somme);            
            v1NormInv = 1.0 / v1Norm;            
            V1(:,1) = V1(:,1)*v1NormInv;
            
        end
        
        % gemm(U1, D, false, V1, false);
        U1=D*V1;
        
        phaseU1 = sum(U1(:));
        
        phaseU1 = phaseU1/ abs(phaseU1);
        
        c = real(phaseU1);
        d = imag(phaseU1);
        
        for cha=1:1:CHA
                       
            v=V1(cha);
            
            a=real(v);
            b=imag(v);
          
            V1(cha,1)=  a*c + b*d +1j*( a*d - b*c);
        end
                  
        coilMap(ro,e1,:)=V1(:,1);
        
    end
end

return

