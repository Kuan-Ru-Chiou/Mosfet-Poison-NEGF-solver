 
clear all 
  
%Constants (all MKS, except energy which is in eV) 
hbar=1.06e-34;
q=1.6e-19;
epsil=10*8.85E-12;
kT=.025; 
m=.25*9.1e-31;
n0=2*m*kT*q/(2*pi*(hbar^2)); 
  
%inputs 
a=3e-10;
t=(hbar^2)/(2*m*(a^2)*q);
beta=q*a*a/epsil; 
Ns=15;Nc=70;Np=Ns+Nc+Ns;XX=a*1e9*[1:1:Np]; 
mu=.318;
Nd=2*((n0/2)^1.5)*FD_int_num(mu/kT,1/2,1e-6,200);
Nd=Nd*[ones(Ns,1);.5*ones(Nc,1);ones(Ns,1)];
  
%d2/dx2 matrix for Poisson solution 
D2=-(2*diag(ones(1,Np)))+(diag(ones(1,Np-1),1))+...
(diag(ones(1,Np-1),-1)); 
D2(1,1)=-1;D2(Np,Np)=-1;%zero field condition 
  
%Hamiltonian matrix FDM 
T=(2*t*diag(ones(1,Np)))-(t*diag(ones(1,Np-1),1))-...
(t*diag(ones(1,Np-1),-1)); 


%energy grid 
NE=501;E=linspace(-.25,.5,NE);dE=E(2)-E(1),zplus=i*1e-12; 
 
  

%gate voltage
% Vgg = -0.3;
% Ug = Vgg*ones(Np,1);



%initial guess for U 
U=[zeros(Ns,1);.2*ones(Nc,1);zeros(Ns,1)]; 


%voltage bias steps 
NV=5;VV=linspace(0,.25,NV);dV=VV(2)-VV(1);
for kV=1:5
V=VV(kV) 
f1=n0*log(1+exp((mu-E)./kT));f2=n0*log(1+exp((mu-V-E)./kT)); 
%   f11=n0*FD_int_num((mu-E(1))/kT,0,1e-6,200);%check numerical and exact
%   f22=n0*FD_int_num((mu-V-E(1))/kT,0,1e-6,200);%check numerical and exact

            in=10; 
            while in>0.01
                        sig1=zeros(Np);sig2=zeros(Np);rho=zeros(Np); II1=0;II2=0; II3=0;
                        
                      
  
                        %%Hamiltonian matrix FEM method
                        
                        H11 = (hbar^2)/(2*m*(a^2)*q);
                        H12 = -(hbar^2)/(2*m*(a^2)*q);
                        S11 = 1/3;
                        S12 = 1/6;
                        S = (2*S11*diag(ones(1,Np)))+(S12*diag(ones(1,Np-1),1))+(S12*diag(ones(1,Np-1),-1));
                        S(1,1) = 1/3;
                        S(Np,Np) = 1/3;
                        HH = (2*H11*diag(ones(1,Np)));
                        HH(1,1) = H11;
                        HH(Np,Np) = H11;
                        HH = HH+ (1/3)*diag([U(1:Np-1); 0]) + (1/3)*diag([0; U(2:Np)]);
                        HH =HH+(1/6)*diag(U(1:Np-1),1) + (1/6)*diag(U(1:Np-1),-1) +diag(H12*ones(Np-1,1),1)+diag(H12*ones(Np-1,1),-1);
    
                     
                        for k=1:NE 
                                    
                                    ck=1-((E(k)+zplus-U(1))/(2*t));ka=acos(ck); 
                                    % sig1(1,1)=-t*exp(i*ka);gam1=i*(sig1-sig1'); 
                                    sig1(1,1)=-1i*t*ka;gam1=i*(sig1-sig1');
                                    
                                    
                                    ck=1-((E(k)+zplus-U(Np))/(2*t));ka=acos(ck); 
%                                     sig2(Np,Np)=-t*exp(i*ka);gam2=i*(sig2-sig2'); 
                                     sig2(Np,Np)=-1i*t*ka;gam2=i*(sig2-sig2'); 
                                    
                                    sig1in = gam1*f1(k);
                                    sig2in = gam2*f2(k);
                                    
           
%                                     G=inv(((E(k)+zplus)*eye(Np))-T-diag(U)-sig1-sig2);%FDM
                                    G=inv(diag(E(k)+zplus).*S-HH-sig1-sig2);%FEM
                                    Gn = G'*(sig1in+sig2in)*G;
                                    A1=G'*gam1*G;A2=G'*gam2*G; 
                                    A = 1i*(G-G');
                                    AAA =A1+A2;
  
                                    rho=rho+(dE*((f1(k)*A1)+(f2(k)*A2))/(2*pi));
                                    II1 = II1+((q*q/hbar/2/pi)*dE)*(trace(sig1in*A-gam1*Gn));
                                    
                                    
                                    II2 = II2+((q*q/hbar/2/pi)*dE)*(trace(sig2in*A-gam2*Gn));
                                    
                                    II3 = II3+((q*q/hbar/2/pi)*dE)*trace(gam1*G*gam2*G')*(f1(k)-f2(k));
                                    
                                    
  
                                    %%%%%check internal                                     
%                                       HHH = T+diag(U);
%                                     for i =1:Np-1
%                                         
%                                         jjj(NE,i,kV) = (1i)*(trace(HHH(i,i+1)*Gn(i+1,i)-HHH(i+1,i)*Gn(i,i+1)));
%                                     end
%                                     
                                    
                                    
                                    
                                    
                                    
                                    
                        end 
                        
                                                n=(1/a)*real(diag(rho)); 
  
                                                
                        % Solve quasi-fermi-level using  bisection method
                        xup=100*ones(Np,1);
                        xlow=-100*ones(Np,1);
                        for kk=1:Np
                            ii_iter=0;
                            while ii_iter<40
%                                 res(kk)= n(kk)./(2*(n0/2)^1.5)-FD_int_num((xup(kk)+xlow(kk))*0.5,1/2,1e-8,3000);
                                res(kk)= n(kk)./(2*(n0/2)^1.5)-FD_int_num(((xup(kk)+xlow(kk))*0.5)-U(kk)/kT,1/2,1e-8,3000);
                                error=max(abs(res(kk)));
                                
                                if res(kk)>0
                                    xlow(kk)=(xup(kk)+xlow(kk))/2;
                                else
                                    xup(kk)=(xup(kk)+xlow(kk))/2;
                                end
                                
                                ii_iter=ii_iter+1;
                            end
                            if (max(abs(res(kk)))>1)
                                disp('anti_dmmuy exceeds iteration limit');
                            end
                            x(kk)=(xup(kk)+xlow(kk))*0.5;
                        end
                        
                                                
                                                
                                    Fn = x; %this term contain Ef/kT  quasi-fermi-level
                                    
                        %correction dU from Poisson 
                 
                        
                        D=zeros(Np,1);
                        for k=1:Np
%                             z=(Fn(k)-U(k))/kT;
                              z =Fn(k)-(U(k)/kT);
%                              D(k,1)=2*((n0/2)^1.5)*((Fhalf(z+.1)-Fhalf(z))/.1)/kT;% OK
                            
                            D(k,1)=2*((n0/2)^1.5)*FD_int_num(z,-1/2,1e-8,3000)/kT;  % later solve unit is 1/ev??
                          
                        end
                        
                                            
                        dN = 0.5*[n(1:Np-1)-Nd(1:Np-1); 0]+ 0.5*[0; n(2:Np)-Nd(2:Np)] + ((1/beta)*D2*U); %FEM
                        D =  0.5*[D(1:Np-1); 0]+ 0.5*[0 ;D(2:Np)];
                        
                        
%                        dN=n-Nd+((1/beta)*D2*U); %FDM
                       %  dU=(-beta)*(inv(D2-(beta*diag(D))))*dN;
                       
                       jjaco = sparse(Np,Np);
                       jjaco =  D2-(beta*diag(D));
                      %  dU = (-beta)*inv(jjaco)*dN;
                        dU = (-beta)*(jjaco\dN);
                                    %Brown Lindsay guauantee find newton root
                                    ddU = dU;
                                    dddU = abs(ddU);
                                    index1 = find(1.0<dddU & dddU< 3.7);
                                    ddU(index1) = sign(ddU(index1)).*(dddU(index1)).^(1/5);
                                    index2 = find(dddU >= 3.7);
                                    ddU(index2) = sign(ddU(index2)).*log(dddU(index2));
                                    dU = ddU;
                                    
                                    
                                    U=U+0.25*dU;
                        
                      
                        % Check for convergence 
%                      in=(max(max(abs(dN))))/(max(max(Nd))) 
                      in=max(abs(dU))
                        
                        
                        
                        
            end
            
            
            UU(:,kV)=U;
            nn(:,kV) = n;
            III1(:,kV) = II1;
            III2(:,kV) = II2;
            III3(:,kV) = II3;
                   
end


% hold on 
% plot(XX,J(:,NV)) 
%plot(XX,Fn)
% plot electron potential in the device
figure(1)
plot(XX,UU(:,NV)) 
% plot source to drain current
figure(2)
plot(VV,III1')
  
