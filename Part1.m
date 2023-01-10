clc
clear all
close all

tic

input = [0.02	2	2	0.3	0.0666666666666667	0.0666666666666667	1.2	10	1.5	9.5	0.5	2];
Nr_exp = 1;

%% System Elements
for iiii = 1:Nr_exp

c = [input(iiii,2) input(iiii,3)]; %Nr. of Machines per Station
k = [input(iiii,12) input(iiii,12)]; %Buffers Capacity

Service = [input(iiii,5) input(iiii,6)]; %Prod. Rate [Secs/Part] of each single Machine (in each Station)
Arrival = 25; %Secs/Part -> Arrival rate to B1
Startup = [input(iiii,1) input(iiii,1)];


mu = 1./Service;
lambda = 1./Arrival;
delta = 1./Startup;

M1 = [0:1:c(1)]; 
M2 = [0:1:c(2)];
B1 = [0:1:k(1)];
B2 = [0:1:k(2)];

nu = lambda + c(1)*(mu(1)+delta(1)) + c(2)*(mu(2)+delta(2)); %Lippman

%% Cost Parameters

w_sup = [input(iiii,10) input(iiii,10)]; %Start-up Cost per each Machine
w_idle = [input(iiii,9) input(iiii,9)]; %Idle Cost per each Machine
w_busy = [input(iiii,8) input(iiii,8)]; %Service Cost per each Machine
w_holding = [input(iiii,11) input(iiii,11)]; %Holding Cost per each Buffer

rew = [0 0]; %Production Reward

rho = 0.80; %Discount Rate

%% System States

System = {B1, B2, M1, M2};
System_Size = length(M1)*length(M2)*length(B1)*length(B2);
States = Allcomb(System{:});
        
%% Matrix P
 
P = zeros(System_Size);
Pdelta = zeros(System_Size);
 
eta = nu + rho;

for ii = 1:System_Size
    
    Prod_rate(ii,1) = mu(1)*States(ii,3);
    Prod_rate(ii,2) = mu(2)*States(ii,4);
    StartUp_rate(ii,1) = delta(1);
    StartUp_rate(ii,2) = delta(2);
    
    % Transition rates for the Arrival/Departures
    
    if States(ii,1) < k(2)
        P(ii, ii + length(B2)*length(M1)*length(M2)) = lambda; %Arrival to Buffer 1
    end
    
    if States(ii,1) > 0 && States(ii,2) < k(2)
        P(ii, ii - length(B2)*length(M1)*length(M2)) = (max(Prod_rate(ii,1), 0)); %Departure from Buffer 1 
        P(ii, ii + length(M1)*length(M2)) = (max(Prod_rate(ii,1), 0)); %Arrival to Buffer 2 
    end
    
    if States(ii,2) > 0
        P(ii, ii - length(M1)*length(M2)) = (max(Prod_rate(ii,2), 0)); %Departure from Buffer 2
    end
    
    % Transition rates for the Startup
    
    if States(ii,3) < c(1)
        P(ii, ii + length(M2)) = StartUp_rate(ii,1); %Startup Completed in Station 2
    end
    
    if States(ii,4) < c(2)
        P(ii, ii + 1) = StartUp_rate(ii,2); %Startup Completed in Station 2
    end

        
end

 for ii = 1:System_Size
    
    for jj = 1:System_Size
        
        if jj ~= ii
            
        Pdelta(ii,jj) = P(ii,jj).*(1/eta);
        
        end
        
    end
    
    Pdelta(ii,ii) = 1 - sum(Pdelta(ii,:)); %Matrix Consistency
     
 end

bbb = [];
 
for ii = 1:System_Size 
    
    bbb(ii) = sum(Pdelta(ii,:)); % Check Matrix Consistency
    
end

if sum(bbb) ~= System_Size 
    
   error('Matrix P NOT consistent')
   
end

% Matrix Pdelta represents P-delta -> P-delta = nu/eta * P' = nu/eta * P/nu = P/eta 

%% Matrix C    
 
C = zeros(System_Size);
eta = nu + rho;

for ii = 1:System_Size
     
    for jj = 1:System_Size
 
        c_holding = w_holding(1)*States(jj,1) + w_holding(2)*States(jj,2);
        
        c_service_1 = min(States(jj,1),States(jj,3))*w_busy(1) + (States(jj,3) - min(States(jj,1),States(jj,3)))*w_idle(1);
        c_service_2 = min(States(jj,2),States(jj,4))*w_busy(2) + (States(jj,4) - min(States(jj,2),States(jj,4)))*w_idle(2);
        
        C(ii,jj) = (c_holding + c_service_1 + c_service_2)/eta;
        
    end
    
end

C = C.';
 
% C Represents the stage cost -> already divided by eta

%% Matrix SU (Startup)

SU = zeros(System_Size);

for ii = 1:System_Size
    
    up_1 = States(ii,3);
    up_2 = States(ii,4);
    
    for jj = 1:System_Size
        
        a_1 = States(jj,3);
        a_2 = States(jj,4);
        
        if a_1 > up_1
            csu_1 = (a_1 - up_1)*w_sup(1);
        else
            csu_1 = 0;
        end
        
        if a_2 > up_2
            csu_2 = (a_2 - up_2)*w_sup(2);
        else
            csu_2 = 0;
        end
        
        SU(ii,jj) = csu_1 + csu_2;
        
    end
    
end

%% TotCost

CTOT = zeros(System_Size);

for ii = 1:System_Size
    
    for jj = 1:System_Size
        
        CTOT(ii,jj) = Pdelta(ii,jj)*(C(ii,jj) + SU(ii,jj));
        
    end
end


%% Send to Excel

alph = 1/System_Size;
alphavec = alph*ones(System_Size,1);
%mkdir iiii
filename = ['data',num2str(iiii),'.xlsx'];

writematrix(System_Size,filename,'Sheet',1)
writematrix(c,filename,'Sheet',2)
writematrix(k,filename,'Sheet',3)

writematrix(alphavec,filename,'Sheet',9)
writematrix(Pdelta,filename,'Sheet',10)
writematrix(CTOT,filename,'Sheet',11)

MachCplex = States(:,3:4);

writematrix(MachCplex,filename,'Sheet',12)

XCplex = [];
OF = [];

writematrix(XCplex,filename,'Sheet',13)
writematrix(OF,filename,'Sheet',15)

%%

COL = System_Size;

if COL <= 26                        % [A..Z]
        CHAR = char(mod(COL-1,26)+1+64);
    elseif COL <= 702                   % [AA..ZZ]
        COL = COL-26;    
        CHAR1 = char(floor((COL-1)/26)+1+64);
        CHAR0 = char(mod(COL-1,26)+1+64);
        CHAR = [CHAR1 CHAR0];
    elseif COL <= 16384                 % [AAA..XFD]
        COL = COL-702; 
        CHAR2 = char(floor((COL-1)/676)+1+64);
        COL=COL-(floor((COL-1)/676))*676;
        CHAR1 = char(floor((COL-1)/26)+1+64);
        CHAR0 = char(mod(COL-1,26)+1+64);
        CHAR = [CHAR2 CHAR1 CHAR0];
    else
        disp('Column does not exist in Excel!');
end




Details_1 = ['Sheet9','!A1:A',num2str(System_Size)];
Details_2 = ['Sheet10','!A1:',CHAR,num2str(System_Size)];
Details_3 = ['Sheet11','!A1:',CHAR,num2str(System_Size)];
Details_4 = ['Sheet12','!A1:B',num2str(System_Size)];
Details_5 = ['Sheet13','!A1:',CHAR,num2str(System_Size)];
Details_6 = ['Sheet15','!A1:A1'];

writematrix(Details_1,filename,'Sheet',4)
writematrix(Details_2,filename,'Sheet',5)
writematrix(Details_3,filename,'Sheet',6)
writematrix(Details_4,filename,'Sheet',7)
writematrix(Details_5,filename,'Sheet',8)
writematrix(Details_6,filename,'Sheet',14)

end

time_elapsed = toc;
