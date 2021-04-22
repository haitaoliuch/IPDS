clear;clc;close all;
Lambda    =[420];
interval = 10;         
strInData = sprintf('../Lambda%d.mat',Lambda(L)); 
load(strInData,'Test_instance')
Tn = 60;
mu = 10;
DT=60; M = 2;  V =11;  LM =1e6;
        for index =1:300
            instance =Test_instance{1,index};
            p=[]; a=[]; r=[];
            for t=1:Tn
                p =[p;cell2mat(instance{4,t})'];
                a =[a;instance{2,t}'];
                r =[r;ones(instance{1,t},1)*t*mu];
            end
            J = length(p); % 100;
            
            d = [1:V]'*DT; % vehicle departure times
            Q =1000*ones(V,1); % Vehicle capacity
            %%---- Decision variables
            X =binvar(J+1,J+1,'full');         % x_{ij}
            Y =binvar(M,J);           % y_{mi}
            D =sdpvar(J,1);           % departure time
            S =sdpvar(J,1);           % STARTING time of processing
            C =sdpvar(J,1);           % completion time of each item
            A =binvar(J,V);         % auxiliary variable
            %% Constraints:
            Con=[];
            for m=1:M
                Con=[Con,sum(Y(m,:))<=1]; 
            end

            for j=1:J
                Con=[Con, sum(X(j,:))-X(j,j)==1];
                Con=[Con,sum(Y(:,j)) + sum(X(1:J,j))-X(j,j) ==1]; % Each job is either the first to be processed on a machine or succeeds another one
                Con = [Con, S(j)>= a(j)]; % starting time greater than release time
                for m=1:M
                    Con = [Con,C(j)>= Y(m,j)*p(j)];
                end
            end

            for j=1:J
                for i=1:J
                    if i~=j
                        Con=[Con,C(j)>=C(i)+p(j)-LM*(1-X(i,j))];
                        Con=[Con,S(j)>=C(i)-LM*(1-X(i,j))]; % NEW
                    end      
                end
            end
            
            for j=1:J
                Con = [Con, C(j)>= S(j)+p(j),D(j)>=C(j)]; % departure time greater than completion time
                Con = [Con, d(V,1)>=C(j)]; 
                
                Con=[Con,sum(A(j,:))==1]; %Departure time:  each order should be assigned to one vehicle
                Con=[Con, A(j,:)*d ==D(j)]; % Assign a departure time
            end

            for v=1:V
                Con=[Con,sum(A(:,v))<=Q(v)];
            end
            %% Objective
            obj = sum(D-a)/J;
            ops = sdpsettings('verbose',1,'solver','gurobi');
            result  = optimize(Con,obj,ops);
            if result.problem== 0
                disp('Processing schedule:') 
                value(obj)
            else
                disp('have error');
            end
            %[result]=Optimization(Test_instance{1,i});
        end





