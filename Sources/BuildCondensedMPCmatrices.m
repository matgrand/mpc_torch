function condensedMatrices=BuildCondensedMPCmatrices(SystemModel,R,Q,F,f,N,S)

%Extracting Matrices from the strucutre (for code readibility only)
A=SystemModel.A;
B=SystemModel.B;
M=SystemModel.M;
C=SystemModel.C;

call_Q=kron(eye(N),Q);
call_R=kron(eye(N),R);

call_F=kron(eye(N),F);
call_f=kron(ones(N,1),f);

call_AC=[];
call_BC=[];   
call_MC=[];
x_term_r=[];

for ind_raw=1:N
    % call_A
    call_AC=[call_AC;
              C*A^(ind_raw)];
    % call_B
    
    call_B_new_col=[];
    for ind_col=1:N
        if ind_col<=ind_raw
            call_B_new_col=[call_B_new_col, C*A^(ind_raw-ind_col)*B];
        else
            call_B_new_col=[call_B_new_col, zeros(size(C*B))];
        end
    end
    
    call_BC=[call_BC;
              call_B_new_col];
          
    % call_M
    call_M_new_col=[];
    for ind_col=1:N
        if ind_col<=ind_raw
            call_M_new_col=[call_M_new_col, C*A^(ind_raw-ind_col)*M];
        else
            call_M_new_col=[call_M_new_col, zeros(size(C*M))];
        end
    end
    
    call_MC=[call_MC;
              call_M_new_col];
end

%% Adding Terminal Cost
if exist('S') && ~isempty(S)
    n=size(A,1);
    m=size(B,2);
    p=size(C,1);
    d=size(M,2);
    
    call_Q=[call_Q(1:end-p,1:end-p),        zeros((N-1)*p,n)
            zeros(n,(N-1)*p),             S];
    
    % call_A    
    call_AC=[call_AC(1:end-p,:);
            A^(N)];  
      
    % call_B    
    call_B_new_col=[];
    for ind_col=1:N
        if ind_col<=N
            call_B_new_col=[call_B_new_col, A^(N-ind_col)*B];
        else
            call_B_new_col=[call_B_new_col, zeros(size(B))];
        end
    end
    
    call_BC=[call_BC(1:end-p,1:end-m), zeros((N-1)*p,m);
            call_B_new_col];
          
    % call_M
    call_M_new_col=[];
    for ind_col=1:N
        if ind_col<=N
            call_M_new_col=[call_M_new_col, A^(N-ind_col)*M];
        else
            call_M_new_col=[call_M_new_col, zeros(size(M))];
        end
    end
    
    call_MC=[call_MC(1:end-p,1:end-d), zeros((N-1)*p,d);
              call_M_new_col];   

    x_term_r=zeros(n,1);
end


   % Building Output Structure:
   condensedMatrices.call_AC=call_AC;
   condensedMatrices.call_BC=call_BC;
   condensedMatrices.call_MC=call_MC;
   condensedMatrices.call_Q=call_Q;
   condensedMatrices.call_R=call_R;
   condensedMatrices.call_F=call_F;
   condensedMatrices.call_f=call_f;
   condensedMatrices.x_term_r=x_term_r;
    
end






