function condensedMatrices=BuildCondensedMPCmatrices(SystemModel,R,Q,F,f,N,S)

%Extracting Matrices from the strucutre (for code readibility only)
A=SystemModel.A;
B=SystemModel.B;
M=SystemModel.M;
C=SystemModel.C;

cQ=kron(eye(N),Q);
cR=kron(eye(N),R);

cF=kron(eye(N),F);
cf=kron(ones(N,1),f);

cAC=[];
cBC=[];   
cMC=[];
x_term_r=[];

for ind_raw=1:N
    % cA
    cAC=[cAC;
              C*A^(ind_raw)];
    % cB
    
    cB_new_col=[];
    for ind_col=1:N
        if ind_col<=ind_raw
            cB_new_col=[cB_new_col, C*A^(ind_raw-ind_col)*B];
        else
            cB_new_col=[cB_new_col, zeros(size(C*B))];
        end
    end
    
    cBC=[cBC;
              cB_new_col];
          
    % cM
    cM_new_col=[];
    for ind_col=1:N
        if ind_col<=ind_raw
            cM_new_col=[cM_new_col, C*A^(ind_raw-ind_col)*M];
        else
            cM_new_col=[cM_new_col, zeros(size(C*M))];
        end
    end
    
    cMC=[cMC;
              cM_new_col];
end

%% Adding Terminal Cost
if exist('S') && ~isempty(S)
    n=size(A,1);
    m=size(B,2);
    p=size(C,1);
    d=size(M,2);
    
    cQ=[cQ(1:end-p,1:end-p),        zeros((N-1)*p,n)
            zeros(n,(N-1)*p),             S];
    
    % cA    
    cAC=[cAC(1:end-p,:);
            A^(N)];  
      
    % cB    
    cB_new_col=[];
    for ind_col=1:N
        if ind_col<=N
            cB_new_col=[cB_new_col, A^(N-ind_col)*B];
        else
            cB_new_col=[cB_new_col, zeros(size(B))];
        end
    end
    
    cBC=[cBC(1:end-p,1:end-m), zeros((N-1)*p,m);
            cB_new_col];
          
    % cM
    cM_new_col=[];
    for ind_col=1:N
        if ind_col<=N
            cM_new_col=[cM_new_col, A^(N-ind_col)*M];
        else
            cM_new_col=[cM_new_col, zeros(size(M))];
        end
    end
    
    cMC=[cMC(1:end-p,1:end-d), zeros((N-1)*p,d);
              cM_new_col];   

    x_term_r=zeros(n,1);
end


   % Building Output Structure:
   condensedMatrices.cAC=cAC;
   condensedMatrices.cBC=cBC;
   condensedMatrices.cMC=cMC;
   condensedMatrices.cQ=cQ;
   condensedMatrices.cR=cR;
   condensedMatrices.cF=cF;
   condensedMatrices.cf=cf;
   condensedMatrices.x_term_r=x_term_r;
    
end






