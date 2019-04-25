function [state_vector_static] = SIS_static(adj,Infect_source,b,u,t)
%individual simulate on static network
%    DataSet: only support adj format now
%    b: infectios probability
%    u: recovery rate 
%    t: total timestep simulated
%    Infect_soure:  primary case��eg.[26,27] reprensets NO.26 and NO.27 are
%                   primary cases in certain simulation

%% ��Ҫ��ʼ���õĲ���
population = length(adj);                                   
timestamps = t;
Infec_initial = zeros(1,population);                   %set each individual's initial state
Infec_initial(1,Infect_source) = 1;
state_vector_static = zeros(timestamps+1,population);  %store each individual's state in every timestep
state_vector_static(1,:) = Infec_initial;

% s = rng;
% rng(s);
for t = 1:1:timestamps
   state_vector_static(t+1,:) = state_vector_static(t,:) ;
   I = find(state_vector_static(t,:) == 1);
   RAN2 = rand(size(I));
   ItoS = RAN2<u;
   state_vector_static(t+1,I(ItoS)) = 0;
    
   S = find(state_vector_static(t,:) == 0);  % find susceptible population
   danger_contact = adj(S,:)* state_vector_static(t,:)';   % 
   P = zeros(size(danger_contact)); %P��ÿ���׸нڵ�ĽӴ���Ⱦ���ʣ�������,sizeͬS
   [r,~] = size(danger_contact);
   for i = 1:1:r
       P(i,1) = 1-(1-b)^danger_contact(i,1);
   end
   RAN = rand(size(danger_contact)); %�������
   StoI = RAN<P; %��Ⱦ��S������S������,size<S
   state_vector_static(t+1,S(StoI)) = 1;
end
end