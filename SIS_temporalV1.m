function [state_vector] = SIS_temporalV1(DataSet,DataFlag,Infect_source,b,u)
%individual scale SIS simulation in temporal network 
%    DataSet: Input temporal network 
%    DataFalg£º a flag for input data format
%               catb = unweighted link triplet(i,j,t)
%               catw =  weighted link triplet(i,j,t,weight)
%               TimA = time adjaceny matrix of size (N^2)*T
%    Infect_soure:  primary case£¬eg.[26,27] reprensets NO.26 and NO.27 are
%                   primary cases in certain simulation
%    b: infectios probability
%    u: recovery rate
 
if DataFlag == 'catb' 
    population = length(union(DataSet(:,1),DataSet(:,2)));  %determin the number of induvidual in the network
    
    Infec_initial = zeros(1,population);                    %set primary case
    Infec_initial(1,Infect_source) = 1;
    timestamps = max(DataSet(:,3));
    
    state_vector = zeros(timestamps+1,population);          %store each individual's state at each timestep
    state_vector(1,:) = Infec_initial;
    
%     s = rng;
%     rng(s);
    for t = 1:1:timestamps
       state_vector(t+1,:) = state_vector(t,:);             %infectious individual recover at rate u
       infSet = find(state_vector(t,:)==1);
       recSet = rand(size(infSet))<u;
       state_vector(t+1,infSet(recSet)) = 0;

       link_row = find(DataSet(:,3) == t);                %epidemics get spread through contacts between individuals
       for i = 1:1:length(link_row)
           p1Index = DataSet(link_row(i),1);
           p2Index = DataSet(link_row(i),2);
           if state_vector(t,p1Index) ~= state_vector(t,p2Index)%infection happens between S and I
               if rand<b                                              
                  if state_vector(t,p1Index) == 1
                      state_vector(t+1,p2Index) = 1;
                  else
                      state_vector(t+1,p1Index) = 1;
                  end
               end
           end
       end 
    end  
elseif DataFlag == 'catw'
    population = length(union(DataSet(:,1),DataSet(:,2)));  %determin the number of induvidual in the network
    
    Infec_initial = zeros(1,population);%set primary case
    Infec_initial(1,Infect_source) = 1;
    timestamps = max(DataSet(:,3));
    
    state_vector = zeros(timestamps+1,population);%store each individual's state at each timestep
    state_vector(1,:) = Infec_initial;
    
%     s = rng;
%     rng(s);
    for t = 1:1:timestamps
       state_vector(t+1,:) = state_vector(t,:);             %infectious individual recover at rate u
       infSet = find(state_vector(t,:)==1);
       recSet = rand(size(infSet))<u;
       state_vector(t+1,infSet(recSet)) = 0;

       link_row = find(DataSet(:,3) == t);                %epidemics get spread through contacts between individuals
       for i = 1:1:length(link_row)
           p1Index = DataSet(link_row(i),1);
           p2Index = DataSet(link_row(i),2);
           weight = DataSet(link_row(i),4);
           if state_vector(t,p1Index) ~= state_vector(t,p2Index)%infection happens between S and I
               if rand<(1-(1-b)^weight)                                              
                  if state_vector(t,p1Index) == 1
                      state_vector(t+1,p2Index) = 1;
                  else
                      state_vector(t+1,p1Index) = 1;
                  end
               end
           end
       end 
    end                 
elseif DataFlag == 'TimA'
    [row, col] = size(DataSet);
    population = sqrt(row);
    timestamps = col;
    
    Infec_initial = zeros(1,population); %set primary case
    Infec_initial(1,Infect_source) = 1;
    state_vector = zeros(timestamps+1,population); %store each individual's state at each timestep
    state_vector(1,:) = Infec_initial;
    
    for t = 1:1:timestamps
        state_vector(t+1,:) = state_vector(t,:); %infectious individual get recovered
        infSet = find(state_vector(t,:)==1);
        recSet = rand(size(infSet))<u;
        state_vector(t+1,infSet(recSet)) = 0;
        
%         s = rng;
%         rng(s);
        snapshot = reshape(DataSet(:,t),[population,population])';
        S = state_vector(t,:) == 0; %S represent susceptible individual in the previous timestamp
        danger_contact = snapshot(S,:)*state_vector(t,:)'; % calculate  number of contacts between each susceptible individual and the infected individual
        P = zeros(size(danger_contact)); %P represent each susceptible indivudual's probability of being infected
        [r,~] = size(danger_contact);
        for i = 1:1:r
           P(i,1) = 1-(1-b)^danger_contact(i,1); 
        end
        RAN = rand(size(danger_contact));
        state_vector(t+1,S(RAN' < P')) = 1;
    end
end     
end