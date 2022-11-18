n_syn=length(dwell);
for j=1:n_syn
    if isempty(dwell{j}) % if cell is not empty
        dwell_time_s{j}=dwell{j};
    else
dwell_m=dwell{j};
ID_particle=unique(dwell_m(:,1)); % find particle IDs
n_particle=length(ID_particle); % find number of different particles
    for k=1:n_particle
s_dwell=sum(dwell_m(find(dwell_m(:,1)==ID_particle(k)),2));
dwell_sum(k,1)=ID_particle(k);
dwell_sum(k,2)=s_dwell;
dwell_time_s{j}=dwell_sum;
    end    
    end
end
dwell_t=cell2mat(dwell_time_s');
figure
hist(dwell_t(:,2),48)
clear ID_particle dwell_m dwell_sum j k n_particle n_syn s_dwell