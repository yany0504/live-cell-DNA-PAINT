function [DarkTime,Darkstd, BrightTime,Brightstd]=calc_tau(intensity,exposure)
time_events=find(intensity~=0);
temp_offtime=diff(time_events);
for i=(1:length(temp_offtime))
    if temp_offtime(i) > 1
        offtime(i)=temp_offtime(i);
    else
        offtime(i)=0;
    end
end
offtime=nonzeros(offtime)*exposure;
DarkTime=mean(offtime);
Darkstd=std(offtime);
% try
%     [DarkTime,Darkstd] = mle(offtime,'distribution','exponential');
% catch
%     fprintf('Data is not enough to calculate DarkTime for statistics.\n')
%     [DarkTime,Darkstd] = [0,0];
% end

%% for tauBright
time_noevents=find(intensity==0);
temp_ontime=diff(time_noevents);
% ontime=zeros(length(temp_ontime),1);
for i=(1:length(temp_ontime))
    if temp_ontime(i) > 1
        ontime(i)=temp_ontime(i);
    else
        ontime(i)=0;
    end
end
ontime=nonzeros(ontime)*exposure;
BrightTime=mean(ontime);
Brightstd=std(ontime);
% try
%     [BrightTime,Brightstd] = mle(ontime,'distribution','exponential');
% catch
%     fprintf('Data is not enough to calculate BrightTime for statistics.\n');
%     [BrightTime,Brightstd] = [0,0];
% end
