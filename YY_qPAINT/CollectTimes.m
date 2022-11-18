function [DarkTime,BrightTime]=CollectTimes(data,exposure)
%% Taking time information from the data
frame=data(:,2);
[count, time]=histc(frame,1:max(frame));
time=((1:length(count))*exposure)';
% figure;
% plot(time,count);
% xlabel('Time (s)');
% ylabel('Events');

%% calculate dark time by fitting CDF with F(tau_Dark)=(1-exp(tau_Dark/<tau_Dark>))a+b
time_events=find(count~=0);
temp_offtime=diff(time_events);
for i=(1:length(temp_offtime))
    if temp_offtime(i) > 1
        offtime(i)=temp_offtime(i);
    else
        offtime(i)=0;
    end
end
DarkTime=nonzeros(offtime)*exposure;
%tauD=mean(offtime);
%stdtauD=std(offtime);
%DarkTime = mle(offtime,'cdf','exp');
% figure;
%hist(offtime);
% h=cdfplot(offtime);
% hold on
% xgrid=linspace(0,h.XData(end-1),100);
% cdf_fit=1-exp(-(xgrid)/BrightTime(1)); % CDF curve for exponential distribution
% plot(xgrid,cdf_fit);
% xlabel('Dark time (s)');
% ylabel('Frequency');
% title('Dark time (s)');


%% for tauBright
time_noevents=find(count==0);
temp_ontime=diff(time_noevents);
for i=(1:length(temp_ontime))
    if temp_ontime(i) > 1
        ontime(i)=temp_ontime(i);
    else
        ontime(i)=0;
    end
end
BrightTime=nonzeros(ontime)*exposure;
%BrightTime = mle(ontime,'cdf','exp');
% figure;
% h2=cdfplot(ontime);
% hold on
% xgrid2=linspace(0,h2.XData(end-1),100);
% cdf_fit2=1-exp(-(xgrid2)/DarkTime(1));
% plot(xgrid2,cdf_fit2);
% xlabel('Bright time (s)');
% ylabel('Frequency');
% title('Bright time (s)');

