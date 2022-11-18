function [DarkTime,DarkTimeCI, BrightTime, BrightTimeCI]=CalculateTauDark(data,exposure,minDarkTime)
%% Taking time information from the data
frame=data(:,2);
[count, time]=histc(frame,min(frame):max(frame));
time=((1:length(count))*exposure)';
% figure;
% plot(time,count);
% xlabel('Time (s)');
% ylabel('Events');

%% calculate dark time by fitting CDF with F(tau_Dark)=(1-exp(tau_Dark/<tau_Dark>))a+b
time_events=find(count~=0);
temp_offtime=diff(time_events);
% offtime=zeros(length(temp_offtime),1);
for i=(1:length(temp_offtime))
    if temp_offtime(i) > minDarkTime;
        offtime(i)=temp_offtime(i);
    else
        offtime(i)=0;
    end
end
offtime=nonzeros(offtime)*exposure;
%tauD=mean(offtime);
%stdtauD=std(offtime);
try
    [DarkTime, DarkTimeCI] = mle(offtime,'dist','exp');
catch
    fprintf('Data is not enough to calculate DarkTime for statistics.\n')
    DarkTime = 0;
    DarkTimeCI = [0,0];
end

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
try
    time_noevents=find(count==0);
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

    try
        [BrightTime, BrightTimeCI] = mle(ontime,'dist','exp');
    catch
        fprintf('Data is not enough to calculate BrightTime for statistics.\n');
        BrightTime = 0;
        BrightTimeCI = [0,0];
    end
catch
    fprintf('Data is not enough to calculate BrightTime for statistics.\n');
    BrightTime = 0;
    BrightTimeCI = [0,0];
end

% figure;
% h2=cdfplot(ontime);
% hold on
% xgrid2=linspace(0,h2.XData(end-1),100);
% cdf_fit2=1-exp(-(xgrid2)/DarkTime(1));
% plot(xgrid2,cdf_fit2);
% xlabel('Bright time (s)');
% ylabel('Frequency');
% title('Bright time (s)');

