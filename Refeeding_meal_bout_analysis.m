%Set up a time vectorx
FP.time = linspace(-10,40,60001);
%Make the first column your corrected data and the second the time vector:

FP.refeeding.B9422_refeeding(:,2) = FP.time;
%%

%   for x = 1:10
%       disp(x)
%   end
% 

%Convert the keydown values into minutes  and separate bout start and bout
%end into temporar variables:

B9422_refeeding_keydown = (B9422_refeeding_keydown ./ 60) - 15;
B9422_refeeding_keydown(2:2:end) = B9422_refeeding_keydown(2:2:end) - 0.16667;

FP.refeeding.B9422_boutstart = downsample(B9422_refeeding_keydown(:),2); %temporary variable 
FP.refeeding.B9422_boutend = downsample(B9422_refeeding_keydown(2:end),2);


%%
%Plot the refeeding session 1data and bout start/end lines:
figure(1)
%plot(FP.refeeding.B9422_refeeding(:,2),FP.refeeding.B9422_refeeding(:,1))
for i = 1:length(FP.refeeding.B9422_boutstart)
    xline(0, "LineWidth",2,"LineStyle","--")
    xline(FP.refeeding.B9422_boutstart(i),"LineWidth",1,"Color",'g')
    xline(FP.refeeding.B9422_boutend(i),"LineWidth",1,"Color",'r')
    xline(30, "LineWidth",2,"LineStyle","--")
end



%% refeeding activity bout quantification analysis:
%The following script will be dedicated to quantifying and analyzing the
%activity of a brain region at the beginning of a bout (meal) and compare
%it to the activity at the end of the bout.

%First, we will want to convert the pre-food access period processed and finalized trace of a
%session to its own zscore:
%Pre-food access recordings are 10 minutes long, because processed traces
%are 20 data points per second we want points up to 12001:

% DELETE THE TIME COLUMN IN YOUR FP.refeeding.ANIMALID_refeeding FILE OR
% YOU WILL GET AN ERROR.


FP.refeeding.zscore_pre.B9422_refeeding = zscore(FP.refeeding.B9422_refeeding(1:12001)); 

%The following will extract the standard deviation and mean of period -7min
%to -2min and we will use that standard deviation and mean of this set to normalize period
%where food is accessible: 

FP.refeeding.zscore_pre.B9422_refeeding_7to2_std = std(FP.refeeding.B9422_refeeding(3601:9601));
FP.refeeding.zscore_pre.B9422_refeeding_7to2_mean = mean(FP.refeeding.B9422_refeeding(3601:9601));

%Next we will convert the period where food is available to a zscore using
%the std and mean we extracted above: 

FP.refeeding.zscore_duringmeal.B9422_refeeding_during = (FP.refeeding.B9422_refeeding(12002:end) - FP.refeeding.zscore_pre.B9422_refeeding_7to2_mean) ./ FP.refeeding.zscore_pre.B9422_refeeding_7to2_std;

%We should have successfully normalized the food available period and post
%to a baseline period (pre food), we will now want to connect the two to
%see what it looks like:

FP.refeeding.zscore_finalized.B9422_refeeding = []; %Initialize the variable as an empty matrix

FP.refeeding.zscore_finalized.B9422_refeeding(1:12001) = FP.refeeding.zscore_pre.B9422_refeeding;
FP.refeeding.zscore_finalized.B9422_refeeding(12002:60001) = FP.refeeding.zscore_duringmeal.B9422_refeeding_during;

figure

plot(FP.time,FP.refeeding.zscore_finalized.B9422_refeeding)

FP.refeeding.zscore_finalized.B9422_refeeding = FP.refeeding.zscore_finalized.B9422_refeeding'
FP.refeeding.zscore_finalized.B9422_refeeding(:,2) = FP.time;

%%
figure(1)
hold on
for i = 1:length(FP.refeeding.B9422_boutstart)
    xline(0, "LineWidth",2,"LineStyle","--")
    xline(FP.refeeding.B9422_boutstart(i),"LineWidth",1,"Color",'g')
    xline(FP.refeeding.B9422_boutend(i),"LineWidth",1,"Color",'r')
    xline(30, "LineWidth",2,"LineStyle","--")

end


%% The following will extract the data 15s after a bout start and average the values, and do the same for 15s before a bout ends:

FP.bl = 300; %number of datapoints to look at before keydown
FP.seg_dur = 300; %number of datapoints to look at after keydown
t = [];
temp = [];
FP.temp_boutstart = FP.refeeding.B9422_boutstart; %temporary variable 
FP.temp_boutend = FP.refeeding.B9422_boutend;
FP.refeeding.afterboutstarts.B9422_refeeding_15safterstart = []; %necessary if you re-run this part of the code changing baseline or seg_dur
FP.refeeding.beforeboutends.B9422_refeeding_15sbeforeend = [];

%% Within bout interval

% The following will extract the data 15s after a bout start and average the values, and do the same for 15s before a bout ends:

FP.bl = 300; %number of datapoints to look at before keydown
FP.seg_dur = 300; %number of datapoints to look at after keydown


FP.refeeding.afterboutstarts.B9422_refeeding_15safterstart = []; %necessary if you re-run this part of the code changing baseline or seg_dur
FP.refeeding.beforeboutends.B9422_refeeding_15sbeforeend = [];


% 15 SECONDS AFTER BOUT STARTS

for i = 1:length(FP.refeeding.B9422_boutstart)
    
       
       t = [];
       temp = [];

       t(1,:)  = find(FP.refeeding.zscore_finalized.B9422_refeeding(:,2) >= FP.refeeding.B9422_boutstart(i)); %find all FP data indices that are greater than or equal to the timestamp of when a bout starts. In this case, the time vector is on the second column because we manually put the FP.time on the variable in the second column.
       temp(1,:) = t(1); %take the first one 
       t = [];

       FP.refeeding.afterboutstarts.B9422_refeeding_15safterstart(i,:) = FP.refeeding.zscore_finalized.B9422_refeeding(temp(1):temp(1)+FP.seg_dur,1);
       temp = [];

       FP.refeeding.afterboutstarts.B9422_15safterstart_avg = mean(FP.refeeding.afterboutstarts.B9422_refeeding_15safterstart, 2)
        
end


%15 SECONDS BEFORE BOUT ENDS

for i = 1:length(FP.refeeding.B9422_boutend)
       
       t = [];
       temp = [];

       t(1,:)  = find(FP.refeeding.zscore_finalized.B9422_refeeding(:,2) >= FP.refeeding.B9422_boutend(i)); %find all FP data indices that are greater than or equal to the timestamp of when a bout starts. In this case, the time vector is on the second column because we manually put the FP.time on the variable in the second column.
       temp(1,:) = t(1); %take the first one 
       t = [];

       FP.refeeding.beforeboutends.B9422_refeeding_15sbeforeend(i,:) = FP.refeeding.zscore_finalized.B9422_refeeding(temp(1)-FP.bl:temp(1),1);
       temp = [];

       FP.refeeding.beforeboutends.B9422_15sbeforeend_avg = mean(FP.refeeding.beforeboutends.B9422_refeeding_15sbeforeend, 2)
        
end

%% Inter bout interval
%15 SECONDS BEFORE BOUT STARTS

FP.refeeding.beforeboutstarts.B9422_refeeding_15sbeforetart = []; %necessary if you re-run this part of the code changing baseline or seg_dur
FP.refeeding.afterboutends.B9422_refeeding_15safterend = [];


for i = 1:length(FP.refeeding.B9422_boutstart)
       
       t = [];
       temp = [];

       t(1,:)  = find(FP.refeeding.zscore_finalized.B9422_refeeding(:,2) >= FP.refeeding.B9422_boutstart(i,1)); %find all FP data indices that are greater than or equal to the timestamp of when a bout starts. In this case, the time vector is on the second column because we manually put the FP.time on the variable in the second column.
       temp(1,:) = t(1); %take the first one 
       t = [];

       FP.refeeding.beforeboutstarts.B9422_refeeding_15sbeforestart(i,:) = FP.refeeding.zscore_finalized.B9422_refeeding(temp(1)-FP.bl:temp(1),1);
       temp = [];

       FP.refeeding.beforeboutstarts.B9422_15sbeforestart_avg = mean(FP.refeeding.beforeboutstarts.B9422_refeeding_15sbeforestart, 2)
        
end


%15 SECONDS AFTER BOUT ENDS

for i = 1:length(FP.refeeding.B9422_boutend)
       
       t = [];
       temp = [];

       t(1,:)  = find(FP.refeeding.zscore_finalized.B9422_refeeding(:,2) >= FP.refeeding.B9422_boutend(i,1)); %find all FP data indices that are greater than or equal to the timestamp of when a bout starts. In this case, the time vector is on the second column because we manually put the FP.time on the variable in the second column.
       temp(1,:) = t(1); %take the first one 
       t = [];

       FP.refeeding.afterboutends.B9422_refeeding_15safterend(i,:) = FP.refeeding.zscore_finalized.B9422_refeeding(temp(1)-FP.bl:temp(1),1);
       temp = [];

       FP.refeeding.afterboutends.B9422_15safterend_avg = mean(FP.refeeding.afterboutends.B9422_refeeding_15safterend, 2)
        
end
%% Extracting traces from 15s before the start of a bout to 15s after the start of a bout:
%This section will be dedicated to extracting data from 15s before a bout
%starts and 15 seconds after it starts to compute an average trace. 

FP.bl = 600; %number of datapoints to look at before keydown
FP.seg_dur = 300; %number of datapoints to look at after keydown
t = [];
temp = [];
FP.temp_boutstart = downsample(B9422_refeeding_keydown(:),2); %temporary variable 
FP.temp_boutend = downsample(B9422_refeeding_keydown(2:end),2);
%FP.refeeding.avg_trace_startbout.B9422_refeeding_15sboutstart = []; %necessary if you re-run this part of the code changing baseline or seg_dur


for i = 1
        
        t(1,:)  = find(FP.refeeding.zscore_finalized.B9422_refeeding(:,2) >= FP.temp_boutstart(2)); %find all FP data indices that are greater than or equal to the timestamp of when a bout starts. In this case, the time vector is on the second column because we manually put the FP.time on the variable in the second column.
        temp(1,:) = t(1); %take the first one        
        
end

for ii = 1 
 
FP.refeeding.avg_trace_startbout.temp(1,:) = FP.refeeding.zscore_finalized.B9422_refeeding(temp(ii)-FP.bl:temp(ii)+FP.seg_dur,1) %Our data is on the first column because we put time on the second column. 

end

%% Extracting 15s of data for AUC analyses (not downsampled)

%FP.refeeding.auc_15s_beforeboutstarts.






%% Area under the curve analysis- Pre food and Post food access

%Timepoints
timepoints_10min = [0:10];
timepoints_5min = [0:5];
timepoints_10min_idx_prefood = [1;1201;2401;3601;4801;6001;7201;8401;9601;10801;11401];
timepoints_10min_idx_postfood = [48601;49201;50401;51601;52801;54001;55201;56401;57601;58801;60001];
timepoints_5min_idx_prefood = [1;1201;2401;3601;4801];
timepoints_5min_idx_postfood = [54001;55201;56401;57601;58801;60001];

for i = 1:length(timepoints_5min_idx_prefood) %Number of seconds being accounted for as shown above.
        
        FP.refeeding.auc.B9422_5min_prefood(i) = FP.refeeding.zscore_finalized.B9422_refeeding((timepoints_5min_idx_prefood(i)),1);
        
end

for i = 1:length(timepoints_5min_idx_postfood) %Number of seconds being accounted for as shown above.
        
        FP.refeeding.auc.B9422_5min_postfood(i) = FP.refeeding.zscore_finalized.B9422_refeeding((timepoints_5min_idx_postfood(i)),1);
        
end

FP.refeeding.auc.B9422_5min_prefood
FP.refeeding.auc.B9422_5min_postfood






%% AUC 5 minutes Post last bout 
%FP.refeeding.B9422_boutend(end)

FP.refeeding.auc.post_last_bout.B9422 = [];
temp = find(FP.refeeding.zscore_finalized.B9422_refeeding(:,2) >= FP.refeeding.B9422_boutend(end));
t = temp(1);

for i = 1:5

FP.refeeding.auc.post_last_bout.B9422(i) = FP.refeeding.zscore_finalized.B9422_refeeding(t,1);

t = t + 1201;

end 

FP.refeeding.auc.post_last_bout.B9422 = FP.refeeding.auc.post_last_bout.B9422'

%% Extracting 10-minutes of data 

FP.refeeding.freq_analysis.B5407_postfood = [];
temp = find(FP.refeeding.zscore_finalized.B5407_refeeding(:,2) >= FP.refeeding.B5407_boutend(end));
t = temp(1);



FP.freq_analysis.B5407_postfood(:,1) = FP.refeeding.zscore_finalized.B9422_refeeding(t:t+6000,1)
FP.freq_analysis.B5407_postfood(:,2) = FP.refeeding.zscore_finalized.B9422_refeeding(t:t+6000,2)





 











