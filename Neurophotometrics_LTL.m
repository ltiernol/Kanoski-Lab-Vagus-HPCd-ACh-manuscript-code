%% Import data:
%Import the data as a numeric matrix where time in on the first column and
%the signals are in all subsequent columns to the right.
%The following code assumes recording sessions should produce 3 csv files (maybe more depending on
%experiment: (1) rawdata; contains photometry data, (2) time; should be a
%vector of time, here we have it as ms in 'time of day', (3) keydown file
%with timestamps.

%The time and keydown files get imported as a numeric matrix.
% The photometry data get imported as a table. We import as:
% table name: AnimalID, and the columns 0G or 1G (whichever region you're
% recording from) and change its name to signal.


D2708.time = D2708_time; %The variable to the right is the raw file that gets imported, it gets imported as numeric matrix,
% and we're assigning it to a table.

%theres usually a weird artifact within the first 1-5 data points so we cut the first ten to be safe.
FP.rawdata(:,1) = D2708.time(20:end);
FP.rawdata(:,2) = D2708.signal(20:end);

%Assign the first time data point of the new table to startpoint.
startpoint = FP.rawdata(1,1);

%Here we normalize the keydown time values to start at zero and convert
%them to seconds, they come in ms originally. 
D2708_sfs_keydown = (D2708_sfs_keydown - startpoint) ./ 1000;

%Here we do the same thing but for the time column of our data. 
FP.rawdata(:,1) = (FP.rawdata(:,1) - startpoint)./1000; 

%quick check to make sure which is isosbestic
plot(downsample(FP.rawdata(2:end,2),2),'b')
hold on
plot(downsample(FP.rawdata(:,2),2),'r')
ylabel({'Raw fluorescence'});
xlabel({'data point #'});
title({'AnimalID Raw de-interleaved isosbestic and Ca dependent traces'});

%% if running behavioral task where keydown corresponds to the start of Anymaze behavioral tracking then run the code below:

FP.rawdata(:,1) = FP.rawdata(:,1) - D2708_sfs_keydown;

% Define time window
tMin = -10;
tMax = 120;

% Logical index based on time column
idx = FP.rawdata(:,1) >= tMin & FP.rawdata(:,1) <= tMax;

% Keep only the desired time range
FP.rawdata = FP.rawdata(idx, :)



%% De interleaving
%Here we de-interleave the signal because the system alternates collection
%between biomolecule-dependent signal and the isosbestic signal.

%FP.calcium_dependent(:,1) = downsample(FP.rawdata(:,1),2);
%FP.calcium_dependent(:,2) = downsample(FP.rawdata(:,2),2);
FP.calcium_dependent(:,1) = downsample(FP.rawdata(2:end,1),2); 
FP.calcium_dependent(:,2) = downsample(FP.rawdata(2:end,2),2);

%FP.isosbestic(:,1) = downsample(FP.rawdata(2:end,1),2); 
%FP.isosbestic(:,2) = downsample(FP.rawdata(2:end,2),2);
FP.isosbestic(:,1) = downsample(FP.rawdata(:,1),2);
FP.isosbestic(:,2) = downsample(FP.rawdata(:,2),2);

%Then we plot them to take another look
f3 = figure;


 subplot(2,1,1)
 plot(FP.calcium_dependent(:,1),FP.calcium_dependent(:,2),'b')
 
 subplot(2,1,2)
 plot(FP.isosbestic(:,1),FP.isosbestic(:,2),'r')

 
 f6 = figure;
 
yyaxis left
plot(FP.calcium_dependent(:,1),FP.calcium_dependent(:,2),'b')
yyaxis right 
plot(FP.isosbestic(:,1),FP.isosbestic(:,2),'r');


 
 %% Methods of neural data processing:
 
 %Method 1: Substracting a running average over a long period of time from
%your data. This works well -- for quick and dirty tests -- but is a bit of
%a hack. It also tends to fail if you have shorter recordings and higher
%SNR. ie Your recording in fiber 2. Let's use that as an example for the
%different methods.


f4 = figure;
subplot(5,1,1)
plot(FP.calcium_dependent(:,1),FP.calcium_dependent(:,2));
title('Deinterleaved Raw Data')
xlabel('Time of day in total ms')


subplot(5,1,2)
plot(FP.calcium_dependent(:,1),FP.calcium_dependent(:,2)-smooth(FP.calcium_dependent(:,2),5455));
title('Linearize by Smoothing')
xlabel('Time of day in total ms')

ylabel('F')


%Another method is to fit the data with a biexponential and subtract that
%fit from the data. This is nicer, since it is based on biology (bleaching
%is monoexponential + bleaching of fiber / heat mediated LED decay gives
%you the second exponential.

%The fiber bleaching + heat mediated decay can be minimized experimentally
%by "pre bleaching" the fiber and pre-heating the LED (you can do this all
%at once -- just set the light power to 100% for ~10 minutes before your
%experient while you are getting everything else set up.)

subplot(5,1,3)
temp_fit = fit(FP.calcium_dependent(:,1),FP.calcium_dependent(:,2),'exp2');
plot(FP.calcium_dependent(:,1),FP.calcium_dependent(:,2)-temp_fit(FP.calcium_dependent(:,1)));
title('Linearize by Fitting with Biexponential')
xlabel('Time of day in total ms')
ylabel('F')


% You can also fit, scale, and subtract the isosbestic signal from the
% calcium dependent signal. This is probably the "best" of the three
% methods desribed here -- as it is less affected by having good signal (ie
% the former 2 methods get wonky when your SNR is super high and your
% recording is short -- as it starts to fit the signal, rather than the
% slow decay, which is bad news bears)

subplot(5,1,4)
%fit isosbestic
temp_fit = fit(FP.calcium_dependent(:,1),FP.isosbestic(:,2),'exp2'); %note, in this case, I am ignoring the first 2000 points where there is this weird fast decay to get a better fit. experimentally, i normally set things up so this isn't an important time in the recording / animal is habituating.
%scale fit to calcium dependent data
fit2 = fitlm(temp_fit(FP.calcium_dependent(:,1)),FP.calcium_dependent(:,2));
%calculate a crude dF/F by subtracting and then dividing by the fit
figure(2)
subplot(3,1,3)
plot(FP.calcium_dependent(:,1),100*(FP.calcium_dependent(:,2)-(fit2.Fitted))./(fit2.Fitted))
xlabel('Time (seconds)')
ylabel('dF/F (%)')
title('A5521 CA1 - processed and traces subtracted')


FP.corrected_data = FP.calcium_dependent(:,1),100*(FP.calcium_dependent(:,2)-(fit2.Fitted))./(fit2.Fitted);
%NOTE: To get a normalized measurement (ie dF/F) you will need to know what
%the background is. For example, if you subtracted this fit from the raw
%data and then divided by the fit -- it would look super wonky -- as you
%are subtracting numbers above and below 0.

%One way to try and get at this, is record for a second without any
%excitation lights on. The signal will never be 0 -- but this will be a
%low-end of what your background signal is. True "background" would be your
%recording with the excitation lights on when no neurons are active (which
%is hard to get). The former works great.

%If your data are clean, the last method is probably the best (imho) -- but
%there is an issue. We now have dF/F values that are less than 0 which
%can't happen. Stupid fix is to add the absolute value of the lowest number
%to the entire vector. Best method is to have some estimate of background
%subtraction.

%For argument's sake, let's say the background for the isosbestic signal
%was 5000.
subplot(5,1,5)
FP.fakebackground =0;% find(min(D2708_bl.signal));
%fit isosbestic
temp_fit = fit(FP.calcium_dependent(:,1),FP.isosbestic(:,2),'exp2'); %note, in this case, I am ignoring the first 2000 points where there is this weird fast decay to get a better fit. experimentally, i normally set things up so this isn't an important time in the recording / animal is habituating.
%scale fit to calcium dependent data
fit2 = fitlm(temp_fit(FP.calcium_dependent(:,1)),FP.calcium_dependent(:,2));
%calculate a crude dF/F by subtracting and then dividing by the fit

plot(FP.calcium_dependent(:,1),100*(FP.calcium_dependent(:,2)-((fit2.Fitted)-FP.fakebackground))./(fit2.Fitted-FP.fakebackground))
xlabel('Time of day in total ms')
ylabel('dF/F (%)')
title('Linearizing + Normalizing Using Isosbestic + Background Subtracting')


%Now, if the ultimate goal is to show short sections of data chopped up
%around an event, then it doesn't make a whole heap of difference which
%method you use. In most cases, even uncorrected data, over 30 second
%segments, looks fine at a shorter time scale (since bleaching occurs over
%many minutes). However, it is good practice + necessary if you want to say
%something like "signal later on in the trial was of a different magnitude
%than signal earlier in the trial."

FP.corrected_data(:,1) = FP.calcium_dependent(:,1); %we'll keep the time vector in the first column to make things easier



temp_fit = fit(FP.calcium_dependent(:,1),FP.isosbestic(:,2),'exp2'); %note, in this case, I am ignoring the first 2000 points where there is this weird fast decay to get a better fit. experimentally, i normally set things up so this isn't an important time in the recording / animal is habituating.
%scale fit to calcium dependent data
    fit2 = fitlm(temp_fit(FP.calcium_dependent(:,1)),FP.calcium_dependent(:,2));
%calculate a crude dF/F by subtracting and then dividing by the fit

    FP.corrected_data(:,2) = 100*(FP.calcium_dependent(:,2)-((fit2.Fitted)-FP.fakebackground))./(fit2.Fitted-FP.fakebackground);
clear temp_fit fit2


f5 = figure;
    subplot(2,1,1)
   
    plot(FP.corrected_data(:,1),FP.corrected_data(:,2))
     title('Corrected Data')
     xlabel('time (seconds)')
     ylabel('dF/F')
     
FP.corrected_data(:,1) = FP.corrected_data(:,1);
FP.corrected_data(:,2) = FP.corrected_data(:,2);


%note, there is still some wonky stuff going on (ie not super duper linear.
%this is partially due to using a fake background (not in subplot 4 of
%previous graph, you don't really see this) -- and partially some artifact
%of the fit function. it tends to work better if you fit with corrected
%time data (ie not time of day in total ms -- but relative time in seconds)
%-- not sure why that is, but didn't want to muddy the waters here




%%
%Referencing behaviors
%1) Indexing: Find datapoint in FP data that is closest to behavioral
%timestamp
%2) Chop up data about these points
%3) Make dope graphs
%4) Profit

%For this, since the FP and keystroke data are in the same "time space" (ie
%timestamped by teh same computer in the same format, we are going to loop
%through each keystroke timestamp and go from there

FP.bl = 12000; %number of datapoints to look at before keydown
FP.seg_dur = 36000; %number of datapoints to look at after keydown
t = [];
temp = [];
FP.saline.D2708 = []; %necessary if you re-run this part of the code changing baseline or seg_dur
for i = 1
        t  = find(FP.corrected_data(:,1) >= D2708_keydown(i)); %find all FP data indices that are greater than or equal to keystroke timestampe
        temp(i,:) = t(1); %take the first one  
        
        %FP.seg_dur = length(FP.corrected_data)-temp;
      
        
        
        %note, these data were already corrected -- so, theoretically, the
        %differences in baseline activity are real. if you'd like to baseline
        %correct these data, you can subtract the mean of the bs period

        %FP.ERF.data(:,i) = FP.ERF.data(:,i) - mean(FP.ERF.data(1:FP.bl,i));

end
new_start = [];
for ii = 1 % 8 is the number of csp and csn trials in a session. 
FP.saline.D2708(1,:) = FP.corrected_data(temp(ii)-FP.bl:temp(ii)+FP.seg_dur,2) 
%new_start(ii) = FP.saline_Day3.data_inj_4014_CA3(ii,1);
%FP.saline_Day3.a4014_CA3_salinex3_final(ii,:) = FP.saline_Day3.data_inj_4014_CA3(ii,:) - new_start(ii);
end

FP.saline.D2708 = FP.saline.D2708';
