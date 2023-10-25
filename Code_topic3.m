clear all
close all
clc
%% load Data Subjects

EEG_rest(:,:,1)=struct2array(load("Data/Subject01_1.mat"));
EEG_rest(:,:,2)=struct2array(load("Data/Subject02_1.mat"));
EEG_rest(:,:,3)=struct2array(load("Data/Subject03_1.mat"));
EEG_rest(:,:,4)=struct2array(load("Data/Subject05_1.mat"));
EEG_rest(:,:,5)=struct2array(load("Data/Subject07_1.mat"));
EEG_rest(:,:,6)=struct2array(load("Data/Subject08_1.mat"));

EEG_calc(:,:,1)=struct2array(load("Data/Subject01_2.mat"));
EEG_calc(:,:,2)=struct2array(load("Data/Subject02_2.mat"));
EEG_calc(:,:,3)=struct2array(load("Data/Subject03_2.mat"));
EEG_calc(:,:,4)=struct2array(load("Data/Subject05_2.mat"));
EEG_calc(:,:,5)=struct2array(load("Data/Subject07_2.mat"));
EEG_calc(:,:,6)=struct2array(load("Data/Subject08_2.mat"));

load helper_code/chanlocs.mat

Fs=500;


%t_rest=[1:length(subj1_rest.C3)]/Fs;
%t_calc=[1:length(subj1_calc.C3)]/Fs;

%Zero-mean signals
for p=1:size(EEG_rest,3)
    for j=1: size(EEG_rest,2)
        EEG_rest_zeromean(:,j,p)= EEG_rest(:,j,p)-mean(EEG_rest(:,j,p));
    end
end

for p=1:size(EEG_calc,3)
    for j=1: size(EEG_calc,2)
        EEG_calc_zeromean(:,j,p)= EEG_calc(:,j,p)-mean(EEG_calc(:,j,p));
    end
end

%plot(EEG_rest_zeromean(:,1,1))
%mean(EEG_rest_zeromean(:,1,2))


%% Power spectral density maps

PSD_rest=[];
f_rest=[];
[rows,columns,pages]=size(EEG_rest_zeromean);
for p=1:pages %for the pages
    for j=1:columns
            current_value=EEG_rest_zeromean(:,j,p);
            k = 39;
            L = ceil((2*length(current_value))/(k+1));
            noverlap=L/2;
            [PSD_rest(:,j,p),f_rest(:,j,p)]=pwelch(current_value,hamming(L),noverlap,[],Fs);
    end
end


PSD_calc=[];
f_calc=[];
[rows,columns,pages]=size(EEG_calc_zeromean);
for p=1:pages %for the pages
    for j=1:columns
            current_value=EEG_calc_zeromean(:,j,p);
            k = 20;
            L = ceil((2*length(current_value))/(k+1));
            noverlap=ceil(L/2);
            [PSD_calc(:,j,p),f_calc(:,j,p)]=pwelch(current_value,hamming(L),noverlap,[],Fs);
    end
end


%% DOWNSAMPLE (rest) 

% LP filter at 60Hz (because we can see that the range is up to about 50 Hz)
for p=1:size(EEG_rest_zeromean,3)
    for j=1:size(EEG_rest_zeromean,2)
        EEG_rest_lp(:,j,p) = lowpass(EEG_rest_zeromean(:,j,p), 60, Fs);
        EEG_rest_ds(:,j,p) = downsample(EEG_rest_lp(:,j,p), 110);
    end
end
% we saw what's the last significant component in freq and chose a
% sampling freq (with the theorem) to downsample the signal keeping only
% the significant samples and shortening the processing time
% the new Fs (with respect to the original signal given) is 
% Fs_new = Fs*fs = 500*110 = 55000 Hz (actually not needed)

%% DFA (rest)
segm_samples=[4 8 16 32 64 128 256]; % numero di sample per segmento (da far ruotare nel ciclo)
fit_order = 1;

for p=1:size(EEG_rest_ds,3)
    for j=1:size(EEG_rest_ds,2)
        current_EEG_rest(1,:) = EEG_rest_ds(:,j,p);
        Fn_rest = DFA(current_EEG_rest, segm_samples, fit_order);
        
        % 6. calcolo la retta di interpolazione di Fn(n) (in scala
        % logaritmica)
        n=log10(segm_samples);
        Fn_log=log10(Fn_rest);
        fn_coeff = polyfit(n,Fn_log, 1);
        fn_fit = polyval(fn_coeff,n);

        % 7. calcolo lo slope delle rette ottenute
        maxy = max(fn_fit);
        miny = min(fn_fit);
        if miny<abs(miny)
            y_slope = maxy + abs(miny);
        else
            y_slope = maxy - miny;
        end 
        x_slope = n(end) - n(1);
        slope_rest(j,p) = y_slope/x_slope;               
    end 
end

% figure
% plot(n,Fn_log,'o')
% hold on
% plot(n,fn_fit)


%% DOWNSAMPLE (calc) 

% LP filter at 60Hz (because we can see that the range is up to about 50 Hz)
for p=1:size(EEG_calc_zeromean,3)
    for j=1:size(EEG_calc_zeromean,2)
        EEG_calc_lp(:,j,p) = lowpass(EEG_calc_zeromean(:,j,p), 60, Fs);
        EEG_calc_ds(:,j,p) = downsample(EEG_calc_lp(:,j,p), 110);
    end
end
% we saw what's the last significant component in freq and chose a
% sampling freq (with the theorem) to downsample the signal keeping only
% the significant samples and shortening the processing time
% the new Fs (with respect to the original signal given) is 
% Fs_new = Fs*fs = 500*110 = 55000 Hz (actually not needed)

%% DFA (calc)

for p=1:size(EEG_calc_ds,3)
    for j=1:size(EEG_calc_ds,2)
        current_EEG_calc(1,:) = EEG_calc_ds(:,j,p);
        Fn_calc = DFA(current_EEG_calc, segm_samples, fit_order);
        
        % 6. calcolo la retta di interpolazione di Fn(n) (in scala
        % logaritmica)
        n=log10(segm_samples);
        Fn_log=log10(Fn_calc);
        fn_coeff = polyfit(n,Fn_log, 1);
        fn_fit = polyval(fn_coeff,n);

        % 7. calcolo lo slope delle rette ottenute
        maxy = max(fn_fit);
        miny = min(fn_fit);
        if miny<abs(miny)
            y_slope = maxy + abs(miny);
        else
            y_slope = maxy - miny;
        end 
        x_slope = n(end) - n(1);
        slope_calc(j,p) = y_slope/x_slope;               
    end 
end

% figure
% plot(n,Fn_log,'o')
% hold on
% plot(n,fn_fit)

topoplot('example')
