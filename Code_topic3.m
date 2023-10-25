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


%% DFA

% Let's consider 1 lead from 1 patient who is resting
current_EEG(1,:) = EEG_rest_zeromean(:,1,1);

segm_samples=ceil(linspace(200,900,15));
Fn=zeros(length(segm_samples),1);
for x=1:length(segm_samples)
    
    N_EEG=length(current_EEG);
    n_segm=floor(N_EEG/segm_samples(x));
    N_new=n_segm*segm_samples(x);
    fit_order=1;
    % 1. calcolo la media del segmento
    mean_segm=mean(current_EEG(1:N_new));

    % 2. tolgo la media da ogni singolo punto e "integreted"
    y=zeros(N_new,1);
    for i=1:N_new
        y(i)=sum(current_EEG(1:i)-mean_segm);
    end

    % 3. considero un segmento per volta e calcolo i coeff del polinomio
    fitcoeff=zeros(n_segm,fit_order+1);
    fitcurve=zeros(N_new,1);
    for j=1:n_segm
    fitcoeff(j,:)=polyfit(1:segm_samples(x),y((j-1)*segm_samples(x)+1:j*segm_samples(x)),fit_order);
     % 4. calcolo la retta per ogni segmento
    fitcurve((j-1)*segm_samples(x)+1:j*segm_samples(x))=polyval(fitcoeff(j,:),1:segm_samples(x));
    end

    
    % 5. calcolo msd e ne prendo la radice quadrata
    msd= sum((y-fitcurve).^2)/N_new;
    Fn(x)=sqrt(msd);
end


    n=log10(segm_samples);
    Fn_log=log10(Fn);

figure
plot(n,Fn_log,'o')


