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

%figure (1)
%subplot 211
%plot(f_rest(:,12,2),PSD_rest(:,12,2))
%subplot 212
%plot(f_calc(:,12,2),PSD_calc(:,12,2))

%figure (2)
%subplot 211
%plot(f_rest(:,12,3),PSD_rest(:,12,3))
%subplot 212
%plot(f_calc(:,12,3),PSD_calc(:,12,3))

%figure(3)
%subplot 211
%plot(f_rest(:,12,6),PSD_rest(:,12,6))
%subplot 212
%plot(f_calc(:,12,6),PSD_calc(:,12,6))

%% DFA

%EEG_rest
% 1.integration of the signal for all the leads
%for p=1:pages
    %for j=1:columns
        %EEG_rest_integ(:,j,p)=cumtrapz(EEG_rest_zeromean(:,j,p));
    %end 
%end

% Let's consider 1 lead from 1 patient who is resting
current_EEG(1,:) = EEG_rest_zeromean(:,1,1);

segm_samples=ceil(linspace(200,900,15));
Fn=zeros(length(segm_samples),1);
for x=1:length(segm_samples)
    x
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




% 2. dividing the signal in sub-sequences
%fitted=[];
%n_sample=ceil(linspace(200,900,100));
%for l=1:length(n_sample)
%n_sample=500;
    %[sub_sequences,n_sequences]=subsequences(current_EEG,n_sample(l));

    %msd_matrix=zeros(n_sequences,1);
    %polyn=[];
    %for i=0:n_sequences-1
        %corrente=sub_sequences(i+1,:);
        
        %coef_poly=polyfit(1:n_sample(l),corrente,1);
       % polyn=[polyn; polyval(coef_poly,1:n_sample(l))];
        
        %msd_matrix(i+1) =mean((corrente-polyn(i+1)).^2,2); % mean-square_deviation
    
   % end
    
    %fitted=[fitted sqrt(mean(msd_matrix,1))];
%end
%n=log10(n_sample);
%Fn=log10(fitted);

%figure
%plot(n,Fn,'o')
%hold on
%coef_fn=polyfit(n,Fn,1);
%polyn_fn=polyval(coef_fn,n);
%plot(polyn_fn)
%hold off

      
   % end
 
%end




