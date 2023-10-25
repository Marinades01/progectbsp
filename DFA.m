function [Fn]=DFA(current, nos, order)
% current = the current lead cosidered
% nos = vector with the different number of samples per segment
% order = order of the polynomial fit
% Fn = vector of F(n) values for one lead with the "nos" vector 

Fn=zeros(length(nos),1);
for x=1:length(nos)
    N_EEG=length(current);
    n_segm=floor(N_EEG/nos(x)); % numero di segmenti 
    N_new=n_segm*nos(x);

    % 1. calcolo la media del segmento
    mean_segm=mean(current(1:N_new)); % media di tutti i segmenti fino a qui (diversa per ogni ciclo)

    % 2. tolgo la media da ogni singolo punto e "integreted"
    y=zeros(N_new,1);
    for i=1:N_new
        y(i)=sum(current(1:i)-mean_segm);
    end

    % 3. considero un segmento per volta e calcolo i coeff del polinomio
    fitcoeff=zeros(n_segm,order+1); % un coeff per ogni esponente più quello per esponente zero
    fitcurve=zeros(N_new,1); % per i punti della curva da inserire ad ogni ciclo (ad ogni diverso num di sample)
    for k=1:n_segm
    fitcoeff(k,:)=polyfit(1:nos(x),y((k-1)*nos(x)+1:k*nos(x)),order);

    % 4. calcolo la retta per ogni segmento
    fitcurve((k-1)*nos(x)+1:k*nos(x))=polyval(fitcoeff(k,:),1:nos(x));
    end        
    
    % 5. calcolo msd e ne prendo la radice quadrata
    msd= sum((y-fitcurve).^2)/N_new;
    Fn(x)=sqrt(msd); 
end


end