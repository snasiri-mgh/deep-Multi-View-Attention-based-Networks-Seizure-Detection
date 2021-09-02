function [denoised] = ekg_artifact(ekg_filtered,input_signal,sfreq)
    P = 15;
    Xc = centre(input_signal) ;
    [AfasticaD] = fasticaamar(Xc,'verbose','on','numOfIC',P, 'lastEig',P);
    SfasticaD = pinv(AfasticaD)*input_signal ;
    cor_ekg = ((SfasticaD ./ (max(abs(SfasticaD')))') * ekg_filtered')/size(SfasticaD,2) ;
    a=find(abs(cor_ekg/sum(abs(cor_ekg)))>0.2);
    
    SfasticaD_ekg=SfasticaD(a,:);
    
    SfasticaD_ekg=SfasticaD_ekg/max(abs(SfasticaD_ekg));
    
    x=find(abs(SfasticaD_ekg)>0.7);
    if length(a)>0
        dist = [];
    for i = 1 : length(x)-1
        dist(i)=x(i+1) - x(i);
    end
    if ~isempty(dist)
    f1 = find(dist > 2*sfreq);
    f2 = find(2 < dist & dist < sfreq/10); 
    
    
    if length(f1) > 0 
        a=[];
    end
    
    if  length(f2) > 0
        a=[];
    end
    end
    end
    
    AfasticaD(:,a)=0;
    denoised = AfasticaD*SfasticaD ;
end

