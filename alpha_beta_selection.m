function [alpha, beta] = alpha_beta_selection(fs,N);


fr = fs/N;

alpha_low = 0.5;

alpha_high = 0.8

alpha = zeros(1,N/2+1);

for fk = 1:1:(N/2) + 1
    
    f = fk*fr;
    
    if f <= 2000;
    
        alpha(fk) =  alpha_low;
    else
        alpha(fk) = (f-2000)*(alpha_high-alpha_low)/((fs/2)-2000) + alpha_low ;
    
    end
end

f = (1:1:257)*fr/1000;




beta_high = 0.2;

beta_low = 1;

A = 16.54

beta = zeros(1,257);

for fk = 1:1:(N/2) + 1
    
    f = fk*fr;
    
    
    beta(fk) = log10((f/A)+1)/log10((fs/2*A)+1)*(beta_high-beta_low) +  beta_low;
    
    
end

alpha_p = alpha'; beta_p = beta';

q=1;
for z = 1:N
    if z<=(N/2) + 1
       alpha(z) =  alpha_p(z);
        
    else
       alpha(z) = alpha_p((N/2) + 1-q);
       q = q + 1;  
    end  
    
end

q=1;
for z = 1:N
    if z<=(N/2) + 1
       beta(z) =  beta_p(z);
        
    else
       beta(z) = beta_p((N/2) + 1-q);
       q = q + 1;  
    end  
    
end

alpha = alpha'; beta = beta';


end