function [sse, p, sd] = pqmodel(x, m, f, pdata, qdata, nrOfHb )
%PQMODEL Numerical solving of windkessel differntial equations

    %assign parameters
    R=x(1);
    C=x(2);
    Z=x(3);
    L=x(4);
    
    % calculate parameters
    N=length(pdata);
    T=1/f;            % cardiacal period [s]
    dt = nrOfHb*T/(N-1);

    %loop over K*nrOfHb cycles
    K=2;
    % Approximate solution of DV
    p=zeros(N,1);
    r = p;
    pa = p;
    pa(N)=pdata(1);
    p(N)=pdata(1);
    pa(1) = pa(N);
    p(1)=p(N);
    for k=1:K
        pa(1)=pa(N);
        p(1)=p(N);


        if(k==K)
            for n=2:1:N
                pa(n)=pa(n-1)*(1-dt/(R*C))+qdata(n-1)*dt/C;
            end
        end

        if m==1
            p = pa;
        elseif m==2
            for n=2:1:N
                p(n) = pa(n)+Z * qdata(n);
            end
        else
            for n=2:1:N
                r(n) = (pa(n)-pa(n-1))/dt + pa(n)*Z/L + (qdata(n)-qdata(n-1))*Z/dt;
                p(n) = (1-Z*dt/L)*p(n-1) + dt*r(n);
            end
        end
    end

    sd = zeros(N,1);
    sse = 0;
    for n=1:1:N
        sdsq = (p(n)-pdata(n))^2;
        sd(n) = sqrt(sdsq);
        sse = sse + sdsq;
    end
end

