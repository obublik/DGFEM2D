function out = myqs(data)
if(isempty(data) ~= 1)
    n = length(data(:,1));
    if(n < 2)
        out = data;
    else
        ind = floor(n/2);
        p = data(ind,:);
        L = zeros(size(data));
        iL = 1;
        R = zeros(size(data));
        iR = 1;
        if(n > 1)
            for i = 1:n
                if(i ~= ind)
                    if(je_mensi(data(i,:),p))
                        L(iL,:) = data(i,:);
                        iL = iL + 1;
                    else
                        R(iR,:) = data(i,:);
                        iR = iR + 1;
                    end
                end
            end
            L = myqs(L(1:iL-1,:));
            R = myqs(R(1:iR-1,:));
            out = [L; p; R];
        end
    end
else
    out = [];
end

function c = je_mensi(a,b)
c = 0;
if(b(1) - a(1) > 1e-5)
    c = 1;
elseif(abs(a(1) -b(1)) < 1e-5)
    if(b(2) - a(2) > 1e-5)
        c = 1;
    elseif(abs(a(2) - b(2)) < 1e-5)
        if(b(3) - a(3) > 1e-5)
            c = 1;
        end
    end
end