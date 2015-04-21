function x = cmMin(fdot, P = [], verbose=0)
    if size(P, 2) < 3
        p = [0.75; 0.75]
        %p
        while any((g = fdot(p)) < 0) && all(finite(p)) 
            g
            p = p .* ((g < 0) + 1)
        end
        
        g = fdot(p);
        P = [0         p'*g/g(2) 0; ...
             p'*g/g(1) 0         0]
    end
    if verbose
        P
    end
    
    maxit = 2000;
    if all(fdot([0; 0]) >= 0)
        x = [0; 0];
        return;
    end
    for it=1:maxit
        numVertices = size(P, 2);
        cm = [0; 0];
        sTotal = 0.0;
        for j=3:numVertices % Center of mass calculation with trivial triangulation
            r1 = P(:, 1);
            r2 = P(:, j - 1);
            r3 = P(:, j);
            cmLocal = (r1 + r2 + r3) / 3;
            r21 = r2 - r1;
            r31 = r3 - r1;
            sLocal = 0.5 * abs(r21(1) * r31(2) - r21(2) * r31(1));

            cm = cm + (cmLocal * sLocal);
            sTotal = sTotal + sLocal;
        end
        if sTotal ~= 0
            cm = cm / sTotal;
        else
            cm = P(:, 1);
        end
        g = fdot(cm);
        d = 0; % diameter of polygon
        for j = 1:numVertices
            for k = j+1:numVertices
                d = max(norm(P(:,k) - P(:,j)), d);
            end
        end
        if verbose
            P
            sTotal
            cm
            g
            d
        end
        if d < 1e-8 || norm(g) < 1e-8
            break;
        end
        newP = [];
        for j = 1:numVertices
            a = P(:, j);
            b = P(:, mod(j, numVertices) + 1);
            dirVec = b - a;
            tnom = (cm - a)' * g;
            tdenom = dirVec' * g;

            if (tnom > 0)
                newP = [newP a];
            end

            if tdenom ~= 0
                t = tnom / tdenom;
                if t >= 0 && t <= 1
                    x = dirVec * t + a;
                    newP = [newP x];                
                end
            end
        end
        P = newP;
    end
    if verbose
        it
    end
    x = cm;
end