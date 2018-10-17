clear;
e = 1.60217662E-19;
speedLight = 299792458;
invLam = importdata('invLam.mat');
itt = 1;
overstep = false;
for N = 380
    if overstep == true
        break;
    end
    sumsq = 0;
    for i = 1:6
        filex = sprintf('x%d.mat', i);
        filey = sprintf('y%d.mat', i);
        x = importdata(filex);
        y = importdata(filey);
        cutoff = max(y)/2;
        found = false;
        it = 1;
        while found == false
            if y(it) > cutoff 
                if it + N/2 > length(x)
                    overstep = true;
                    break;
                elseif it - N/2 < 0
                    overstep = true;
                    break;
                end
                outX = x((it-N/2):(it+N/2));
                outY = y((it-N/2):(it+N/2));
                found = true;
            end
            it = it + 1;
        end
        [FO, goodness] = fit(outX, outY, 'poly1');
        %disp(goodness.rmse)
        sumsq = sumsq + goodness.adjrsquare^2;
        VD = -FO.p2/FO.p1;
        lineY = outX * FO.p1 + FO.p2;
        c = confint(FO);
        errIntercept = c(:,2);
        errIntercept = (errIntercept(2) - errIntercept(1))/2;
        errSlope = c(:,1);
        errSlope = (errSlope(2) - errSlope(1))/2;
        errVD = sqrt((errIntercept/FO.p1)^2 + (FO.p2*errSlope/(FO.p1^2))^2);
        %disp(VD); disp(errVD); disp(errIntercept); disp(errSlope);
        %plot(outX,outY,'.', outX, lineY, '-');
        %title(N);
        %drawnow;
        VDs(i,1) = VD; VDs(i,2) = errVD;
    end
    store(itt) = sumsq/6;
    if store(itt) == max(store)
        disp(N);
    end
    plot(store); drawnow;
    V = VDs(:,1);
    W = VDs(:,2);
    [finalFit, finalgoodness] = fit(invLam,V,'poly1', 'Weight', W);
    %plot(finalFit);
    h = finalFit.p1*e/speedLight;
    sigSlope = confint(finalFit);
    sigSlope = (sigSlope(2,1) - sigSlope(1,1))/2;
    errH = sqrt((e*sigSlope/speedLight)^2);
    disp(VDs);
    disp(h*10^-9); disp(errH*10^-9);
    outH(itt,1) = h; outH(itt,2) = errH;
    itt = itt + 1;
end
%plot(linspace(10,N,size(outH,1)),outH)