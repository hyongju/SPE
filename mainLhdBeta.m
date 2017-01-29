% Spatial distribution estimation (SPE) by mobile sensor networks (MSNs)
% using Bayesian Approach (MMLE + FIR filter)
% by Hyongju Park
clear all;close all;clc

%% 
nRuns = 2;  % total number of runs
n = 10;         % number of robots
orderK = n;    % order k, available options, 1, 2, n
fixed = 0;      % robots: will move (1), will not move (0) 
d = 2;          % dimension of the space  
nGrid = 20;  % number of grids
errTol = 0.003;
% generate halton sequence 
hSet1 = haltonset(d,'Skip',1e3,'Leap',1e2);
hScrambled1 = scramble(hSet1,'RR2');
%pos = net(hScrambled1,n);

% generate some initial configuration of n robots
pos = 1/4 *(net(hScrambled1,n) -0.25 * ones(n,d)) + 0.25 * ones(n,d);

% random seed
rng('shuffle')
varPos = 0.2^2;       % detection function (e.g., 1: large, 0.1: small, relat
                    % ive to the workspace size)
varInfo = 0.5;     % how noisy the sensor is (e.g., 0: perfect)
% boundary of the workspace
squareQ = [0 0;1 0;1 1;0 1];   % square worksapce
bndIdx = convhull(squareQ);   
bndPnts = squareQ(bndIdx,:);
fault = [1];                 % faulty robots' indice 
nsampPos = 1000;                  % number of samples (for target locations)
% number of samples for target information state
nsampInfo = 100;
sampPosTmp = net(hScrambled1,nsampPos);   % generate quasi-random samples
% consider points inside the boundary
sampPos = sampPosTmp(inhull(sampPosTmp,bndPnts),:);  

%% generate ground truth distribution
% mixture of gausians
posPeaks = [0.2785 0.9649;
        0.5469 0.1576;
        0.9575 0.9706];
mix1 = mvnpdf(sampPos,posPeaks(1,:),0.05*eye(2));
mix2 = mvnpdf(sampPos,posPeaks(2,:),0.05*eye(2));
mix3 = mvnpdf(sampPos,posPeaks(3,:),0.05*eye(2));
gTruthTmp = mix1 + mix2 + mix3;
gTruth = gTruthTmp / sum(gTruthTmp);

fig0 = figure('position',[100 100 600 600],'Color',[1 1 1]);
plot3(sampPos(:,1),sampPos(:,2),gTruth(:),'Marker','.','MarkerSize',1,'LineStyle','none');
hold on;
bdp = convhull(bndPnts);
% plot(posPeaks(:,1),posPeaks(:,2),'x');hold on;
plot(bndPnts(bdp,1),bndPnts(bdp,2),'b-');
hold on;
plot(pos(:,1),pos(:,2),'o');axis([0 1 0 1]);
set(gca,'xtick',[0 1]);
set(gca,'ytick',[0 1]);
xlabel('X');ylabel('Y');zlabel('Target distribution');

normExpWgt = gridDisplay([sampPos(:,1) sampPos(:,2) gTruth(:)],nGrid,max(gTruth));
fig1 = figure('position',[100 100 600 600],'Color',[1 1 1]);
imshow(normExpWgt','InitialMagnification','fit')
hold on;
plot(posPeaks(:,1)*nGrid,(posPeaks(:,2))*nGrid,'s','MarkerFaceColor','w','MarkerEdgeColor','r','MarkerSize',10);hold on;

hSet2 = haltonset(1,'Skip',1e3,'Leap',1e2);
hScrambled2 = scramble(hSet2,'RR2');
sampInfo = net(hScrambled2,nsampInfo)*max(gTruth);

% initialize weights...
prvWgt = 1/nsampInfo  * ones(nsampInfo ,size(sampPos,1));

% save samples...
particleWgt = ones(1,size(sampPos,1))*mean(gTruth);

normExpWgt = gridDisplay([sampPos particleWgt'],nGrid,max(gTruth));
fig3 = figure('position',[100 100 600 600],'Color',[1 1 1]);
imshow(1-normExpWgt','InitialMagnification','fit');hold on;
plot(pos(:,1)*nGrid,pos(:,2)*nGrid,'p','MarkerFaceColor','w','MarkerEdgeColor','c','MarkerSize',10);hold on;

savData{1,1} = pos;
savData{1,3} = sampPos;
savData{1,2} = particleWgt;
cstPrv = 0;
diff = 1/eps;
for curStep = 2:nRuns
    curStep
    if fixed == 1
        cst = 0;
    else
        if norm(diff) > errTol
            switch(orderK)
                case 1
                    [cst,optPos] = ctrlOrder1(pos, sampPos,particleWgt',varPos);
                case 2
                    [cst,optPos] = ctrlOrder2(pos, sampPos,particleWgt',varPos);
                case size(pos,1)
                    [cst,optPos] = ctrlOrdern(pos, sampPos,particleWgt',varPos);
                otherwise
                    error('cannot do that');
            end
            pos = optPos;
            diff = cstPrv - cst(size(cst,2));
            cstPrv = cst(size(cst,2));
        end
    end
    for i = 1:size(sampPos,1)
        mvnMax = 1*mvnpdf([0 0],[0 0],eye(2)*varPos);
        idxl = findMinDist(sampPos(i,:),pos,size(pos,1));
        switch(orderK)
            case 1
                if ismember(idxl(1), fault)
                    detectLhd(i)=0;
                else
                    detectLhd(i) = mvnpdf(sampPos(i,:),pos(idxl(1),:),eye(2)*varPos)/mvnpdf([0 0],[0 0],eye(2)*varPos);
                end
            case 2
                if ismember(idxl(1), fault)
                    detectLhd(i)= 1 - (1-mvnpdf(sampPos(i,:),pos(idxl(2),:),eye(2)*varPos)/mvnpdf([0 0],[0 0],eye(2)*varPos));
                elseif ismember(idxl(2), fault)
                    detectLhd(i)= 1 - (1-mvnpdf(sampPos(i,:),pos(idxl(1),:),eye(2)*varPos)/mvnpdf([0 0],[0 0],eye(2)*varPos));
                else                
                    detectLhd(i)=1-(1-mvnpdf(sampPos(i,:),pos(idxl(2),:),eye(2)*varPos)/mvnpdf([0 0],[0 0],eye(2)*varPos))*...
                        (1-mvnpdf(sampPos(i,:),pos(idxl(1),:),eye(2)*varPos)/mvnpdf([0 0],[0 0],eye(2)*varPos));
                end
            case size(pos,1)
                prdTerm = 1;
                for j = 1:size(pos,1)
                    if ismember(idxl(j),fault)
                        prdTerm = prdTerm * 1;
                    else
                        prdTerm = prdTerm * ...
                            (1-mvnpdf(sampPos(i,:),pos(idxl(j),:),eye(2)*varPos)/mvnpdf([0 0],[0 0],eye(2)*varPos));
                    end
                end
                detectLhd(i)=1-prdTerm;
        end
    end
    % clear previous weights (robots' previous location has no meaning...)
    wgt1 = detectLhd/sum(detectLhd);

    % findMinDist
    for i = 1:size(sampPos,1)
        for j = 1:length(sampInfo)
              infoLhd{i}(j)=normpdf(gTruth(i,:),sampInfo(j),max(sampInfo)*varInfo);
        end
    end

    for i = 1:size(sampPos,1)
        wgt2(:,i) = sirFilter(infoLhd{i}',prvWgt(:,i));
        prvWgt(:,i) = wgt2(:,i);
    end

    for i = 1:size(sampPos,1)
        wgt2Avg(i) = sampInfo' * wgt2(:,i) * wgt1(i);
    end

    particleWgt = wgt2Avg/sum(wgt2Avg);

    normExpWgt = gridDisplay([sampPos particleWgt'],nGrid,max(particleWgt));
    fig4 = figure('position',[100 100 600 600],'Color',[1 1 1]);
    imshow(normExpWgt','InitialMagnification','fit');hold on;
    plot(pos(:,1)*nGrid,pos(:,2)*nGrid,'p','MarkerFaceColor','w','MarkerEdgeColor','c','MarkerSize',10);hold on;

    savData{curStep,1} = pos;
    savData{curStep,3} = sampPos;
    savData{curStep,2} = particleWgt;
    savData{curStep,4} = cst;
end

% K-L Divergence 
for i = 1:curStep
    KLDiv(i) = 0;
    for j = 1:size(particleWgt,2)
        if gTruth(j) ~=0 && savData{i,2}(j) ~=0
        KLDiv(i) = KLDiv(i) + savData{i,2}(j) * log(savData{i,2}(j) / gTruth(j));
        end
    end
end
figure,
plot(1:curStep,KLDiv(1,1:curStep))
xlabel('time step')
ylabel('KL divergence')

% plot positions changes
for j = 1:size(pos,1)
    for i = 1:curStep
        posPlt{j}(i,:) = savData{i,1}(j,:);
    end
end
figure,
for j = 1:size(pos,1)
    plot(posPlt{j}(:,1),posPlt{j}(:,2),'o-');
    hold on;
end
axis([0 1 0 1]);
axis('square');

figure,
for j = 1:size(pos,1)
    plot(1:curStep,posPlt{j}(:,1),'-');
    hold on;
end
xlabel('time step');
ylabel('X');
figure,
for j = 1:size(pos,1)
    plot(1:curStep,posPlt{j}(:,2),'-');
    hold on;
end
xlabel('time step');
ylabel('Y');


