clear; close all;

load('mouse_average.mat');
facial = load('mouse_facial.mat');

clr = lines(5);

%% average mouse
[fig, ax] = initFigure();
for i = 1:30
    plot3(whiskers{i}(:,1), whiskers{i}(:,2), whiskers{i}(:,3),...
        '-', 'LineWidth', 1.5, 'Color', clr(names(i,1)-'A'+1, :));
end
plotFacial(ax, facial, 0);
finishFigure(fig, ax, ...
    'Average mouse, aligned with average row plane',...
    'Avg_RowPlane');

%% average mouse - bregma lambda
theta = cart2pol(facial.bregma(2)-facial.lambda(2), ...
                 facial.bregma(3)-facial.lambda(3));
[fig, ax] = initFigure();
for i = 1:30
    whisker = whiskers{i}*rotx(rad2deg(theta));
    plot3(whisker(:,1), whisker(:,2), whisker(:,3),...
        '-', 'LineWidth', 1.5, 'Color', clr(names(i,1)-'A'+1, :));
end
plotFacial(ax, facial, theta);
finishFigure(fig, ax, ...
    'Average mouse, aligned with bregma-lambda',...
    'Avg_BregmaLambda');

%% equation based mouse
[fig, ax] = initFigure();
basepoint = cell(1, 3);
whiskers = cell(30, 1);
for i = 1:30
    row = names(i,1)-'A'+1;
    col = names(i,3)-'0';
    eq_ThetaBP = 11.7 * col - 41.2;                                 % eq.1
    eq_PhiBP = -16.8 * row + 56.3;                                  % eq.2
    eq_RBP = -0.0152 * eq_ThetaBP - 0.0165 * eq_PhiBP + 4.44;       % eq.3b
    [basepoint{:}] = sph2cart(deg2rad(eq_ThetaBP), deg2rad(eq_PhiBP), eq_RBP);
    
    eq_S = 11.7 * exp(-0.018 * eq_ThetaBP);                         % eq.4b
    eq_A = 0.00037 * eq_ThetaBP + 0.017;                            % eq.5
    eq_ThetaW = 0.922 * eq_ThetaBP + 85.1;                          % eq.7b
    eq_PhiW = 1.00 * eq_PhiBP + 10.2;                               % eq.8
    eq_ZetaW = 0.620 * eq_ThetaBP + 1.00 * eq_PhiBP +43.0;          % eq.10
    
    [y, z] = as2xy(eq_A, eq_S, 100);
    whisker2D = [zeros(size(y)), -y, -z];
    whiskers{i} = [basepoint{:}] + ...
        whisker2D * roty(eq_ZetaW) * rotx(eq_PhiW) * rotz(eq_ThetaW)';
    
    plot3(whiskers{i}(:,1), whiskers{i}(:,2), whiskers{i}(:,3),...
        '-', 'LineWidth', 1.5, 'Color', clr(names(i,1)-'A'+1, :));
end
plotFacial(ax, facial, 0);
finishFigure(fig, ax, ...
    'Equation-based mouse, aligned with average row plane',...
    'Eq_RowPlane');

%% equation based mouse - bregma lambda
[fig, ax] = initFigure();
for i = 1:30
    whisker = whiskers{i}*rotx(rad2deg(theta));
    plot3(whisker(:,1), whisker(:,2), whisker(:,3),...
        '-', 'LineWidth', 1.5, 'Color', clr(names(i,1)-'A'+1, :));
end
plotFacial(ax, facial, theta);
finishFigure(fig, ax, ...
    'Equation-based mouse, aligned with bregma-lambda',...
    'Eq_BregmaLambda');


%% helper functions
function [fig, ax] = initFigure()
    fig = figure('Color', 'w'); hold on;
    ax = gca;
end

function finishFigure(fig, ax, title, savename)
%     view(90, 0);
    axis equal;
    plot3(ax, [0, 0], ax.YLim, [0, 0], 'k:');
    plot3(ax, ax.XLim, [0, 0], [0, 0], 'k:');
    ax.Title.String = title;
    savefig(fig, savename);
end

function plotFacial(ax, facial, angle)
    data = cell2mat(struct2cell(facial))*rotx(rad2deg(angle));
    plot3(ax, data(:,1), data(:,2), data(:,3), 'ko');
end

function R = rotx(t)
    ct = cosd(t);
    st = sind(t);
    R = [1   0    0
         0   ct  -st
         0   st   ct];
end

function R = roty(t)
    ct = cosd(t);
    st = sind(t);
    R = [ct  0   st
         0   1   0
        -st  0   ct];
end

function R = rotz(t)
    ct = cosd(t);
    st = sind(t);
    R = [ct  -st  0
         st   ct  0
         0    0   1];
end

function [x,y] = as2xy(a,s,numPts)
    xMax = fminsearch(@LOCAL_FitS, s, ...
        optimset('Algorithm','active-set','TolFun',1e-2,'display','off'),...
        a, s, numPts);
    x = (0:xMax/(numPts-1):xMax)';
    y = a.*(x.^2);
end


function e = LOCAL_FitS(x,a,s,numPts)
    xVals = 0:x/(numPts-1):x;
    yVals = a.*(xVals.^2);
    sGuess = sum(sqrt(...
        (xVals(:,2:end)-xVals(:,1:end-1)).^2 +...
        (yVals(:,2:end)-yVals(:,1:end-1)).^2));
    e = abs(s - sGuess);
end