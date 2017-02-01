function alpha = dfa(datetimeArray,dataArray,order,timeScaleRange,varargin)
%DFA Detrended fluctuation analysis (DFA) is a scaling analysis method used
% to estimate long-range power-law correlation exponents in noise signals.
%   
%   Originally introduced by:
%   Peng, C.K. et al. (1994). "Mosaic organization of DNA nucleotides". 
%   Phys. Rev. E 49: 1685–1689. doi:10.1103/physreve.49.1685.
%   
%   Method adapted from:
%   Hu, K. et al (2001). "Effect of trends on detrended fluctuation 
%   analysis". Phys. Rev. E 64 (1): 011114. doi:10.1103/physreve.64.011114.
%   
%   Input:
%       - datetimeArray is a vector of datetimes
%       - dataArray is the vector of the time series signal that is to be 
%         analyzed. It is assumed to have a constant sampling rate.
%       - order is the order ploynomial to fit for the local trend. Default
%         is 1.
%       - timeScaleRange is a duration pair [min_duration,max_duration]
%         specifying the range of time scales to evaluate. Default 1.5
%         hours to 8 hours.
%       - Optional: figure handle for loglog plot of F(n)
%   
%   Interpretting alpha
%       alpha < 0.5	: anti-correlated
%       alpha = 0.5	: uncorrelated, white noise
%       alpha > 0.5	: correlated
%       alpha = 1   : 1/f-noise, pink noise
%       alpha > 1   : non-stationary, unbounded
%       alpha = 1.5 : Brownian noise

% % Open the pool
% if isempty(gcp('nocreate'))
%     parpool;
% end

% Rename the variables to match Hu, et al. 2001
u = dataArray(:);
l = order;

% Integrate the time series u
N = numel(u);
y = zeros(N,1);
for j1 = 1:N
    y(j1) = sum(u(1:j1) - mean(u));
end

% Calculate time scales to test
epoch = mode(diff(datetimeArray));
n = (ceil(timeScaleRange(1)/epoch):floor(timeScaleRange(2)/epoch))';

F = zeros(numel(n),1);
for i1 = 1:numel(n)
    N_max = floor(N/n(i1))*n(i1);
    y2 = y(1:N_max);
    
    y_fit = zeros(size(y2)); % preallocate y_fit
    for i2 = 1:floor(N/n(i1))
        boxIdx = (i2-1)*n(i1)+1:i2*n(i1);
        x = boxIdx';
        % Calculate the local trend (y_fit)
        p = polyfit(x,y2(boxIdx),l);
        y_fit(boxIdx) = polyval(p,x);
    end
    % Calculate the detrended fluctuation function (Y)
    Y = y2 - y_fit;
    
    F(i1) = sqrt(sum(Y.^2)/N_max);
    
%     if i1 == floor(numel(n)/2)
%         localtrendplot(datetimeArray,dataArray,y2,y_fit,Y,N_max)
%     end
end

logF = log10(F);
logn = log10(n);
p_F = polyfit(logn,logF,1);

alpha = p_F(1);

% Create plot if a figure handle is given
if nargin >= 5
    h = varargin{1};
    figure(h);
    
    n_hrs = hours(n*epoch);
    hLog = loglog(n_hrs,F,n_hrs,10.^polyval(p_F,logn));
    ylabel('Detrended fluctuation F(n)');
    xlabel('Time scale n (hours)');
    hAxes = gca;
    xText = .75*10^mean(log10(hAxes.XLim));
    yText = 1.25*10^mean(log10(hAxes.YLim));
    textString = ['\alpha = ',num2str(alpha)];
    text(xText,yText,textString);
end

if nargin >= 6
    figureTitle = varargin{2};
    title(figureTitle);
end

end


function localtrendplot(datetimeArray,dataArray,y2,y_fit,Y,N_max)
figure
subplot(3,1,1)
plot(datetimeArray(1:N_max),dataArray(1:N_max))
subplot(3,1,2)
plot(datetimeArray(1:N_max),[y2,y_fit])
subplot(3,1,3)
plot(datetimeArray(1:N_max),Y)

end
