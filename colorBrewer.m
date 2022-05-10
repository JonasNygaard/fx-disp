function brewedColors = colorBrewer(numColors,colorSet)

%% colorBrewer.m
% ########################################################################### %
% function  brewedColors = colorBrewer(numColors,colorSet)
% Purpose:  Create a vector of RGB colors based on a pre-defined set
%
% Input:    numColors   = Scalar indicating the number of colors needed
%           colorSet    = The pre-defined set to be used
%
% Output:   brewedColor = Matrix of RGB codes for colors for plots
%               
% Author:
% Jonas N. Eriksen
% Department of Economics and Business Economics
% Aarhus University and CREATES
%
% Encoding: UTF8
% Last modified: April, 2018
% ########################################################################### %

if (nargin > 3)
    error('colorBrewer.m: Too many input arguments');
end

if (nargin < 1)
    error('colorBrewer.m: Not enough input arguments');
end

if (numColors > 9)
    error('colorBrewer.m: A most five colors are currently supported');
end

if (numColors < 0)
    error('colorBrewer.m: Negative number of colors not a valid input');
end

if (nargin < 2)
    colorSet = 6;
end

%% Setting RGB values for color palette
% ########################################################################### %
%{
    We consider different set of colors. The first three are based heavily 
    on colorbrewer2.org, whereas the remaining schemes are mixes of different 
    colors that simply looks nice on print. 
%}
% ########################################################################### %

% Setting RGB values
switch colorSet

    case 1

        rgbValues = [
            0.21569     0.49412     0.72157             % Blue
            0.89412     0.10196     0.1098              % Red
            0.30196     0.68627     0.2902              % Green
            1           0.49804     0                   % Orange
            0.59608     0.30588     0.63922             % Purple
            1           1           0.2                 % Yellow
            0.65098     0.33725     0.15686             % Brown 
            0.96863     0.50588     0.74902             % Pink
        ];

    case 2

        rgbValues = [
            0.55294     0.62745     0.79608             % Blue
            0.98824     0.55294     0.38431             % Orange
            0.4         0.76078     0.64706             % Teal
            0.90588     0.54118     0.76471             % Purple
            0.65098     0.84706     0.32941             % Green
            1           0.85098     0.18431             % Yellow
            0.89804     0.76863     0.58039             % Brown 
            0.70196     0.70196     0.70196             % Gray
        ];

    case 3

        rgbValues = [
            0.50196     0.69412     0.82745             % Blue
            0.98431     0.50196     0.44706             % Red
            0.55294     0.82745     0.78039             % Teal
            0.99216     0.70588     0.38431             % Orange
            0.7451      0.72941     0.8549              % Purple
            1           1           0.70196             % Yellow
            0.70196     0.87059     0.41176             % Green
            0.98824     0.80392     0.89804             % Pink
        ];

    case 4

        rgbValues = [
            0.000   0.447   0.741
            0.850   0.325   0.098
            0.929   0.694   0.125
            0.494   0.184   0.556
            0.466   0.674   0.188
            0.301   0.745   0.933
            0.635   0.078   0.184
            0.000   0.447   0.741
            0.16    0.16    0.16      
        ];

    case 5

        rgbValues = [
            42      52      122     % Blue 
            220     76      51      % Red
            44      162     95      % Green
            166     118     29      % Brown
            179     179     179     % Gray
            139     0       139     % Purple        
            253     180     98      % Orange
            252     205     229     % Pink
            40      40      40      % Black
        ]./255;

    case 6

        rgbValues = [
            97      158     207
            189     104     133
            88      170     100
            153     126     63
            124     108     203
            173     174     61
            208     109     55
            212     70      93
            197     87      178
        ]./255;

end

% Setting output
if ismember(numColors,1:9)

    brewedColors  = rgbValues(numColors,:);

else

    brewedColors  = rgbValues;

end

end

% ########################################################################### %
% [EOF]
% ########################################################################### %