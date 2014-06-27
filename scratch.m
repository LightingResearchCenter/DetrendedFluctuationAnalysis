% Reset
close all
clear
clc

% Enable CDF Toolkit
[parentDir,~,~] = fileparts(pwd);
cdfToolkitDir = fullfile(parentDir,'LRC-CDFtoolkit');
addpath(cdfToolkitDir);

% Import test data
testFilePath = 'test.cdf';
TestData = ProcessCDF(testFilePath);

% Set valid test data range
TestData.Variables.logicalArray = ...
    TestData.Variables.time >= datenum(2014,06,09,06,37,30) & ...
    TestData.Variables.time <= datenum(2014,06,14,11,03,00);

% Reassign variables
logicalArray     = TestData.Variables.logicalArray;
datenumArray     = TestData.Variables.time(logicalArray);
datetimeArray    = datetime(datenumArray,'ConvertFrom','datenum');
illuminanceArray = TestData.Variables.illuminance(logicalArray);
claArray         = TestData.Variables.CLA(logicalArray);
csArray          = TestData.Variables.CS(logicalArray);
activityArray    = TestData.Variables.activity(logicalArray);

timeScaleRange   = [duration(1,30,0),duration(8,0,0)];

h = figure;

alpha = dfa(datetimeArray,activityArray,1,timeScaleRange,h);

display(alpha);