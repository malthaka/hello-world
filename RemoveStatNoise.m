% function [outClean] = RemoveStatNoise(fNoisy, fs)
% This function would take a noisy waveform, estimatie its
% stationary noise floor and then attempt to reduce the noise
% Input: fNoisy (the noisy waveform), 
%        param (structure holding parameters):
%           fs (sampling frequency of waveform)
%           frameSec (frame time in seconds)
%           binIncre (estimator adjustment rate)
%           gainAlpha (signal level smoother rate)
%           win (the windowing function to be applied)
% output: outClean (the cleaner waveform
%
% Copyright © 2020 The MathWorks, Inc.  
% Francis Tiong (ftiong@mathworks.com)
%
function [outClean] = RemoveStatNoise(fNoisy, param) 

%% setting parameters
fs = param.fs;
frameSec = param.frameSec;
binIncre = param.binIncre;
gainAlpha = param.gainAlpha;
hh = param.win;


%coder.varsize('fs',[8000 48000],0);

%frameSec = 0.01;       % frame time in seconds
                        % The overlap size is fixed at 50%
%binIncre = 0.05;       % estimator adjustment rate
%gainAlpha = 0.95;      % signal level smoother rate
%hh = hanning(960);     % hanning(fftsize)
%% 
ll = length(fNoisy);
frameSize = fs * frameSec;
coder.varsize('frameSize',[1 1]);
fftsize = frameSize*2;  % change size to double for 50% overlap   
if (fftsize > 192000)
    fftsize = 192000
end
% windowing
%hh = hanning(960);% hanning(fftsize)
hh = hh';

totNumFrame = single(floor(ll/frameSize));

outClean=zeros(1,ll);

binMem = ones(1,fftsize);               % noise level estimate per bin
previousFrameLvl = zeros(1,fftsize);    % frame level of the previous frame
startFrame = single(2);
overlapFrame = single(2);               % overlapFrame = 1 + numFrameOverlap

coder.varsize('totNumFrame',[1 1]);
    for ii= startFrame:totNumFrame
        
       %% obtain the amplitude of the signal frame in frequency domain       
       inn = fNoisy((ii-startFrame) *frameSize +1: (ii-startFrame) *frameSize +frameSize *overlapFrame);       
       fftin = single(inn);
       fff = fft(inn);
       anglefff = angle(fff);           % remmeber the phase angle, to be applied back after removing noise
       absfff = abs(fff);               % the amplitude will be used to track and estimate the noise floor
       meanbefore = mean(absfff);
       
       %% estimate the noise floor
       idx = absfff > binMem;           % binMem is the estimate of the noise floor
       binMem(idx) = binMem(idx) + binIncre;
       binMem(~idx) = binMem(~idx) - binIncre;
       binMem = movmean(binMem,5);      % smoothing across bin

       temp = absfff - binMem;          % remove the noise floor estimated from the signal
       idxx = temp < 0;                 
       temp(idxx) = binIncre;

       %% calculate and apply an appropriate adjustment to the signal amplitude
       temp(1:10) = absfff(1:10);
       lvlDiff = absfff - temp;         % the temporary clean signal is used to obtain an esimtated level change measure

       % LPF gain change
       meanafter = mean(temp);
        if ii==startFrame               
            previousFrameLvl = zeros(1,fftsize);  
        end
       previousFrameLvl = previousFrameLvl*(1-gainAlpha) + lvlDiff*gainAlpha; % filter the level change measure to obtain a more controlled measure
       temp = absfff - previousFrameLvl;    % obtain an estimate of a clean spectrum, temp

       %% put the phase back and obtain the time domain signal
       tempout = temp.* exp(j*(anglefff));  % put the phase back in
       tempout(fftsize/2+2:fftsize) = conj(tempout(fftsize/2:-1:2));
       ifftin = single(tempout);
       outt = real(ifft( tempout));

       %% windowing and overlap the output
       if ii==startFrame
           outClean((ii-startFrame)*frameSize+1: (ii-startFrame)*frameSize+frameSize*2) = outt(1:frameSize*2);
       else
           oldd = outClean((ii-startFrame)*frameSize+1: (ii-startFrame)*frameSize+frameSize); 
           outClean((ii-startFrame)*frameSize+1: (ii-startFrame)*frameSize+frameSize) = outt(1:frameSize).*hh(1:frameSize) + oldd.*hh(frameSize+1:end);

           outClean((ii-startFrame)*frameSize+1+frameSize: (ii-startFrame)*frameSize+frameSize*2) = outt(frameSize+1:end);
       end

    end

end