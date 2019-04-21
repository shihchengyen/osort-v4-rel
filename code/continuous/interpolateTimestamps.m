%
%returns a timestamp for every data sample in dataSamples, interpolated
%based on the timestamps provided in timestamps (one stamp for every 512
%samples). Used to convert from neuralynx data structures to internal data
%structures.
%
%
%urut/aug10/MP initial version
%urut/dec11 new version: need sampling rate to interpolate
%correctly(required)
%
%
function [yi] = interpolateTimestamps(dataSamples, timestamps, Fs)
blockSize=512;

timestamps = [ timestamps timestamps(end)+512*Fs/1e6 ]; % need to add one timepoint

x = 0:blockSize:length(dataSamples);
xi = 1:length(dataSamples);
yi = interp1q(x',timestamps', xi' );