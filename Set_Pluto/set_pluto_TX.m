function   [tx_radio]=set_pluto_TX(fc,fs,sn,Tx_Gain)

sdrdev('Pluto');
% fc=2400e6;
% fs=61.44e6
% samples_per_frame = 2560;
samples_per_frame=length(sn);
%Initialize SDR Receiver Functionality
% rx_radio = sdrrx('Pluto');
% rx_radio.RadioID =              'ip:192.168.2.3';
% rx_radio.CenterFrequency =      fc;
% rx_radio.BasebandSampleRate =   fs;
% rx_radio.GainSource =           'Manual';
% rx_radio.Gain =                 Rx_Gain;
% rx_radio.OutputDataType =       'double';
% rx_radio.SamplesPerFrame =      samples_per_frame*2;


% load('tx_low_pass_filter.mat');

% Initialize SDR Transmit Functionality
tx_radio = sdrtx('Pluto');
tx_radio.RadioID =                  'ip:192.168.2.5';
tx_radio.CenterFrequency =          fc;
tx_radio.BasebandSampleRate =       fs;
tx_radio.Gain =                     Tx_Gain;
tx_radio.ShowAdvancedProperties =   true;
tx_radio.DataSourceSelect =         'Input Port';
constDiagram = comm.ConstellationDiagram('SamplesPerSymbol',samples_per_frame, ...
    'SymbolsToDisplaySource','Property','SymbolsToDisplay',10000, 'ShowTrajectory', true);
end