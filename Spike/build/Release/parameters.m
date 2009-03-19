DT=0.000010;
TotalMS=4800;
rfpath=Results//Network/Servers/mac0.cns.ox.ac.uk/Volumes/Data/Users/evans/Code/Spike/diagnostic.prm/;

% Diagnostic simulation
loops = 2;
nRecordsPL = 24;

nExcit = 24;
	nSynEfE = 24;
	nSynIE = 6;
nInhib = 12;
	nSynEI = 24;
	nSynII = 6;

transP_Train = 0.4;
transP_Test = 0.4;

nStimuli = 3;
nTransPS = 4;
shift = 2;

current 	= 1.25e-9;
gMax = 48.0e-9;