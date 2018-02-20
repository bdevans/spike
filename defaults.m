%# Default parameters for Spiking Neural Network Simulation

%### Simulation ###
DT          = 0.0001;   % Numerical time step
loops 		= 5;    % Number of training epochs         [int >= 1]
pretrain 	= 1;    % Test the network before training  {0,1}
train 		= 1;    % Train the network                 {0,1}
priorPhases = 0;    % Preliminary (testing &) training  {0,1}
isolateEfE  = 0;    % If PP, train ElE then EfE         {0,1}
trainPause	= 0;    % Reset neurons between stimuli     {0,1}
noise		= 0;    % Noise model {NONE, UNIFORM, GAUSSIAN}
noiseScale	= 0.015;	% Masquelier et al 2009
normalise   = 0;    % Normalise weights                 {0,1,2}
nRecordsPL 	= 0;    % Neurons to record per layer       [0 <= i <= nExcit]
printConnections = 1;	% Output connectivity files		{0,1}
probConnect	= 0;
loadWeights = 0;    % Skip training and load weights    {0,1}

%### Stimuli ###
randStimOrder 	= 1;    % Randomize stimuli during training     {0,1}
randTransOrder	= 0;    % Randomize transforms during training  {0,1}
randTransDirection  = 0;% Random direction of transforms during training {0,1}
interleaveTrans	= 0;    % Interleave transforms during training {0,1}
localRep 	= 1;        % Distributed or local stimuli          {0,1}
current 	= 1.25e-9;	% (float) 0.00000000125 # Amps
currentSpread	= 0.25e-9;
stimGroups  = 0;
%nGroups     = 1; % Internal
nStimuli 	= 6;        % Number of stimuli                     [int >= 1]
nTransPS 	= 4;        % Number of transforms per stimulus     [int >= 1]
newTestSet	= 0;
K           = 1;
M           = 2;
%nTestGroups = 1;
transP_Train 	= 1.00;	% Transform presentation period (training)  # s
transP_Test 	= 1.00;	% Transform presentation period (testing)   # s
shift 		= 2;		% Number of neurons to move each transform by
a		= 0.0;          % Sparsity of input layer -> nFiringNeurons [0,1]
useFilteredImages = 0;
    gabor	= 1;

%### Network ###
nLayers = 2;            % Number of layers      [int >= 2]
inputInhib = 1;         % Presence of input inhibitory connections  {0,1}
	nExcit = 120;       % Number of Excitatory neurons per layer    (int)
		nSynEfE = 60;   % Number of feedforward E -> E synapes per E neuron
		nSynElE = 0;    % Number of lateral E -> E synapses per E neuron
		nSynIE = 40;    % Number of (lateral) I -> E synapses per E neuron
	nInhib = 40;        % Number of Inhibitory neurons per layer (int)
		nSynEI = 60;    % Number of (lateral) E -> I synapses per I neuron
		nSynII = 5;     % Number of (lateral) I -> I synapses per I neuron
initEfE = 1;            % Initialisation of EfE weights: {Constant,Uniform,Gaussian,SOM}
iEfE = 0.0;             % Initialisation value when initEfE==Constant
%axonDelay = 0;          % Now an internal variable. Axonal delay model {0,1,2,3} {NONE, CONST, UNIFORM, GAUSSIAN}
delayEfE = 0;           % {0:NONE, 1:CONST, 2:UNIFORM, 3:GAUSSIAN}
delayElE = 0;           % {0:NONE, 1:CONST, 2:UNIFORM, 3:GAUSSIAN}
delayEI = 0;            % {0:NONE, 1:CONST, 2:UNIFORM, 3:GAUSSIAN}
% Constant delay (1)
	d_const = 0.005;	% Constant delay (float)        # s
% delay = UNIFORM; (2)
	d_min 	= 0.002;    % Minimum delay (float)         # s
	d_max 	= 0.008;    % Maximum delay (float)         # s
% delay = GAUSSIAN; (3)
	d_mean 	= 0.005;    % Mean delay (float)            # s
	d_sd 	= 0.002;    % Standard deviation (float)    # s
% delay = SOM; (4)
    spatialScale = 16.0;% Neurons
    condSpeed = 1.0;    % Neurons/s
    maxDelay = 0.025;   % Delay (s) across 1/2 longest input dimension
%SOM         = 0;        % Use SOM architecture (bool)
    SOMinput    = 1;    % Use SOM architecture in input layer (bool)
    SOMsigE     = 4;    % S.D. of SOM for excitatory neurons
    SOMsigI     = 4;    % S.D. of SOM for inhibitory neurons
    SOMclip     = 3;    % S.D.s at which to clip connections i.e. 3sigma
%trainEfE    = 1;        % Make EfE weights plastic (added for priorPhases)
trainElE    = 0;        % Make ElE weights plastic {true,false}
initElE     = 1;        % Initialisation of ElE weights: {Constant,Uniform,Gaussian,SOM}
iElE = 0.0;             % Initialisation value when initElE==Constant
vSquare     = [0,1];

%### Cell bodies ###
%# See Song, Miller & Abbot 2000
%# WB uses IF parameters from this paper, which cites mainly cat v1 literature.
%# Troyer TW, Krukowski AE, Priebe NJ, Miller KD (1998) Contrast-invariant orientation tuning in cat visual cortex: thalamocortical input tuning and correlation-based intracortical connectivity. Journal of Neuroscience 18:5908--5927.

%# Two groups that do a lot of this kind of modelling are:  Ken D Miller and McLaughlin & Shapley
%# McLaughlin D, Shapley R, Shelley M, Wielaard DJ (2000) A neuronal network model of macaque primary visual cortex (V1): orientation selectivity and dynamics in the input layer 4Calpha. Proc Natl Acad Sci U S A 97:8087--8092.

capE 	= 0.1e-9;	%# 0.0000000002 	# Farads
capI 	= 5.0e-12;	%# 0.00000000001	# Farads
gLeakE 	= 5.0e-9;	%# 0.00000001		# Siemens
gLeakI 	= 5.0e-9;	%# 0.000000005		# Siemens

VrestE 	= -0.070;	%# Volts (-70mV)
VrestI 	= -0.070;	%# Volts
VhyperE	= -0.085;	%# Volts
VhyperI = -0.085;	%# Volts
ThreshE = -0.038;	%# Volts
ThreshI = -0.045;	%# Volts
VrevE 	= -0.005;	%# Volts
VrevI 	= -0.090;	%# Volts

refract = 0.001;	%# Seconds

adaptation = 0;         %# Liu+Wang2001
    alphaCa = 200e-9;   %# 0.00000002; %# DgK/gK
	tauCa   = 0.050;    %# Seconds
	gAHP    = 15e-6;    %# 0.000015; %# S/cm^2 %# gK
	VK      = -0.080;   %# Volts

%### Synapses ###
%# Perrinet+++01
noSTDPdelay = 0;    %# {0: modified form; 1: original form}
alphaC 	= 0.5;		%# (Dimensionless)
tauC 	= 0.005;	%# s
alphaD	= 0.5;		%# (Dimensionless)
tauD 	= 0.005;	%# s
learnR 	= 0.1;		%# (Dimensionless)
gMax 	= 50.0e-9;	%# 0.0000000048
DgEfE   = 1.0;      %modEf   = 1.0;
DgElE   = 1.0;      %Dg_ElE  = 1.0;
DgIE 	= 0.1;      %Dg_IE 	= 0.1;
DgEI 	= 1.0;      %Dg_EI 	= 1.0;		
DgII 	= 0.1;      %Dg_II 	= 0.1;
tauEE 	= 0.010;	%# s
tauIE 	= 0.001;	%# s
tauEI 	= 0.001;	%# s
tauII 	= 0.001;	%# s
