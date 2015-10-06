% Battery of cell analyses to confirm new model details work properly
loops       = 2;
pretrain    = 1;
train       = 1;

probConnect = 1;
printConnections    = 1;

current     = 0.25e-9;
currentSpread   = 0.25e-9;

noise 	    = 2;
adaptation  = 1;

nStimuli    = 2;
nTransPS    = 4;
transP_Train    = 0.25;
transP_Test     = 0.25;

shift   = 2;
a       = 0.125;

nRecordsPL  = 2;
nLayers     = 3;
vExcit  = [64,64,64];
rInhib  = 0.25;
pCnxEfE = [0, 1.0, 0.5];
pCnxElE = [0.5, 0.5, 0.5];
pCnxEI  = [1.0, 1.0, 1.0];
pCnxIE  = [1.0, 1.0, 1.0];
pCnxII  = [0, 0, 0];

tauC = 0.015;
tauD = 0.025;

capI = 50e-12;

modEf = 1.5;
Dg_ElE = 0.1;
Dg_IE = 1.0;
Dg_EI = 1.0;

gMax = 5e-9;
