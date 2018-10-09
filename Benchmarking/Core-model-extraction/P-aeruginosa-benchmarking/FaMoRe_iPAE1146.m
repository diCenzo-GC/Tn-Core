function subnetworks = FaMoRe_meliloti()
%% Initial stuff

% Original Network
model = load('iPAE1146_consistent.mat');
load('iPAE1146_consistent.mat')

% Get maximal growth rate
sol = optimizeCbModel(iPAE1146_consistent, 'max');

% Number of functionalities
n_functionalities = 1;

%% Protected metabolites

prot_mets = [];

%% Protected Reactions

prot_rxns = [{'rJB00259'};
{'rJB00264'};
{'rJB00270'};
{'rJB00271'};
{'rPY00164'};
{'rPY00165'};
{'rPY00169'};
{'rPY00170'};
{'rxn00001'};
{'rxn00003'};
{'rxn00029'};
{'rxn00060'};
{'rxn00083'};
{'rxn00097'};
{'rxn00117'};
{'rxn00119'};
{'rxn00126'};
{'rxn00141'};
{'rxn00154'};
{'rxn00193'};
{'rxn00237'};
{'rxn00239'};
{'rxn00285'};
{'rxn00288'};
{'rxn00293'};
{'rxn00313'};
{'rxn00337'};
{'rxn00364'};
{'rxn00409'};
{'rxn00410'};
{'rxn00412'};
{'rxn00459'};
{'rxn00461'};
{'rxn00533'};
{'rxn00549'};
{'rxn00555'};
{'rxn00611'};
{'rxn00612'};
{'rxn00686'};
{'rxn00704'};
{'rxn00710'};
{'rxn00747'};
{'rxn00770'};
{'rxn00777'};
{'rxn00785'};
{'rxn00786'};
{'rxn00789'};
{'rxn00790'};
{'rxn00800'};
{'rxn00832'};
{'rxn00834'};
{'rxn00838'};
{'rxn00839'};
{'rxn00898'};
{'rxn00907'};
{'rxn00917'};
{'rxn00931'};
{'rxn01018'};
{'rxn01069'};
{'rxn01100'};
{'rxn01101'};
{'rxn01106'};
{'rxn01116'};
{'rxn01127'};
{'rxn01200'};
{'rxn01211'};
{'rxn01219'};
{'rxn01241'};
{'rxn01255'};
{'rxn01256'};
{'rxn01353'};
{'rxn01360'};
{'rxn01362'};
{'rxn01434'};
{'rxn01465'};
{'rxn01466'};
{'rxn01477'};
{'rxn01485'};
{'rxn01509'};
{'rxn01512'};
{'rxn01517'};
{'rxn01519'};
{'rxn01520'};
{'rxn01603'};
{'rxn01643'};
{'rxn01644'};
{'rxn01672'};
{'rxn01673'};
{'rxn01678'};
{'rxn01732'};
{'rxn01739'};
{'rxn01740'};
{'rxn01917'};
{'rxn01964'};
{'rxn01973'};
{'rxn01974'};
{'rxn02008'};
{'rxn02011'};
{'rxn02201'};
{'rxn02212'};
{'rxn02264'};
{'rxn02271'};
{'rxn02285'};
{'rxn02286'};
{'rxn02288'};
{'rxn02303'};
{'rxn02373'};
{'rxn02504'};
{'rxn02789'};
{'rxn02811'};
{'rxn02867'};
{'rxn02895'};
{'rxn02914'};
{'rxn02929'};
{'rxn02933'};
{'rxn03031'};
{'rxn03084'};
{'rxn03136'};
{'rxn03137'};
{'rxn03147'};
{'rxn03164'};
{'rxn03175'};
{'rxn03408'};
{'rxn03638'};
{'rxn03841'};
{'rxn03904'};
{'rxn03907'};
{'rxn03908'};
{'rxn03909'};
{'rxn03910'};
{'rxn03958'};
{'rxn04954'};
{'rxn05003'};
{'rxn05028'};
{'rxn05034'};
{'rxn05229'};
{'rxn05231'};
{'rxn05233'};
{'rxn05293'};
{'rxn05329'};
{'rxn05330'};
{'rxn05331'};
{'rxn05332'};
{'rxn05333'};
{'rxn05334'};
{'rxn05335'};
{'rxn05336'};
{'rxn05337'};
{'rxn05338'};
{'rxn05339'};
{'rxn05340'};
{'rxn05341'};
{'rxn05342'};
{'rxn05343'};
{'rxn05344'};
{'rxn05345'};
{'rxn05346'};
{'rxn05347'};
{'rxn05348'};
{'rxn05350'};
{'rxn05457'};
{'rxn05458'};
{'rxn05459'};
{'rxn05460'};
{'rxn05461'};
{'rxn05462'};
{'rxn05465'};
{'rxn06075'};
{'rxn06076'};
{'rxn06493'};
{'rxn06591'};
{'rxn06672'};
{'rxn06673'};
{'rxn06937'};
{'rxn07980'};
{'rxn07993'};
{'rxn08043'};
{'rxn08066'};
{'rxn08094'};
{'rxn08309'};
{'rxn08335'};
{'rxn08352'};
{'rxn08756'};
{'rxn09037'};
{'rxn09111'};
{'rxn09112'};
{'rxn09113'};
{'rxn09114'};
{'rxn09179'};
{'rxn09200'};
{'rxn09201'};
{'rxn09202'};
{'rxn09203'};
{'rxn09208'};
{'rxn09209'};
{'rxn09210'};
{'rxn09211'};
{'rxn09889'};
{'rxn10042'};
{'rxn10053'};
{'rxn10126'};
{'rxn10199'};
{'rxn11245'};
{'rxn12218'};
{'rxn12879'};
{'rxn12880'};
{'rxn12892'};
{'rxn13044'};
{'rxn13108'};
{'rxn13689'};
{'rxn13705'};
{'rxn13726'};
{'rxn13826'};
{'rxn13840'};
{'rxn13843'};
{'rxn13844'};
{'rxn13846'};
{'rxn13851'};
{'rxn13874'};
{'rxn13880'};
{'rxn13881'};
{'rxn13882'};
{'rxn13883'};
{'rxn13884'};
{'rxn13889'};];

% Minimal DoF
dof = 1;

%% Functionalities
f{1}.names = [{'PAO1_Biomass'}];
f{1}.type = ['>'];
f{1}.rhsValues = [0.5 * sol.f];

%% Tolerance
tol = 1e-06;

% Number of wanted minimal Networks:
number = 1;

% fast version 
fast = 0;

% delete blocked reactions:
use_F2C2 = 0;

% BigM formulation:
use_bigM = 1;

% Reduce the network:
requirements.n_functionalities = n_functionalities;
requirements.prot_mets = prot_mets;
requirements.prot_rxns = prot_rxns;
requirements.dof = dof;
requirements.f = f;
requirements.tol = tol;
requirements.number = number;
requirements.fast = fast;
requirements.use_F2C2 = use_F2C2;
requirements.use_bigM = use_bigM;

% get the model from the loaded struct
fieldname = fieldnames(model);
model = model.(fieldname{1});

% compute the subnetworks
subnetworks = reduce_model(model, requirements);

