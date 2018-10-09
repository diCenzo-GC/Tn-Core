function subnetworks = FaMoRe_meliloti()
%% Initial stuff

% Original Network
model = load('iGD1575_consistent.mat');
load('iGD1575_consistent.mat')

% Get maximal growth rate
sol = optimizeCbModel(iGD1575_consistent, 'max');

% Number of functionalities
n_functionalities = 1;


%% Protected metabolites

prot_mets = [];

%% Protected Reactions

prot_rxns = [{'rxn00001_c0'};
{'rxn00011_c0'};
{'rxn00060_c0'};
{'rxn00077_c0'};
{'rxn00085_c0'};
{'rxn00086_c0'};
{'rxn00097_c0'};
{'rxn00100_c0'};
{'rxn00117_c0'};
{'rxn00119_c0'};
{'rxn00122_c0'};
{'rxn00126_c0'};
{'rxn00141_c0'};
{'rxn00190_c0'};
{'rxn00199_c0'};
{'rxn00209_c0'};
{'rxn00237_c0'};
{'rxn00248_c0'};
{'rxn00250_c0'};
{'rxn00256_c0'};
{'rxn00262_c0'};
{'rxn00285_c0'};
{'rxn00288_c0'};
{'rxn00293_c0'};
{'rxn00299_c0'};
{'rxn00302_c0'};
{'rxn00313_c0'};
{'rxn00337_c0'};
{'rxn00338_c0'};
{'rxn00392_c0'};
{'rxn00407_c0'};
{'rxn00409_c0'};
{'rxn00410_c0'};
{'rxn00412_c0'};
{'rxn00414_c0'};
{'rxn00420_c0'};
{'rxn00423_c0'};
{'rxn00441_c0'};
{'rxn00459_c0'};
{'rxn00461_c0'};
{'rxn00474_c0'};
{'rxn00515_c0'};
{'rxn00599_c0'};
{'rxn00604_c0'};
{'rxn00623_c0'};
{'rxn00646_c0'};
{'rxn00704_c0'};
{'rxn00710_c0'};
{'rxn00726_c0'};
{'rxn00727_c0'};
{'rxn00770_c0'};
{'rxn00781_c0'};
{'rxn00782_c0'};
{'rxn00789_c0'};
{'rxn00790_c0'};
{'rxn00791_c0'};
{'rxn00792_c0'};
{'rxn00799_c0'};
{'rxn00800_c0'};
{'rxn00832_c0'};
{'rxn00834_c0'};
{'rxn00838_c0'};
{'rxn00839_c0'};
{'rxn00851_c0'};
{'rxn00902_c0'};
{'rxn00917_c0'};
{'rxn00973_c0'};
{'rxn00974_c0'};
{'rxn01000_c0'};
{'rxn01018_c0'};
{'rxn01100_c0'};
{'rxn01127_c0'};
{'rxn01255_c0'};
{'rxn01256_c0'};
{'rxn01268_c0'};
{'rxn01269_c0'};
{'rxn01301_c0'};
{'rxn01302_c0'};
{'rxn01304_c0'};
{'rxn01331_c0'};
{'rxn01332_c0'};
{'rxn01353_c0'};
{'rxn01360_c0'};
{'rxn01362_c0'};
{'rxn01387_c0'};
{'rxn01388_c0'};
{'rxn01434_c0'};
{'rxn01465_c0'};
{'rxn01476_c0'};
{'rxn01512_c0'};
{'rxn01517_c0'};
{'rxn01520_c0'};
{'rxn01603_c0'};
{'rxn01643_c0'};
{'rxn01672_c0'};
{'rxn01673_c0'};
{'rxn01678_c0'};
{'rxn01682_c0'};
{'rxn01739_c0'};
{'rxn01790_c0'};
{'rxn01791_c0'};
{'rxn01871_c0'};
{'rxn01872_c0'};
{'rxn01964_c0'};
{'rxn01973_c0'};
{'rxn01974_c0'};
{'rxn01975_c0'};
{'rxn01977_c0'};
{'rxn02008_c0'};
{'rxn02155_c0'};
{'rxn02175_c0'};
{'rxn02185_c0'};
{'rxn02186_c0'};
{'rxn02187_c0'};
{'rxn02200_c0'};
{'rxn02201_c0'};
{'rxn02212_c0'};
{'rxn02213_c0'};
{'rxn02264_c0'};
{'rxn02285_c0'};
{'rxn02286_c0'};
{'rxn02331_c0'};
{'rxn02341_c0'};
{'rxn02342_c0'};
{'rxn02373_c0'};
{'rxn02376_c0'};
{'rxn02380_c0'};
{'rxn02402_c0'};
{'rxn02405_c0'};
{'rxn02465_c0'};
{'rxn02473_c0'};
{'rxn02474_c0'};
{'rxn02476_c0'};
{'rxn02503_c0'};
{'rxn02504_c0'};
{'rxn02507_c0'};
{'rxn02508_c0'};
{'rxn02775_c0'};
{'rxn02789_c0'};
{'rxn02811_c0'};
{'rxn02834_c0'};
{'rxn02835_c0'};
{'rxn02895_c0'};
{'rxn02897_c0'};
{'rxn02914_c0'};
{'rxn02929_c0'};
{'rxn02937_c0'};
{'rxn02938_c0'};
{'rxn02988_c0'};
{'rxn03062_c0'};
{'rxn03068_c0'};
{'rxn03084_c0'};
{'rxn03135_c0'};
{'rxn03136_c0'};
{'rxn03137_c0'};
{'rxn03146_c0'};
{'rxn03147_c0'};
{'rxn03159_c0'};
{'rxn03164_c0'};
{'rxn03174_c0'};
{'rxn03175_c0'};
{'rxn03181_c0'};
{'rxn03194_c0'};
{'rxn03394_c0'};
{'rxn03395_c0'};
{'rxn03397_c0'};
{'rxn03408_c0'};
{'rxn03419_c0'};
{'rxn03421_c0'};
{'rxn03435_c0'};
{'rxn03436_c0'};
{'rxn03437_c0'};
{'rxn03445_c0'};
{'rxn03491_c0'};
{'rxn03512_c0'};
{'rxn03514_c0'};
{'rxn03534_c0'};
{'rxn03536_c0'};
{'rxn03537_c0'};
{'rxn03538_c0'};
{'rxn03540_c0'};
{'rxn03638_c0'};
{'rxn03893_c0'};
{'rxn03901_c0'};
{'rxn03904_c0'};
{'rxn03907_c0'};
{'rxn03908_c0'};
{'rxn03909_c0'};
{'rxn03910_c0'};
{'rxn03958_c0'};
{'rxn04070_c0'};
{'rxn04113_c0'};
{'rxn04384_c0'};
{'rxn04385_c0'};
{'rxn04413_c0'};
{'rxn04996_c0'};
{'rxn05144_c0'};
{'rxn05149_c0'};
{'rxn05229_c0'};
{'rxn05231'};
{'rxn05233'};
{'rxn05239_c0'};
{'rxn05256_c0'};
{'rxn05293_c0'};
{'rxn05294_c0'};
{'rxn05322_c0'};
{'rxn05323_c0'};
{'rxn05324_c0'};
{'rxn05325_c0'};
{'rxn05326_c0'};
{'rxn05327_c0'};
{'rxn05328_c0'};
{'rxn05329_c0'};
{'rxn05330_c0'};
{'rxn05331_c0'};
{'rxn05332_c0'};
{'rxn05333_c0'};
{'rxn05334_c0'};
{'rxn05335_c0'};
{'rxn05343_c0'};
{'rxn05344_c0'};
{'rxn05345_c0'};
{'rxn05346_c0'};
{'rxn05347_c0'};
{'rxn05348_c0'};
{'rxn05349_c0'};
{'rxn05350_c0'};
{'rxn05460_c0'};
{'rxn05462_c0'};
{'rxn05465'};
{'rxn05744_c0'};
{'rxn06075'};
{'rxn06076'};
{'rxn06280'};
{'rxn06300'};
{'rxn06432'};
{'rxn06435'};
{'rxn06438'};
{'rxn06439'};
{'rxn06440'};
{'rxn06441'};
{'rxn06442'};
{'rxn06443'};
{'rxn06444'};
{'rxn06445'};
{'rxn06446'};
{'rxn06447'};
{'rxn06448'};
{'rxn06449'};
{'rxn06538'};
{'rxn06672_c0'};
{'rxn06673_c0'};
{'rxn06723_c0'};
{'rxn06729_c0'};
{'rxn06864'};
{'rxn06865_c0'};
{'rxn06936'};
{'rxn06937'};
{'rxn08335_c0'};
{'rxn08336_c0'};
{'rxn08352_c0'};
{'rxn08547_c0'};
{'rxn08548_c0'};
{'rxn08551_c0'};
{'rxn08552_c0'};
{'rxn08756_c0'};
{'rxn08885_c0'};
{'rxn09177_c0'};
{'rxn09203_c0'};
{'rxn09211_c0'};
{'rxn10030_c0'};
{'rxn10042_c0'};
{'rxn10052_c0'};
{'rxn10060_c0'};
{'rxn10202_c0'};
{'rxn10203_c0'};
{'rxn10204_c0'};
{'rxn10476_c0'};
{'rxn10534'};
{'rxn11544_c0'};
{'rxn11545_c0'};
{'rxn11595'};
{'rxn11946_c0'};
{'rxn12192'};
{'rxn12224_c0'};
{'rxn12510_c0'};
{'rxn12512_c0'};
{'rxn13784_c0'};
{'rxnDMPE'};
{'rxnMMPE'};
{'rxnPC'}];

% Minimal DoF
dof = 1;

%% Functionalities
f{1}.names = [{'biomass_bulk_c0'}];
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

