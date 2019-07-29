# Tn-Core

Tn-Core is Matlab toolbox for generation of gene-centric, context-specific core metabolic network reconstructions through integration of experimental Tn-seq and RNA-seq datasets. The last named version is version 2.3. As the toolbox will continue to be updated, and to avoid countless version numbers, each change will no longer be associated with a new version number. Instead, changes will be made in Tn-Core Toolbox version stored in the directory entitled *Most-Recent-Version*. When using this version of the toolbox, the commit number can be used to document the version of the toolbox used.

Sample output data can be found in the Sample-Output-Data directory. Scripts used in the initial benchmarking of Tn-Core, along with the output data, are provided in the Benchmarking directory.

The manuscript describing Tn-Core can be found at:
diCenzo GC, Mengoni A, Fondi M (2019) Tn-Core: a toolbox for integrating Tn-seq gene essentiality data and constraint-based metabolic modelling. ACS Synthetic Biology. 8(1):158-169. doi:10.1021/acssynbio.8b00432.

## Dependencies

Tn-Core is dependent on the following seven softwares/toolboxes. Please note that in the case of FASTCORE, it is necessary to download the original version of FASTCORE as Tn-Core is not compatible with the version provided in the COBRA Toolbox.

1.    [MATLAB](https://mathworks.com)
2.    [COBRA Toolbox](https://opencobra.github.io/cobratoolbox/stable)
3.    [TIGER Toolbox version 1.2-beta](https://csbl.bitbucket.io/tiger/download.html)
4.    [FASTCORE version 1.0](https://uni.lu/forschung/fstc/life_sciences_research_unit/research_areas/systems_biology/software/fastcore)
5.    [libSBML](https://sourceforge.net/projects/sbml/files/libsbml)
6.    [SBML Toolbox](https://sourceforge.net/projects/sbml/files/SBMLToolbox)
7.    [iLOG CPLEX Studio](https://ibm.com/products/ilog-cplex-optimization-studio) or the [Gurobi solver](https://www.gurobi.com)
