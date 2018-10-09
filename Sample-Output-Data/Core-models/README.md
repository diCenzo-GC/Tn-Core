# Core and refined models

The Excel_models directory contains each of the models described below formatted as tables in Excel workbooks. The COBRA_models directory contains each of the models described below formatted as COBRA models in MATLAB files.

The following models are included in this data set:
- iGD1575_model: The starting version of the S. meliloti iGD1575 genome-scale metabolic model
- Draft_meliloti_model: A draft, fully automated S. meliloti genome-scale metabolic model prepared with KBase
- iGD726_model: The S. meliloti manually prepared core metabolic model iGD726
- TnCore_model_1: The core model extracted from iGD1575_model using Tn-Core with Tn-seq data but without RNA-seq data
- TnCore_model_2: The core model extracted from iGD1575_model using Tn-Core with Tn-seq data and with RNA-seq data, without re-introduction of highly expressed genes
- TnCore_model_3: The core model extracted from iGD1575_model using Tn-Core with Tn-seq data and with RNA-seq data, with re-introduction of highly expressed genes
- FASTCORE_model: The core model extracted from iGD1575_model using a FASTCORE-based pipeline adapted to deal with Tn-seq data
- minNW_model: The core model extracted from iGD1575_model using a minNW-based pipeline adapted to deal with Tn-seq data
- GIMME_model_1: The core model extracted from iGD1575_model using a GIMME-based pipeline adapted to deal with Tn-seq data
- GIMME_model_2: The core model extracted form iGD1575_model using a GIMME-based pipeline using RNA-seq data
