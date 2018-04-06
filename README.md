8 Photosynthesis Models. Each model has a two-level parameter hierarchy (genotype & individual level). 
Some parameters estimated only at genotype level of hierarchy others are predicted for each individual.
Models are based on 3 combinations of the following 3 assumptions:

CaCc_Jm - Assumes infinite mesophyll conductance, electron transport based on empirical formula using a maximum rate (Jmax), quantum yield and a curvature parameters there is no temperature dependency on parameters.

CaCc_Jf - Assumes infinite mesophyll conductance, electron transport based on chlorophyll fluoresces there is no temperature dependency on parameters.

CaCc_Jm_Temp - Assumes infinite mesophyll conductance, electron transport based on empirical formula using a maximum rate (Jmax), quantum yield and a curvature parameters there is an Arrhenius style temperature dependency on parameters.

CaCc_Jf_Temp - Assumes infinite mesophyll conductance, electron transport based on empirical formula using a maximum rate (Jmax), quantum yield and a curvature parameters there is an Arrhenius style temperature dependency on parameters.

CiCc_Jm - Assumes a mesophyll conductance limitation, electron transport based on empirical formula using a maximum rate (Jmax), quantum yield and a curvature parameters there is no temperature dependency on parameters.

CiCc_Jf - Assumes a mesophyll conductance limitation, electron transport based on chlorophyll fluoresces there is no temperature dependency on parameters.

CiCc_Jm_Temp - Assumes a mesophyll conductance limitation, electron transport based on empirical formula using a maximum rate (Jmax), quantum yield and a curvature parameters there is an Arrhenius style temperature dependency on parameters.

CiCc_Jf_Temp - Assumes a mesophyll conductance limitation, electron transport based on empirical formula using a maximum rate (Jmax), quantum yield and a curvature parameters there is an Arrhenius style temperature dependency on parameters.

_Model provides a model text for implementation in rjags

_Script provides model implementation on A/Ci data, diagnostics and plotting.

Key References

Original Model:

Farquhar, G.D., Caemmerer, S. V. and Berry, J. A. (1980) A biochemical model of photosynthetic CO2 assimilation in leaves of C3 species. Planta, 149, 78-90.

On incorporation of mesophyll conductance:

Ethier G.J. & Livingston N.J. (2004) On the need to incorporate sensitivity to CO2 transfer conductance into the Farquhar–von Caemmerer–Berry leaf photosynthesis model. Plant, Cell & Environment, 27, 137-153.

Similar Bayesian Implementation:

Patrick L.D., Ogle K. & Tissue D.T. (2009) A hierarchical Bayesian approach for estimation of photosynthetic parameters of C3 plants. Plant, Cell & Environment, 32, 1695-1709.

If you use any of these models, the compiled prior or posterior data please cite  publication:

Pleban JR 1, Mackay DS 1, Aston TL 2,
Ewers BE 2,3 and Weinig C 2,3,4 (2018)
Phenotypic Trait Identification Using a
Multimodel Bayesian Method: A Case
Study Using Photosynthesis in
Brassica rapa Genotypes.
Front. Plant Sci. 9:448.
doi: 10.3389/fpls.2018.00448

1-Geography, University at Buffalo, Buffalo, NY, USA 2-Botany, University of Wyoming, Laramie, WY, USA
