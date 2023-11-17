
# Introduction
  
...TEXT...

Utilizes the isopattern function of the R package Envipat (Loos et al, Analytical Chemistry, 2015, 87, 5738-5744).



## Instructions    
Choose the parameters in the _Initial settings_ tab. Press submit and wait for calculation to finish. A table will then be generated with all ions that conform with the initial setting parameters. The table can be exported to excel by clicking on the "Excel" button at the top.  
The _Interfering ions_ tab can be used to check for ions that interfer with each other at the set resolution of the mass spectrometer. Default is set to R=60,000. The plots and tables are interactive and the user can filter the _"Interference at MS res?"_ by clicking on _FALSE_ on the plot legend (and thereby keeping all _TRUE_ ions, which will remove all ions that can be resolved by the set MS resolution).


## Initial settings tab
  
__C atoms min__ and __C atoms max__: is the range of number of carbon atoms to generate. Minimum value is C=3 and maximum is C=30. 
  
__Cl atoms min__ and __Cl atoms max__: is the range of number of chlorine atoms to generate. Maximum number of Cl has been limited to 15.  
  
__Add adducts/fragments__: refers to the formula of adducts and/or fragments to generate from a set list of available options. Multiple selections are possible.  
  
__[CP]__ and __[CO]__: refers to either chlorinated paraffins [CP] or chlorinated (mono)olefins [CO]. 
  
__[xx-yy]__ or __[xx-yy-zz]__: where _-yy_ refers to the adduct/fragment ions or _-yy-zz_ which are the consecutive loss fragments. Currently, a limited selection is available but more can be added upon request.  
  
__[xx]-__ and __[xx]+__ refers to the charge of the ion (limited to single charged species, +1 or -1).  
  
_Note that [M+Cl-HCl]- can also be written as [M-H]-_  
  
__Isotope rel ab threshold (%)__: is the threshold for relative abundance for isotopologues for each chemical formula of the adduct/fragment ion. Ions below this threshold will not be included into the generated ion table.
  
### Output table  
  
__Parent_Formula__: the chemical formula of the molecular ion.  

__Cl_perc__: the chlorination degree (mass percentage of Cl).  
  
__Charge__: The charge of the ion.  
  
__Fragment__: The fragment and isotopic type of the ion species.  
  
__Frag_MonoIso_Formula__: the chemical formula of the adduct/fragment ion.  
  
__Isotopologue__: the isotopologue in relation to the monoisotopic ion.  
  
__Isotope_Formula__: the exact isotopic formula of the adduct/fragment ion.  
  
__m/z__: the mass-over-charge of the adduct/fragment ion.  
  
__Rel_ab__: the relative abundance of the different isotopologues of each adduct/fragment ion.  
  
__12C__, __13C__, __1H__, __2H__, __35Cl__, __37Cl__: the number of atoms for each element
  
## Interfering ions tab  
  
__difflag__, __difflead__: internal calculations for the difference in m/z values between the two nearest ions. If "interference at MS res?" filter has been used, then the previous/next ions might not be shown. 
  
__reslag__, __reslead__: internal calculations for the MS resolution needed to separate the two nearest ions. If "interference at MS res?" filter has been used, then the previous/next ions might not be shown.  
  
__interference__: indicate whether or not the m/z two nearest ions can interfere with each other at the set MS resolution value. _"false"_ means no interference and _"true"_ means there is interference (and therefore the MS resolution cannot resolve these peaks).  
  

### Plot outputs  
  
Note: Some bars can exceed 100% in relative abundance in the y-axis, and this indicates that some isotopologues have the exact same ion formula. You can hover over the different segments of these bars to check the overlapping mass ions.

__Hover text__:   
__difflag & difflead (prev and next)__: the difference between the m/z with previous and next ions (axis ordered from lowest to highest m/z).  
__reslag & reslead (prev and next)__: the MS resolution needed to resolve previous and next ion (axis ordered from lowest to highest m/z).  
  
  
## Skyline tab  

__TO BE ADDED__
  
  



