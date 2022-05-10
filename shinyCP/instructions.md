
# Introduction
  
...TEXT...

## Instructions    
Choose the parameters in the _Initial settings_ tab. Press submit and wait for calculation to finish. A table will then be generated with all ions that conform with the initial setting parameters. The table can be exported to excel by clicking on the "Excel" button at the top.  
The _Interfering ions_ tab can be used to check for ions that interfer with each other at the set resolution of the mass spectrometer. Default is set to R=60,000. The plots and tables are interactive and the user can filter the _Interference at MS res?_ by clicking on _FALSE_ on the plot legend (and thereby keeping all _TRUE_ ions, which will remove all ions that can be resolved by the set MS resolution).


## Initial settings tab
  
_C atoms min_ and _C atoms max_: is the range of number of carbon atoms to generate. Minimum value is C=3 and maximum is C=30. 
  
_Cl atoms min_ and _Cl atoms max_: is the range of number of chlorine atoms to generate. Minimum value is C=3 and maximum is C=30.  
  
_Add adducts/fragments_: refers to the formula of adducts and/or fragments to generate from a set list of available options. Multiple selections are possible.  
  
_[CP]_ and _[CO]_: refers to either chlorinated paraffins [CP] or chlorinated (mono)olefins [CO]. 
  
_[xx-yy]_ or _[xx-yy-zz]_: where _-yy_ or _-yy-zz refers to the adduct/fragment ions. Currently, a limited selection is available but more can be added later if requested.  
  
_[xx]-_ and _[xx]+_ refers to the charge of the ion (limited to single charged species, +1 or -1).  
  
[M+Cl-HCl]- can be written as [M-H]-  
  
_Isotope rel ab threshold (%)_: is the threshold for relative abundance for isotopologues for each chemical formula of the adduct/fragment ion. 
  
### Output table  
  
_Parent_Formula_: the chemical formula of the molecular ion.  
  
_Charge_: The charge of the ion.  
  
_Fragment_: The fragment and isotopic type of the ion species.  
  
_Frag_MonoIso_Formula_: the chemical formula of the adduct/fragment ion.  
  
_Isotopologue_: the isotopologue in relation to the monoisotopic ion.  
  
_Isotope_Formula_: the exact isotopic formula of the adduct/fragment ion.  
  
_m/z_: the mass-over-charge of the adduct/fragment ion.  
  
_Rel_ab_: the relative abundance of the different isotopologues of each adduct/fragment ion.  
  
_12C, 13C, 1H, 2H, 35Cl, 37Cl_: the number of atoms for each element
  
## Interfering ions tab  
  
_difflag, difflead_: internal calculations for the difference in m/z values between the two nearest ions.  
  
_reslag, reslead_: internal calculations for the MS resolution needed to separate the two nearest ions.  
  
_interference_: indicate whether or not the m/z two nearest ions can interfere with each other at the set MS resolution value. _"false"_ means no interference and _"true"_ means there is interference (and therefore the MS resolution cannot resolve these peaks).  
  

### Hover text
_difflag & difflead (prev and next)_: the difference between the m/z with previous and next ions (axis ordered from lowest to highest m/z).  
_reslag & reslead (prev and next)_: the MS resolution needed to resolve previous and next ion (axis ordered from lowest to highest m/z).  






