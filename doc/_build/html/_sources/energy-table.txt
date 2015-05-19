=====================   ==================================================  ===================================================================
Attribute               Output file                                         Description
=====================   ==================================================  ===================================================================
energy                  ``HF Done`` (last occurence)                        E :sub:`el`
energies                ``HF Done`` (all occurences)                        list of all E :sub:`el`
zpe                     ``Zero-point correction``                           ZPE
thermalcorrection       ``Thermal correction to Energy``                    E :sub:`corr` = ZPE + E :sub:`trans` + E :sub:`rot` + E :sub:`vib`
enthalpycorrection      ``Thermal correction to Enthalpy``                  H :sub:`corr` = E :sub:`corr` + k :sub:`b` T
gibbscorrection         ``Thermal correction to Gibbs Free Energy``         G :sub:`corr` = H :sub:`corr` - T S :sub:`tot`
zpesum                  ``Sum of electronic and zero-point Energies``       E :sub:`el` + ZPE
thermal                 ``Sum of electronic and thermal Energies``          E :sub:`tot` = E :sub:`el` + E :sub:`corr`
enthalpy                ``Sum of electronic and thermal Enthalpies``        H = E :sub:`el` + H :sub:`corr`
gibbs                   ``Sum of electronic and thermal Free Energies``     G = E :sub:`el` + G :sub:`corr`
=====================   ==================================================  ===================================================================
