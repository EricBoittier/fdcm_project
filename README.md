#  fMDCM - A Geometrically Dependent Charge Model for CHARMM

##  How to Make an fMDCM Model
1. Create a topology file for CHARMM, atom ordering is important.
2. Calculate the ESP using Gaussian, taking care to preserve ordering.
3. Fit an MDCM model, generation of the frames.txt file should also preserve ordering.
4. For fMDCM, an atom-charge-dictionary linking the charges to their parent atom is required. 
This must also preserve ordering.

## Example scripts
### Exporting a standard MDCM model to CHARMM
$BINDIR/comb-xyz-to-dcm.pl $XYZFILE $DCUBE $FRAMEFILE ${NAME}.dcm


