from molvs import charge, fragment, standardize, tautomer
from rdkit import Chem
from rdkit.Chem.Descriptors import ExactMolWt
from loguru import logger
from dimorphite_dl import DimorphiteDL


class SMILESStandardizer:
    """
    A class for standardizing SMILES (Simplified Molecular Input Line Entry System) strings.

    This class provides functionality to standardize chemical structures represented as SMILES strings
    through various chemical standardization steps including protonation, tautomer canonicalization,
    and fragment selection.

    Attributes
    ----------
    ALLOWED_ATOMIC_NUMBERS : set
        Set of allowed atomic numbers in the molecules.
        Contains: {H(1), Li(3), C(6), N(7), O(8), F(9), Na(11), Mg(12), Si(14),
                  P(15), S(16), Cl(17), K(19), Ca(20), Br(35), I(53)}
    """
    # Class-level constants
    ALLOWED_ATOMIC_NUMBERS = {1, 3, 6, 7, 8, 9, 11, 12, 14, 15, 16, 17, 19, 20, 35, 53}
    # ALLOWED ATOMS :: {H, Li, C, N, O, F, Na, Mg, Si, P, S, Cl, K, Ca, Br, I}

    def __init__(self):
        """
        Initialize the SMILESStandardizer with various chemical standardization tools.

        Sets up standardizer, fragment chooser, uncharger, and tautomer canonicalizer
        from the molvs package.
        """
        # Initialize standardization tools
        self.molvs_standardizer = standardize.Standardizer()
        self.fragment_chooser = fragment.LargestFragmentChooser()
        self.uncharger = charge.Uncharger()
        self.tautomer_canonicalizer = tautomer.TautomerCanonicalizer()

    def protonate_smiles(self, smiles, min_ph=7.4, max_ph=7.4):
        """
        Protonate a given SMILES string using DimorphiteDL in the specified pH range.
        
        Parameters
        ----------
        smiles : str
            The SMILES string to protonate
        min_ph : float, optional
            The minimum pH value for protonation, by default 7.4
        max_ph : float, optional
            The maximum pH value for protonation, by default 7.4
            
        Returns
        -------
        str
            The protonated SMILES string. If protonation fails, returns the original SMILES.

        Notes
        -----
        Uses DimorphiteDL for pH-dependent protonation state prediction.
        Logs the protonation attempt and result using loguru logger.
        """
        try:
            #logger.debug(f"Attempting to protonate SMILES: {smiles} at pH {min_ph} to {max_ph}")
            dimorphite = DimorphiteDL(min_ph=min_ph, max_ph=max_ph, pka_precision=0)
            protonated_smiles = dimorphite.protonate(smiles)
            result = protonated_smiles[0] if protonated_smiles else smiles
            #logger.debug(f"Protonated SMILES: {result}")
            return result
        except Exception as e:
            #logger.error(f"Error protonating SMILES {smiles}: {e}")
            return smiles

    def standardize(self, smi, pH=7.4):
        """
        Standardize a SMILES string using various chemical standardization steps.
        
        Parameters
        ----------
        smi : str
            The SMILES string to standardize
        pH : float, optional
            The pH value for protonation state prediction, by default 7.4
            
        Returns
        -------
        str or None
            Standardized SMILES string if successful, None if standardization fails
            
        Notes
        -----
        The standardization process includes:
        1. Checking for allowed atomic numbers
        2. Basic structure standardization
        3. Selecting the largest fragment
        4. Removing charges
        5. Canonicalizing tautomers
        6. Protonating at specified pH
        
        Additional requirements for successful standardization:
        - Molecular weight must be less than 1500
        - Must contain more than 3 carbon atoms
        """
        try:
            mol = Chem.MolFromSmiles(smi)
            atoms_in_mol = set([i.GetAtomicNum() for i in mol.GetAtoms()])
            
            if atoms_in_mol - self.ALLOWED_ATOMIC_NUMBERS:
                logger.exception(f"Failed to Standardize :: {smi}")
                return None
                
            mol = self.molvs_standardizer(mol)
            #logger.debug(f"Molvs Standardized :: {Chem.MolToSmiles(mol)}")
            mol = self.fragment_chooser(mol)
            #logger.debug(f"Fragment Chosen :: {Chem.MolToSmiles(mol)}")
            mol = self.uncharger.uncharge(mol)
            #logger.debug(f"Uncharged :: {Chem.MolToSmiles(mol)}")
            mol = self.tautomer_canonicalizer.canonicalize(mol)
            #logger.debug(f"Tautomer Canonicalized :: {Chem.MolToSmiles(mol)}")
            
            if (
                ExactMolWt(mol) < 1500
                and len(mol.GetSubstructMatches(Chem.MolFromSmarts("[#6]"))) > 3
            ):
                _smi = Chem.MolToSmiles(mol)
                _smi = self.protonate_smiles(_smi, min_ph=pH, max_ph=pH)
                #logger.debug(f"Protonated at pH {pH} :: {_smi}")
                _smi = _smi.replace('@','')
                #logger.debug(f"Removed Stereochemistry :: {_smi}")
                return _smi
                
            #logger.exception(f"Failed to Standardize :: {smi}")
            return None
        except:
            #logger.exception(f"Failed to Standardize :: {smi}")
            return None
