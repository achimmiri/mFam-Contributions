# mFam-Contributions
Workflow to process tandem MS files and build MassBank records.Given the input metadata file along with raw MSP file(generated by MS-DIAL from raw abf files) .mFam Contributons functions include automated extraction of tandem MS spectra, spectrum cleanup,formula assignment to tandem MS fragments,SMILES,INCHI,INCHIKEY assignment from public database,recalibration of tandem MS spectra with assigned fragments,and extraction MassBank records.

Given your Input Project folder; should have following sub folders Information to be supplied


the Sheet for the meta data of the sample, compound, spectrum, and authors. The table contains one row for each spectrum


    File MSP file containing the spectrum; located under path raw data/exported as raw msp in the project folder (mandatory, e.g. Rutin.msp)
    Name Preferred compound name (mandatory, e.g. Rutin)
    Synonyms Common compound name synonyms (optional, e.g. Rutinoside; Eldrin)
    Structure / Unique ID The compound structure - one of the following three options (mandatory, e.g. 5280805)
        InChI
        SMILES
        PubChem CID
    Formula Molecular formula of the compound (optional, e.g. C27H30O16)
    Exact mass The exact mass of the compound (optional, e.g. 610.15339)
    CHROMATOGRAPHY The name of the chromatography instrumentation (mandatory, e.g. Symmetry C18 Column, Waters)
    RT (min) The retention time of the parent ion (mandatory, e.g. 9.5)
    INSTRUMENT The name of the mass spectrometry instrumentation (mandatory, e.g. LTQ Orbitrap XL, Thermo Scientfic; HP-1100 HPLC, Agilent)
    INSTRUMENT_TYPE The type of the mass spectrometry instrumentation (mandatory, e.g. LC-ESI-ITFT)
    IONIZATION The kind of ionization (mandatory, e.g. ESI) (mandatory, e.g. ESI)
    Collision energy Collision energy for fragmentation (mandatory, e.g. 15eV)
    Adduct The parent ion species (mandatory, e.g. [M-H]-)
    MS/MS Acquisition mode The acquisition mode (mandatory, e.g. DDA)
    Ionization mode The ionization mode (mandatory, e.g. Negative)
    Confidence The confidence of compound identification (mandatory, Pure standard or Predicted)
    Compound class The compound class (optional, e.g. Flavonoid)
    Found in The sample origin (optional, e.g. Tomato leaf (Solanum habrochaites LA1777))
    Database links Links to compound databases (optional, e.g. KEGG C05625; ChEBI 28527)
    Authors All authors who substantially contributed to the spectrum (mandatory, e.g. Ales Svatos, Ravi Kumar Maddula, MPI for Chemical Ecology, Jena, Germany)
    Publication Relevant pblication related to the spectrum (optional, e.g. F. Rasche, A. Svatos, R.K. Maddula, C. Boettcher and S. Boecker. Computing
    fragmentation trees from tandem mass spectrometry data. Anal. Chem., 2011, 83, 1243-1251)                                                                                  

Folder Strcture for the input data 

raw data/exported as raw msp[Exported as Raw MSP file generated by MS_DIAL which has the following information]

NAME
PRECURSORMZ: 
PRECURSORTYPE
SCANNUMBER
RETENTIONTIME
INTENSITY
ISOTOPE
DECONVOLUTED_PEAKHEIGHT
DECONVOLUTED_PEAKAREA
UNIQUEMS
Num Peaks

converted to msp
Results of the Filtered MSP files in relation to meta data 

converted to MassBank
MassBank Records for the passed MSP files


