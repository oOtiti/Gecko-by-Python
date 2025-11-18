import math, keyparameter, keyflag

# BA: Ugly subroutine ... waiting for a nicer one to write rate constants
# for C1 species with OH. Rate consistent with singlec_mech.dat file. 
def wrt_kc1():
    koh = None
    kno3 = None
    kT = None
    keyparameter.kohu.write(f"TREF TMECH  {keyflag.Tref:5.1f} {keyflag.TK:5.1f}\n")
    keyparameter.kno3u.write(f"TREF TMECH  {keyflag.Tref:5.1f} {keyflag.TK:5.1f}\n")

    # GCH4 + GHO => GCH3O2                          : 2.15E-12  0.   1735 ;
    koh = 2.15e-12 * math.exp(-1735. / keyflag.Tref)
    kT = 2.15e-12 * math.exp(-1735. / keyflag.TK)
    keyparameter.kohu.write(f"CH4    {koh:10.2e} {kT:10.2e} DATA, DATA\n")

    # GCH3O2 + GNO3 => GCH3O + GNO2                 : 1.30E-12  0.      0 ;
    kno3 = 1.30e-12
    kT = kno3
    keyparameter.kno3u.write(f"CH3O2  {kno3:10.2e} {kT:10.2e} DATA, COPY\n")

    # GCH3OH + GHO => 0.85 GCH2O + 0.85 GHO2 + 0.15 GCH3O : 3.10E-12  0. 360 ;
    koh = 3.10e-12 * math.exp(-360. / keyflag.Tref)
    kT = 3.10e-12 * math.exp(-360. / keyflag.TK)
    keyparameter.kohu.write(f"CH3OH  {koh:10.2e} {kT:10.2e} DATA, DATA\n")

    # GCH3OH + GNO3 => GCH2O + GHO2 + GHNO3         : 9.4E-13  0.   2650 ;
    kno3 = 9.40e-13
    kT = kno3
    keyparameter.kno3u.write(f"CH3OH  {kno3:10.2e} {kT:10.2e} DATA, COPY\n")

    # GCH3OOH + GHO => GHO + GCH2O                  : 1.00E-12  0.   -190 ;
    # GCH3OOH + GHO => GCH3O2                       : 1.90E-12  0.   -190 ;
    koh = 2.90e-12 * math.exp(190. / keyflag.Tref)
    kT = 2.90e-12 * math.exp(190. / keyflag.TK)
    keyparameter.kohu.write(f"CH3OOH {koh:10.2e} {kT:10.2e} DATA, DATA\n")

    # GN01003 + GHO => GH2O + GCH2O + GNO           : 0.301E-12  0.    0. ;
    koh = 3.01e-13
    kT = koh
    keyparameter.kohu.write(f"!N01003 {koh:10.2e} {kT:10.2e} DATA, COPY\n")  # commented - not produced or negligible

    # GN01001 + GHO => GCH2O + GNO2                 : 4.00E-13  0.    845 ;
    koh = 4.00e-13 * math.exp(-845. / keyflag.Tref)
    kT = 4.00e-13 * math.exp(-845. / keyflag.TK)
    keyparameter.kohu.write(f"!N01001 {koh:10.2e} {kT:10.2e} DATA, DATA\n")  # commented - not produced or negligible

    # GCH2O + GHO => GCO + GHO2                     : 8.60E-12  0.    -20 ;
    koh = 8.60e-12 * math.exp(20. / keyflag.Tref)
    kT = 8.60e-12 * math.exp(20. / keyflag.TK)
    keyparameter.kohu.write(f"CH2O   {koh:10.2e} {kT:10.2e} DATA, DATA\n")

    # GCH2O + GNO3 => GHNO3 + GCO + GHO2            : 5.80E-16  0.      0 ;
    kno3 = 5.80e-16
    kT = kno3
    keyparameter.kno3u.write(f"CH2O   {kno3:10.2e} {kT:10.2e} DATA, COPY\n")

    # GNO1001 + GHO => GHCOOH + GNO2                : 3.10E-12  0.      0 ;
    koh = 3.10e-12
    kT = koh
    keyparameter.kohu.write(f"!NO1001  {koh:10.2e} {kT:10.2e} DATA, COPY\n")  # commented - not produced or negligible

    # GHCOOH + GHO => GHO2 + GCO2                   : 4.50E-13  0.      0 ;
    kT = 4.50e-13
    kT = koh
    keyparameter.kohu.write(f"HCOOH  {koh:10.2e} {kT:10.2e} DATA, COPY\n")

    # GHCOO2H + GHO => GHCOO2                       : 1.90E-12  0.   -190 ;
    koh = 1.90e-12 * math.exp(190. / keyflag.Tref)
    kT = 1.90e-12 * math.exp(190. / keyflag.TK)
    keyparameter.kohu.write(f"HCOO2H  {koh:10.2e} {kT:10.2e} DATA, DATA\n")

    # GHCOO2 + GNO3 => 0.5 GCH3O2 + 0.5 GCO2 + GNO2   : 5.00E-12  0.      0 ;
    kno3 = 5.00e-12
    kT = kno3
    keyparameter.kno3u.write(f"HCOO2  {kno3:10.2e} {kT:10.2e} DATA, COPY\n")