from keyparameter import logu  # 对应 USE keyparameter, ONLY: logu
from keyflag import (  # 对应 USE keyflag，导入所需变量
    maxgen, TK, Tref, Mref,
    g2pfg, g2wfg, dhffg,
    pvap_sar, kohadd_sar, kisom_sar, kdiss_sar,
    rx_ro2_oh, rx_ro2_no2, rx_ro_no2, C_NO2, bimolecrx4criegee,
    critvp, isomerfg, brcut, highnoxfg, rx_ro2_multiclass
)

# ======================================================================
# Purpose: write operating conditions and flags for memo to a log file.
# ======================================================================
def wrtlog():
    # write mechanism parameters to stdout, for logging
    logu.write(' ---------------------------- \n')
    logu.write(' GECKO VERSION=v0.0.0 \n')
    logu.write(' ---------------------------- \n')
    logu.write(f' maxgen = {maxgen} \n')
    logu.write(f'  mechanism temperature (K) used for branching ratio: {TK:6.1f} \n')
    logu.write(f'  reference temperature (K): {Tref:6.1f} \n')
    logu.write(f'  default 3rd body M (molec cm-3): {Mref:1pe10.2} \n')
    logu.write('           \n')

    # -----

    logu.write(' ---------------------------- \n')
    logu.write('    ----  MECHANISM  ---- \n')
    logu.write(' ---------------------------- \n')

    if g2pfg:
        logu.write(' gas/particle mass transfer: included \n')
    else:
        logu.write(' gas/particle mass transfer: not included \n')

    if g2wfg:
        logu.write(' gas/wall mass transfer: included \n')
    else:
        logu.write(' gas/wall mass transfer: not included \n')

    if dhffg:
        logu.write(' DHF isomerisation: yes \n')
    else:
        logu.write(' DHF isomerisation: no \n')
    logu.write('           \n')

    # ----

    logu.write(' ---------------------------- \n')
    logu.write('    --------- SAR -------- \n')
    logu.write(' ---------------------------- \n')
    if pvap_sar == 1:
        logu.write(' vapor pressure scheme=M&Y \n')
    elif pvap_sar == 2:
        logu.write(' vapor pressure scheme=Nannoolal 2008 \n')
    elif pvap_sar == 3:
        logu.write(' vapor pressure scheme=SIMPOL-1 \n')

    if kohadd_sar == 1:
        logu.write(' OH addition SAR =Peeters 1997 \n')
    elif kohadd_sar == 2:
        logu.write(' OH addition SAR =Ziemann 2009 \n')
    elif kohadd_sar == 3:
        logu.write(' OH addition SAR =Jenkin 2018  \n')

    if kisom_sar == 1:
        logu.write(' Alkoxy isomerisation: Atkinson 2007 \n')
    elif kisom_sar == 2:
        logu.write(' Alkoxy isomerisation: Vereecken 2009 \n')

    if kdiss_sar == 1:
        logu.write(' Alkoxy dissociation: Atkinson 2007 \n')
    elif kdiss_sar == 2:
        logu.write(' Alkoxy dissociation: Vereecken 2009 \n')

    logu.write('           \n')
    logu.write(' ---------------------------- \n')
    logu.write('   -- OPTIONAL CHEMISTRY --   \n')
    logu.write(' ---------------------------- \n')
    if rx_ro2_oh:
        logu.write(' RO2 + OH reactions activated \n')
    if rx_ro2_oh:
        logu.write('    >> Warning : chemistry for (OOOH) is simple here and needs to be improved \n')
    if rx_ro2_no2:
        logu.write(' RO2 + NO2 reactions activated \n')
    if rx_ro2_no2:
        logu.write(f'    >> Warning : species bearing (OONO2) only react by ROONO2 -> RO2 + NO2 decomposition \n')

    if rx_ro_no2:
        logu.write(' RO + NO2 -> R(ONO2) reactions activated \n')
    if rx_ro_no2:
        logu.write(f'    >> Warning : rate calculated with [NO2]={C_NO2:1pe10.2} molec.cm-3 to eliminate minor pathways of RO reactions \n')
    if bimolecrx4criegee:
        logu.write(' Bimolecular reaction for stabilized criegee allowed (keep H2O) \n')
    logu.write('           \n')


    # ----

    logu.write(' ---------------------------- \n')
    logu.write('    ----- Reductions ----- \n')
    logu.write(' ---------------------------- \n')
    logu.write(f' critical vapor pressure (atm) = {critvp} \n')
    if isomerfg:
        logu.write(' isomerisation allowed \n')
    else:
        logu.write(' no isomerisation \n')
    logu.write(f' Cut off branching ratio below which a reaction pathway is ignored: {brcut} \n')
    logu.write(f' high-NOx flag: {highnoxfg} \n')
    logu.write(f' All classes of RO2 in RO2+RO2: {rx_ro2_multiclass} \n')
    logu.write('           \n')
    logu.write(' ---------------------------- \n')