# 模块级变量
nrx_all = 0
nrx_n = 0
nrx_hv = 0
nrx_extra = 0
nrx_fo = 0
nrx_tb = 0
nrx_o2 = 0
nrx_meo2 = 0
nrx_ro2 = 0
nrx_isom = 0
nrx_tabcf = 0
nrx_ain = 0
nrx_aou = 0
nrx_win = 0
nrx_wou = 0


def rxinit(r, s, p, arrh, idreac, nlabel, xlabel, folow, fotroe):
    for i in range(len(r)):
        r[i] = ' '
    for i in range(len(s)):
        s[i] = 0.0
    for i in range(len(p)):
        p[i] = ' '
    for i in range(len(folow)):
        folow[i] = 0.0
    for i in range(len(fotroe)):
        fotroe[i] = 0.0
    for i in range(len(arrh)):
        arrh[i] = 0.0
    idreac = 0
    nlabel = 0
    xlabel = 0.0

    # return idreac, nlabel, xlabel so callers can receive updated values
    return idreac, nlabel, xlabel


def rxwrit(lout, r, s, p, arrh, idreac, nlab, xlabel, folow, fotroe, com=None, phase=None):
    from keyparameter import refu, mxnp, mxlco
    from keyflag import wrtref
    from references import mxlcod
    from tempflag import iflost, xxc
    from toolbox import addref

    mnp = len(s)
    if phase is not None:
        if len(phase) != 1:
            print("in rxwrit,unexpected phase string")
            raise Exception("in rxwrit")
    if phase is not None:
        rphase = phase
    else:
        rphase = 'G'

    if com is not None:
        ref = com[:]
        nref = 0
        for i in range(len(ref)):
            if ref[i] != ' ':
                nref = nref + 1
    else:
        ref = [' ']
        nref = 0
    ref1 = [' '] * len(ref)

    rc1 = arrh[0]
    signc = [' '] * len(s)
    charstoi = [' '] * len(s)
    pg = [' '] * len(s)
    s1 = s[:]
    p1 = p[:]
    xlab = xlabel
    s2 = [0.0] * len(s)
    p2 = [' '] * len(s)

    if iflost == 1:
        for i in range(len(s1)):
            s1[i] = 0.0
        for i in range(len(p1)):
            p1[i] = ' '
        s1[0] = float(xxc)
        p1[0] = 'XCLOST'

    for i in range(len(s1)):
        if s1[i] == 0.0:
            p1[i] = ' '
    for i in range(len(p1)):
        if p1[i] == ' ':
            s1[i] = 0.0

    for i in range(mnp - 1):
        if p1[i][0] != ' ':
            for j in range(i + 1, mnp):
                if p1[j] == p1[i]:
                    s1[i] = s1[i] + s1[j]
                    s1[j] = 0.0
                    p1[j] = ' '

    np = 0
    for i in range(mnp):
        if p1[i] != ' ':
            np = np + 1
            s2[np-1] = s1[i]
            p2[np-1] = p1[i]

    for i in range(np):
        if s2[i] == 0.0:
            continue
        if s2[i] == 1.0:
            continue
        charstoi[i] = f"{abs(s2[i]):6.3f}"

    rg = [' '] * len(r)
    idrx = 0
    for i in range(len(r)):
        if r[i][0] != ' ':
            if r[i][:3] == 'HV ':
                idrx = 1
                rg[i] = 'HV'
            elif r[i][:6] == 'EXTRA ':
                idrx = 2
                rg[i] = 'EXTRA'
            elif r[i][:4] == '(+M)':
                idrx = 3
                rg[i] = 'FALLOFF'
            elif r[i][:2] == 'M ':
                idrx = 4
                rg[i] = 'TBODY'
            elif r[i][:6] == 'OXYGEN':
                idrx = 5
                rg[i] = 'OXYGEN'
            elif r[i][:4] == 'PERO':
                idrx = 6
                rg[i] = r[i]
            elif r[i][:6] == 'MEPERO':
                idrx = 7
                rg[i] = 'MEPERO'
            elif r[i][:6] == 'ISOM  ':
                idrx = 8
                rg[i] = 'ISOM'
            elif r[i][:6] == 'TABCF ':
                idrx = 9
                rg[i] = 'TABCF'
            else:
                rg[i] = rphase + r[i]

    for i in range(np):
        if p2[i] != ' ':
            locheck = False
            if p2[i][:4] == '(+M)':
                locheck = True
                p2[i] = ' '
            if p2[i] == 'EMPTY ':
                locheck = True
            if locheck:
                pg[i] = p2[i]
                if p2[i] == 'EMPTY ':
                    pg[i] = 'NOTHING'
            else:
                pg[i] = rphase + p2[i]

    c1 = ' '
    c2 = ' '
    if rg[1] != ' ':
        c1 = '+'
    if rg[2] != ' ':
        c2 = '+'

    for i in range(1, np):
        j = i - 1
        if p2[i] != ' ':
            if s2[i] < 0.0:
                signc[j] = '-'
            else:
                signc[j] = '+'

    if wrtref:
        for i in range(nref - 1):
            for j in range(i + 1, nref):
                if ref[i] == ref[j]:
                    ref[j] = ' '
        nref1 = 0
        for i in range(nref):
            if ref[i] == ' ':
                continue
            nref1 = nref1 + 1
            ref1[nref1-1] = ref[i]

        refline = ' '
        for i in range(nref1):
            ipos = 13 * (i) + 1
            refline = refline + " ; " + ref1[i]

    count4rxn(idrx)

    nsplit = 1
    if np > mxnp:
        nsplit = (np - 1) // mxnp + 1

    if nsplit == 1:
        j1 = 1
        j2 = 4
        line = f"{rg[0]:7s} {c1} {rg[1]:7s} {c2} {rg[2]:7s}  =>"
        for jj in range(j1-1, j2):
            if jj < np:
                line += f" {charstoi[jj]:6s} {pg[jj]:7s} {signc[jj]}"
            else:
                line += " " * 16
        line += f"    {rc1:10.3e} {arrh[1]:4.1f} {arrh[2]:7.0f}"
        print(line, file=lout)
        if wrtref:
            refline_out = f"{rg[0]:7s} {c1} {rg[1]:7s} {c2} {rg[2]:7s}  =>"
            for jj in range(j1-1, j2):
                if jj < np:
                    refline_out += f" {charstoi[jj]:6s} {pg[jj]:7s} {signc[jj]}"
                else:
                    refline_out += " " * 16
            refline_out += f"    {rc1:10.3e} {arrh[1]:4.1f} {arrh[2]:7.0f}  {refline}"
            print(refline_out, file=refu)
    else:
        for k in range(1, nsplit + 1):
            j1 = 1 + mxnp * (k - 1)
            j2 = j1 + mxnp - 1

            if k == 1:
                line = f"{rg[0]:7s} {c1} {rg[1]:7s} {c2} {rg[2]:7s}  =>"
                for jj in range(j1-1, min(j2, np)):
                    line += f" {charstoi[jj]:6s} {pg[jj]:7s} {signc[jj]}"
                print(line, file=lout)
                if wrtref:
                    refline_out = f"{rg[0]:7s} {c1} {rg[1]:7s} {c2} {rg[2]:7s}  =>"
                    for jj in range(j1-1, min(j2, np)):
                        refline_out += f" {charstoi[jj]:6s} {pg[jj]:7s} {signc[jj]}"
                    print(refline_out, file=refu)

            elif k == nsplit:
                line = " " * 32
                for jj in range(j1-1, min(j2, np)):
                    line += f" {charstoi[jj]:6s} {pg[jj]:7s} {signc[jj]}"
                line += f"    {rc1:10.3e} {arrh[1]:4.1f} {arrh[2]:7.0f}"
                print(line, file=lout)
                if wrtref:
                    refline_out = " " * 32
                    for jj in range(j1-1, min(j2, np)):
                        refline_out += f" {charstoi[jj]:6s} {pg[jj]:7s} {signc[jj]}"
                    refline_out += f"    {rc1:10.3e} {arrh[1]:4.1f} {arrh[2]:7.0f}  {refline}"
                    print(refline_out, file=refu)

            else:
                line = " " * 32
                for jj in range(j1-1, min(j2, np)):
                    line += f" {charstoi[jj]:6s} {pg[jj]:7s} {signc[jj]}"
                print(line, file=lout)
                if wrtref:
                    refline_out = " " * 32
                    for jj in range(j1-1, min(j2, np)):
                        refline_out += f" {charstoi[jj]:6s} {pg[jj]:7s} {signc[jj]}"
                    print(refline_out, file=refu)

    if idreac > 0:
        if idreac == 1:
            print(f"  HV / {nlab:5d}  {xlab:6.3f} /", file=lout)
            if wrtref:
                print(f"  HV / {nlab:5d}  {xlab:6.3f} /", file=refu)
        elif idreac == 3:
            line = f"  FALLOFF /{folow[0]:9.3e}{folow[1]:5.1f}{folow[2]:7.0f}"
            for i in range(4):
                line += f" {fotroe[i]:6.1f}"
            line += "/"
            print(line, file=lout)
            if wrtref:
                refline_out = f"  FALLOFF /{folow[0]:9.2e}{folow[1]:5.1f}{folow[2]:7.0f}"
                for i in range(4):
                    refline_out += f" {fotroe[i]:6.1f}"
                refline_out += "/"
                print(refline_out, file=refu)
        elif idreac == 2:
            print(f"  EXTRA / {nlab:5d} /", file=lout)
            if wrtref:
                print(f"  EXTRA / {nlab:5d} /", file=refu)
        else:
            print('--error--, in rxwrit. Idreac is > 3, in reaction :')
            print(rg[0], c1, rg[1], c2, rg[2], ' =>')
            raise Exception("in rxwrit")


def rxwrit_dyn(lout, r, s, p, arrh, idreac, auxinfo, charfrom, charto):
    from keyparameter import mxnp

    mnp = len(s)

    if r[0] != p[0]:
        print('--error--, in rxwrit_dyn. distinct reactant and product:')
        print(r[0], '=>', p[0])
        raise Exception("in rxwrit_dyn")

    if s[0] != 1.0:
        print('--error--, in rxwrit_dyn. Expect stoe. coef. equal 1')
        print(s[0])
        raise Exception("in rxwrit_dyn")

    if p[0][:2] == '  ':
        print('--error--, in rxwrit_dyn. No product found')
        print(p[0])
        raise Exception("in rxwrit_dyn")

    for i in range(1, mnp):
        if p[i][:2] != ' ':
            print('--error--, in rxwrit_dyn. More than one product found')
            print(p[i])
            raise Exception("in rxwrit_dyn")

    if idreac < 1 or idreac > 5:
        print('--error--, in rxwrit_dyn. Idreac out of bound')
        print(p[i])
        raise Exception("in rxwrit_dyn")

    locheck = False
    if r[1][:4] == 'AIN ':
        locheck = True
    if r[1][:4] == 'AOU ':
        locheck = True
    if r[1][:4] == 'WIN ':
        locheck = True
    if r[1][:4] == 'WOU ':
        locheck = True
    if not locheck:
        print('--error--, in rxwrit_dyn. Unexpected keyword :')
        print(r[1])
        raise Exception("in rxwrit_dyn")

    charstoi = [' '] * len(s)
    signc = [' '] * len(s)
    pg = [' '] * len(s)
    rg = [' '] * len(r)

    rg[0] = charfrom + r[0]
    rg[1] = r[1]
    pg[0] = charto + p[0]

    j1 = 1
    j2 = 4
    c1 = '+'
    c2 = ' '

    count4rxn(20 + idreac)
    line = f"{rg[0]:7s}{c1}{rg[1]:7s}{c2}{rg[2]:7s} =>"
    for jj in range(j1-1, j2):
        if jj < 1:  # 只有1个产物
            line += f" {charstoi[jj]:5s} {pg[jj]:7s} {signc[jj]}"
        else:
            line += " " * 16
    line += f"    {arrh[0]:10.3e} {arrh[1]:4.1f} {arrh[2]:7.0f}"
    print(line, file=lout)

    if idreac == 4:
        print(f'  WOU/ {auxinfo[0]:10.2e}/', file=lout)
    elif idreac == 3:
        print(f'  WIN/ {auxinfo[0]:10.2e}/', file=lout)


def count4rxn(idrx):
    global nrx_all, nrx_n, nrx_hv, nrx_extra, nrx_fo, nrx_tb, nrx_o2
    global nrx_meo2, nrx_ro2, nrx_isom, nrx_tabcf, nrx_ain, nrx_aou
    global nrx_win, nrx_wou

    nrx_all = nrx_all + 1
    if idrx == 0:
        nrx_n = nrx_n + 1
    elif idrx == 1:
        nrx_hv = nrx_hv + 1
    elif idrx == 2:
        nrx_extra = nrx_extra + 1
    elif idrx == 3:
        nrx_fo = nrx_fo + 1
    elif idrx == 4:
        nrx_tb = nrx_tb + 1
    elif idrx == 5:
        nrx_o2 = nrx_o2 + 1
    elif idrx == 6:
        nrx_ro2 = nrx_ro2 + 1
    elif idrx == 7:
        nrx_meo2 = nrx_meo2 + 1
    elif idrx == 8:
        nrx_isom = nrx_isom + 1
    elif idrx == 9:
        nrx_tabcf = nrx_tabcf + 1
    elif idrx == 21:
        nrx_ain = nrx_ain + 1
    elif idrx == 22:
        nrx_aou = nrx_aou + 1
    elif idrx == 23:
        nrx_win = nrx_win + 1
    elif idrx == 24:
        nrx_wou = nrx_wou + 1
    else:
        raise Exception("reaction ID not identified")