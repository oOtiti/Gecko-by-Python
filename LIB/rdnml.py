def rd_nml():
    from keyparameter import dirgecko, dirout, nmlu
    from keyflag import critvp, brcut, yldcut, rxloss, maxgen, TK, rx_ro_no2, rx_ro2_no2, C_NO2, pvap_sar, g2pfg, g2wfg, isomerfg, highnoxfg, dhffg, chafg, rx_ro2_multiclass, bimolecrx4criegee

    filein = 'gecko.nml'
    ierr = 0

    try:
        with open(filein, 'r') as f:
            lines = f.readlines()
    except:
        print('--error--, problem reading namelist file')
        raise

    current_section = None
    for line in lines:
        line = line.strip()
        if line.startswith('&dir'):
            current_section = 'dir'
        elif line.startswith('&thresholds'):
            current_section = 'thresholds'
        elif line.startswith('&env_cond'):
            current_section = 'env_cond'
        elif line.startswith('&reductions'):
            current_section = 'reductions'
        elif line.startswith('&sar'):
            current_section = 'sar'
        elif line.startswith('&process'):
            current_section = 'process'
        elif line.startswith('/'):
            current_section = None
        elif current_section and '=' in line:
            key, value = line.split('=', 1)
            key = key.strip()
            value = value.strip().rstrip(',')

            if current_section == 'dir':
                if key == 'dirgecko':
                    dirgecko = value.strip("'")
                elif key == 'dirout':
                    dirout = value.strip("'")
            elif current_section == 'thresholds':
                if key == 'critvp':
                    critvp = float(value)
                elif key == 'brcut':
                    brcut = float(value)
                elif key == 'yldcut':
                    yldcut = float(value)
                elif key == 'rxloss':
                    rxloss = float(value)
            elif current_section == 'env_cond':
                if key == 'TK':
                    TK = float(value)
                elif key == 'C_NO2':
                    C_NO2 = float(value)
            elif current_section == 'reductions':
                if key == 'maxgen':
                    maxgen = int(value)
                elif key == 'isomerfg':
                    isomerfg = int(value)
                elif key == 'highnoxfg':
                    highnoxfg = int(value)
                elif key == 'rx_ro2_multiclass':
                    rx_ro2_multiclass = int(value)
            elif current_section == 'sar':
                if key == 'pvap_sar':
                    pvap_sar = int(value)
            elif current_section == 'process':
                if key == 'g2pfg':
                    g2pfg = int(value)
                elif key == 'g2wfg':
                    g2wfg = int(value)
                elif key == 'dhffg':
                    dhffg = int(value)
                elif key == 'chafg':
                    chafg = int(value)
                elif key == 'rx_ro2_no2':
                    rx_ro2_no2 = int(value)
                elif key == 'rx_ro_no2':
                    rx_ro_no2 = int(value)
                elif key == 'bimolecrx4criegee':
                    bimolecrx4criegee = int(value)