def tryfloat(x):
    try:
        return float(x)
    except ValueError:
        return float(0)

def read_radex(file,flow,fupp,bw=0.01,debug=False):
    """ 
    less hack-ey way to read radex.out files
    """
    linenum = 0
    line = file.readline()
    linenum+=1
    if debug: print line
    if line == '':
        return 0
    words = line.split()
    if words[1] == '--':
        freq = tryfloat(words[4])
    else:
        freq = 0
    while not(freq*(1-bw) < flow < freq*(1+bw)):
        if words[1] == 'T(kin)':
            tkin = tryfloat(words[3])
        elif line.find("Density of H2") != -1:
            dens = tryfloat(words[5])
        elif line.find("Density of pH2") != -1:
            pdens = tryfloat(words[5])
        elif line.find("Density of oH2") != -1:
            odens = tryfloat(words[5])
        elif line.find("Column density") != -1:
            col = tryfloat(words[4])
        line = file.readline(); linenum+=1
        words = line.split()
        if words[1] == '--':
            freq = tryfloat(words[4])
    TexLow   = tryfloat(words[6])
    TauLow   = tryfloat(words[7])
    TrotLow  = tryfloat(words[8])
    FluxLow  = tryfloat(words[11])
    line = file.readline(); linenum+=1
    words = line.split()
    if words[1] == '--':
        freq = tryfloat(words[4])
    while  not(freq*(1-bw) < fupp < freq*(1+bw)):
        line = file.readline(); linenum+=1
        if debug: print freq,flow,line
        words = line.split()
        if words[1] == '--':
            freq = tryfloat(words[4])
    TexUpp   = tryfloat(words[6])
    TauUpp   = tryfloat(words[7])
    TrotUpp  = tryfloat(words[8])
    FluxUpp  = tryfloat(words[11])
    while len(words) > 0 and words[1] == '--':
        line = file.readline(); linenum+=1
        words = line.split()
    return tkin,dens,col,TexLow,TexUpp,TauLow,TauUpp,TrotLow,TrotUpp,FluxLow,FluxUpp
