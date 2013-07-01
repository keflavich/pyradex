"""
This code is exclusively for use of testing read_radex.  Don't use it in any
production anything until you're sure it works.
"""

bw    = 0.01   # "bandwidth": free spectral range around line (used to say which lines get printer)

def read_radex_old(outfile,flow,fupp,bw=bw):
    """
    A hack-ey means of reading a radex.out file.  
    Cycles through it based on fixed format
    """
    line  = outfile.readline()
    words = line.split()
    while (words[-1] != "FLUX"):
        line  = outfile.readline()
        words = line.split()
        if words[1] == "T(kin)":
            temp  = float(words[-1])
        if words[1] == "Density" and words[3] == "H2":
            dens  = float(words[-1])
        if words[1] == "Column":
            col  = float(words[-1])
    line  = outfile.readline()
    line  = outfile.readline()
    words = line.split()
    ftmp  = float(words[4])
    while ((ftmp < flow*(1-bw)) or (ftmp > flow/(1-bw))):
        line  = outfile.readline()
        words = line.split()
        ftmp  = float(words[4])
    low   = float(words[-2])
    TexLow   = float(words[6])
    TauLow   = float(words[7])
    TrotLow  = float(words[8])
    FluxLow  = float(words[11])
    line  = outfile.readline()
    words = line.split()
    ftmp  = float(words[4])
    while ((ftmp < fupp*(1-bw)) or (ftmp > fupp/(1-bw))):
        line  = outfile.readline()
        words = line.split()
        ftmp  = float(words[4])
    upp   = float(words[-2])
    TexUpp   = float(words[6])
    TauUpp   = float(words[7])
    TrotUpp  = float(words[8])
    FluxUpp  = float(words[11])
    return temp,dens,col,TexLow,TexUpp,TauLow,TauUpp,TrotLow,TrotUpp,FluxLow,FluxUpp

def tryfloat(x):
    try:
        return float(x)
    except ValueError:
        return float(0)

def read_radex(file,flow,fupp,bw=bw,debug=False):
    """ 
    less hack-ey way to read radex.out files
    """
    linenum = 0
    line = file.readline(); linenum+=1
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
