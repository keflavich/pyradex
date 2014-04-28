import numpy as np
import pyradex
from astropy.utils.console import ProgressBar

def grid_wrapper(molecule,
                 temperatures=[],
                 densities=[],
                 columns=[],
                 h2columns=[],
                 abundances=[],
                 transition_indices=[],
                 orthopararatios=[],
                 observable_parameters=['tex','source_line_surfbrightness','tau'],
                 ):


    ntemp = len(temperatures)
    ndens = len(densities)
    coltype = 'h2' if h2columns else 'mol'
    columns = columns or h2columns
    ncols = len(columns) or len(h2columns)
    nabund = len(abundances)
    nopr = len(orthopararatios)

    grids = {tid:
             {par: np.empty([nopr, ncols, nabund, ntemp, ndens])
              for par in observable_parameters}
             for tid in transition_indices}

    # Just a quick first run to get things initialized
    if coltype == 'mol':
        R = pyradex.Radex(species=molecule, column=columns[0], abundance=abundances[0])
    else:
        R = pyradex.Radex(species=molecule, h2column=h2columns[0], abundance=abundances[0])
    R.run_radex()

    # get the table so we can look at the frequency grid
    # table = R.get_table()

    # Target frequencies:
    # frequencies = table[np.array(transition_indices)]

    pb = ProgressBar(ntemp*ndens*ncols*nabund*nopr)

    for opri, opr in enumerate(orthopararatios):
        fortho = opr/(1+opr)
        for coli, col in enumerate(columns):
            for abundi, abund in enumerate(abundances):
                for temi,tem in enumerate(temperatures):
                    R.temperature = tem
                    for densi,dens in enumerate(densities):
                        R.density = {'oH2':dens*fortho,'pH2':dens*(1-fortho)}
                        if coltype == 'mol':
                            R.column = col
                        else:
                            R.h2column = col
                        R.abundance = abund # reset column to the appropriate value
                        R.run_radex(reuse_last=False, reload_molfile=True)

                        # type is important: MUST be tuple!
                        ind = (opri, coli, abundi, temi, densi)

                        pb.update()

                        for tid in transition_indices:
                            for par in observable_parameters:
                                val = getattr(R, par)[tid]
                                if hasattr(val,'value'):
                                    val = val.value
                                grids[tid][par][ind] = val

    return grids
