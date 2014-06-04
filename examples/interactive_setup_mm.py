styleargs = {'linewidth': 2, 'alpha': 0.5, 'color':'#5A228B'}

def setup(tem=5,dens=5,taugrid=taugrid_303,texgrid=texgrid_303):
    fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2, sharex='col', sharey='row', squeeze=True, figsize=(10,7))
    plt.subplots_adjust(hspace=0,wspace=0)
    lines1, = ax1.plot(temperatures, taugrid[dens,:], **styleargs)
    ax1.set_ylim(-0.2,0.2)
    p1, = ax1.plot(temperatures[tem],taugrid[dens,tem], 'o',alpha=0.5, markeredgecolor='none')
    lines2, = ax3.plot(temperatures, texgrid[dens,:], **styleargs)
    ax3.set_ylim(0,20)
    p3, =ax3.plot(temperatures[tem],texgrid[dens,tem], 'o',alpha=0.5, markeredgecolor='none')

    lines3, = ax2.semilogx(densities, taugrid[:,tem], **styleargs)
    ax2.set_ylim(-0.2,0.2)
    p2, = ax2.plot(densities[dens],taugrid[dens,tem], 'o',alpha=0.5, markeredgecolor='none')
    lines4, = ax4.semilogx(densities, texgrid[:,tem], **styleargs)
    p4, = ax4.plot(densities[dens],texgrid[dens,tem], 'o',alpha=0.5, markeredgecolor='none')
    ax4.set_ylim(0,20)
    plt.suptitle("$T=%i$ K, $n=10^{%0.1f}$ cm$^{-3}$" % (temperatures[tem],np.log10(densities[dens])),
                 fontsize=20)
    ax4.set_xlabel('$n(H_2)$',fontsize=20)
    ax3.set_xlabel("T",fontsize=20)
    ax1.set_ylabel(r"$\tau$",fontsize=20)
    ax3.set_ylabel("$T_{ex}$",fontsize=20)
    plt.show()
    return fig,lines1,lines2,lines3,lines4,p1,p2,p3,p4

def run_plot_temden(taugrid=taugrid_303, texgrid=texgrid_303):
    fig,lines1,lines2,lines3,lines4,p1,p2,p3,p4 = setup(taugrid=taugrid,texgrid=texgrid)
    
    @interact(tem=(0,temperatures.size-1),dens=(0,densities.size-1))
    def plot_temden(tem,dens):
        lines1.set_data(temperatures, taugrid[dens,:])
        lines2.set_data(temperatures, texgrid[dens,:])
        lines3.set_data(densities, taugrid[:,tem])
        lines4.set_data(densities, texgrid[:,tem])
        p1.set_data(temperatures[tem],taugrid[dens,tem])
        p2.set_data(densities[dens],taugrid[dens,tem])
        p3.set_data(temperatures[tem],texgrid[dens,tem])
        p4.set_data(densities[dens],texgrid[dens,tem])

        display(fig)
