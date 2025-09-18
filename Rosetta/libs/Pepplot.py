import re
import numpy as np # Importa numpy, que é uma biblioteca voltada a cálculos matemáticos
import seaborn as sns  # Importa o seaborn, outra biblioteca para visualização de dados
import pandas as pd # Importa pandas, que é uma biblioteca para trabalhar com dados tabelados
import matplotlib.pyplot as plt # Importa o matplotlib.pyplot, que é uma biblioteca para visualização de dados

# Esta função vai ordenar os textos corretamente como um humano ordenaria: 0, 1, 2, 3, 4, ..., 10, 11, 12, 13... e não 0, 1, 10, 11,12,...., 2, 3, 4, 5...
def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)

# PepPlot é aquele plot que fizemos com as "colunas de peptídeos"
def PepPlot(data, ligcolumn, reccolumn, hue, palette, ax=None):
    receptors = list(natural_sort(pd.unique(data[reccolumn])))
    ligands = list(natural_sort(pd.unique(data[ligcolumn])))
    huetypes = list(pd.unique(data[hue]))
    figwidth = plt.gcf().get_figwidth()
    figheight = plt.gcf().get_figheight()
    dpi = plt.gcf().get_dpi()
    maxstack = data.groupby([reccolumn,ligcolumn]).count().groupby([reccolumn]).count().max()[0]
    maxwidth = figwidth*dpi/(3*len(receptors)*0.7)
    maxheight = figheight*dpi/maxstack
    if maxwidth < maxheight:
        maxsize = maxwidth
        top = maxwidth/maxheight
    else:
        top = 1
    maxsize = min([maxwidth,maxheight])
    ticks = np.linspace(0,1,len(receptors))
    yticks = np.linspace(0,top,maxstack)
    if ax==None:
        fig,ax = plt.subplots(figsize=(figwidth, figheight))
        ax.set_xlim(-0.2,1.2)
        ax.set_ylim(-0.2,1.2)
        ax.set_xticks(ticks,rotation=90)
        ax.set_xticklabels(receptors,rotation=90)
    else:
        ax.set_xticks(ticks)
        ax.set_xlim(-0.01,1.03)
        ax.set_ylim(-0.005,1.03)
        ax.set_xticklabels(receptors,rotation=90)
        i = 0
        for rec in receptors:
            j = 0
            used = []
            for lig in reversed(natural_sort(list(data.loc[data[reccolumn]==rec,ligcolumn]))):
                if lig not in used:
                    h = data.loc[(data[reccolumn]==rec)&(data[ligcolumn]==lig),hue].values[0]
                    k = huetypes.index(h)
                    c = palette[k]
                    ax.text(ticks[i]-0.01,yticks[j], lig.replace('_',''), fontsize=maxsize*0.9, color=c, fontname='monospace', fontweight='bold')
                    j += 1
                    used.append(lig)
            i += 1
            
        for i in range(len(huetypes)):
            ax.text(0.30*i+0.13, 1.05, huetypes[i],fontsize=maxsize*0.9, color=palette[i])