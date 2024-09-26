import matplotlib.pyplot as plt
import os

from rich.table import Table
from rich.console import Console

class GraphPalette:
    def __init__(self, colors: list):
        self.colors = colors
        self.N = len(colors)

    def __repr__(self): # take str(Palette)
        str_color = ""
        for k in range(self.N):
            str_color += str(self.colors[k])
            if k!= self.N - 1:
                str_color += " | "
        return str_color
    
    def __getitem__(self, index): # take Palette[k]
        return self.colors[index % self.N]
    
    def __iter__(self): # take list(Palette)
        return iter(self.colors)
    
class GraphData:
    def __init__(self, X, Y, tex_mode=False, ascii_label = "line", latex_label = "$line$"):
        self.X = X
        self.Y = Y

        self.func_ascii = ascii_label
        if tex_mode:
            self.func_latex = latex_label
        else:
            self.func_latex = ascii_label

class GraphDataBase:
    def __init__(self, filename = "graph", tex_mode=False, ascii_labels = ("Title", "X", "Y"), latex_labels = ("$Title$", "$X$", "$Y$")):
        self.DataBase = []

        self.filename = filename
        self.title_ascii= ascii_labels[0]
        self.x_ascii = ascii_labels[1]
        self.y_ascii = ascii_labels[2]

        if tex_mode: 
            self.title_latex = latex_labels[0]
            self.x_latex = latex_labels[1]
            self.y_latex = latex_labels[2]
        else:
            self.title_latex = ascii_labels[0]
            self.x_latex = ascii_labels[1]
            self.y_latex = ascii_labels[2]
        
    def __len__(self):
        return len(self.DataBase)
        
    def __getitem__(self, index): # take DataBase[k]
        return self.DataBase[index]
    
    def __iter__(self): # take list(DataBase)
        return iter(self.DataBase)

class GraphSize:
    def __init__(self, figsize = (10, 10), x_range = None, y_range = None):
        self.figsize = figsize
        self.x_range = x_range
        self.y_range = y_range

class GraphStyle:
    def __init__(self, palette, linestyle='-', marker=None, legend=False, loglog=False, linewidth=1, markersize=1):
        self.palette = palette
        self.linestyle = linestyle
        self.marker = marker
        self.legend = legend
        self.loglog = loglog
        self.linewidth = linewidth
        self.markersize = markersize

class Graph:
    def __init__(self, DataBase, Size, Style):
        self.DataBase = DataBase
        self.N_graph = len(DataBase)
        self.Size = Size
        self.Style = Style
    def draw_picture(self, draw=True, save=False):
        plt.style.use('seaborn-v0_8-whitegrid')
        fig, ax = plt.subplots(figsize=self.Size.figsize, layout='constrained') 

        if self.Size.x_range != None:
            a, b = self.Size.x_range
            ax.set_xlim(a, b)
        if self.Size.y_range != None:
            a, b = self.Size.y_range
            ax.set_ylim(a, b)
        
        ax.set_title(self.DataBase.title_latex, loc='center', fontsize=20)
        ax.set_xlabel(self.DataBase.x_latex)
        ax.set_ylabel(self.DataBase.y_latex)
        
        for k in range(self.N_graph):
            Data = self.DataBase[k]
            if self.Style.loglog:
                ax.loglog(Data.X, Data.Y, color=self.Style.palette[k],
                                          linestyle=self.Style.linestyle,
                                          marker=self.Style.marker, 
                                          linewidth=self.Style.linewidth, 
                                          markersize=self.Style.markersize,
                                          label=Data.func_latex)
            else:
                ax.plot(Data.X, Data.Y, color=self.Style.palette[k],
                                        linestyle=self.Style.linestyle,
                                        marker=self.Style.marker, 
                                        linewidth=self.Style.linewidth, 
                                        markersize=self.Style.markersize,
                                        label=Data.func_latex)
        if self.Style.legend:
            ax.legend()

        if save:
            if not os.path.isdir('graph'):
                os.mkdir('graph')
            plt.savefig('graph/' + self.DataBase.filename + '.pdf')
        
        if draw:
            plt.show()

    def draw_table(self, save=False):
        table = Table(title = self.DataBase.title_ascii, title_style="bold cyan", header_style="cyan")
        X = self.DataBase[0].X
        Y_arr = []
        for n in range(self.N_graph):
            Y = self.DataBase[n].Y
            Y_arr.append(Y)

        table.add_column("N", justify = "right", style="cyan")
        table.add_column(self.DataBase.x_ascii, justify = "center", style=self.Style.palette[self.N_graph])
        
        for n in range(self.N_graph):
            Data = self.DataBase[n]
            table.add_column(f"{self.DataBase.y_ascii} | {Data.func_ascii}", justify = "center", style=self.Style.palette[n])

        for k in range(len(X)):
            line = [str(k+1), str(X[k])]
            for n in range(self.N_graph):
                line.append(str(Y_arr[n][k]))
            table.add_row(*line)
        
        console = Console()
        console.print('')
        console.print(table)

