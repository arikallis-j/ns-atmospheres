import matplotlib.pyplot as plt
import numpy as np
fig, ax = plt.subplots(figsize=(10, 10), layout='constrained') 

def init_graph():
    fig, ax = plt.subplots(figsize=(10, 10), layout='constrained') 

def show_graph():
    plt.show() 

def draw_graph(X, Y, color='green', title="Simple Plot", label='test', xlabel = 'x label', ylabel='y label', xrange=None, yrange=None, loglog=False):
    if loglog:
        ax.loglog(X, Y, label=label, color=color, linewidth=1, linestyle='-')   
    else:
        ax.plot(X, Y, label=label, color=color, linewidth=1, linestyle='-')   
    if xrange != None:
        left, right = xrange
        ax.set_xlim(left, right)   
    if yrange != None:
        left, right = yrange
        ax.set_ylim(left, right)   
    
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.get_legend()