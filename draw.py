import networkx as nx
from matplotlib import pyplot as plt
import numpy as np

def pos_from_helix_classes(G):
    pos = {}
    for helix in G.nodes():

        pos[helix] = (helix[0] + .5*helix[2],
                    helix[1] - .5*helix[2])

    return pos

from matplotlib.lines import Line2D
class LineDataUnits(Line2D):
    def __init__(self, *args, **kwargs):
        _lw_data = kwargs.pop("linewidth", 1)
        super().__init__(*args, **kwargs)
        self._lw_data = _lw_data

    def _get_lw(self):
        if self.axes is not None:
            ppd = 72./self.axes.figure.dpi
            trans = self.axes.transData.transform
            return ((trans((1, self._lw_data))-trans((0, 0)))*ppd)[1]
        else:
            return 1

    def _set_lw(self, lw):
        self._lw_data = lw

    _linewidth = property(_get_lw, _set_lw)

import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.patheffects as PathEffects
def plot_helix_class_diagram(
        ax,
        graph,
        selected_helix_classes = None,
        sequence_length = None,
        component_labels = None,
        filename=None,
        edge_color_key=None):
    if selected_helix_classes is None:
        selected_helix_classes = ()

    for helix in graph.nodes():
        i,j,k = helix
        if helix in selected_helix_classes:
            color = "lightgreen"
        else:
            color = "lightsteelblue"
        helix_line = LineDataUnits(
            [i, i + k],
            [j, j - k],
            linewidth=.7,
            solid_capstyle="round",
            color=color,
            zorder=-10)
        ax.add_line(helix_line)

    node_positions = pos_from_helix_classes(graph)

    edges = graph.edges()

    ax.grid(True)

    if edge_color_key is not None:
        edge_color_dict = nx.get_edge_attributes(graph, edge_color_key)
        edge_colors = [edge_color_dict[edge] for edge in edges]

        mcl = nx.draw_networkx_edges(
            graph,
            node_positions,
            edgelist=edges,
            width=1.,
            alpha=0.6,
            ax=ax,
            edge_color=edge_colors,
            arrows=False)

    else:
        mcl = nx.draw_networkx_edges(
            graph,
            node_positions,
            edgelist=edges,
            width=1.,
            alpha=0.6,
            ax=ax,
            edge_color="black",
            arrows=False)

    if component_labels is not None:
        for component in nx.connected_components(graph):
            max_i = 0
            min_j = 10000000000
            for helix in component:
                max_i = max(helix[0] + helix[2] - 1, max_i)
                min_j = min(helix[1] - helix[2] + 1, min_j)

            txt = ax.text(max_i + 0.5, min_j + 0.5, 
                    component_labels[tuple(sorted(component))],
                    color='darkred')

            txt.set_path_effects([PathEffects.withStroke(linewidth=1., foreground='w')])

    if sequence_length is None:
        ax.autoscale(True)
    else:
        ax.set_xlim(0, sequence_length + 1)
        ax.set_ylim(0, sequence_length + 1)

    if filename is None:
        plt.show()
    else:
        plt.savefig(filename,dpi=300)

    plt.close()

def plot_radial_diagram(
        helix_structure,
        sequence,
        filename = None,
        label = None):
    import RNA
    import data

    dot_bracket_string, skipped_pairs = data.To_Dot_Bracket(
        helix_structure, len(sequence))

    #print(helix_structure)
    #print(dot_bracket_string)

    fig, ax = plt.subplots(figsize=(8.,6.))

    coords = RNA.get_xy_coordinates(dot_bracket_string)

    coords = np.array(
        [(coords.get(idx).X, coords.get(idx).Y) 
        for idx in range(len(dot_bracket_string))])

    ax.plot(coords[:,0], coords[:,1], zorder=0, color="black")
    ax.scatter(coords[:,0], coords[:,1], zorder=1, color="orange", s=50, edgecolor="black")

    basepair_list = []
    for helix in helix_structure:
        i,j,k = helix

        for idx in range(k):
            pair = (i + idx, j - idx)
            if pair not in skipped_pairs:
                basepair_list.append(pair)

    bp_coords = np.array([
        [coords[pair[0]], coords[pair[1]]] for pair in basepair_list])

    ax.plot(bp_coords[:,:,0].T, bp_coords[:,:,1].T, zorder=0, color="black", linewidth=3)

    ax.set_aspect(1)
    ax.tick_params(left=False, labelleft=False, bottom=False, labelbottom=False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["bottom"].set_visible(False)

    if label is not None:
        txt = ax.text(0.9,0.9, label, transform = ax.transAxes, size=30)
        txt.set_path_effects([PathEffects.withStroke(linewidth=4.,foreground='w')])

    fig.tight_layout(pad=0.05)
    
    if filename is None:
        plt.show()
    else:
        plt.savefig(filename,dpi=300)

    plt.close()



def plot_clustered_stems(
        ax,
        graph,
        stem_clusters,
        sequence_length = None,
        component_labels = None,
        filename = None):

    helix_clusters = []
    for cluster in stem_clusters:
        helix_cluster = []
        for stem in cluster:
            helix_cluster += stem
        helix_clusters.append(helix_cluster)

    plot_clustered_helix_classes(
        ax,
        graph,
        helix_clusters,
        sequence_length,
        component_labels,
        filename)
    

def plot_clustered_helix_classes(
        ax,
        graph,
        helix_clusters,
        sequence_length = None,
        component_labels = None,
        filename = None):

    import numpy as np
    helix_class_cluster_dict = {}
    for idx, cluster in enumerate(helix_clusters):
        for helix_class in cluster:
            helix_class_cluster_dict[helix_class] = idx

    colormap = plt.cm.tab10
    for helix in graph.nodes():
        i,j,k = helix
        color = colormap(helix_class_cluster_dict[helix])
        helix_line = LineDataUnits(
            [i, i + k],
            [j, j - k],
            linewidth=.7,
            solid_capstyle="round",
            color=color,
            zorder=-10)
        ax.add_line(helix_line)

    node_positions = pos_from_helix_classes(graph)

    edges = graph.edges()
    cmap = plt.cm.get_cmap('gist_heat')

    ax.grid(True)

    mcl = nx.draw_networkx_edges(
        graph,
        node_positions,
        edgelist=edges,
        width=1.,
        alpha=0.6,
        ax=ax,
        arrows=False)

    if component_labels is not None:
        for component in nx.connected_components(graph):
            max_i = 0
            min_j = 10000000000
            for helix in component:
                max_i = max(helix[0] + helix[2] - 1, max_i)
                min_j = min(helix[1] - helix[2] + 1, min_j)

            txt = ax.text(max_i + 0.5, min_j + 0.5, 
                    component_labels[tuple(sorted(component))],
                    color='darkred')

            txt.set_path_effects([PathEffects.withStroke(linewidth=1., foreground='w')])

    if sequence_length is None:
        ax.autoscale(True)
    else:
        ax.set_xlim(0, sequence_length + 1)
        ax.set_ylim(0, sequence_length + 1)

    if filename is None:
        plt.show()
    else:
        plt.savefig(filename,dpi=300)

    plt.close()

def plot_clustered_dotplot(
        ax,
        clustered_basepairs,
        sequence_length = None,
        noise_index = None,
        filename = None):

    import numpy as np
    colormap = plt.cm.tab20

    for idx, cluster in enumerate(clustered_basepairs):
        x = [basepair[0] for basepair in cluster]
        y = [basepair[1] for basepair in cluster]

        color = np.array([colormap(idx)])
        ax.scatter(x, y, c = color)

    if sequence_length is None:
        ax.autoscale(True)
    else:
        ax.set_xlim(0, sequence_length + 1)
        ax.set_ylim(0, sequence_length + 1)

    if filename is None:
        plt.show()
    else:
        plt.savefig(filename,dpi=300)

    plt.close()

def generate_arc_diagram(
        filename,
        sequence_length,
        helix_structure_list,
        important_helices = None):

    print(important_helices)
    seen_helices = set()
    seen_important_helices = set()

    for structure in helix_structure_list:
        for helix in structure:

            if helix in important_helices:
                seen_important_helices.add(helix)
            else:
                seen_helices.add(helix)
    if filename.endswith(".png"):
        filename = filename[:4] + ".txt"

    with open(filename, "w") as f:

        f.write("# {}\n".format(sequence_length))
        f.write("i\tj\tlength\tvalue\n")

        for helix in seen_helices:
            i,j,k = helix
            f.write("{}\t{}\t{}\t0.0\n".format(i,j,k))
        
        for helix in seen_important_helices:
            i,j,k = helix
            f.write("{}\t{}\t{}\t1.0\n".format(i,j,k))
    
        f.write("\n")
    
    import os
    import subprocess
    with open(os.devnull, "w") as FNULL:
        retcode = subprocess.call(
            ['Rscript', 'arc_diagram_emph.R',filename], 
            stdout=FNULL, 
            stderr=subprocess.STDOUT)

def generate_arc_diagram_mpl(
        sequence_length,
        helix_classes,
        important_classes = None,
        max_diameter = None,
        filename = None,
        label_dict = None):

    from matplotlib.patches import Arc
    from matplotlib.collections import PatchCollection
    import matplotlib.patheffects as PathEffects

    patches = []
    if max_diameter is None:
        max_diameter = 1

    for helix_class in helix_classes:

        i,j,k = helix_class

        center = (i+j)/2

        if helix_class in important_classes:
            edge_color = (0.,0.,0.)
        else:
            edge_color = (0.5,0.5,0.5)

        for idx in range(k):
            diameter = (j - i) - (2 * idx)

            patch = Arc(
                xy=(center, 0), 
                width=diameter, 
                height=diameter, 
                theta1=0,
                theta2=180,
                linewidth=1.5,
                edgecolor=edge_color,
                alpha=0.9)

            patches.append(patch)

            if diameter > max_diameter:
                max_diameter = diameter

    figure_ratio = max_diameter / sequence_length * 0.6
    fig, ax = plt.subplots(figsize=(8., 8. * figure_ratio))

    for patch in patches:
        ax.add_patch(patch)

    if label_dict is not None:
        seen_labels = set()

        for helix_class in sorted(helix_classes):
            if label_dict[helix_class] in seen_labels:
                continue
            seen_labels.add(label_dict[helix_class])

            i, j, k = helix_class
            diameter = j - i
            root_2_over_2 = 0.7071

            x = (j + i) / 2 - diameter * root_2_over_2
            y = diameter * root_2_over_2

            txt = ax.text(x, y, 
                label_dict[helix_class],
                color="black")

            txt.set_path_effects([PathEffects.withStroke(linewidth=1.,foreground='w')])

    ax.set_xlim(0,sequence_length)
    ax.set_ylim(0, max_diameter / 2)
    ax.set_aspect(1)
    ax.tick_params(left=False, labelleft=False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["top"].set_visible(False)

    fig.tight_layout(pad=0.05)

    if filename is None:
        plt.show()
    else:
        plt.savefig(filename,dpi=300)

    plt.close()

def make_contiguous(approximate_class_sets):
    return approximate_class_sets

def find_outer_arc(class_set):
    return (6, 80)

def find_inner_arcs(class_set):
    return [(9, 27), (33, 76)]

def shift(helix):
    return (helix[0] + 1, helix[1] + 1, helix[2])
    
from operator import sub
def get_aspect(ax):
    # Total figure size
    figW, figH = ax.get_figure().get_size_inches()
    # Axis size on figure
    _, _, w, h = ax.get_position().bounds
    # Ratio of display units
    disp_ratio = (figH * h) / (figW * w)
    # Ratio of data units
    # Negative over negative because of the order of subtraction
    data_ratio = sub(*ax.get_ylim()) / sub(*ax.get_xlim())

    return disp_ratio / data_ratio

def generate_region_arc_diagram(
        sequence_length,
        helix_classes,
        approximate_class_sets = None,
        important_classes = None,
        negative_classes = None,
        figure_ratio = 0.6,
        filename = None,
        label = None,
        label_dict = None):

    flat_approximate_set = set()
    [flat_approximate_set.update(class_set) for class_set in approximate_class_sets]

    #contiguous_approximate_classes = make_contiguous(approximate_class_sets)
    #outer_arc_list =  [find_outer_arc( class_set) for class_set in contiguous_approximate_classes]
    #inner_arcs_list = [find_inner_arcs(class_set) for class_set in contiguous_approximate_classes]
    #outer_arc_list = []
    #inner_arcs_list = []

    from matplotlib.patches import Arc, PathPatch
    from matplotlib.path import Path
    from matplotlib.collections import PatchCollection

    patches = []
    max_height = 0
    min_height = 0

    top_max_height = figure_ratio * (0.65) * sequence_length
    bottom_max_height = figure_ratio * (0.35) * sequence_length

    top_height = 1
    bottom_height = 1
    for helix_class in helix_classes:
        i,j,k = helix_class
        diameter = (j - i)
        height = diameter * 0.5

        if helix_class in negative_classes:
            if height > bottom_height:
                bottom_height = height
        else:
            if height > top_height:
                top_height = height

    top_ratio = top_max_height / top_height
    bottom_ratio = bottom_max_height / bottom_height

    top_ratio = (1 if top_ratio > 1 else top_ratio)
    bottom_ratio = (1 if bottom_ratio > 1 else bottom_ratio)

    for helix_class in helix_classes:

        i,j,k = helix_class

        center = (i+j)/2

        if helix_class in important_classes:
            edge_color = (0.,0.,0.)
            if helix_class in negative_classes:
                edge_color = (0.5, 0., 0.)
        else:
            edge_color = (0.5,0.5,0.5)

        if helix_class in flat_approximate_set:
            line_style = ":"
        else:
            line_style = "-"

        for idx in range(k):
            diameter = (j - i) - (2 * idx)

            if helix_class in negative_classes:
                theta1, theta2 = 180, 360
                height = diameter * bottom_ratio
            else:
                theta1, theta2 = 0, 180
                height = diameter * top_ratio

            patch = Arc(
                xy=(center, 0), 
                width=diameter, 
                height=height, 
                theta1=theta1,
                theta2=theta2,
                linewidth=2,
                edgecolor=edge_color,
                linestyle=line_style,
                alpha=0.9)

            patches.append(patch)
    
    fig, ax = plt.subplots(figsize=(8., 8. * figure_ratio))

    for patch in patches:
        ax.add_patch(patch)

    ax.set_xlim(0, sequence_length)
    ax.set_ylim(-bottom_max_height-1, top_max_height+1)
    ax.set_aspect(1)

    txt_list = []
    if label_dict is not None:
        seen_labels = set()

        for helix_class in sorted(helix_classes):
            if label_dict[shift(helix_class)] in seen_labels:
                continue
            seen_labels.add(label_dict[shift(helix_class)])

            i, j, k = helix_class
            radius = (j - i) * 0.5
            root_2_over_2 = 0.7071

            x = (j + i) / 2 - radius * root_2_over_2
            y = radius * root_2_over_2

            if helix_class in negative_classes:
                y = (-y * bottom_ratio) - 2
            else:
                y = y * top_ratio

            txt = ax.text(x, y, 
                label_dict[shift(helix_class)],
                color="black",
                size=20)

            txt.set_path_effects([PathEffects.withStroke(linewidth=4.,foreground='w')])
            txt_list.append(txt)

    #ax.set_ylim(0, max_diameter / 2 + 2)

    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["top"].set_visible(False)

    #if (negative_classes is not None 
    #        and len(negative_classes) < len(helix_classes)
    #        and len(negative_classes) > 0):
    #    ax.set_ylim(min_diameter / 2 - 2, max_diameter / 2 + 2)
    #    #ax.spines["bottom"].set_position("center")
    #    ax.set_aspect("auto")
    #elif (negative_classes is not None and 
    #        len(negative_classes) > 0):
    #    ax.set_ylim(min_diameter / 2 - 2, 0)
    #    #ax.spines["bottom"].set_position(("data",0))
    #    ax.spines["top"].set_visible(True)
    #    ax.spines["bottom"].set_visible(False)

    #if get_aspect(ax) > 1:
    #    ax.set_aspect(1)

    ax.tick_params(left=False, labelleft=False)
    ax.spines["bottom"].set_position(("data",0))
    #ax.spines["top"].set_position(("data",0))

    fig.tight_layout(pad=0.05)

    if label is not None:
        txt = ax.text(0.9,0.9, label, transform = ax.transAxes, size=30)
        txt.set_path_effects([PathEffects.withStroke(linewidth=2.,foreground='w')])

    if filename is None:
        plt.show()
    else:
        plt.savefig(filename,dpi=300)

    plt.close()


if __name__ == "__main__":

    helix_list = [(1, 100, 5), (6, 30, 4), (10, 23, 4), (29, 80, 5), (40, 55, 5), (58, 70, 5)]
    sequence_length = 100
    approximate_class_sets = [[(6, 30, 4), (29, 80, 5)]]
    important_classes = []
    
    generate_region_arc_diagram(
        sequence_length,
        helix_list,
        approximate_class_sets,
        important_classes)

