import numpy as np
import pyscal.core as pc
import pyscal.crystal_structures as pcs
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib import colors as mcolors
from matplotlib.patches import Polygon
import ipywidgets as widgets
import random


def get_plotvectors(atom):
    """
    When an atom object is passed, find the plot vectors
    """
    vecs = []
    allvertexnos = []
    unq_vnos = np.unique(atom.vertex_numbers)
    for vno in unq_vnos:
        #print(vno)
        ipos = atom.vertex_vectors[vno*3:vno*3+3]
        #print(len(ipos))
        vecs.append(ipos)
    modvecs = []
    st = 1
    for vno in atom.face_vertices:
        vphase = atom.vertex_numbers[st:st+vno]
        #print(vphase)
        dummymodvecs = []
        for v in vphase:
            dummymodvecs.append(vecs[v])
        modvecs.append(dummymodvecs)
        st += (vno+1)
    return modvecs

def plot_3d(modveclist):
    """
    Make a 3d plot
    """
    fig = plt.figure(figsize=(18, 18))
    ax = fig.add_subplot(1,1,1,projection='3d')

    colors = ['#F57F17', '#455A64']

    for count, modvec in enumerate(modveclist):
        for pp in modvec:
            ax.scatter(np.array(pp)[:,0], np.array(pp)[:,1], np.array(pp)[:,2], color=colors[count], s=40)
            ax.plot(np.array(pp)[:,0], np.array(pp)[:,1], np.array(pp)[:,2], color=colors[count], linewidth=2, alpha=0.3)

    fig.canvas.layout.width = '400px'
    fig.canvas.layout.height = '400px'
    fig.canvas.layout.margin = '0px'
    #axis fills
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False

    # Now set color to white (or whatever is "invisible")
    ax.xaxis.pane.set_edgecolor('w')
    ax.yaxis.pane.set_edgecolor('w')
    ax.zaxis.pane.set_edgecolor('w')

    # Bonus: To get rid of the grid as well:
    ax.grid(False)
    ax.set_axis_off()
    #get rid of axis


def create_dropdown():
    w1 = widgets.Dropdown(
        options=['liquid', 'bcc', 'fcc', 'distorted-bcc', 'distorted-fcc'],
        value='bcc',
        description='Structure:',
        disabled=False,
    )
    return w1


def return_possible_polys(struct):
    if struct=='bcc':
        atoms, box = pcs.make_crystal('bcc', lattice_constant=4, repetitions=[3,3,3])
        sys = pc.System()
        sys.atoms = atoms
        sys.box = box
        sys.find_neighbors(method='voronoi')
        sys.calculate_vorovector()
        voros = [" ".join(np.array(atom.vorovector).astype(str)) for atom in sys.atoms]
        #take unique
        unique_voros = np.unique(voros)
    elif struct=='fcc':
        atoms, box = pcs.make_crystal('fcc', lattice_constant=4, repetitions=[3,3,3])
        sys = pc.System()
        sys.atoms = atoms
        sys.box = box
        sys.find_neighbors(method='voronoi')
        sys.calculate_vorovector()
        voros = [" ".join(np.array(atom.vorovector).astype(str)) for atom in sys.atoms]
        #take unique
        unique_voros = np.unique(voros)
    elif struct=='liquid':
        sys = pc.System()
        sys.read_inputfile('datafiles/lqd.dat')
        sys.find_neighbors(method='voronoi')
        sys.calculate_vorovector()
        voros = [" ".join(np.array(atom.vorovector).astype(str)) for atom in sys.atoms]
        #take unique
        unique_voros = np.unique(voros)
    elif struct=='distorted-bcc':
        sys = pc.System()
        sys.read_inputfile('datafiles/mixed_bcc.dat')
        sys.find_neighbors(method='voronoi')
        sys.calculate_vorovector(edge_cutoff=0.0, area_cutoff=0.0)
        voros = [" ".join(np.array(atom.vorovector).astype(str)) for atom in sys.atoms if atom.vorovector[1] >= 2]
        #take unique
        unique_voros = np.unique(voros)

    elif struct=='distorted-fcc':
        sys = pc.System()
        sys.read_inputfile('datafiles/mixed_fcc.dat')
        sys.find_neighbors(method='voronoi')
        sys.calculate_vorovector()
        voros = [" ".join(np.array(atom.vorovector).astype(str)) for atom in sys.atoms]
        #take unique
        unique_voros = np.unique(voros)

    #now create a dropdown widget with this option
    poly_widget = widgets.Dropdown(
        options=unique_voros,
        value=unique_voros[0],
        description='polyhedra:',
        disabled=False,
    )
    return widgets.interact(visualise_poly, poly=poly_widget, atoms=widgets.fixed(sys.atoms))

def visualise_poly(poly, atoms):
    refined_atoms = [atom for atom in atoms if " ".join(np.array(atom.vorovector).astype(str)) == poly]
    satom = random.choice(refined_atoms)
    pvecs = get_plotvectors(satom)
    plot_3d([pvecs])

def create_double_widgets():
    w1 = create_dropdown()
    w1 = widgets.interactive(return_possible_polys, struct=w1)
    w2 = create_dropdown()
    w2 = widgets.interactive(return_possible_polys, struct=w2)
    return widgets.HBox([w1, w2])
