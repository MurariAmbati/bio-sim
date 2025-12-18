"""
pathway diagram visualization
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle, Rectangle, Wedge
import numpy as np
from typing import Optional, Dict, Tuple

from ..models.ar_pathway import ArPathwayModel, PathwayState


def plot_pathway_schematic(model: Optional[ArPathwayModel] = None,
                          state: Optional[PathwayState] = None,
                          show_values: bool = True,
                          figsize: tuple = (18, 14),
                          save: Optional[str] = None,
                          show: bool = True) -> plt.Figure:
    """
    create detailed ar signaling pathway schematic
    
    args:
        model: ar pathway model
        state: current state (for showing molecule counts)
        show_values: whether to show numerical values
        figsize: figure size
        save: path to save
        show: whether to display
        
    returns:
        matplotlib figure
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    # color scheme
    colors = {
        'ligand': '#F18F01',
        'receptor': '#2E86AB',
        'nucleus': '#E8E8E8',
        'cytoplasm': '#F5F5F5',
        'dna': '#6A4C93',
        'protein': '#E63946',
        'coregulator': '#52B788'
    }
    
    # --- compartments ---
    
    # cell membrane
    cell_boundary = Rectangle((0, 0), 20, 18, 
                             linewidth=4, edgecolor='black', 
                             facecolor=colors['cytoplasm'], alpha=0.3)
    ax.add_patch(cell_boundary)
    ax.text(1, 17, 'CELL', fontsize=14, fontweight='bold')
    
    # nuclear membrane
    nuclear_boundary = Rectangle((8, 5), 11, 11,
                                linewidth=3, edgecolor='#6A4C93',
                                facecolor=colors['nucleus'], alpha=0.5)
    ax.add_patch(nuclear_boundary)
    ax.text(9, 15.2, 'NUCLEUS', fontsize=13, fontweight='bold', color='#6A4C93')
    
    # --- extracellular ligands ---
    
    # testosterone
    testosterone_circle = Circle((2, 14), 0.6, color=colors['ligand'], 
                                ec='black', linewidth=2, alpha=0.8)
    ax.add_patch(testosterone_circle)
    ax.text(2, 14, 'T', ha='center', va='center', 
           fontsize=12, fontweight='bold', color='white')
    ax.text(2, 12.8, 'Testosterone', ha='center', fontsize=10, fontweight='bold')
    
    # 5-alpha reductase conversion
    arrow_5ar = FancyArrowPatch((2.6, 14), (4.4, 12),
                               arrowstyle='->,head_width=0.3,head_length=0.5',
                               color='black', linewidth=2.5)
    ax.add_patch(arrow_5ar)
    ax.text(3.5, 13.5, '5α-reductase', ha='center', fontsize=8,
           bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.7))
    
    # dht
    dht_circle = Circle((5, 11), 0.6, color='#C73E1D',
                       ec='black', linewidth=2, alpha=0.8)
    ax.add_patch(dht_circle)
    ax.text(5, 11, 'DHT', ha='center', va='center',
           fontsize=11, fontweight='bold', color='white')
    ax.text(5, 9.8, 'Dihydrotestosterone', ha='center', fontsize=9, fontweight='bold')
    
    # --- cytoplasmic ar ---
    
    # free ar with hsp90
    ar_cyto_x, ar_cyto_y = 3.5, 7
    ar_box = FancyBboxPatch((ar_cyto_x-1, ar_cyto_y-0.7), 2, 1.4,
                           boxstyle="round,pad=0.15",
                           facecolor=colors['receptor'],
                           edgecolor='black', linewidth=2, alpha=0.7)
    ax.add_patch(ar_box)
    ax.text(ar_cyto_x, ar_cyto_y+0.3, 'AR', fontsize=13, 
           fontweight='bold', ha='center', color='white')
    ax.text(ar_cyto_x, ar_cyto_y-0.3, 'HSP90', fontsize=9,
           ha='center', color='white', style='italic')
    ax.text(ar_cyto_x, ar_cyto_y-1.5, 'Cytoplasmic AR', 
           ha='center', fontsize=9, fontweight='bold')
    
    if show_values and state:
        ax.text(ar_cyto_x, ar_cyto_y-2, f'n={state.ar_free_cytoplasm:.0f}',
               ha='center', fontsize=8,
               bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # --- ligand binding ---
    
    # testosterone binding
    arrow_t_bind = FancyArrowPatch((2.6, 13.5), (3, 8),
                                  arrowstyle='->,head_width=0.3,head_length=0.5',
                                  color='#F18F01', linewidth=2.5, linestyle='--')
    ax.add_patch(arrow_t_bind)
    
    # dht binding  
    arrow_dht_bind = FancyArrowPatch((4.5, 10.5), (3.8, 8),
                                    arrowstyle='->,head_width=0.3,head_length=0.5',
                                    color='#C73E1D', linewidth=2.5, linestyle='--')
    ax.add_patch(arrow_dht_bind)
    
    # ar:ligand complex
    ar_ligand_x, ar_ligand_y = 3.5, 4
    ar_ligand_box = FancyBboxPatch((ar_ligand_x-1, ar_ligand_y-0.7), 2, 1.4,
                                  boxstyle="round,pad=0.15",
                                  facecolor='#A23B72',
                                  edgecolor='black', linewidth=2, alpha=0.8)
    ax.add_patch(ar_ligand_box)
    ax.text(ar_ligand_x, ar_ligand_y, 'AR:DHT', fontsize=12,
           fontweight='bold', ha='center', color='white')
    ax.text(ar_ligand_x, ar_ligand_y-1.3, 'Activated AR',
           ha='center', fontsize=9, fontweight='bold')
    
    # conformational change arrow
    arrow_conform = FancyArrowPatch((ar_cyto_x, ar_cyto_y-0.8), 
                                   (ar_ligand_x, ar_ligand_y+0.8),
                                   arrowstyle='->,head_width=0.4,head_length=0.6',
                                   color='red', linewidth=3)
    ax.add_patch(arrow_conform)
    ax.text(ar_cyto_x+0.5, 5.5, 'conformational\nchange', ha='center', fontsize=7,
           bbox=dict(boxstyle='round', facecolor='pink', alpha=0.7))
    
    # --- nuclear translocation ---
    
    # nuclear pore
    pore_x, pore_y = 7.5, 9
    pore = Rectangle((pore_x, pore_y-0.5), 0.5, 1,
                    facecolor='gray', edgecolor='black', linewidth=2)
    ax.add_patch(pore)
    ax.text(pore_x+0.25, pore_y+1.3, 'NPC', ha='center', fontsize=8,
           fontweight='bold')
    
    # translocation arrow
    arrow_import = FancyArrowPatch((ar_ligand_x+1, ar_ligand_y),
                                  (9.5, ar_ligand_y),
                                  arrowstyle='->,head_width=0.4,head_length=0.6',
                                  color='green', linewidth=3.5)
    ax.add_patch(arrow_import)
    ax.text(6.5, ar_ligand_y+0.8, 'nuclear import\n(importins)', ha='center', fontsize=8,
           bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.8))
    
    # --- nuclear ar ---
    
    # ar monomer in nucleus
    ar_nuc_x, ar_nuc_y = 11, 9
    ar_nuc_box = FancyBboxPatch((ar_nuc_x-0.8, ar_nuc_y-0.5), 1.6, 1,
                               boxstyle="round,pad=0.1",
                               facecolor='#A23B72',
                               edgecolor='black', linewidth=2, alpha=0.8)
    ax.add_patch(ar_nuc_box)
    ax.text(ar_nuc_x, ar_nuc_y, 'AR:DHT', fontsize=11,
           fontweight='bold', ha='center', color='white')
    
    if show_values and state:
        ax.text(ar_nuc_x, ar_nuc_y-1.2, f'n={state.ar_dht_nucleus:.0f}',
               ha='center', fontsize=8,
               bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # --- dimerization ---
    
    # dimer
    dimer_x, dimer_y = 13.5, 12
    
    # two connected boxes for dimer
    dimer_box1 = FancyBboxPatch((dimer_x-1.2, dimer_y-0.4), 1, 0.8,
                               boxstyle="round,pad=0.08",
                               facecolor='#C73E1D',
                               edgecolor='black', linewidth=2, alpha=0.8)
    ax.add_patch(dimer_box1)
    
    dimer_box2 = FancyBboxPatch((dimer_x-0.2, dimer_y-0.4), 1, 0.8,
                               boxstyle="round,pad=0.08",
                               facecolor='#C73E1D',
                               edgecolor='black', linewidth=2, alpha=0.8)
    ax.add_patch(dimer_box2)
    
    ax.text(dimer_x-0.7, dimer_y, 'AR', fontsize=10,
           fontweight='bold', ha='center', color='white')
    ax.text(dimer_x+0.3, dimer_y, 'AR', fontsize=10,
           fontweight='bold', ha='center', color='white')
    ax.text(dimer_x-0.2, dimer_y-1.1, 'AR Homodimer',
           ha='center', fontsize=9, fontweight='bold')
    
    # dimerization arrow
    arrow_dimer = FancyArrowPatch((ar_nuc_x+0.8, ar_nuc_y+0.3),
                                 (dimer_x-1.2, dimer_y+0.2),
                                 arrowstyle='->,head_width=0.3,head_length=0.5',
                                 color='purple', linewidth=2.5)
    ax.add_patch(arrow_dimer)
    ax.text(12.2, 10.5, 'dimerization', ha='center', fontsize=7,
           bbox=dict(boxstyle='round', facecolor='plum', alpha=0.7))
    
    # --- dna binding ---
    
    # dna (are)
    dna_x, dna_y = 16, 13
    
    # draw dna helix (simplified)
    for i in range(3):
        dna_seg = Rectangle((dna_x-0.2+i*0.3, dna_y-0.1), 0.25, 0.2,
                           facecolor=colors['dna'], edgecolor='black', linewidth=1)
        ax.add_patch(dna_seg)
    
    ax.plot([dna_x, dna_x+0.8], [dna_y, dna_y], 'k-', linewidth=3)
    ax.text(dna_x+0.4, dna_y-0.7, 'ARE', ha='center', fontsize=11,
           fontweight='bold', color=colors['dna'])
    ax.text(dna_x+0.4, dna_y-1.3, '(Androgen Response\nElement)', ha='center', fontsize=7)
    
    # dna binding arrow
    arrow_dna = FancyArrowPatch((dimer_x+0.8, dimer_y+0.3),
                               (dna_x, dna_y+0.3),
                               arrowstyle='->,head_width=0.4,head_length=0.6',
                               color='darkviolet', linewidth=3)
    ax.add_patch(arrow_dna)
    ax.text(14.8, dimer_y+1, 'DNA binding', ha='center', fontsize=8,
           bbox=dict(boxstyle='round', facecolor='lavender', alpha=0.8))
    
    # ar:dna complex
    ax.text(dimer_x-0.2, dimer_y+1.5, 'AR₂:DNA', ha='center', fontsize=11,
           fontweight='bold',
           bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.6))
    
    if show_values and state:
        total_dna_bound = state.ar_testosterone_dna + state.ar_dht_dna
        ax.text(dimer_x-0.2, dimer_y+2.3, f'n={total_dna_bound:.0f}',
               ha='center', fontsize=8,
               bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # --- coregulator recruitment ---
    
    # coactivators
    coact_x, coact_y = 15, 10
    coact_circle = Circle((coact_x, coact_y), 0.5, color=colors['coregulator'],
                         ec='black', linewidth=2, alpha=0.7)
    ax.add_patch(coact_circle)
    ax.text(coact_x, coact_y, 'SRC', ha='center', va='center',
           fontsize=9, fontweight='bold', color='white')
    ax.text(coact_x, coact_y-1, 'Coactivators', ha='center', fontsize=8,
           fontweight='bold')
    
    # recruitment arrow
    arrow_coact = FancyArrowPatch((coact_x, coact_y+0.5), (dimer_x, dimer_y-0.5),
                                 arrowstyle='->,head_width=0.2,head_length=0.4',
                                 color=colors['coregulator'], linewidth=2, linestyle='--')
    ax.add_patch(arrow_coact)
    
    # mediator complex
    mediator_circle = Circle((16.5, 10.5), 0.4, color='#FFB703',
                            ec='black', linewidth=1.5, alpha=0.7)
    ax.add_patch(mediator_circle)
    ax.text(16.5, 10.5, 'MED', ha='center', va='center',
           fontsize=7, fontweight='bold')
    
    # --- transcription ---
    
    # rna polymerase ii
    pol_x, pol_y = 17, 11.5
    pol_box = FancyBboxPatch((pol_x-0.5, pol_y-0.3), 1, 0.6,
                            boxstyle="round,pad=0.08",
                            facecolor='orange',
                            edgecolor='black', linewidth=2, alpha=0.7)
    ax.add_patch(pol_box)
    ax.text(pol_x, pol_y, 'Pol II', ha='center', va='center',
           fontsize=9, fontweight='bold', color='white')
    
    # transcription arrow
    arrow_tx = FancyArrowPatch((pol_x, pol_y-0.4), (pol_x, 8.5),
                              arrowstyle='->,head_width=0.3,head_length=0.5',
                              color='red', linewidth=3, linestyle='-')
    ax.add_patch(arrow_tx)
    ax.text(pol_x+1, 10, 'transcription', ha='center', fontsize=8,
           bbox=dict(boxstyle='round', facecolor='mistyrose', alpha=0.8))
    
    # mrna
    mrna_x, mrna_y = 17, 8
    for i in range(4):
        mrna_circle = Circle((mrna_x-0.15+i*0.1, mrna_y), 0.08,
                            color=colors['protein'], alpha=0.8)
        ax.add_patch(mrna_circle)
    
    ax.text(mrna_x, mrna_y-0.6, 'mRNA', ha='center', fontsize=10,
           fontweight='bold', color=colors['protein'])
    
    # target genes
    targets = ['PSA', 'KLK2', 'TMPRSS2', 'NKX3.1']
    for i, target in enumerate(targets):
        ax.text(mrna_x, mrna_y-1.2-i*0.4, target, ha='center', fontsize=7,
               bbox=dict(boxstyle='round', facecolor='lightcoral', alpha=0.5))
    
    if show_values and state:
        ax.text(mrna_x, mrna_y-3, f'PSA mRNA\nn={state.mrna_psa:.0f}',
               ha='center', fontsize=7,
               bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # --- translation ---
    
    # ribosome and protein
    protein_x, protein_y = 13, 7
    
    # nuclear export
    arrow_export = FancyArrowPatch((mrna_x-0.5, mrna_y), (9, 7),
                                  arrowstyle='->,head_width=0.3,head_length=0.5',
                                  color='blue', linewidth=2, linestyle=':')
    ax.add_patch(arrow_export)
    ax.text(13, 8.2, 'nuclear\nexport', ha='center', fontsize=7,
           bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.7))
    
    # ribosome
    ribosome = Wedge((protein_x, protein_y), 0.4, 0, 180,
                    facecolor='brown', edgecolor='black', linewidth=1.5, alpha=0.7)
    ax.add_patch(ribosome)
    ax.text(protein_x, protein_y-0.7, 'Ribosome', ha='center', fontsize=7)
    
    # protein
    arrow_tl = FancyArrowPatch((protein_x, protein_y-0.5), (protein_x, 4),
                              arrowstyle='->,head_width=0.3,head_length=0.5',
                              color='green', linewidth=2.5)
    ax.add_patch(arrow_tl)
    ax.text(protein_x+1, 5.5, 'translation', ha='center', fontsize=7,
           bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.7))
    
    protein_box = FancyBboxPatch((protein_x-0.6, 3-0.4), 1.2, 0.8,
                                boxstyle="round,pad=0.1",
                                facecolor='#06D6A0',
                                edgecolor='black', linewidth=2, alpha=0.8)
    ax.add_patch(protein_box)
    ax.text(protein_x, 3, 'PSA', ha='center', va='center',
           fontsize=11, fontweight='bold', color='white')
    ax.text(protein_x, 2, 'Protein', ha='center', fontsize=9, fontweight='bold')
    
    if show_values and state:
        ax.text(protein_x, 1.3, f'n={state.protein_psa:.0f}',
               ha='center', fontsize=8,
               bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # --- feedback regulation ---
    
    # fkbp5 feedback
    feedback_x, feedback_y = 10, 13.5
    feedback_circle = Circle((feedback_x, feedback_y), 0.4,
                            color='red', ec='black', linewidth=1.5, alpha=0.6)
    ax.add_patch(feedback_circle)
    ax.text(feedback_x, feedback_y, '⊖', ha='center', va='center',
           fontsize=18, fontweight='bold', color='white')
    ax.text(feedback_x, feedback_y-0.9, 'FKBP5\nfeedback', ha='center', fontsize=7,
           fontweight='bold')
    
    arrow_feedback = FancyArrowPatch((feedback_x, feedback_y-0.5), (ar_nuc_x, ar_nuc_y+0.5),
                                    arrowstyle='->,head_width=0.2,head_length=0.3',
                                    color='red', linewidth=2, linestyle='--')
    ax.add_patch(arrow_feedback)
    
    # --- legend ---
    
    legend_elements = [
        mpatches.Patch(color=colors['ligand'], label='Androgen'),
        mpatches.Patch(color=colors['receptor'], label='AR (inactive)'),
        mpatches.Patch(color='#A23B72', label='AR (activated)'),
        mpatches.Patch(color='#C73E1D', label='AR dimer'),
        mpatches.Patch(color=colors['dna'], label='DNA/ARE'),
        mpatches.Patch(color=colors['protein'], label='mRNA'),
        mpatches.Patch(color='#06D6A0', label='Protein'),
        mpatches.Patch(color=colors['coregulator'], label='Coregulator')
    ]
    
    ax.legend(handles=legend_elements, loc='lower left', 
             frameon=True, shadow=True, fontsize=9)
    
    # set limits and style
    ax.set_xlim(-1, 21)
    ax.set_ylim(-0.5, 18.5)
    ax.set_aspect('equal')
    ax.axis('off')
    
    title_text = 'Androgen Receptor Signaling Pathway'
    if state and show_values:
        title_text += ' (with molecule counts)'
    ax.set_title(title_text, fontsize=20, fontweight='bold', pad=20)
    
    plt.tight_layout()
    
    if save:
        plt.savefig(save, dpi=300, bbox_inches='tight')
    
    if show:
        plt.show()
    else:
        plt.close()
    
    return fig
