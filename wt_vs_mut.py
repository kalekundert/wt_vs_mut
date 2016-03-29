#!/usr/bin/env python2
# encoding: utf-8

import sys, re, collections
from pymol import cmd
from pymol.wizard import Wizard
from pprint import pprint

amino_acids = {
        'ALA': 'A',
        'CYS': 'C',
        'ASP': 'D',
        'GLU': 'E',
        'GLY': 'G',
        'HIS': 'H',
        'ILE': 'I',
        'LYS': 'K',
        'LEU': 'L',
        'MET': 'M',
        'ASN': 'N',
        'PHE': 'F',
        'PRO': 'P',
        'GLN': 'Q',
        'ARG': 'R',
        'SER': 'S',
        'THR': 'T',
        'VAL': 'V',
        'TRP': 'W',
        'TYR': 'Y',
}

blosum_62 = {
        ('A', 'A') :  4,
        ('A', 'R') : -1,
        ('A', 'N') : -2,
        ('A', 'D') : -2,
        ('A', 'C') :  0,
        ('A', 'Q') : -1,
        ('A', 'E') : -1,
        ('A', 'G') :  0,
        ('A', 'H') : -2,
        ('A', 'I') : -1,
        ('A', 'L') : -1,
        ('A', 'K') : -1,
        ('A', 'M') : -1,
        ('A', 'F') : -2,
        ('A', 'P') : -1,
        ('A', 'S') :  1,
        ('A', 'T') :  0,
        ('A', 'W') : -3,
        ('A', 'Y') : -2,
        ('A', 'V') :  0,

        ('R', 'A') : -1,
        ('R', 'R') :  5,
        ('R', 'N') :  0,
        ('R', 'D') : -2,
        ('R', 'C') : -3,
        ('R', 'Q') :  1,
        ('R', 'E') :  0,
        ('R', 'G') : -2,
        ('R', 'H') :  0,
        ('R', 'I') : -3,
        ('R', 'L') : -2,
        ('R', 'K') :  2,
        ('R', 'M') : -1,
        ('R', 'F') : -3,
        ('R', 'P') : -2,
        ('R', 'S') : -1,
        ('R', 'T') : -1,
        ('R', 'W') : -3,
        ('R', 'Y') : -2,
        ('R', 'V') : -3,

        ('N', 'A') : -2,
        ('N', 'R') :  0,
        ('N', 'N') :  6,
        ('N', 'D') :  1,
        ('N', 'C') : -3,
        ('N', 'Q') :  0,
        ('N', 'E') :  0,
        ('N', 'G') :  0,
        ('N', 'H') :  1,
        ('N', 'I') : -3,
        ('N', 'L') : -3,
        ('N', 'K') :  0,
        ('N', 'M') : -2,
        ('N', 'F') : -3,
        ('N', 'P') : -2,
        ('N', 'S') :  1,
        ('N', 'T') :  0,
        ('N', 'W') : -4,
        ('N', 'Y') : -2,
        ('N', 'V') : -3,

        ('D', 'A') : -2,
        ('D', 'R') : -2,
        ('D', 'N') :  1,
        ('D', 'D') :  6,
        ('D', 'C') : -3,
        ('D', 'Q') :  0,
        ('D', 'E') :  2,
        ('D', 'G') : -1,
        ('D', 'H') : -1,
        ('D', 'I') : -3,
        ('D', 'L') : -4,
        ('D', 'K') : -1,
        ('D', 'M') : -3,
        ('D', 'F') : -3,
        ('D', 'P') : -1,
        ('D', 'S') :  0,
        ('D', 'T') : -1,
        ('D', 'W') : -4,
        ('D', 'Y') : -3,
        ('D', 'V') : -3,

        ('C', 'A') :  0,
        ('C', 'R') : -3,
        ('C', 'N') : -3,
        ('C', 'D') : -3,
        ('C', 'C') :  9,
        ('C', 'Q') : -3,
        ('C', 'E') : -4,
        ('C', 'G') : -3,
        ('C', 'H') : -3,
        ('C', 'I') : -1,
        ('C', 'L') : -1,
        ('C', 'K') : -3,
        ('C', 'M') : -1,
        ('C', 'F') : -2,
        ('C', 'P') : -3,
        ('C', 'S') : -1,
        ('C', 'T') : -1,
        ('C', 'W') : -2,
        ('C', 'Y') : -2,
        ('C', 'V') : -1,

        ('Q', 'A') : -1,
        ('Q', 'R') :  1,
        ('Q', 'N') :  0,
        ('Q', 'D') :  0,
        ('Q', 'C') : -3,
        ('Q', 'Q') :  5,
        ('Q', 'E') :  2,
        ('Q', 'G') : -2,
        ('Q', 'H') :  0,
        ('Q', 'I') : -3,
        ('Q', 'L') : -2,
        ('Q', 'K') :  1,
        ('Q', 'M') :  0,
        ('Q', 'F') : -3,
        ('Q', 'P') : -1,
        ('Q', 'S') :  0,
        ('Q', 'T') : -1,
        ('Q', 'W') : -2,
        ('Q', 'Y') : -1,
        ('Q', 'V') : -2,

        ('E', 'A') : -1,
        ('E', 'R') :  0,
        ('E', 'N') :  0,
        ('E', 'D') :  2,
        ('E', 'C') : -4,
        ('E', 'Q') :  2,
        ('E', 'E') :  5,
        ('E', 'G') : -2,
        ('E', 'H') :  0,
        ('E', 'I') : -3,
        ('E', 'L') : -3,
        ('E', 'K') :  1,
        ('E', 'M') : -2,
        ('E', 'F') : -3,
        ('E', 'P') : -1,
        ('E', 'S') :  0,
        ('E', 'T') : -1,
        ('E', 'W') : -3,
        ('E', 'Y') : -2,
        ('E', 'V') : -2,

        ('G', 'A') :  0,
        ('G', 'R') : -2,
        ('G', 'N') :  0,
        ('G', 'D') : -1,
        ('G', 'C') : -3,
        ('G', 'Q') : -2,
        ('G', 'E') : -2,
        ('G', 'G') :  6,
        ('G', 'H') : -2,
        ('G', 'I') : -4,
        ('G', 'L') : -4,
        ('G', 'K') : -2,
        ('G', 'M') : -3,
        ('G', 'F') : -3,
        ('G', 'P') : -2,
        ('G', 'S') :  0,
        ('G', 'T') : -2,
        ('G', 'W') : -2,
        ('G', 'Y') : -3,
        ('G', 'V') : -3,

        ('H', 'A') : -2,
        ('H', 'R') :  0,
        ('H', 'N') :  1,
        ('H', 'D') : -1,
        ('H', 'C') : -3,
        ('H', 'Q') :  0,
        ('H', 'E') :  0,
        ('H', 'G') : -2,
        ('H', 'H') :  8,
        ('H', 'I') : -3,
        ('H', 'L') : -3,
        ('H', 'K') : -1,
        ('H', 'M') : -2,
        ('H', 'F') : -1,
        ('H', 'P') : -2,
        ('H', 'S') : -1,
        ('H', 'T') : -2,
        ('H', 'W') : -2,
        ('H', 'Y') :  2,
        ('H', 'V') : -3,

        ('I', 'A') : -1,
        ('I', 'R') : -3,
        ('I', 'N') : -3,
        ('I', 'D') : -3,
        ('I', 'C') : -1,
        ('I', 'Q') : -3,
        ('I', 'E') : -3,
        ('I', 'G') : -4,
        ('I', 'H') : -3,
        ('I', 'I') :  4,
        ('I', 'L') :  2,
        ('I', 'K') : -3,
        ('I', 'M') :  1,
        ('I', 'F') :  0,
        ('I', 'P') : -3,
        ('I', 'S') : -2,
        ('I', 'T') : -1,
        ('I', 'W') : -3,
        ('I', 'Y') : -1,
        ('I', 'V') :  3,

        ('L', 'A') : -1,
        ('L', 'R') : -2,
        ('L', 'N') : -3,
        ('L', 'D') : -4,
        ('L', 'C') : -1,
        ('L', 'Q') : -2,
        ('L', 'E') : -3,
        ('L', 'G') : -4,
        ('L', 'H') : -3,
        ('L', 'I') :  2,
        ('L', 'L') :  4,
        ('L', 'K') : -2,
        ('L', 'M') :  2,
        ('L', 'F') :  0,
        ('L', 'P') : -3,
        ('L', 'S') : -2,
        ('L', 'T') : -1,
        ('L', 'W') : -2,
        ('L', 'Y') : -1,
        ('L', 'V') :  1,

        ('K', 'A') : -1,
        ('K', 'R') :  2,
        ('K', 'N') :  0,
        ('K', 'D') : -1,
        ('K', 'C') : -3,
        ('K', 'Q') :  1,
        ('K', 'E') :  1,
        ('K', 'G') : -2,
        ('K', 'H') : -1,
        ('K', 'I') : -3,
        ('K', 'L') : -2,
        ('K', 'K') :  5,
        ('K', 'M') : -1,
        ('K', 'F') : -3,
        ('K', 'P') : -1,
        ('K', 'S') :  0,
        ('K', 'T') : -1,
        ('K', 'W') : -3,
        ('K', 'Y') : -2,
        ('K', 'V') : -2,

        ('M', 'A') : -1,
        ('M', 'R') : -1,
        ('M', 'N') : -2,
        ('M', 'D') : -3,
        ('M', 'C') : -1,
        ('M', 'Q') :  0,
        ('M', 'E') : -2,
        ('M', 'G') : -3,
        ('M', 'H') : -2,
        ('M', 'I') :  1,
        ('M', 'L') :  2,
        ('M', 'K') : -1,
        ('M', 'M') :  5,
        ('M', 'F') :  0,
        ('M', 'P') : -2,
        ('M', 'S') : -1,
        ('M', 'T') : -1,
        ('M', 'W') : -1,
        ('M', 'Y') : -1,
        ('M', 'V') :  1,

        ('F', 'A') : -2,
        ('F', 'R') : -3,
        ('F', 'N') : -3,
        ('F', 'D') : -3,
        ('F', 'C') : -2,
        ('F', 'Q') : -3,
        ('F', 'E') : -3,
        ('F', 'G') : -3,
        ('F', 'H') : -1,
        ('F', 'I') :  0,
        ('F', 'L') :  0,
        ('F', 'K') : -3,
        ('F', 'M') :  0,
        ('F', 'F') :  6,
        ('F', 'P') : -4,
        ('F', 'S') : -2,
        ('F', 'T') : -2,
        ('F', 'W') :  1,
        ('F', 'Y') :  3,
        ('F', 'V') : -1,

        ('P', 'A') : -1,
        ('P', 'R') : -2,
        ('P', 'N') : -2,
        ('P', 'D') : -1,
        ('P', 'C') : -3,
        ('P', 'Q') : -1,
        ('P', 'E') : -1,
        ('P', 'G') : -2,
        ('P', 'H') : -2,
        ('P', 'I') : -3,
        ('P', 'L') : -3,
        ('P', 'K') : -1,
        ('P', 'M') : -2,
        ('P', 'F') : -4,
        ('P', 'P') :  7,
        ('P', 'S') : -1,
        ('P', 'T') : -1,
        ('P', 'W') : -4,
        ('P', 'Y') : -3,
        ('P', 'V') : -2,

        ('S', 'A') :  1,
        ('S', 'R') : -1,
        ('S', 'N') :  1,
        ('S', 'D') :  0,
        ('S', 'C') : -1,
        ('S', 'Q') :  0,
        ('S', 'E') :  0,
        ('S', 'G') :  0,
        ('S', 'H') : -1,
        ('S', 'I') : -2,
        ('S', 'L') : -2,
        ('S', 'K') :  0,
        ('S', 'M') : -1,
        ('S', 'F') : -2,
        ('S', 'P') : -1,
        ('S', 'S') :  4,
        ('S', 'T') :  1,
        ('S', 'W') : -3,
        ('S', 'Y') : -2,
        ('S', 'V') : -2,

        ('T', 'A') :  0,
        ('T', 'R') : -1,
        ('T', 'N') :  0,
        ('T', 'D') : -1,
        ('T', 'C') : -1,
        ('T', 'Q') : -1,
        ('T', 'E') : -1,
        ('T', 'G') : -2,
        ('T', 'H') : -2,
        ('T', 'I') : -1,
        ('T', 'L') : -1,
        ('T', 'K') : -1,
        ('T', 'M') : -1,
        ('T', 'F') : -2,
        ('T', 'P') : -1,
        ('T', 'S') :  1,
        ('T', 'T') :  5,
        ('T', 'W') : -2,
        ('T', 'Y') : -2,
        ('T', 'V') :  0,

        ('W', 'A') : -3,
        ('W', 'R') : -3,
        ('W', 'N') : -4,
        ('W', 'D') : -4,
        ('W', 'C') : -2,
        ('W', 'Q') : -2,
        ('W', 'E') : -3,
        ('W', 'G') : -2,
        ('W', 'H') : -2,
        ('W', 'I') : -3,
        ('W', 'L') : -2,
        ('W', 'K') : -3,
        ('W', 'M') : -1,
        ('W', 'F') :  1,
        ('W', 'P') : -4,
        ('W', 'S') : -3,
        ('W', 'T') : -2,
        ('W', 'W') : 11,
        ('W', 'Y') :  2,
        ('W', 'V') : -3,

        ('Y', 'A') : -2,
        ('Y', 'R') : -2,
        ('Y', 'N') : -2,
        ('Y', 'D') : -3,
        ('Y', 'C') : -2,
        ('Y', 'Q') : -1,
        ('Y', 'E') : -2,
        ('Y', 'G') : -3,
        ('Y', 'H') :  2,
        ('Y', 'I') : -1,
        ('Y', 'L') : -1,
        ('Y', 'K') : -2,
        ('Y', 'M') : -1,
        ('Y', 'F') :  3,
        ('Y', 'P') : -3,
        ('Y', 'S') : -2,
        ('Y', 'T') : -2,
        ('Y', 'W') :  2,
        ('Y', 'Y') :  7,
        ('Y', 'V') : -1,

        ('V', 'A') :  0,
        ('V', 'R') : -3,
        ('V', 'N') : -3,
        ('V', 'D') : -3,
        ('V', 'C') : -1,
        ('V', 'Q') : -2,
        ('V', 'E') : -2,
        ('V', 'G') : -3,
        ('V', 'H') : -3,
        ('V', 'I') :  3,
        ('V', 'L') :  1,
        ('V', 'K') : -2,
        ('V', 'M') :  1,
        ('V', 'F') : -1,
        ('V', 'P') : -2,
        ('V', 'S') : -2,
        ('V', 'T') :  0,
        ('V', 'W') : -3,
        ('V', 'Y') : -1,
        ('V', 'V') :  4,
}


class WildtypeVsMutant (Wizard):
    """
    Find the residues that differ between two selections and allow the user to
    view each difference one at a time.
    """

    def __init__(self, wildtype_obj='', mutant_obj='', focus_sele=''):
        self.mutations = []
        self.active_environments = []
        self.active_mutations = []
        self.aligned_seqs = '', ''
        self.aligned_resis = [], []
        self.neighbor_radius = 4
        self.zoom_padding = 5
        self.wildtype_hilite = 'white'
        self.mutant_hilite = 'yellow'
        self.show_polar_h = False
        self.active_prompt = ''
        self.original_view = cmd.get_view()
        self.wildtype_obj = ''
        self.mutant_obj = ''
        self.focus_sele = ''

        self.set_focus_sele(focus_sele)
        self.set_wildtype_object(wildtype_obj)
        self.set_mutant_object(mutant_obj)

    def set_wildtype_object(self, wildtype_obj):
        self.wildtype_obj = wildtype_obj
        self.active_prompt = ''; self.redraw()
        self.update_mutation_list()

    def set_mutant_object(self, mutant_obj):
        self.mutant_obj = mutant_obj
        self.active_prompt = ''; self.redraw()
        self.update_mutation_list()

    def set_focus_sele(self, focus_sele):
        self.focus_sele = focus_sele
        self.active_prompt = ''; self.redraw()
        self.update_mutation_list()

    def update_mutation_list(self):
        """
        Find sequence differences between the wildtype and mutant structures to
        highlight.  Insertions and deletions are handled smoothly because the
        sequences are aligned (using the Nedleman-Wunsch algorithm) before they
        being compared.
        """

        if not self.wildtype_obj or not self.mutant_obj:
            return

        # unaligned_seqmaps: A tuple of two ordered dictionaries mapping
        # (residue id, chain id) pairs to one-letter residue names.  These two
        # mappings encode the sequences on the wildtype and mutant objects.

        unaligned_seqmaps = (
                get_sequence(self.wildtype_obj),
                get_sequence(self.mutant_obj),
        )

        # focus_resis: A list of all the (residue id, chain id) pairs for all
        # the residues the user is interested in, or undefined if the user is
        # interested in every residue.  If defined, the list will also contain
        # None tacked onto the end, indicating that the user is also interested
        # in gaps.  Mutations not in this list should not be shown.

        if self.focus_sele:
            focus_resis = get_sequence(self.focus_sele).keys() + [None]

        # unaligned_seqs: A tuple of two strings representing the sequences of
        # the wildtype and mutant objects.  The mapping to residue and chain
        # ids has been removed.

        unaligned_seqs = tuple(
                ''.join(seqmap.values())
                for seqmap in unaligned_seqmaps
        )

        # self.aligned_seqs: A tuple of two strings representing the aligned
        # wildtype and mutant sequences.  These sequences may differ from the
        # unaligned sequences in that they may have dashes inserted where the
        # two sequences don't align well.

        self.aligned_seqs = get_alignment(*unaligned_seqs)

        # self.aligned_resis: A tuple of two lists which associate (residue id,
        # chain id) pairs with each residue in the aligned sequences.  Missing
        # residues (i.e. dashes) are represented in these lists by None.

        self.aligned_resis = [], []

        for i in range(2):
            for res in self.aligned_seqs[i]:
                self.aligned_resis[i].append(
                        unaligned_seqmaps[i].popitem(last=False)[0]
                        if res != '-' else (None, None))

        # self.mutations: A list of indices of positions that differ between
        # the two aligned sequences, excluding terminal gaps.  These indices
        # can be used with both self.aligned_seqs and self.aligned_resis.

        def in_focus_sele(i):
            if not self.focus_sele:
                return True
            else:
                return (self.aligned_resis[0][i] in focus_resis or
                        self.aligned_resis[1][i] in focus_resis)

        self.mutations = [
                i for i in find_mutations(self.aligned_seqs)
                if in_focus_sele(i)
        ]

        # Automatically zoom in on the first mutation.

        self.cycle()

    def set_active_mutation(self, index):
        """
        Set the mutation being currently highlighted.
        """
        self.active_mutations = [index]
        self.redraw()

    def set_neighbor_radius(self, radius):
        """
        Change how many sidechains are included in the zoomed view.
        """
        self.neighbor_radius = radius
        self.redraw()

    def set_zoom_padding(self, padding):
        """
        Change the zoom level on the displayed sidechains.
        """
        self.zoom_padding = padding
        self.redraw()

    def set_wildtype_hilite(self, color):
        """
        Change the color used to highlight the wildtype residue.
        """
        self.wildtype_hilite = color
        self.redraw()

    def set_mutant_hilite(self, color):
        """
        Change the color used to highlight the mutant residue.
        """
        self.mutant_hilite = color
        self.redraw()

    def set_show_polar_h(self, option):
        """
        Change whether or not polar hydrogens are displayed.
        """
        self.show_polar_h = option
        self.redraw()

    def get_next_mutation(self):
        """
        Return the index of the next mutation.  The index is relative to the
        sequence alignment. Resets to first mutation if multiple mutations are active.
        """
        if not self.mutations:
            return None

        if len(self.active_mutations) == 1:
            index = self.mutations.index(self.active_mutations[0]) + 1
            return self.mutations[index]
        else:
            return self.mutations[0]

    def get_mutation_name(self, muti=None):
        """
        Return the name of given mutation.  The `muti' parameter is an index
        into the sequence alignment that defaults to the currently active
        mutation.  The names are formatted like so:

        Mutation:  D38E
        Insertion: -38E
        Deletion:  D38-
        """
        if muti is None:
            mutis = self.active_mutations
        else:
            mutis = [muti]

        return_string = ''
        for i, muti in enumerate(mutis):
            # Figure out the one-letter residue names from the alignment.

            wildtype_res = self.aligned_seqs[0][muti]
            mutant_res = self.aligned_seqs[1][muti]

            # Figure out the residue number from the wildtype (residue id, chain
            # id) list.  Insertions are handled specially.  For this case there is
            # no wildtype residue number, to we search back for the nearest non-gap
            # in the wildtype sequence and use its number instead.

            resi = None
            while resi is None:
                resi, chain = self.aligned_resis[0][muti]; muti -= 1

            if i + 1 >= len(mutis):
                format_string = '{}{}{}'
            else:
                format_string = '{}{}{}, '

            return_string += format_string.format(wildtype_res, resi, mutant_res)

        # Return the concatenated name.

        return return_string

    def get_panel(self):
        if not self.mutant_obj or not self.wildtype_obj:
            return [
                [1, 'Wildtype vs Mutant Wizard', ''],
                [2, 'Cancel', 'cmd.get_wizard().cleanup()'],
            ]

        buttons = [
            [1, 'Wildtype vs Mutant Wizard', ''],
            [2, 'View next mutation', 'cmd.get_wizard().cycle()'],
            [2, 'Show all mutations', 'cmd.get_wizard().show_all_mutations()'],
            [3, 'Wildtype highlight: {}'.format(self.wildtype_hilite), 'wt_hilite'],
            [3, 'Mutant highlight: {}'.format(self.mutant_hilite), 'mut_hilite'],
            [3, 'Neighbor radius: {0:.1f}A'.format(self.neighbor_radius), 'radius'],
            [3, 'Zoom padding: {0:.1f}A'.format(self.zoom_padding), 'padding'],
            [3, 'Polar hydrogens: {}'.format('show' if self.show_polar_h else 'hide'), 'hydrogen'],
        ]

        for muti in self.mutations:
            command = 'cmd.get_wizard().set_active_mutation({})'
            name = self.get_mutation_name(muti)
            if muti in self.active_mutations: name += ' <--'
            buttons += [[2, name, command.format(muti)]]

        buttons += [[2, 'Done', 'cmd.get_wizard().cleanup()']]
        return buttons

    def get_menu(self, tag):
        menus = {
            'wt_hilite': [[2, 'Highlight Color', '']],
            'mut_hilite': [[2, 'Highlight Color', '']],
            'radius': [[2, 'Neighbor Radius', '']],
            'padding': [[2, 'Zoom Padding', '']],
            'hydrogen': [[2, 'Polar Hydrogens', '']],
        }

        # Define the highlight color menu.
        colors = (
                '\\900red',
                '\\090green',
                '\\009blue',
                '\\990yellow',
                '\\909magenta',
                '\\099cyan',
                '\\950orange',
                '\\999white',
                '\\999none',
        )
        for color in colors:
            menus['wt_hilite'] += [[
                1, color,
                'cmd.get_wizard().set_wildtype_hilite("{}")'.format(color[4:])]]
            menus['mut_hilite'] += [[
                1, color,
                'cmd.get_wizard().set_mutant_hilite("{}")'.format(color[4:])]]

        # Define the neighbor radius menu.
        for radius in range(1, 11):
            menus['radius'] += [[
                1, '{0:0.1f}A'.format(radius),
                'cmd.get_wizard().set_neighbor_radius({})'.format(radius)]]

        # Define the zoom padding menu.
        for padding in range(15 + 1):
            menus['padding'] += [[
                1, '{0}A'.format(padding),
                'cmd.get_wizard().set_zoom_padding({})'.format(padding)]]

        # Define the polar hydrogen menu.
        toggled_option = not self.show_polar_h
        menus['hydrogen'] += [[
            1, 'show' if toggled_option else 'hide',
            'cmd.get_wizard().set_show_polar_h({})'.format(toggled_option)]]

        # Return the right menu.
        return menus[tag]

    def get_prompt(self):
        # The \999 code changes the text color to white.
        if not self.wildtype_obj:
            return ["Select the wildtype object: \\999{}".format(self.active_prompt)]
        elif not self.mutant_obj:
            return ["Select the mutant object: \\999{}".format(self.active_prompt)]
        else:
            return ["Press <Ctrl-Space> to view the next mutation..."]

    def do_key(self, key, x, y, mod):
        """
        Take responsibility for handling key presses.
        """
        # If the user presses Ctrl-Space (key, mod == 0, 2) cycle to the next
        # mutation.  Otherwise let pymol handle the key press.

        if self.wildtype_obj and self.mutant_obj:
            if (key, mod) == (0, 2):
                self.cycle()
            else:
                return 0

        # If the user hasn't specified a wildtype and a mutant object, prompt
        # for them.  Letters (key >= 32), Backspace (key in 8, 127), and Enter
        # (key in 10, 13) are handled specially.  Everything else is passed on
        # to pymol.

        else:
            if key in (8, 127):
                self.active_prompt = self.active_prompt[:-1]
            elif key >= 32:
                self.active_prompt += chr(key)
            elif key in (10, 13):
                if not self.wildtype_obj: self.set_wildtype_object(self.active_prompt)
                elif not self.mutant_obj: self.set_mutant_object(self.active_prompt)
            else:
                return 0

        cmd.refresh_wizard()
        return 1

    def get_event_mask(self):
        return Wizard.event_mask_key

    def cycle(self):
        """
        Focus on the next mutant, or close the wizard if the last mutant is
        currently active.
        """
        try: self.set_active_mutation(self.get_next_mutation())
        except IndexError: self.cleanup()

    def redraw(self):
        """
        Highlight the sidechains around the current mutation.  Everything this
        method adds to the scene is put in its own object, so that it can be
        easily undone by the next call to redraw() or cleanup().
        """
        cmd.refresh_wizard()

        self.delete_active_environments()

        for muti in self.active_mutations:
            wt_obj = self.wildtype_obj
            wt_resi, wt_chain = self.aligned_resis[0][muti]
            wt_sele = 'none' if wt_resi is None else \
                      'resi {} and chain {}'.format(wt_resi, wt_chain)
            mut_obj = self.mutant_obj
            mut_resi, mut_chain = self.aligned_resis[1][muti]
            mut_sele = 'none' if mut_resi is None else \
                      'resi {} and chain {}'.format(mut_resi, mut_chain)
            h_sele = (
                    '(elem H and (neighbor elem C))' if self.show_polar_h else
                    '(elem H)')
            env_sele = (
                    '(byres {{}} within {self.neighbor_radius} of '
                    '(({wt_obj} and {wt_sele}) or ({mut_obj} and {mut_sele}))) '
                    'and not {h_sele}'.format(**locals()))

            initial_view = cmd.get_view()

            if len(self.active_mutations) > 1:
                wt_env_name = 'wt_env_' + self.get_mutation_name(muti)
                mut_env_name = 'mut_env_' + self.get_mutation_name(muti)
            else:
                wt_env_name = 'wt_env'
                mut_env_name = 'wt_env'

            self.create_environment(wt_env_name, env_sele.format(wt_obj))
            self.create_environment(mut_env_name, env_sele.format(mut_obj))
            self.draw_self_self_hbonds(mut_env_name)
            if self.wildtype_hilite != 'none':
                cmd.color(self.wildtype_hilite, '{wt_env_name} and {wt_sele} and elem C'.format(**locals()))
            if self.mutant_hilite != 'none':
                cmd.color(self.mutant_hilite, '{mut_env_name} and {mut_sele} and elem C'.format(**locals()))
            cmd.show_as('sticks', wt_env_name)
            cmd.show_as('sticks', mut_env_name)

        if len(self.active_mutations) == 1:
            cmd.set_view(initial_view)
            cmd.zoom(wt_env_name, buffer=self.zoom_padding, animate=-1)

    def create_environment(self, env_name, env):
        """
        Create a new environment in PyMOL, and keep track of it in member variable.
        """
        cmd.create(env_name, env)
        self.active_environments.append(env_name)

    def draw_self_self_hbonds(self, env_name):
        hbonds_env_name = env_name + '_hbonds'
        cmd.dist(hbonds_env_name, env_name, env_name, mode=2)
        self.active_environments.append(hbonds_env_name)

    def delete_active_environments(self):
        """
        Delete any created environments
        """
        for env in self.active_environments:
            cmd.delete(env)
        self.active_environments = []

    def cleanup(self):
        """
        Remove any changes this wizard has made to the scene.
        """
        self.delete_active_environments()
        cmd.set_view(self.original_view)
        cmd.set_wizard()

    def show_all_mutations(self):
        """
        Draws all mutations simultaneously
        """
        self.delete_active_environments()
        self.active_mutations = self.mutations
        self.redraw()
        cmd.set_view(self.original_view)


def wt_vs_mut(mutant_obj=None, wildtype_obj=None, focus_sele=None):
    """
    Provide a convenient way to launch the wizard from the command-line.
    """
    wizard = WildtypeVsMutant(mutant_obj, wildtype_obj, focus_sele)
    cmd.set_wizard(wizard)

def get_sequence(selection):
    """
    Return the protein sequence associated with the given object or selection.
    """
    sequence = collections.OrderedDict()
    str_to_int = int; aa_table = amino_acids
    sequence_builder = \
            'sequence[str_to_int(resi), chain] = aa_table.get(resn, "X")'
    cmd.iterate(selection, sequence_builder, space=locals())
    return sequence

def get_alignment(seq1, seq2, score_matrix=blosum_62, gap_penalty=-20):
    """
    Return the alignment between the two given sequences.  The alignment is
    performed using the Needleman-Wunsch algorithm, as implemented by Evan
    Dempsey and published on his blog:

    https://evandempsey.wordpress.com/2013/01/08/needleman-wunsch-algorithm-for-dna-sequence-alignment/

    The default score function is BLOSUM62 and the default gap penalty is -20.
    This is a very large gap penalty because I want to discourage gaps as much
    as possible.  Each gap is taken as a mutation, so the more gaps there are
    the more things the user has to scroll through.  The point of the alignment
    is just to allow slightly different sequences to be compared, not to really
    figure out the best way to match up each residue.
    """
    n = len(seq1)
    m = len(seq2)

    # Make two-dimensional list for subproblem solutions.

    subproblems = [[0 for x in range(m+1)] for x in range(n+1)]

    # Fill in zeros on both dimensions with gap penalties.

    for i in range(n+1):
        subproblems[i][0] = i * gap_penalty

    for j in range(m+1):
        subproblems[0][j] = j * gap_penalty

    # Calculate subproblem solutions.

    for i in range(1, n+1):
        for j in range(1, m+1):
            case1 = subproblems[i-1][j-1] + get_score_from_matrix(seq1[i-1], seq2[j-1], score_matrix=score_matrix)
            case2 = subproblems[i-1][j] + gap_penalty
            case3 = subproblems[i][j-1] + gap_penalty
            subproblems[i][j] = max(case1, case2, case3)

    # Backtrace to reconstruct optimal alignment.

    alignment1 = ''
    alignment2 = ''

    i = n
    j = m
    while i > 0 or j > 0:
        pos = subproblems[i][j]
        case1 = subproblems[i-1][j-1] + get_score_from_matrix(seq1[i-1], seq2[j-1], score_matrix=score_matrix)
        case2 = subproblems[i-1][j] + gap_penalty
        case3 = subproblems[i][j-1] + gap_penalty

        if i > 0 and pos == case1:
            alignment1 = seq1[i-1] + alignment1
            alignment2 = seq2[j-1] + alignment2
            i -= 1
            j -= 1
        elif i > 0 and pos == case2:
            alignment1 = seq1[i-1] + alignment1
            alignment2 = '-' + alignment2
            i -= 1
        elif j > 0 and pos == case3:
            alignment1 = '-' + alignment1
            alignment2 = seq2[j-1] + alignment2
            j -= 1

    return alignment1, alignment2

def find_mutations(alignment):
    """
    Return a list of any indices that differ between the aligned sequences.
    Internal gaps are counted, but terminal gaps are not.
    """
    first_non_gap = max(
            re.search('^-*', seq).end()
            for seq in alignment
    )
    last_non_gap = min(
            re.search('-*$', seq).start()
            for seq in alignment
    )
    return [
            i for i in range(first_non_gap, last_non_gap)
            if alignment[0][i] != alignment[1][i]
    ]

def does_sele_exist(sele):
    from types import ListType
    session = cmd.get_session()
    for i in session["names"]:
        if type(i) is ListType:
            if sele == i[0]:
                return True
    return False

def get_score_from_matrix(pos1, pos2, score_matrix=blosum_62):
    """
    Return score from score_matrix if it exists. Otherwise (as for non-canonicals),
    return generic score values for matches and mismatches.
    """
    try:
        return score_matrix[pos1, pos2]
    except KeyError:
        if pos1 == pos2:
            return 5
        else:
            return -4

## Add "wt_vs_mut" as pymol command
cmd.extend('wt_vs_mut', wt_vs_mut)
cmd.auto_arg[0]['wt_vs_mut'] = cmd.auto_arg[0]['zoom']
cmd.auto_arg[1]['wt_vs_mut'] = cmd.auto_arg[0]['zoom']

## Trick to get "wizard wt_vs_mut" working
sys.modules['pymol.wizard.wt_vs_mut'] = sys.modules[__name__]

## Add item to plugin menu
try:
    from pymol.plugins import addmenuitem
    def __init_plugin__(self): # (no fold)
        addmenuitem('Wildtype vs. Mutant', lambda s=self: wt_vs_mut())
except:
    def __init__(self): # (no fold)
        self.menuBar.addmenuitem(
                'Plugin', 'command', 'Wildtype vs. Mutant',
                label='Wildtype vs. Mutant', command=lambda s=self: wt_vs_mut())
