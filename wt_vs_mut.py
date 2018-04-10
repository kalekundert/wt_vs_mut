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
        ('A', 'A') :  4, ('R', 'A') : -1, ('N', 'A') : -2, ('D', 'A') : -2, ('C', 'A') :  0, ('Q', 'A') : -1, ('E', 'A') : -1, ('G', 'A') :  0, ('H', 'A') : -2, ('I', 'A') : -1, ('L', 'A') : -1, ('K', 'A') : -1, ('M', 'A') : -1, ('F', 'A') : -2, ('P', 'A') : -1, ('S', 'A') :  1, ('T', 'A') :  0, ('W', 'A') : -3, ('Y', 'A') : -2, ('V', 'A') :  0,
        ('A', 'R') : -1, ('R', 'R') :  5, ('N', 'R') :  0, ('D', 'R') : -2, ('C', 'R') : -3, ('Q', 'R') :  1, ('E', 'R') :  0, ('G', 'R') : -2, ('H', 'R') :  0, ('I', 'R') : -3, ('L', 'R') : -2, ('K', 'R') :  2, ('M', 'R') : -1, ('F', 'R') : -3, ('P', 'R') : -2, ('S', 'R') : -1, ('T', 'R') : -1, ('W', 'R') : -3, ('Y', 'R') : -2, ('V', 'R') : -3,
        ('A', 'N') : -2, ('R', 'N') :  0, ('N', 'N') :  6, ('D', 'N') :  1, ('C', 'N') : -3, ('Q', 'N') :  0, ('E', 'N') :  0, ('G', 'N') :  0, ('H', 'N') :  1, ('I', 'N') : -3, ('L', 'N') : -3, ('K', 'N') :  0, ('M', 'N') : -2, ('F', 'N') : -3, ('P', 'N') : -2, ('S', 'N') :  1, ('T', 'N') :  0, ('W', 'N') : -4, ('Y', 'N') : -2, ('V', 'N') : -3,
        ('A', 'D') : -2, ('R', 'D') : -2, ('N', 'D') :  1, ('D', 'D') :  6, ('C', 'D') : -3, ('Q', 'D') :  0, ('E', 'D') :  2, ('G', 'D') : -1, ('H', 'D') : -1, ('I', 'D') : -3, ('L', 'D') : -4, ('K', 'D') : -1, ('M', 'D') : -3, ('F', 'D') : -3, ('P', 'D') : -1, ('S', 'D') :  0, ('T', 'D') : -1, ('W', 'D') : -4, ('Y', 'D') : -3, ('V', 'D') : -3,
        ('A', 'C') :  0, ('R', 'C') : -3, ('N', 'C') : -3, ('D', 'C') : -3, ('C', 'C') :  9, ('Q', 'C') : -3, ('E', 'C') : -4, ('G', 'C') : -3, ('H', 'C') : -3, ('I', 'C') : -1, ('L', 'C') : -1, ('K', 'C') : -3, ('M', 'C') : -1, ('F', 'C') : -2, ('P', 'C') : -3, ('S', 'C') : -1, ('T', 'C') : -1, ('W', 'C') : -2, ('Y', 'C') : -2, ('V', 'C') : -1,
        ('A', 'Q') : -1, ('R', 'Q') :  1, ('N', 'Q') :  0, ('D', 'Q') :  0, ('C', 'Q') : -3, ('Q', 'Q') :  5, ('E', 'Q') :  2, ('G', 'Q') : -2, ('H', 'Q') :  0, ('I', 'Q') : -3, ('L', 'Q') : -2, ('K', 'Q') :  1, ('M', 'Q') :  0, ('F', 'Q') : -3, ('P', 'Q') : -1, ('S', 'Q') :  0, ('T', 'Q') : -1, ('W', 'Q') : -2, ('Y', 'Q') : -1, ('V', 'Q') : -2,
        ('A', 'E') : -1, ('R', 'E') :  0, ('N', 'E') :  0, ('D', 'E') :  2, ('C', 'E') : -4, ('Q', 'E') :  2, ('E', 'E') :  5, ('G', 'E') : -2, ('H', 'E') :  0, ('I', 'E') : -3, ('L', 'E') : -3, ('K', 'E') :  1, ('M', 'E') : -2, ('F', 'E') : -3, ('P', 'E') : -1, ('S', 'E') :  0, ('T', 'E') : -1, ('W', 'E') : -3, ('Y', 'E') : -2, ('V', 'E') : -2,
        ('A', 'G') :  0, ('R', 'G') : -2, ('N', 'G') :  0, ('D', 'G') : -1, ('C', 'G') : -3, ('Q', 'G') : -2, ('E', 'G') : -2, ('G', 'G') :  6, ('H', 'G') : -2, ('I', 'G') : -4, ('L', 'G') : -4, ('K', 'G') : -2, ('M', 'G') : -3, ('F', 'G') : -3, ('P', 'G') : -2, ('S', 'G') :  0, ('T', 'G') : -2, ('W', 'G') : -2, ('Y', 'G') : -3, ('V', 'G') : -3,
        ('A', 'H') : -2, ('R', 'H') :  0, ('N', 'H') :  1, ('D', 'H') : -1, ('C', 'H') : -3, ('Q', 'H') :  0, ('E', 'H') :  0, ('G', 'H') : -2, ('H', 'H') :  8, ('I', 'H') : -3, ('L', 'H') : -3, ('K', 'H') : -1, ('M', 'H') : -2, ('F', 'H') : -1, ('P', 'H') : -2, ('S', 'H') : -1, ('T', 'H') : -2, ('W', 'H') : -2, ('Y', 'H') :  2, ('V', 'H') : -3,
        ('A', 'I') : -1, ('R', 'I') : -3, ('N', 'I') : -3, ('D', 'I') : -3, ('C', 'I') : -1, ('Q', 'I') : -3, ('E', 'I') : -3, ('G', 'I') : -4, ('H', 'I') : -3, ('I', 'I') :  4, ('L', 'I') :  2, ('K', 'I') : -3, ('M', 'I') :  1, ('F', 'I') :  0, ('P', 'I') : -3, ('S', 'I') : -2, ('T', 'I') : -1, ('W', 'I') : -3, ('Y', 'I') : -1, ('V', 'I') :  3,
        ('A', 'L') : -1, ('R', 'L') : -2, ('N', 'L') : -3, ('D', 'L') : -4, ('C', 'L') : -1, ('Q', 'L') : -2, ('E', 'L') : -3, ('G', 'L') : -4, ('H', 'L') : -3, ('I', 'L') :  2, ('L', 'L') :  4, ('K', 'L') : -2, ('M', 'L') :  2, ('F', 'L') :  0, ('P', 'L') : -3, ('S', 'L') : -2, ('T', 'L') : -1, ('W', 'L') : -2, ('Y', 'L') : -1, ('V', 'L') :  1,
        ('A', 'K') : -1, ('R', 'K') :  2, ('N', 'K') :  0, ('D', 'K') : -1, ('C', 'K') : -3, ('Q', 'K') :  1, ('E', 'K') :  1, ('G', 'K') : -2, ('H', 'K') : -1, ('I', 'K') : -3, ('L', 'K') : -2, ('K', 'K') :  5, ('M', 'K') : -1, ('F', 'K') : -3, ('P', 'K') : -1, ('S', 'K') :  0, ('T', 'K') : -1, ('W', 'K') : -3, ('Y', 'K') : -2, ('V', 'K') : -2,
        ('A', 'M') : -1, ('R', 'M') : -1, ('N', 'M') : -2, ('D', 'M') : -3, ('C', 'M') : -1, ('Q', 'M') :  0, ('E', 'M') : -2, ('G', 'M') : -3, ('H', 'M') : -2, ('I', 'M') :  1, ('L', 'M') :  2, ('K', 'M') : -1, ('M', 'M') :  5, ('F', 'M') :  0, ('P', 'M') : -2, ('S', 'M') : -1, ('T', 'M') : -1, ('W', 'M') : -1, ('Y', 'M') : -1, ('V', 'M') :  1,
        ('A', 'F') : -2, ('R', 'F') : -3, ('N', 'F') : -3, ('D', 'F') : -3, ('C', 'F') : -2, ('Q', 'F') : -3, ('E', 'F') : -3, ('G', 'F') : -3, ('H', 'F') : -1, ('I', 'F') :  0, ('L', 'F') :  0, ('K', 'F') : -3, ('M', 'F') :  0, ('F', 'F') :  6, ('P', 'F') : -4, ('S', 'F') : -2, ('T', 'F') : -2, ('W', 'F') :  1, ('Y', 'F') :  3, ('V', 'F') : -1,
        ('A', 'P') : -1, ('R', 'P') : -2, ('N', 'P') : -2, ('D', 'P') : -1, ('C', 'P') : -3, ('Q', 'P') : -1, ('E', 'P') : -1, ('G', 'P') : -2, ('H', 'P') : -2, ('I', 'P') : -3, ('L', 'P') : -3, ('K', 'P') : -1, ('M', 'P') : -2, ('F', 'P') : -4, ('P', 'P') :  7, ('S', 'P') : -1, ('T', 'P') : -1, ('W', 'P') : -4, ('Y', 'P') : -3, ('V', 'P') : -2,
        ('A', 'S') :  1, ('R', 'S') : -1, ('N', 'S') :  1, ('D', 'S') :  0, ('C', 'S') : -1, ('Q', 'S') :  0, ('E', 'S') :  0, ('G', 'S') :  0, ('H', 'S') : -1, ('I', 'S') : -2, ('L', 'S') : -2, ('K', 'S') :  0, ('M', 'S') : -1, ('F', 'S') : -2, ('P', 'S') : -1, ('S', 'S') :  4, ('T', 'S') :  1, ('W', 'S') : -3, ('Y', 'S') : -2, ('V', 'S') : -2,
        ('A', 'T') :  0, ('R', 'T') : -1, ('N', 'T') :  0, ('D', 'T') : -1, ('C', 'T') : -1, ('Q', 'T') : -1, ('E', 'T') : -1, ('G', 'T') : -2, ('H', 'T') : -2, ('I', 'T') : -1, ('L', 'T') : -1, ('K', 'T') : -1, ('M', 'T') : -1, ('F', 'T') : -2, ('P', 'T') : -1, ('S', 'T') :  1, ('T', 'T') :  5, ('W', 'T') : -2, ('Y', 'T') : -2, ('V', 'T') :  0,
        ('A', 'W') : -3, ('R', 'W') : -3, ('N', 'W') : -4, ('D', 'W') : -4, ('C', 'W') : -2, ('Q', 'W') : -2, ('E', 'W') : -3, ('G', 'W') : -2, ('H', 'W') : -2, ('I', 'W') : -3, ('L', 'W') : -2, ('K', 'W') : -3, ('M', 'W') : -1, ('F', 'W') :  1, ('P', 'W') : -4, ('S', 'W') : -3, ('T', 'W') : -2, ('W', 'W') : 11, ('Y', 'W') :  2, ('V', 'W') : -3,
        ('A', 'Y') : -2, ('R', 'Y') : -2, ('N', 'Y') : -2, ('D', 'Y') : -3, ('C', 'Y') : -2, ('Q', 'Y') : -1, ('E', 'Y') : -2, ('G', 'Y') : -3, ('H', 'Y') :  2, ('I', 'Y') : -1, ('L', 'Y') : -1, ('K', 'Y') : -2, ('M', 'Y') : -1, ('F', 'Y') :  3, ('P', 'Y') : -3, ('S', 'Y') : -2, ('T', 'Y') : -2, ('W', 'Y') :  2, ('Y', 'Y') :  7, ('V', 'Y') : -1,
        ('A', 'V') :  0, ('R', 'V') : -3, ('N', 'V') : -3, ('D', 'V') : -3, ('C', 'V') : -1, ('Q', 'V') : -2, ('E', 'V') : -2, ('G', 'V') : -3, ('H', 'V') : -3, ('I', 'V') :  3, ('L', 'V') :  1, ('K', 'V') : -2, ('M', 'V') :  1, ('F', 'V') : -1, ('P', 'V') : -2, ('S', 'V') : -2, ('T', 'V') :  0, ('W', 'V') : -3, ('Y', 'V') : -1, ('V', 'V') :  4,
}


class WildtypeVsMutant (Wizard):
    """
    Find the residues that differ between two selections and allow the user to
    view each difference one at a time.
    """

    def __init__(self, wildtype_obj='', mutant_obj='', focus_sele=''):
        self.mutations = []
        self.active_mutations = set()
        self.aligned_seqs = '', ''
        self.aligned_resis = [], []
        self.extra_positions = set()
        self.neighbor_radius = 4
        self.zoom_padding = 5
        self.flip_dist_cutoff = 4.0  # 2.0 counts rotating aromatic rings 180°
        self.wildtype_hilite = 'white'
        self.mutant_hilite = 'yellow'
        self.show_polar_h = False
        self.mutation_mode = 'mutations'  # also allowed: 'muts+flips'
        self.active_prompt = ''
        self.original_view = cmd.get_view()
        self.wildtype_obj = ''
        self.mutant_obj = ''
        self.focus_sele = ''

        self.set_focus_sele(focus_sele)
        self.set_wildtype_object(wildtype_obj)
        self.set_mutant_object(mutant_obj)

    def set_wildtype_object(self, wildtype_obj):
        """
        Specify which object should be considered the "wildtype".
        """
        self.wildtype_obj = wildtype_obj
        self.active_prompt = ''
        self.update_mutation_list()

    def set_mutant_object(self, mutant_obj):
        """
        Specify which object should be considered the "mutant".
        """
        self.mutant_obj = mutant_obj
        self.active_prompt = ''
        self.update_mutation_list()

    def set_focus_sele(self, focus_sele):
        """
        Specify which atoms the user is interested in comparing.  For example, 
        if the given selection comprises the interface between two proteins, 
        the rest of the wizard will only consider mutations in that region.
        """
        self.focus_sele = focus_sele
        self.active_prompt = ''
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

        # mutations: A set of indices of positions that differ between
        # the two aligned sequences, excluding terminal gaps.  These indices
        # can be used with both self.aligned_seqs and self.aligned_resis.

        def in_focus_sele(i): #
            if not self.focus_sele:
                return True
            else:
                return (self.aligned_resis[0][i] in focus_resis or
                        self.aligned_resis[1][i] in focus_resis)

        mutations = {
                i for i in find_mutations(self.aligned_seqs)
                if in_focus_sele(i)
        }

        # If the user requested include sidechains that are packed very 
        # differently in the list of mutations (i.e. "Include: muts+flips"), 
        # find those positions and include them (if they are within the focus 
        # selection).

        if self.mutation_mode == 'muts+flips':
            mutations |= {
                    i for i in self.find_flips()
                    if in_focus_sele(i)
            }

        # self.mutations: A list of indices of positions that should be 
        # displayed for the user to inspect.  This could include any positions 
        # that have mutated, have undergone sidechain flips, or have just been 
        # selected by the user.  These indices can be used with both 
        # self.aligned_seqs and self.aligned_resis.

        old_mutations = self.mutations
        self.mutations = sorted(mutations | self.extra_positions)

        # Make sure that all the residues in self.active_residues are still in 
        # self.mutations.  If are are not, replace them with the next residue 
        # that is still active. 

        old_active_mutations = self.active_mutations
        self.active_mutations = set()

        for muti in old_active_mutations:
            if muti in self.mutations:
                self.active_mutations.add(muti)
            else:
                for mutj in self.mutations:
                    if mutj > muti:
                        self.active_mutations.add(mutj)
                        break

        # If there aren't any active mutations, make the first one active.

        if not self.active_mutations and self.mutations:
            self.active_mutations = {self.mutations[0]}

        self.redraw()

    def find_flips(self):
        """
        Yield the indices of any positions where the wildtype and mut 
        sidechains have at least one atom further apart than the distance 
        threshold.
        """
        from math import sqrt

        for i in range(len(self.aligned_resis[0])):
            wt_seq = self.aligned_seqs[0][i]
            mut_seq = self.aligned_seqs[1][i]

            if wt_seq != mut_seq: continue
            if wt_seq in 'X-' or mut_seq in 'X-': continue

            # I wanted to use `cmd.rms_cur()` to calculate this, but it refuses 
            # to work on residues with different numbers, which I often have.

            wt_xyzs = {}
            mut_xyzs = {}
            sele = '({}) and resi {} and chain {} and not hydro'

            cmd.iterate_state(1,
                    sele.format(self.wildtype_obj, *self.aligned_resis[0][i]),
                    'wt_xyzs[name] = x,y,z',
                    space=locals(),
            )
            cmd.iterate_state(1,
                    sele.format(self.mutant_obj, *self.aligned_resis[1][i]),
                    'mut_xyzs[name] = x,y,z',
                    space=locals(),
            )

            # If we have two of the same residues, the heavy-atoms should have 
            # the same names.
            assert wt_xyzs.keys() == mut_xyzs.keys()

            max_dist = 0
            for k in wt_xyzs:
                dist = sqrt(sum(
                    (wt_xyzs[k][j] - mut_xyzs[k][j])**2
                    for j in range(3)
                ))
                max_dist = max(dist, max_dist)

            if max_dist > self.flip_dist_cutoff:
                yield i

    def set_active_mutation(self, index):
        """
        Set the mutation being currently highlighted.
        """
        self.active_mutations = {index}
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

    def set_mutation_mode(self, mode):
        """
        Change which residues are included in the mutation list (e.g. just 
        mutations or mutations and rotamer flips).
        """
        self.mutation_mode = mode
        self.update_mutation_list()

    def get_next_mutation(self):
        """
        Return the index of the next mutation.  The index is relative to the
        sequence alignment. Resets to first mutation if multiple mutations are 
        active.
        """
        if not self.mutations:
            raise IndexError

        if len(self.active_mutations) == 1:
            i = self.mutations.index(next(iter(self.active_mutations)))
            return self.mutations[i + 1]
        else:
            return self.mutations[0]

    def get_mutation_name(self, muti=None):
        """
        Return the name of given mutation.  The `muti' parameter is an index
        into the sequence alignment that defaults to the currently active
        mutation.  The returned name will be formatted like so:

        Mutation:    D38E
        Insertion:   -38E
        Deletion:    D38-
        Unchanged:   D38
        """
        if muti is None:
            mutis = sorted(self.active_mutations)
        else:
            mutis = [muti]

        names = []
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

            if wildtype_res != mutant_res:
                name = '{}{}{}'.format(wildtype_res, resi, mutant_res)
            else:
                name = '{}{}'.format(wildtype_res, resi)

            names.append(name)

        # Return the concatenated name.

        return ', '.join(names)

    def get_panel(self):
        """
        Return a list of lists describing the entries that should appear in the 
        right-hand panel of the GUI (below all the selections) when the wizard 
        is running.  The first few entries control common actions and basic 
        settings.  The remaining entries are buttons that the user can press to 
        see specific mutations.
        """

        # Each entry is described by a list with 3 elements.  The first element 
        # is a number that specifies the type of entry: 1 for text, 2 for a 
        # button, 3 for a menu.  The second element is the text that will be 
        # seen by the user.  The third element is an argument that means 
        # different things for each type of entry.  For buttons, it's a piece 
        # of code (given as a string) that should be executed when the button 
        # is pressed.  For menus, it's a "tag" that can be passed to get_menu()
        # to get a description of the menu in question.

        # If the user hasn't provided wildtype and mutant objects yet, then 
        # get_prompt() should be leading them through the process of doing 
        # that.  Until then, just display the name of the wizard and the option 
        # to quit.

        if not self.mutant_obj or not self.wildtype_obj:
            return [
                [1, "Wildtype vs Mutant Wizard", ''],
                [2, "Cancel", 'cmd.get_wizard().cleanup()'],
            ]

        # Make the buttons and menus that control the basic actions and 
        # settings that don't depend on the specific system being studied.

        buttons = [
            [1, "Wildtype vs Mutant Wizard", ''],
            [2, "Show next mutation", 'cmd.get_wizard().show_next_mutation()'],
            [2, "Show all mutations", 'cmd.get_wizard().show_all_mutations()'],
            [2, "Add (sele) positions", 'cmd.get_wizard().add_selected_positions()'],
            [3, "Wildtype highlight: {}".format(self.wildtype_hilite), 'wt_hilite'],
            [3, "Mutant highlight: {}".format(self.mutant_hilite), 'mut_hilite'],
            [3, "Neighbor radius: {0:.1f}A".format(self.neighbor_radius), 'radius'],
            [3, "Zoom padding: {0:.1f}A".format(self.zoom_padding), 'padding'],
            [3, "Polar hydrogens: {}".format('show' if self.show_polar_h else 'hide'), 'hydrogen'],
            [3, "Include: {}".format(self.mutation_mode), 'mutation_mode'],
        ]

        # Make a button for each mutation.  The user can click on these buttons 
        # to view specific mutations.  The button text turns green if the user 
        # is currently viewing the associated mutation.

        for muti in self.mutations:
            command = 'cmd.get_wizard().set_active_mutation({})'
            name = self.get_mutation_name(muti)
            if muti in self.active_mutations:
                name = '\\090' + name
            buttons += [[2, name, command.format(muti)]]

        if not self.mutations:
            from textwrap import wrap
            buttons += [
                [1, '\\900' + x, '']
                for x in wrap("No mutations found.  Do `{0}' and `{1}' have the same sequence?".format(self.wildtype_obj, self.mutant_obj), 26)
            ]

        # Make a quit button and return the buttons data structure.

        buttons += [[2, "Done", 'cmd.get_wizard().cleanup()']]
        return buttons

    def get_menu(self, tag):
        """
        Return a dictionary describing the entries in the menu associated with 
        the given tag.  These tags come from get_panel(), so each menu is 
        associated with one buttons on the right-hand side of the GUI.
        """
        menus = {
            'wt_hilite': [[2, 'Highlight Color', '']],
            'mut_hilite': [[2, 'Highlight Color', '']],
            'radius': [[2, 'Neighbor Radius', '']],
            'padding': [[2, 'Zoom Padding', '']],
            'hydrogen': [[2, 'Polar Hydrogens', '']],
            'mutation_mode': [[2, 'Include', '']],
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
        for radius in range(11):
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

        # Define the mutation mode menu.
        menus['mutation_mode'] += [[1, 'mutations', 'cmd.get_wizard().set_mutation_mode("mutations")']]
        menus['mutation_mode'] += [[1, 'muts+flips', 'cmd.get_wizard().set_mutation_mode("muts+flips")']]

        # Return the right menu.
        return menus[tag]

    def get_prompt(self):
        """
        Return text to be displayed in the top left corner of the view area.
        If the user has not yet provided wildtype and mutant objects, prompt 
        for that.  Otherwise, tell the user about the <Ctrl-Space> hotkey.
        """
        # The \999 code changes the text color to white.
        if not self.wildtype_obj:
            return ["Select the wildtype object: \\999{}".format(self.active_prompt)]
        elif not self.mutant_obj:
            return ["Select the mutant object: \\999{}".format(self.active_prompt)]
        elif not self.mutations:
            return
        else:
            return ["Press <Ctrl-Space> to view the next mutation..."]

    def do_key(self, key, x, y, mod):
        """
        Take responsibility for handling key presses.
        """
        # If the user presses Ctrl-Space (key, mod == 0, 2), show the next
        # mutation.  Otherwise let pymol handle the key press.

        if self.wildtype_obj and self.mutant_obj:
            if (key, mod) == (0, 2):
                self.show_next_mutation()
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

    def show_next_mutation(self):
        """
        Focus on the next mutant, or close the wizard if the last mutant is
        currently active.
        """
        try:
            self.set_active_mutation(self.get_next_mutation())
        except IndexError:
            self.cleanup() if self.mutations else self.redraw()

    def show_all_mutations(self):
        """
        Show all the mutations simultaneously.
        """
        self.active_mutations = set(self.mutations)
        self.redraw()

    def add_selected_positions(self, selection='sele'):
        """
        Add any positions currently in the given selection to the list of 
        mutations to show (even if those positions are not in fact mutations).
        """
        wt_obj = self.wildtype_obj
        mut_obj = self.mutant_obj

        if 'wt_env' in cmd.get_names("selections"):
            wt_obj = '({}) or wt_env'.format(wt_obj)
        if 'mut_env' in cmd.get_names("selections"):
            mut_obj = '({}) or mut_env'.format(mut_obj)

        wt_resis = get_residues('({}) and ({})'.format(selection, wt_obj))
        mut_resis = get_residues('({}) and ({})'.format(selection, mut_obj))
        
        for resi in wt_resis:
            muti = self.aligned_resis[0].index(resi)
            self.extra_positions.add(muti)

        for resi in mut_resis:
            muti = self.aligned_resis[1].index(resi)
            self.extra_positions.add(muti)

        self.update_mutation_list()
        self.redraw()

    def redraw(self):
        """
        Highlight the sidechains around the active mutation(s).  Everything 
        this method adds to the scene is put in its own object, so that it can 
        be easily undone by the next call to redraw() or cleanup().
        """
        cmd.refresh_wizard()
        if not self.active_mutations:
            return

        wt_obj = self.wildtype_obj
        mut_obj = self.mutant_obj

        # wt_seles and mut_seles: lists containing individual pymol selection 
        # expressions for each residue that the user wants to see.

        wt_seles = []
        mut_seles = []

        for muti in self.active_mutations:
            wt_resi, wt_chain = self.aligned_resis[0][muti]
            mut_resi, mut_chain = self.aligned_resis[1][muti]
            if wt_resi is not None:
                wt_seles.append('(resi {wt_resi} and chain {wt_chain})'.format(**locals()))
            if mut_resi is not None:
                mut_seles.append('(resi {mut_resi} and chain {mut_chain})'.format(**locals()))

        # wt_sele and mut_sele: Pymol selection expressions that combine all 
        # the individual selections in wt_seles and mut_seles.

        wt_sele = '({})'.format(' or '.join(wt_seles)) if wt_seles else 'none'
        mut_sele = '({})'.format(' or '.join(mut_seles)) if mut_seles else 'none'

        # env_sele: A pymol selection expression that specifies all which 
        # residues are close enough to the mutation(s) to be considered part of 
        # the "environment".  This selection is actually a template, and needs 
        # to be formatted with either wt_obj or mut_obj before use.

        h_sele = (
                '(elem H and (neighbor elem C))' if self.show_polar_h else
                '(elem H)')
        env_sele = (
                '(byres {{}} within {self.neighbor_radius} of '
                  '(({wt_obj} and {wt_sele}) or ({mut_obj} and {mut_sele})))'
                'and not {h_sele}'.format(**locals()))

        # Render the scene into the wt_env, mut_env, wt_hbonds, and mut_hbonds 
        # objects.  Limiting ourselves to these objects makes it easy to redraw 
        # our comparisons without getting in the user's way, or vice versa.

        initial_view = cmd.get_view()

        cmd.delete('wt_env')
        cmd.delete('mut_env')
        cmd.delete('wt_hbonds')
        cmd.delete('mut_hbonds')

        cmd.create('wt_env', env_sele.format(wt_obj))
        cmd.create('mut_env', env_sele.format(mut_obj))

        # Hide the wildtype polar contacts by default, because I'm more often 
        # interested in the contacts being made in the design.

        cmd.distance('wt_hbonds', 'wt_env', 'wt_env', mode=2)
        cmd.distance('mut_hbonds', 'mut_env', 'mut_env', mode=2)
        cmd.disable('wt_hbonds')

        cmd.show_as('sticks', 'wt_env')
        cmd.show_as('sticks', 'mut_env')

        if self.wildtype_hilite != 'none':
            cmd.color(self.wildtype_hilite, 'wt_env and {wt_sele} and elem C'.format(**locals()))
        if self.mutant_hilite != 'none':
            cmd.color(self.mutant_hilite, 'mut_env and {mut_sele} and elem C'.format(**locals()))

        # Pan smoothly to the new scene.

        cmd.set_view(initial_view)
        cmd.zoom('mut_env', buffer=self.zoom_padding, animate=-1)

    def cleanup(self):
        """
        Remove any changes this wizard made to the user's session.
        """
        cmd.delete('wt_env')
        cmd.delete('mut_env')
        cmd.delete('wt_hbonds')
        cmd.delete('mut_hbonds')

        cmd.set_view(self.original_view)
        cmd.set_wizard()



def wt_vs_mut(mutant_obj=None, wildtype_obj=None, focus_sele=None):
    """
    Provide a convenient way to launch the wizard from the command-line.
    """
    wizard = WildtypeVsMutant(mutant_obj, wildtype_obj, focus_sele)
    cmd.set_wizard(wizard)

def get_residues(selection):
    """
    Return a set of (residue number, chain id) tuples representing all the 
    residues contained in the given selection.
    """
    residues = set()
    int_from_str = int
    cmd.iterate(selection, 'residues.add((int_from_str(resi), chain))', space=locals())
    return residues

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

def get_score_from_matrix(pos1, pos2, score_matrix=blosum_62):
    """
    Return score from score_matrix if it exists. Otherwise (e.g. for 
    non-canonicals), return generic score values for matches and 
    mismatches.
    """
    try:
        return score_matrix[pos1, pos2]
    except KeyError:
        return 5 if pos1 == pos2 else -4


## Add "wt_vs_mut" as pymol command
cmd.extend('wt_vs_mut', wt_vs_mut)
cmd.auto_arg[0]['wt_vs_mut'] = cmd.auto_arg[0]['zoom']
cmd.auto_arg[1]['wt_vs_mut'] = cmd.auto_arg[0]['zoom']

## Trick to get "wizard wt_vs_mut" working
sys.modules['pymol.wizard.wt_vs_mut'] = sys.modules[__name__]

## Add item to plugin menu
try:
    from pymol.plugins import addmenuitem
    def __init_plugin__(self): #
        addmenuitem('Wildtype vs. Mutant', lambda s=self: wt_vs_mut())
except:
    def __init__(self): #
        self.menuBar.addmenuitem(
                'Plugin', 'command', 'Wildtype vs. Mutant',
                label='Wildtype vs. Mutant', command=lambda s=self: wt_vs_mut())
