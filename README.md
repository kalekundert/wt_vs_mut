Wildtype vs. Mutant Pymol Plugin
================================
The wt_vs_mut plugin makes it easy to visualize the difference between two 
similar protein structures.  The plugin automatically finds all the residues 
that differ between the two structures, then it lets the user zoom in on the 
interactions being made by those residues one at a time.  The purpose of this 
is to let the user judge the effect that each mutation might have on one 
structure or the other.  This plugin was originally written to help with the 
manual validation step of a protein design project, but is quite general and 
should be useful in many other contexts as well.

Installation and Usage
----------------------
Start by downloading the `wt_vs_mut.py` file somewhere on your computer.  The 
entire plugin is contained in this file to make installation easier.  Then, in 
pymol, open the two structures you want to compare.  I'll assume one is called 
`wildtype` and the other is called `mutant`.  Enter the following two commands:

    run /path/to/wt_vs_mut.py
    wt_vs_mut wildtype, mutant

Of course, `path/to/wt_vs_mut.py` is the path to wherever you downloaded the 
plugin, relative to wherever pymol is running.  The first command installs the 
plugin.  You'll have to type it every time you run pymol unless you put it in 
your `~/.pymolrc` file.  The second command starts the plugin.  Once the plugin 
is started, it will zoom in on the first difference between the two structures 
and show all the sidechains within a few angstroms.  Press [Ctrl-Space] and it 
will zoom to the next difference and show the same thing.

Note that this plugin requires at least version 1.7 of pymol.  Older versions 
don't have a recent enough distribution of python.  If you run into other 
troubles or would like to suggest a new feature, feel free to submit an issue 
through the GitHub GUI.

