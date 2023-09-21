#!/bin/bash

# Makes permanent dipoles
~/Multiwfn_3.7_bin_Linux_noGUI/Multiwfn << EOF
geometry.fchk
18
5
geometry.out
2
EOF

# Makes transition dipoles
~/Multiwfn_3.7_bin_Linux_noGUI/Multiwfn << EOF
geometry.fchk
18
5
geometry.out
4
EOF

