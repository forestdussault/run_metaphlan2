from __future__ import print_function
import subprocess


def create_cladogram(overall_abundance_profile, workdir):

    # Step 1: Create GraPhlAn input files
    print("\nCreating GraPhlAn input files...")

    annotated_file = overall_abundance_profile.replace('profile', 'profile.annot')

    p = subprocess.Popen('python2 '
                         'export2graphlan.py '
                         '--skip_rows 1,2 '
                         '-i {} '
                         '--tree merged_abundance.tree.txt '
                         '--annotation {} '
                         '--most_abundant 100 '
                         '--abundance_threshold 1 '
                         '--least_biomarkers 10 '
                         '--annotations 5,6 '
                         '--external_annotations 7 '
                         '--min_clade_size 1'.format(overall_abundance_profile, annotated_file),
                         shell=True,
                         cwd=workdir,
                         executable='/bin/bash')
    p.wait()

    # Step 2: Create cladogram pieces
    print("\nCreating cladogram input files...")

    p = subprocess.Popen('python2 '
                         'graphlan_annotate.py '
                         '--annot {} '
                         'merged_abundance.tree.txt '
                         'merged_abundance.xml'.format(annotated_file),
                         shell=True,
                         cwd=workdir,
                         executable='/bin/bash')
    p.wait()

    # Step 3: Visualize cladogram
    print("\nVisualizing cladogram...")
    p = subprocess.Popen('python2 '
                         'graphlan.py '
                         '--dpi 300 '
                         '--pad 2 '
                         'merged_abundance.xml '
                         'merged_abundance.png ',
                         shell=True,
                         cwd=workdir,
                         executable='/bin/bash')
    p.wait()
