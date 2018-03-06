## Wrapper for MetaPhlAn2

This is a simple wrapper for MetaPhlAn2 to quickly profile a metagenome
and generate some preliminary visualizations with GraPhlAn.

#### Requirements:
- [MetaPhlAn2](https://bitbucket.org/biobakery/metaphlan2#markdown-header-installation)
- [GraPhlAn](https://huttenhower.sph.harvard.edu/graphlan)
- [BBMap](https://sourceforge.net/projects/bbmap/)

#### Notes:
The following must be available in your *PATH*:
- bbmerge.sh
- metaphlan2.py
- export2graphlan.py
- graphlan_annotate.py
- graphlan.py

#### Usage:
```bash
run_metaphlan2.py [your_file.fastq.gz]
```