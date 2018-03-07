## Wrappers for MetaPhlAn2 + HUMAnN2

These scripts are simple wrappers for MetaPhlAn2 and HUMAnN2.

run_metaphlan2.py will quickly profile a metagenome and generate
preliminary visualizations with GraPhlAn.

run_humann2.py will run the full HUMAnN2 pipeline on a given sample.

#### Requirements:
- Python==2.7
- MetaPhlAn2==2.6.0
- [GraPhlAn](https://huttenhower.sph.harvard.edu/graphlan)
- [BBMap](https://sourceforge.net/projects/bbmap/)
- [HUMAnN2](http://huttenhower.sph.harvard.edu/humann2)

#### Notes:
The following must be available in your *PATH*:
- bbmerge.sh
- metaphlan2.py
- export2graphlan.py
- graphlan_annotate.py
- graphlan.py
- humann2

#### Usage:
```bash
run_metaphlan2.py [your_file.fastq.gz]
```