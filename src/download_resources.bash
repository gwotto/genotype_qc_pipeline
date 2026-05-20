#!/usr/bin/env bash
set -euo pipefail

curl -L --fail -o high-LD-regions.txt https://raw.githubusercontent.com/meyer-lab-cshl/plinkQC/master/inst/extdata/high-LD-regions-hg19-GRCh37.txt

curl -L --fail -o HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz

echo "Downloads complete"
