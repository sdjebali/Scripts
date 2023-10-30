#!/bin/bash
set -Eexo pipefail

# 3way_venn.sh accepts 8 inputs
# - file for 1st set of the venn (string)
# - file for 2nd set of the venn (string)
# - file for 3rd set of the venn (string)
# - name of 1st set (string)
# - name of 2nd set (string)
# - name of 3rd set (string)
# - title for the venn (string)
# - file for venn plot (string)

# has been tested on genologin with R 3.3.3 and fragencode configuration
# module use /work/project/fragencode/.local/privatemodules/
# module load fragencode-2017.10.06

# pgm=/work/project/fragencode/tools/multi/Scripts/3way_venn.sh
# set1=/work/project/seqoccin/svdetection/nanopore/bos_taurus/compare.vcfs/vcfs.diff.inds/connected.components/sniffles.minimap.DEL.allrun.cc.of.indoffspring.txt 
# set2=/work/project/seqoccin/svdetection/nanopore/bos_taurus/compare.vcfs/vcfs.diff.inds/connected.components/sniffles.minimap.DEL.allrun.cc.of.indfather.txt 
# set3=/work/project/seqoccin/svdetection/nanopore/bos_taurus/compare.vcfs/vcfs.diff.inds/connected.components/sniffles.minimap.DEL.allrun.cc.of.indmother.txt 
# output=/work/project/seqoccin/svdetection/nanopore/bos_taurus/compare.vcfs/vcfs.diff.runs/venn.diagrams/sniffles.minimap.DEL.allrun.cc.3indiv.venn.png
# time $pgm $set1 $set2 $set3 offspring father mother "trio1-minimap+sniffles-DEL-allrun" $output

# on April 14th 2023 added 3 params for the titel

echo '
     library("futile.logger")
     library("VennDiagram")
     png("'$8'")
     run1=unlist(read.delim("'$1'", sep="\t", h=FALSE))
     run2=unlist(read.delim("'$2'", sep="\t", h=FALSE))
     run3=unlist(read.delim("'$3'", sep="\t", h=FALSE))
     venn.plot <- venn.diagram(
     x=list( 
     A=run1,
     B=run2,
     C=run3
     ),
	main="'$7'",   
	filename = NULL,
	fill=c("#FBB4AE", "#B3CDE3", "#CCEBC5"), 
	alpha=0.50,
	category = c("'$4'","'$5'","'$6'"), #legend
	fontfamily="arial",
	main.cex=2,
	main.fontfamily="arial",
	sub.cex=1.3,
	sub.fontfamily="arial",
	cat.fontfamily=c("arial","arial","arial") 
	)
	grid.draw(venn.plot)
	dev.off()
' | R --vanilla
