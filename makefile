all: doc test

doc: doc/manual.six

doc/manual.six: makedoc.g \
		PackageInfo.g \
		doc/PrimaryDecomposition.bib \
		gap/*.gd gap/*.gi examples/*.g examples/doc/*.g
	        gap makedoc.g

clean:
	(cd doc ; ./clean)

test:	doc
	gap maketest.g

archive: test
	(mkdir -p ../tar; cd ..; tar czvf tar/PrimaryDecomposition.tar.gz --exclude ".DS_Store" --exclude "*~" PrimaryDecomposition/doc/*.* PrimaryDecomposition/doc/clean PrimaryDecomposition/gap/*.{gi,gd} PrimaryDecomposition/{PackageInfo.g,README,COPYING,VERSION,init.g,read.g,makedoc.g,makefile,maketest.g} PrimaryDecomposition/examples/*.g PrimaryDecomposition/examples/doc/*.g)

WEBPOS=public_html
WEBPOS_FINAL=~/Sites/homalg-project/PrimaryDecomposition

towww: archive
	echo '<?xml version="1.0" encoding="UTF-8"?>' >${WEBPOS}.version
	echo '<mixer>' >>${WEBPOS}.version
	cat VERSION >>${WEBPOS}.version
	echo '</mixer>' >>${WEBPOS}.version
	cp PackageInfo.g ${WEBPOS}
	cp README ${WEBPOS}/README.PrimaryDecomposition
	cp doc/manual.pdf ${WEBPOS}/PrimaryDecomposition.pdf
	cp doc/*.{css,html} ${WEBPOS}
	rm -f ${WEBPOS}/*.tar.gz
	mv ../tar/PrimaryDecomposition.tar.gz ${WEBPOS}/PrimaryDecomposition-`cat VERSION`.tar.gz
	rm -f ${WEBPOS_FINAL}/*.tar.gz
	cp ${WEBPOS}/* ${WEBPOS_FINAL}
	ln -s PrimaryDecomposition-`cat VERSION`.tar.gz ${WEBPOS_FINAL}/PrimaryDecomposition.tar.gz

