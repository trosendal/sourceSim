# Determine package name and version from DESCRIPTION file
PKG_VERSION=$(shell grep -i ^version DESCRIPTION | cut -d : -d \  -f 2)
PKG_NAME=$(shell grep -i ^package DESCRIPTION | cut -d : -d \  -f 2)

# Name of built package
PKG_TAR=$(PKG_NAME)_$(PKG_VERSION).tar.gz

all: install

# Install package
install:
	cd .. && R CMD build $(PKG_NAME)
	cd .. && R CMD INSTALL $(PKG_TAR)

# Build documentation with roxygen
# 1) Remove old doc
# 2) Generate documentation
roxygen:
	rm -f man/*.Rd
	cd .. && Rscript -e "library(roxygen2); roxygenize('$(PKG_NAME)')"

# Generate PDF output from the Rd sources
# 1) Rebuild documentation with roxygen
# 2) Generate pdf, overwrites output file if it exists
pdf: roxygen
	cd .. && R CMD Rd2pdf --force $(PKG_NAME)

# Build package
build:
	cd .. && R CMD build --compact-vignettes=both $(PKG_NAME)

# Check package
check: build
	cd .. && OMP_THREAD_LIMIT=2 _R_CHECK_CRAN_INCOMING_=FALSE R CMD check \
        --no-stop-on-test-error --as-cran --run-dontrun $(PKG_TAR)

# Check package (without manual and vignettes)
check_quick:
	cd .. && R CMD build --no-build-vignettes --no-manual $(PKG_NAME)
	cd .. && \
        _R_CHECK_CRAN_INCOMING_=FALSE \
        R CMD check \
        --no-stop-on-test-error --no-vignettes --no-manual --as-cran $(PKG_TAR)

.PHONY: install roxygen pdf check check_quick build all
