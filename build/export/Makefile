# Locations to save compressed tar file for distribution

CHARLTON:=/var/anonymous/ftp/charlton
FTP_CHARLTON:=ftp://brrcrftp.cr.usgs.gov/pub/charlton
GEOCHEM:=/var/anonymous/ftp/dlpark/geochem
FTP_GEOCHEM:=ftp://brrcrftp.cr.usgs.gov/geochem
HTML_HOME:=/z/linarcolkr/home/www/projects/GWC_coupled
HTML_PHREEQC:=$(HTML_HOME)/phreeqc
HTML_PHREEQCI:=$(HTML_HOME)/phreeqci

#PhreeqcI version

PHREEQCI_DIR:=$(CHARLTON)/phreeqci
PHREEQCI_TARGET_DIR:=$(GEOCHEM)/pc/phreeqc
PHREEQCI_PUB_DIR:=$(FTP_GEOCHEM)/pc/phreeqc
PHREEQCI_FILE:=$(shell cd "$(PHREEQCI_DIR)"; ls -t phreeqci*.exe | awk '{if (NR == 1) {print}}')
PHREEQCI_VERSION:=$(shell echo $(PHREEQCI_FILE) | sed "s/phreeqci-//" | sed "s/.exe//" )
PHREEQCI_SIZE:=$(shell ls -s -h "$(PHREEQCI_DIR)/$(PHREEQCI_FILE)" | awk '{print $$1}')
PHREEQCI_DATE:=$(shell ls -l --time-style "+%B %e, %G" "$(PHREEQCI_DIR)/$(PHREEQCI_FILE)" | awk '{print $$6,$$7,$$8}')

#Phreeqc For Windows version

PFW_DIR:=$(GEOCHEM)/pc/phreeqc
PFW_PUB_DIR:=$(FTP_GEOCHEM)/pc/phreeqc
PFW_FILE:=$(shell cd "$(PFW_DIR)"; ls -t psetup*.exe | awk '{if (NR == 1) {print}}')
PFW_VERSION:=$(shell echo $(PFW_FILE) | sed "s/psetup//" | sed "s/.exe//" | awk '{print substr($$1,1,1)"."substr($$1,2,2)"."substr($$1,4)}')
PFW_SIZE:=$(shell ls -s -h "$(PFW_DIR)/$(PFW_FILE)" | awk '{print $$1}')
PFW_DATE:=$(shell ls -l --time-style "+%B %e, %G" "$(PFW_DIR)/$(PFW_FILE)" | awk '{print $$6,$$7,$$8}')

#Batch PHREEQC versions



#Batch Windows

PHREEQC_WIN_DIR:=$(CHARLTON)/phreeqc
PHREEQC_WIN_TARGET_DIR:=$(GEOCHEM)/pc/phreeqc
PHREEQC_WIN_PUB_DIR:=$(FTP_GEOCHEM)/pc/phreeqc
PHREEQC_WIN_FILE:=$(shell cd "$(PHREEQC_WIN_DIR)"; ls -t phreeqc*.exe | awk '{if (NR == 1) {print}}')
PHREEQC_WIN_VERSION:=$(shell echo $(PHREEQC_WIN_FILE) | sed "s/phreeqc-//" | sed "s/.exe//" )
PHREEQC_WIN_SIZE:=$(shell ls -s -h "$(PHREEQC_WIN_DIR)/$(PHREEQC_WIN_FILE)" | awk '{print $$1}')
PHREEQC_WIN_DATE:=$(shell ls -l --time-style "+%B %e, %G" "$(PHREEQC_WIN_DIR)/$(PHREEQC_WIN_FILE)" | awk '{print $$6,$$7,$$8}')

#Batch Linux

PHREEQC_LINUX_DIR:=$(CHARLTON)/phreeqc
PHREEQC_LINUX_TARGET_DIR:=$(GEOCHEM)/unix/phreeqc
PHREEQC_LINUX_PUB_DIR:=$(FTP_GEOCHEM)/unix/phreeqc
PHREEQC_LINUX_FILE:=$(shell cd "$(PHREEQC_LINUX_DIR)"; ls -t phreeqc-*Linux.tar.gz | awk '{if (NR == 1) {print}}')
PHREEQC_LINUX_VERSION:=$(shell echo $(PHREEQC_LINUX_FILE) | sed "s/phreeqc-//" | sed "s/.Linux.tar.gz//" )
PHREEQC_LINUX_SIZE:=$(shell ls -s -h "$(PHREEQC_LINUX_DIR)/$(PHREEQC_LINUX_FILE)" | awk '{print $$1}')
PHREEQC_LINUX_DATE:=$(shell ls -l --time-style "+%B %e, %G" "$(PHREEQC_LINUX_DIR)/$(PHREEQC_LINUX_FILE)" | awk '{print $$6,$$7,$$8}')

#Batch source

PHREEQC_SOURCE_DIR:=$(CHARLTON)/phreeqc
PHREEQC_SOURCE_TARGET_DIR:=$(GEOCHEM)/unix/phreeqc
PHREEQC_SOURCE_PUB_DIR:=$(FTP_GEOCHEM)/unix/phreeqc
PHREEQC_SOURCE_FILE:=$(shell cd "$(PHREEQC_SOURCE_DIR)"; ls -t phreeqc-*.source.tar.gz | awk '{if (NR == 1) {print}}')
PHREEQC_SOURCE_VERSION:=$(shell echo $(PHREEQC_SOURCE_FILE) | sed "s/phreeqc-//" | sed "s/.source.tar.gz//" )
PHREEQC_SOURCE_SIZE:=$(shell ls -s -h "$(PHREEQC_SOURCE_DIR)/$(PHREEQC_SOURCE_FILE)" | awk '{print $$1}')
PHREEQC_SOURCE_DATE:=$(shell ls -l --time-style "+%B %e, %G" "$(PHREEQC_SOURCE_DIR)/$(PHREEQC_SOURCE_FILE)" | awk '{print $$6,$$7,$$8}')

#Batch Mac

PHREEQC_MAC_DIR:=$(GEOCHEM)/mac/phreeqc
PHREEQC_MAC_PUB_DIR:=$(FTP_GEOCHEM)/mac/phreeqc
PHREEQC_MAC_FILE:=$(shell cd "$(PHREEQC_MAC_DIR)"; ls -t phreeqc-*.dmg.gz | awk '{if (NR == 1) {print}}')
PHREEQC_MAC_VERSION:=$(shell echo $(PHREEQC_MAC_FILE) | sed "s/phreeqc-//" | sed "s/.dmg.gz//" )
PHREEQC_MAC_SIZE:=$(shell ls -s -h "$(PHREEQC_MAC_DIR)/$(PHREEQC_MAC_FILE)" | awk '{print $$1}')
PHREEQC_MAC_DATE:=$(shell ls -l --time-style "+%B %e, %G" "$(PHREEQC_MAC_DIR)/$(PHREEQC_MAC_FILE)" | awk '{print $$6,$$7,$$8}')

SED_ARGS= \
	-e "s^@PHREEQCI_PUB_DIR@^$(PHREEQCI_PUB_DIR)^"  \
	-e "s^@PHREEQCI_FILE@^$(PHREEQCI_FILE)^" \
	-e "s^@PHREEQCI_VERSION@^$(PHREEQCI_VERSION)^" \
	-e "s^@PHREEQCI_DATE@^$(PHREEQCI_DATE)^" \
	-e "s^@PHREEQCI_SIZE@^$(PHREEQCI_SIZE)^" \
	-e "s^@PFW_VERSION@^$(PFW_VERSION)^" \
	-e "s^@PFW_DATE@^$(PFW_DATE)^" \
	-e "s^@PFW_SIZE@^$(PFW_SIZE)^" \
	-e "s^@PFW_PUB_DIR@^$(PFW_PUB_DIR)^" \
	-e "s^@PFW_FILE@^$(PFW_FILE)^" \
	-e "s^@PHREEQC_WIN_PUB_DIR@^$(PHREEQC_WIN_PUB_DIR)^" \
	-e "s^@PHREEQC_WIN_FILE@^$(PHREEQC_WIN_FILE)^" \
	-e "s^@PHREEQC_WIN_VERSION@^$(PHREEQC_WIN_VERSION)^" \
	-e "s^@PHREEQC_WIN_DATE@^$(PHREEQC_WIN_DATE)^" \
	-e "s^@PHREEQC_WIN_SIZE@^$(PHREEQC_WIN_SIZE)^" \
	-e "s^@PHREEQC_LINUX_PUB_DIR@^$(PHREEQC_LINUX_PUB_DIR)^" \
	-e "s^@PHREEQC_LINUX_FILE@^$(PHREEQC_LINUX_FILE)^" \
	-e "s^@PHREEQC_LINUX_VERSION@^$(PHREEQC_LINUX_VERSION)^" \
	-e "s^@PHREEQC_LINUX_DATE@^$(PHREEQC_LINUX_DATE)^" \
	-e "s^@PHREEQC_LINUX_SIZE@^$(PHREEQC_LINUX_SIZE)^" \
	-e "s^@PHREEQC_SOURCE_PUB_DIR@^$(PHREEQC_SOURCE_PUB_DIR)^" \
	-e "s^@PHREEQC_SOURCE_FILE@^$(PHREEQC_SOURCE_FILE)^" \
	-e "s^@PHREEQC_SOURCE_VERSION@^$(PHREEQC_SOURCE_VERSION)^" \
	-e "s^@PHREEQC_SOURCE_DATE@^$(PHREEQC_SOURCE_DATE)^" \
	-e "s^@PHREEQC_SOURCE_SIZE@^$(PHREEQC_SOURCE_SIZE)^" \
	-e "s^@PHREEQC_MAC_PUB_DIR@^$(PHREEQC_MAC_PUB_DIR)^" \
	-e "s^@PHREEQC_MAC_FILE@^$(PHREEQC_MAC_FILE)^" \
	-e "s^@PHREEQC_MAC_VERSION@^$(PHREEQC_MAC_VERSION)^" \
	-e "s^@PHREEQC_MAC_DATE@^$(PHREEQC_MAC_DATE)^" \
	-e "s^@PHREEQC_MAC_SIZE@^$(PHREEQC_MAC_SIZE)^"



web: 
	sed $(SED_ARGS) < template_phreeqc_index.html > web_phreeqc_index.html
	sed $(SED_ARGS) < template_phreeqci_index.html > web_phreeqci_index.html

	cp -f web_phreeqc_index.html $(HTML_PHREEQC)/index.html 
	cp -f web_phreeqci_index.html $(HTML_PHREEQCI)/index.html 
	cd $(HTML_PHREEQC); rcs -l index.html; ci -l -m"Latest" index.html
	cd $(HTML_PHREEQCI); rcs -l index.html; ci -l -m"Latest" index.html

	cp -f $(PHREEQCI_DIR)/$(PHREEQCI_FILE) $(PHREEQCI_TARGET_DIR)/$(PHREEQCI_FILE)
	cp -f $(PHREEQC_WIN_DIR)/$(PHREEQC_WIN_FILE) $(PHREEQC_WIN_TARGET_DIR)/$(PHREEQC_WIN_FILE)
	cp -f $(PHREEQC_LINUX_DIR)/$(PHREEQC_LINUX_FILE) $(PHREEQC_LINUX_TARGET_DIR)/$(PHREEQC_LINUX_FILE)
	cp -f $(PHREEQC_SOURCE_DIR)/$(PHREEQC_SOURCE_FILE) $(PHREEQC_SOURCE_TARGET_DIR)/$(PHREEQC_SOURCE_FILE)
	cp -f $(PHREEQC_WIN_DIR)/README.Win.TXT $(HTML_PHREEQC)
	cp -f $(PHREEQC_LINUX_DIR)/README.Unix.TXT $(HTML_PHREEQC)
	cp -f $(PHREEQC_MAC_DIR)/README.TXT $(HTML_PHREEQC)/README.Mac.TXT

	mkdir -p temp
	cp -f $(PHREEQC_LINUX_DIR)/$(PHREEQC_LINUX_FILE) temp
	cd temp; tar -xzf $(PHREEQC_LINUX_FILE); cp -f ./*/doc/RELEASE.TXT $(HTML_PHREEQC); cd ..; rm -rf temp


