# Locations to save compressed tar file for distribution
ANONYMOUS:=/var/anonymous/ftp/pub
CHARLTON:=$(ANONYMOUS)/charlton
FTP_CHARLTON:=ftp://brrftp.cr.usgs.gov/pub/charlton
GEOCHEM:=$(ANONYMOUS)/dlpark/geochem
FTP_GEOCHEM:=ftp://brrftp.cr.usgs.gov/geochem
HTML_HOME:=/z/linarcolkr/home/www/projects/GWC_coupled
HTML_PHREEQC:=$(HTML_HOME)/phreeqc
HTML_PHREEQCI:=$(HTML_HOME)/phreeqci
HTML_PHAST:=$(HTML_HOME)/phast
HTML_NETPATH:=$(HTML_HOME)/netpath

#PhreeqcI version

PHREEQCI_DIR:=$(CHARLTON)/phreeqci
PHREEQCI_TARGET_DIR:=$(GEOCHEM)/pc/phreeqc
PHREEQCI_PUB_DIR:=$(FTP_CHARLTON)/phreeqci
PHREEQCI_FILE:=$(shell cd "$(PHREEQCI_DIR)"; ls -t phreeqci*.msi | awk '{if (NR == 1) {print}}')
PHREEQCI_VERSION:=$(shell echo $(PHREEQCI_FILE) | sed "s/phreeqci-//" | sed "s/.msi//" )
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
#PHREEQC_WIN_DIR:=$(GEOCHEM)/pc/phreeqc
PHREEQC_WIN_TARGET_DIR:=$(GEOCHEM)/pc/phreeqc
PHREEQC_WIN_PUB_DIR:=$(FTP_CHARLTON)/phreeqc
#PHREEQC_WIN_PUB_DIR:=$(FTP_GEOCHEM)/pc/phreeqc
PHREEQC_WIN_FILE:=$(shell cd "$(PHREEQC_WIN_DIR)"; ls -t phreeqc*.msi | awk '{if (NR == 1) {print}}')
PHREEQC_WIN_VERSION:=$(shell echo $(PHREEQC_WIN_FILE) | sed "s/phreeqc-//" | sed "s/.msi//" )
PHREEQC_WIN_SVN:=$(shell echo $(PHREEQC_WIN_FILE) | sed "s/phreeqc-//" | sed "s/.msi//" | sed "s/^.*-//" )
PHREEQC_WIN_SIZE:=$(shell ls -s -h "$(PHREEQC_WIN_DIR)/$(PHREEQC_WIN_FILE)" | awk '{print $$1}')
PHREEQC_WIN_DATE:=$(shell ls -l --time-style "+%B %e, %G" "$(PHREEQC_WIN_DIR)/$(PHREEQC_WIN_FILE)" | awk '{print $$6,$$7,$$8}')

#Batch Linux

PHREEQC_LINUX_DIR:=$(CHARLTON)/phreeqc
PHREEQC_LINUX_TARGET_DIR:=$(GEOCHEM)/unix/phreeqc
PHREEQC_LINUX_PUB_DIR:=$(FTP_CHARLTON)/phreeqc
PHREEQC_LINUX_FILE:=$(shell cd "$(PHREEQC_LINUX_DIR)"; ls -t phreeqc-*Linux.tar.gz | awk '{if (NR == 1) {print}}')
PHREEQC_LINUX_VERSION:=$(shell echo $(PHREEQC_LINUX_FILE) | sed "s/phreeqc-//" | sed "s/.Linux.tar.gz//" )
PHREEQC_LINUX_SVN:=$(shell echo $(PHREEQC_LINUX_FILE) | sed "s/phreeqc-//" | sed "s/.Linux.tar.gz//" | sed "s/^.*-//" )
PHREEQC_LINUX_SIZE:=$(shell ls -s -h "$(PHREEQC_LINUX_DIR)/$(PHREEQC_LINUX_FILE)" | awk '{print $$1}')
PHREEQC_LINUX_DATE:=$(shell ls -l --time-style "+%B %e, %G" "$(PHREEQC_LINUX_DIR)/$(PHREEQC_LINUX_FILE)" | awk '{print $$6,$$7,$$8}')

#Batch source

PHREEQC_SOURCE_DIR:=$(CHARLTON)/phreeqc
PHREEQC_SOURCE_TARGET_DIR:=$(GEOCHEM)/unix/phreeqc
PHREEQC_SOURCE_PUB_DIR:=$(FTP_CHARLTON)/phreeqc
PHREEQC_SOURCE_FILE:=$(shell cd "$(PHREEQC_SOURCE_DIR)"; ls -t phreeqc-*.source.tar.gz | awk '{if (NR == 1) {print}}')
PHREEQC_SOURCE_VERSION:=$(shell echo $(PHREEQC_SOURCE_FILE) | sed "s/phreeqc-//" | sed "s/.source.tar.gz//" )
PHREEQC_SOURCE_SVN:=$(shell echo $(PHREEQC_SOURCE_FILE) | sed "s/phreeqc-//" | sed "s/.source.tar.gz//" | sed "s/^.*-//")
PHREEQC_SOURCE_SIZE:=$(shell ls -s -h "$(PHREEQC_SOURCE_DIR)/$(PHREEQC_SOURCE_FILE)" | awk '{print $$1}')
PHREEQC_SOURCE_DATE:=$(shell ls -l --time-style "+%B %e, %G" "$(PHREEQC_SOURCE_DIR)/$(PHREEQC_SOURCE_FILE)" | awk '{print $$6,$$7,$$8}')

#Batch Mac

PHREEQC_MAC_DIR:=$(GEOCHEM)/mac/phreeqc
PHREEQC_MAC_PUB_DIR:=$(FTP_GEOCHEM)/mac/phreeqc
PHREEQC_MAC_FILE:=$(shell cd "$(PHREEQC_MAC_DIR)"; ls -t phreeqc-*.dmg.gz | awk '{if (NR == 1) {print}}')
PHREEQC_MAC_VERSION:=$(shell echo $(PHREEQC_MAC_FILE) | sed "s/phreeqc-//" | sed "s/.dmg.gz//" )
PHREEQC_MAC_SVN:=$(shell echo $(PHREEQC_MAC_FILE) | sed "s/phreeqc-//" | sed "s/.dmg.gz//" | sed "s/^.*-//" )
PHREEQC_MAC_SIZE:=$(shell ls -s -h "$(PHREEQC_MAC_DIR)/$(PHREEQC_MAC_FILE)" | awk '{print $$1}')
PHREEQC_MAC_DATE:=$(shell ls -l --time-style "+%B %e, %G" "$(PHREEQC_MAC_DIR)/$(PHREEQC_MAC_FILE)" | awk '{print $$6,$$7,$$8}')

# Batch PHAST

#Batch Windows

PHAST_WIN_DIR:=$(CHARLTON)/phast
PHAST_WIN_TARGET_DIR:=$(GEOCHEM)/pc/phast
PHAST_WIN_PUB_DIR:=$(FTP_CHARLTON)/phast
PHAST_WIN_X64_FILE:=$(shell cd "$(PHAST_WIN_DIR)"; ls -t phast*-x64.msi | awk '{if (NR == 1) {print}}')
PHAST_WIN_FILE:=$(shell echo $(PHAST_WIN_X64_FILE) | sed "s/-x64//" )
PHAST_WIN_VERSION:=$(shell echo $(PHAST_WIN_FILE) | sed "s/phast-//" | sed "s/.exe//" )
PHAST_WIN_SVN:=$(shell echo $(PHAST_WIN_FILE) | sed "s/phast-//" | sed "s/.exe//" | sed "s/^.*-//" )
PHAST_WIN_SIZE:=$(shell ls -s -h "$(PHAST_WIN_DIR)/$(PHAST_WIN_FILE)" | awk '{print $$1}')
PHAST_WIN_X64_SIZE:=$(shell ls -s -h "$(PHAST_WIN_DIR)/$(PHAST_WIN_X64_FILE)" | awk '{print $$1}')
PHAST_WIN_DATE:=$(shell ls -l --time-style "+%B %e, %G" "$(PHAST_WIN_DIR)/$(PHAST_WIN_FILE)" | awk '{print $$6,$$7,$$8}')
PHAST_WIN_MD5:=$(shell md5sum "$(PHAST_WIN_DIR)/$(PHAST_WIN_FILE)" | awk '{print $$1}')
PHAST_WIN_X64_MD5:=$(shell md5sum "$(PHAST_WIN_DIR)/$(PHAST_WIN_X64_FILE)" | awk '{print $$1}')

#Batch Linux

PHAST_LINUX_DIR:=$(CHARLTON)/phast
PHAST_LINUX_TARGET_DIR:=$(GEOCHEM)/unix/phast
PHAST_LINUX_PUB_DIR:=$(FTP_CHARLTON)/phast
PHAST_LINUX_FILE:=$(shell cd "$(PHAST_LINUX_DIR)"; ls -t phast-*i586.tar.gz | awk '{if (NR == 1) {print}}')
PHAST_LINUX_VERSION:=$(shell echo $(PHAST_LINUX_FILE) | sed "s/phast-//" | sed "s/.i586.tar.gz//" )
PHAST_LINUX_SVN:=$(shell echo $(PHAST_LINUX_FILE) | sed "s/phast-//" | sed "s/.i586.tar.gz//" | sed "s/^.*-//")
PHAST_LINUX_SIZE:=$(shell ls -s -h "$(PHAST_LINUX_DIR)/$(PHAST_LINUX_FILE)" | awk '{print $$1}')
PHAST_LINUX_DATE:=$(shell ls -l --time-style "+%B %e, %G" "$(PHAST_LINUX_DIR)/$(PHAST_LINUX_FILE)" | awk '{print $$6,$$7,$$8}')
PHAST_LINUX_MD5:=$(shell md5sum "$(PHAST_LINUX_DIR)/$(PHAST_LINUX_FILE)" | awk '{print $$1}')

PHAST_LINUX_X64_FILE:=$(shell cd "$(PHAST_LINUX_DIR)"; ls -t phast-*x86_64.tar.gz | awk '{if (NR == 1) {print}}')
PHAST_LINUX_X64_SIZE:=$(shell ls -s -h "$(PHAST_LINUX_DIR)/$(PHAST_LINUX_X64_FILE)" | awk '{print $$1}')
PHAST_LINUX_X64_MD5:=$(shell md5sum "$(PHAST_LINUX_DIR)/$(PHAST_LINUX_X64_FILE)" | awk '{print $$1}')

PHAST_CONFIGURE_FILE:=$(shell echo $(PHAST_LINUX_X64_FILE) | sed "s/x86_64\.//")
PHAST_CONFIGURE_SIZE:=$(shell ls -s -h "$(PHAST_LINUX_DIR)/$(PHAST_CONFIGURE_FILE)" | awk '{print $$1}')
PHAST_CONFIGURE_MD5:=$(shell md5sum "$(PHAST_LINUX_DIR)/$(PHAST_CONFIGURE_FILE)" | awk '{print $$1}')


#RPM Linux

PHAST_LINUX_RPM_DIR:=$(CHARLTON)/phast
PHAST_LINUX_TARGET_RPM_DIR:=$(GEOCHEM)/unix/phast
PHAST_LINUX_PUB_RPM_DIR:=$(FTP_CHARLTON)/phast
PHAST_LINUX_RPM_FILE:=$(shell cd "$(PHAST_LINUX_RPM_DIR)"; ls -t phast-*i586.rpm | awk '{if (NR == 1) {print}}')
PHAST_LINUX_RPM_VERSION:=$(shell echo $(PHAST_LINUX_RPM_FILE) | sed "s/phast-//" | sed "s/.i586.rpm//" )
PHAST_LINUX_RPM_SVN:=$(shell echo $(PHAST_LINUX_RPM_FILE) | sed "s/phast-//" | sed "s/.i586.rpm//" | sed "s/^.*-//" )
PHAST_LINUX_RPM_SIZE:=$(shell ls -s -h "$(PHAST_LINUX_RPM_DIR)/$(PHAST_LINUX_RPM_FILE)" | awk '{print $$1}')
PHAST_LINUX_RPM_DATE:=$(shell ls -l --time-style "+%B %e, %G" "$(PHAST_LINUX_RPM_DIR)/$(PHAST_LINUX_RPM_FILE)" | awk '{print $$6,$$7,$$8}')
PHAST_LINUX_RPM_MD5:=$(shell md5sum "$(PHAST_LINUX_RPM_DIR)/$(PHAST_LINUX_RPM_FILE)" | awk '{print $$1}')

PHAST_LINUX_X64_RPM_FILE:=$(shell cd "$(PHAST_LINUX_RPM_DIR)"; ls -t phast-*x86_64.rpm | awk '{if (NR == 1) {print}}')
PHAST_LINUX_X64_RPM_VERSION:=$(shell echo $(PHAST_LINUX_X64_RPM_FILE) | sed "s/phast-//" | sed "s/.x86_64.rpm//" )
PHAST_LINUX_X64_RPM_SVN:=$(shell echo $(PHAST_LINUX_X64_RPM_FILE) | sed "s/phast-//" | sed "s/.x86_64.rpm//" | sed "s/^.*-//" )
PHAST_LINUX_X64_RPM_SIZE:=$(shell ls -s -h "$(PHAST_LINUX_RPM_DIR)/$(PHAST_LINUX_X64_RPM_FILE)" | awk '{print $$1}')
PHAST_LINUX_X64_RPM_DATE:=$(shell ls -l --time-style "+%B %e, %G" "$(PHAST_LINUX_RPM_DIR)/$(PHAST_LINUX_X64_RPM_FILE)" | awk '{print $$6,$$7,$$8}')
PHAST_LINUX_X64_RPM_MD5:=$(shell md5sum "$(PHAST_LINUX_RPM_DIR)/$(PHAST_LINUX_X64_RPM_FILE)" | awk '{print $$1}')


#WPHAST version

#WPHAST_DIR:=$(CHARLTON)/wphast
#WPHAST_TARGET_DIR:=$(GEOCHEM)/pc/wphast
#WPHAST_PUB_DIR:=$(FTP_CHARLTON)/wphast
#WPHAST_FILE:=$(shell cd "$(WPHAST_DIR)"; ls -t wphast*.msi | awk '{if (NR == 1) {print}}')
#WPHAST_VERSION:=$(shell echo $(WPHAST_FILE) | sed "s/wphast-//" | sed "s/.msi//" )
#WPHAST_SVN:=$(shell echo $(WPHAST_FILE) | sed "s/wphast-//" | sed "s/.msi//" | sed "s/^.*-//")
#WPHAST_SIZE:=$(shell ls -s -h "$(WPHAST_DIR)/$(WPHAST_FILE)" | awk '{print $$1}')
#WPHAST_DATE:=$(shell ls -l --time-style "+%B %e, %G" "$(WPHAST_DIR)/$(WPHAST_FILE)" | awk '{print $$6,$$7,$$8}')

#NetpathXL version

NETPATHXL_DIR:=$(CHARLTON)/netpathxl
NETPATHXL_PUB_DIR:=$(FTP_CHARLTON)/netpathxl
NETPATHXL_FILE:=$(shell cd "$(NETPATHXL_DIR)"; ls -t NetpathXL*.msi | awk '{if (NR == 1) {print}}')
NETPATHXL_VERSION:=$(shell echo $(NETPATHXL_FILE) | sed "s/NetpathXL-//" | sed "s/.msi//" )
NETPATHXL_SVN:=$(shell echo $(NETPATHXL_FILE) | sed "s/NetpathXL-//" | sed "s/.msi//" | sed "s/^.*-//")
NETPATHXL_SIZE:=$(shell ls -s -h "$(NETPATHXL_DIR)/$(NETPATHXL_FILE)" | awk '{print $$1}')
NETPATHXL_DATE:=$(shell ls -l --time-style "+%B %e, %G" "$(NETPATHXL_DIR)/$(NETPATHXL_FILE)" | awk '{print $$6,$$7,$$8}')


#Netpath-win version

NETPATHWIN_DIR:=$(GEOCHEM)/pc/netpath
NETPATHWIN_PUB_DIR:=$(FTP_GEOCHEM)/pc/netpath
NETPATHWIN_FILE:=NPWsetup.zip
#NETPATHWIN_FILE:=$(shell cd "$(NETPATHWIN_DIR)"; ls -t NPWsetup.zip | awk '{if (NR == 1) {print}}')
#NETPATHWIN_VERSION:=$(shell echo $(NETPATHWIN_FILE) | sed "s/NetpathXL-//" | sed "s/.msi//" )
#NETPATHWIN_SVN:=$(shell echo $(NETPATHWIN_FILE) | sed "s/NetpathXL-//" | sed "s/.msi//" | sed "s/^.*-//")
NETPATHWIN_SIZE:=$(shell ls -s -h "$(NETPATHWIN_DIR)/$(NETPATHWIN_FILE)" | awk '{print $$1}')
NETPATHWIN_DATE:=$(shell ls -l --time-style "+%B %e, %G" "$(NETPATHWIN_DIR)/$(NETPATHWIN_FILE)" | awk '{print $$6,$$7,$$8}')

SED_ARGS= \
	-e "s^@PHREEQCI_PUB_DIR@^$(PHREEQCI_PUB_DIR)^"  \
	-e "s^@PHREEQCI_FILE@^$(PHREEQCI_FILE)^" \
	-e "s^@PHREEQCI_VERSION@^$(PHREEQCI_VERSION)^" \
	-e "s^@PHREEQCI_SVN@^$(PHREEQCI_SVN)^" \
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
	-e "s^@PHREEQC_WIN_SVN@^$(PHREEQC_WIN_SVN)^" \
	-e "s^@PHREEQC_WIN_DATE@^$(PHREEQC_WIN_DATE)^" \
	-e "s^@PHREEQC_WIN_SIZE@^$(PHREEQC_WIN_SIZE)^" \
	-e "s^@PHREEQC_LINUX_PUB_DIR@^$(PHREEQC_LINUX_PUB_DIR)^" \
	-e "s^@PHREEQC_LINUX_FILE@^$(PHREEQC_LINUX_FILE)^" \
	-e "s^@PHREEQC_LINUX_VERSION@^$(PHREEQC_LINUX_VERSION)^" \
	-e "s^@PHREEQC_LINUX_SVN@^$(PHREEQC_LINUX_SVN)^" \
	-e "s^@PHREEQC_LINUX_DATE@^$(PHREEQC_LINUX_DATE)^" \
	-e "s^@PHREEQC_LINUX_SIZE@^$(PHREEQC_LINUX_SIZE)^" \
	-e "s^@PHREEQC_SOURCE_PUB_DIR@^$(PHREEQC_SOURCE_PUB_DIR)^" \
	-e "s^@PHREEQC_SOURCE_FILE@^$(PHREEQC_SOURCE_FILE)^" \
	-e "s^@PHREEQC_SOURCE_VERSION@^$(PHREEQC_SOURCE_VERSION)^" \
	-e "s^@PHREEQC_SOURCE_SVN@^$(PHREEQC_SOURCE_SVN)^" \
	-e "s^@PHREEQC_SOURCE_DATE@^$(PHREEQC_SOURCE_DATE)^" \
	-e "s^@PHREEQC_SOURCE_SIZE@^$(PHREEQC_SOURCE_SIZE)^" \
	-e "s^@PHREEQC_MAC_PUB_DIR@^$(PHREEQC_MAC_PUB_DIR)^" \
	-e "s^@PHREEQC_MAC_FILE@^$(PHREEQC_MAC_FILE)^" \
	-e "s^@PHREEQC_MAC_VERSION@^$(PHREEQC_MAC_VERSION)^" \
	-e "s^@PHREEQC_MAC_SVN@^$(PHREEQC_MAC_SVN)^" \
	-e "s^@PHREEQC_MAC_DATE@^$(PHREEQC_MAC_DATE)^" \
	-e "s^@PHREEQC_MAC_SIZE@^$(PHREEQC_MAC_SIZE)^" \
	-e "s^@PHAST_WIN_DIR@^$(PHAST_WIN_DIR)^g" \
	-e "s^@PHAST_WIN_PUB_DIR@^$(PHAST_WIN_PUB_DIR)^g" \
	-e "s^@PHAST_WIN_FILE@^$(PHAST_WIN_FILE)^g" \
	-e "s^@PHAST_WIN_X64_FILE@^$(PHAST_WIN_X64_FILE)^g" \
	-e "s^@PHAST_WIN_VERSION@^$(PHAST_WIN_VERSION)^g" \
	-e "s^@PHAST_WIN_SVN@^$(PHAST_WIN_SVN)^g" \
	-e "s^@PHAST_WIN_DATE@^$(PHAST_WIN_DATE)^g" \
	-e "s^@PHAST_WIN_SIZE@^$(PHAST_WIN_SIZE)^g" \
	-e "s^@PHAST_WIN_X64_SIZE@^$(PHAST_WIN_X64_SIZE)^g" \
	-e "s^@PHAST_WIN_MD5@^$(PHAST_WIN_MD5)^g" \
	-e "s^@PHAST_WIN_X64_MD5@^$(PHAST_WIN_X64_MD5)^g" \
	-e "s^@PHAST_LINUX_DIR@^$(PHAST_LINUX_DIR)^g" \
	-e "s^@PHAST_LINUX_PUB_DIR@^$(PHAST_LINUX_PUB_DIR)^g" \
	-e "s^@PHAST_LINUX_FILE@^$(PHAST_LINUX_FILE)^g" \
	-e "s^@PHAST_LINUX_X64_FILE@^$(PHAST_LINUX_X64_FILE)^g" \
	-e "s^@PHAST_LINUX_VERSION@^$(PHAST_LINUX_VERSION)^g" \
	-e "s^@PHAST_LINUX_SVN@^$(PHAST_LINUX_SVN)^g" \
	-e "s^@PHAST_LINUX_DATE@^$(PHAST_LINUX_DATE)^g" \
	-e "s^@PHAST_LINUX_SIZE@^$(PHAST_LINUX_SIZE)^g" \
	-e "s^@PHAST_LINUX_MD5@^$(PHAST_LINUX_MD5)^g" \
	-e "s^@PHAST_LINUX_X64_SIZE@^$(PHAST_LINUX_X64_SIZE)^g" \
	-e "s^@PHAST_LINUX_X64_MD5@^$(PHAST_LINUX_X64_MD5)^g" \
	-e "s^@PHAST_LINUX_RPM_DIR@^$(PHAST_LINUX_RPM_DIR)^g" \
	-e "s^@PHAST_LINUX_PUB_RPM_DIR@^$(PHAST_LINUX_PUB_RPM_DIR)^g" \
	-e "s^@PHAST_LINUX_RPM_FILE@^$(PHAST_LINUX_RPM_FILE)^g" \
	-e "s^@PHAST_LINUX_RPM_VERSION@^$(PHAST_LINUX_RPM_VERSION)^g" \
	-e "s^@PHAST_LINUX_RPM_SVN@^$(PHAST_LINUX_RPM_SVN)^g" \
	-e "s^@PHAST_LINUX_RPM_DATE@^$(PHAST_LINUX_RPM_DATE)^g" \
	-e "s^@PHAST_LINUX_RPM_SIZE@^$(PHAST_LINUX_RPM_SIZE)^g" \
	-e "s^@PHAST_LINUX_RPM_MD5@^$(PHAST_LINUX_RPM_MD5)^g" \
	-e "s^@PHAST_LINUX_X64_RPM_FILE@^$(PHAST_LINUX_X64_RPM_FILE)^g" \
	-e "s^@PHAST_LINUX_X64_RPM_VERSION@^$(PHAST_LINUX_X64_RPM_VERSION)^g" \
	-e "s^@PHAST_LINUX_X64_RPM_SVN@^$(PHAST_LINUX_X64_RPM_SVN)^g" \
	-e "s^@PHAST_LINUX_X64_RPM_DATE@^$(PHAST_LINUX_X64_RPM_DATE)^g" \
	-e "s^@PHAST_LINUX_X64_RPM_SIZE@^$(PHAST_LINUX_X64_RPM_SIZE)^g" \
	-e "s^@PHAST_LINUX_X64_RPM_MD5@^$(PHAST_LINUX_X64_RPM_MD5)^g" \
	-e "s^@PHAST_CONFIGURE_FILE@^$(PHAST_CONFIGURE_FILE)^g" \
	-e "s^@PHAST_CONFIGURE_SIZE@^$(PHAST_CONFIGURE_SIZE)^g" \
	-e "s^@PHAST_CONFIGURE_MD5@^$(PHAST_CONFIGURE_MD5)^g" \
	-e "s^@WPHAST_PUB_DIR@^$(WPHAST_PUB_DIR)^"  \
	-e "s^@WPHAST_FILE@^$(WPHAST_FILE)^" \
	-e "s^@WPHAST_VERSION@^$(WPHAST_VERSION)^" \
	-e "s^@WPHAST_SVN@^$(WPHAST_SVN)^" \
	-e "s^@WPHAST_DATE@^$(WPHAST_DATE)^" \
	-e "s^@WPHAST_SIZE@^$(WPHAST_SIZE)^" \
	-e "s^@NETPATHXL_PUB_DIR@^$(NETPATHXL_PUB_DIR)^"  \
	-e "s^@NETPATHXL_FILE@^$(NETPATHXL_FILE)^" \
	-e "s^@NETPATHXL_VERSION@^$(NETPATHXL_VERSION)^" \
	-e "s^@NETPATHXL_SVN@^$(NETPATHXL_SVN)^" \
	-e "s^@NETPATHXL_DATE@^$(NETPATHXL_DATE)^" \
	-e "s^@NETPATHXL_SIZE@^$(NETPATHXL_SIZE)^" \
	-e "s^@NETPATHWIN_PUB_DIR@^$(NETPATHWIN_PUB_DIR)^"  \
	-e "s^@NETPATHWIN_FILE@^$(NETPATHWIN_FILE)^" \
	-e "s^@NETPATHWIN_DATE@^$(NETPATHWIN_DATE)^" \
	-e "s^@NETPATHWIN_SIZE@^$(NETPATHWIN_SIZE)^" 


trial: Makefile template_phreeqc_index.html template_phreeqci_index.html template_phast_index.html template_netpath_index.html
	sed $(SED_ARGS) < template_phreeqc_index.html > web_phreeqc_index.html
	sed $(SED_ARGS) < template_phreeqci_index.html > web_phreeqci_index.html
	sed $(SED_ARGS) < template_phast_index.html > web_phast_index.html
	sed $(SED_ARGS) < template_netpath_index.html > web_netpath_index.html
	
test:
	@echo "***PhreeqcI version ***"
	@echo "PHREEQCI_DIR                = $(PHREEQCI_DIR)"
	@echo "PHREEQCI_TARGET_DIR         = $(PHREEQCI_TARGET_DIR)"
	@echo "PHREEQCI_PUB_DIR            = $(PHREEQCI_PUB_DIR)"
	@echo "PHREEQCI_FILE               = $(PHREEQCI_FILE)"
	@echo "PHREEQCI_VERSION            = $(PHREEQCI_VERSION)"
	@echo "PHREEQCI_SIZE               = $(PHREEQCI_SIZE)"
	@echo "PHREEQCI_DATE               = $(PHREEQCI_DATE)"
	@echo "*** Phreeqc For Windows version ***"
	@echo "PFW_DIR                     = $(PFW_DIR)"
	@echo "PFW_PUB_DIR                 = $(PFW_PUB_DIR)"
	@echo "PFW_FILE                    = $(PFW_FILE)"
	@echo "PFW_VERSION                 = $(PFW_VERSION)"
	@echo "PFW_SIZE                    = $(PFW_SIZE)"
	@echo "PFW_DATE                    = $(PFW_DATE)"
	@echo "*** Batch phreeqc Windows ***"
	@echo "PHREEQC_WIN_DIR             = $(PHREEQC_WIN_DIR)"
	@echo "PHREEQC_WIN_TARGET_DIR      = $(PHREEQC_WIN_TARGET_DIR)"
	@echo "PHREEQC_WIN_PUB_DIR         = $(PHREEQC_WIN_PUB_DIR)"
	@echo "PHREEQC_WIN_FILE            = $(PHREEQC_WIN_FILE)"
	@echo "PHREEQC_WIN_VERSION         = $(PHREEQC_WIN_VERSION)"
	@echo "PHREEQC_WIN_SVN             = $(PHREEQC_WIN_SVN)"
	@echo "PHREEQC_WIN_SIZE            = $(PHREEQC_WIN_SIZE)"
	@echo "PHREEQC_WIN_DATE            = $(PHREEQC_WIN_DATE)"
	@echo "*** Batch phreeqc Linux ***"
	@echo "PHREEQC_LINUX_DIR           = $(PHREEQC_LINUX_DIR)"
	@echo "PHREEQC_LINUX_TARGET_DIR    = $(PHREEQC_LINUX_TARGET_DIR)"
	@echo "PHREEQC_LINUX_PUB_DIR       = $(PHREEQC_LINUX_PUB_DIR)"
	@echo "PHREEQC_LINUX_FILE          = $(PHREEQC_LINUX_FILE)"
	@echo "PHREEQC_LINUX_VERSION       = $(PHREEQC_LINUX_VERSION)"
	@echo "PHREEQC_LINUX_SVN           = $(PHREEQC_LINUX_SVN)"
	@echo "PHREEQC_LINUX_SIZE          = $(PHREEQC_LINUX_SIZE)"
	@echo "PHREEQC_LINUX_DATE          = $(PHREEQC_LINUX_DATE)"
	@echo "*** Batch phreeqc source ***"
	@echo "PHREEQC_SOURCE_DIR          = $(PHREEQC_SOURCE_DIR)"
	@echo "PHREEQC_SOURCE_TARGET_DIR   = $(PHREEQC_SOURCE_TARGET_DIR)"
	@echo "PHREEQC_SOURCE_PUB_DIR      = $(PHREEQC_SOURCE_PUB_DIR)"
	@echo "PHREEQC_SOURCE_FILE         = $(PHREEQC_SOURCE_FILE)"
	@echo "PHREEQC_SOURCE_VERSION      = $(PHREEQC_SOURCE_VERSION)"
	@echo "PHREEQC_SOURCE_SVN          = $(PHREEQC_SOURCE_SVN)"
	@echo "PHREEQC_SOURCE_SIZE         = $(PHREEQC_SOURCE_SIZE)"
	@echo "PHREEQC_SOURCE_DATE         = $(PHREEQC_SOURCE_DATE)"
	@echo "*** Batch Mac ***"
	@echo "PHREEQC_MAC_DIR             = $(PHREEQC_MAC_DIR)"
	@echo "PHREEQC_MAC_PUB_DIR         = $(PHREEQC_MAC_PUB_DIR)"
	@echo "PHREEQC_MAC_FILE            = $(PHREEQC_MAC_FILE)"
	@echo "PHREEQC_MAC_VERSION         = $(PHREEQC_MAC_VERSION)"
	@echo "PHREEQC_MAC_SVN             = $(PHREEQC_MAC_SVN)"
	@echo "PHREEQC_MAC_SIZE            = $(PHREEQC_MAC_SIZE)"
	@echo "PHREEQC_MAC_DATE            = $(PHREEQC_MAC_DATE)"
	@echo "*** Batch PHAST Windows ***"
	@echo "PHAST_WIN_DIR               = $(PHAST_WIN_DIR)"
	@echo "PHAST_WIN_TARGET_DIR        = $(PHAST_WIN_TARGET_DIR)"
	@echo "PHAST_WIN_PUB_DIR           = $(PHAST_WIN_PUB_DIR)"
	@echo "PHAST_WIN_X64_FILE          = $(PHAST_WIN_X64_FILE)"
	@echo "PHAST_WIN_FILE              = $(PHAST_WIN_FILE)"
	@echo "PHAST_WIN_VERSION           = $(PHAST_WIN_VERSION)"
	@echo "PHAST_WIN_SVN               = $(PHAST_WIN_SVN)"
	@echo "PHAST_WIN_SIZE              = $(PHAST_WIN_SIZE)"
	@echo "PHAST_WIN_X64_SIZE          = $(PHAST_WIN_X64_SIZE)"
	@echo "PHAST_WIN_DATE              = $(PHAST_WIN_DATE)"
	@echo "PHAST_WIN_MD5               = $(PHAST_WIN_MD5)"
	@echo "PHAST_WIN_X64_MD5           = $(PHAST_WIN_X64_MD5)"
	@echo "*** Batch Linux ***"
	@echo "PHAST_LINUX_DIR             = $(PHAST_LINUX_DIR)"
	@echo "PHAST_LINUX_TARGET_DIR      = $(PHAST_LINUX_TARGET_DIR)"
	@echo "PHAST_LINUX_PUB_DIR         = $(PHAST_LINUX_PUB_DIR)"
	@echo "PHAST_LINUX_FILE            = $(PHAST_LINUX_FILE)"
	@echo "PHAST_LINUX_VERSION         = $(PHAST_LINUX_VERSION)"
	@echo "PHAST_LINUX_SVN             = $(PHAST_LINUX_SVN)"
	@echo "PHAST_LINUX_SIZE            = $(PHAST_LINUX_SIZE)"
	@echo "PHAST_LINUX_DATE            = $(PHAST_LINUX_DATE)"
	@echo "PHAST_LINUX_MD5             = $(PHAST_LINUX_MD5)"
	@echo "PHAST_LINUX_X64_FILE        = $(PHAST_LINUX_X64_FILE)"
	@echo "PHAST_LINUX_X64_SIZE        = $(PHAST_LINUX_X64_SIZE)"
	@echo "PHAST_LINUX_X64_MD5         = $(PHAST_LINUX_X64_MD5)"
	@echo "PHAST_CONFIGURE_FILE        = $(PHAST_CONFIGURE_FILE)"
	@echo "PHAST_CONFIGURE_SIZE        = $(PHAST_CONFIGURE_SIZE)"
	@echo "PHAST_CONFIGURE_MD5         = $(PHAST_CONFIGURE_MD5)"
	@echo "*** RPM Linux ***"
	@echo "PHAST_LINUX_RPM_DIR         = $(PHAST_LINUX_RPM_DIR)"
	@echo "PHAST_LINUX_TARGET_RPM_DIR  = $(PHAST_LINUX_TARGET_RPM_DIR)"
	@echo "PHAST_LINUX_PUB_RPM_DIR     = $(PHAST_LINUX_PUB_RPM_DIR)"
	@echo "PHAST_LINUX_RPM_FILE        = $(PHAST_LINUX_RPM_FILE)"
	@echo "PHAST_LINUX_RPM_VERSION     = $(PHAST_LINUX_RPM_VERSION)"
	@echo "PHAST_LINUX_RPM_SVN         = $(PHAST_LINUX_RPM_SVN)"
	@echo "PHAST_LINUX_RPM_SIZE        = $(PHAST_LINUX_RPM_SIZE)"
	@echo "PHAST_LINUX_RPM_DATE        = $(PHAST_LINUX_RPM_DATE)"
	@echo "PHAST_LINUX_RPM_MD5         = $(PHAST_LINUX_RPM_MD5)"
	@echo "*** RPM Linux x64 ***"
	@echo "PHAST_LINUX_X64_RPM_FILE    = $(PHAST_LINUX_X64_RPM_FILE)"
	@echo "PHAST_LINUX_X64_RPM_VERSION = $(PHAST_LINUX_X64_RPM_VERSION)"
	@echo "PHAST_LINUX_X64_RPM_SVN     = $(PHAST_LINUX_X64_RPM_SVN)"
	@echo "PHAST_LINUX_X64_RPM_SIZE    = $(PHAST_LINUX_X64_RPM_SIZE)"
	@echo "PHAST_LINUX_X64_RPM_DATE    = $(PHAST_LINUX_X64_RPM_DATE)"
	@echo "PHAST_LINUX_X64_RPM_MD5     = $(PHAST_LINUX_X64_RPM_MD5)"




web: 
# PHREEQC
	sed $(SED_ARGS) < template_phreeqc_index.html > web_phreeqc_index.html
	cp -f web_phreeqc_index.html $(HTML_PHREEQC)/index.html 
	cd $(HTML_PHREEQC); rcs -l index.html; ci -l -m"Latest" index.html
#	cp -f $(PHREEQC_WIN_DIR)/$(PHREEQC_WIN_FILE) $(PHREEQC_WIN_TARGET_DIR)/$(PHREEQC_WIN_FILE)
#	cp -f $(PHREEQC_LINUX_DIR)/$(PHREEQC_LINUX_FILE) $(PHREEQC_LINUX_TARGET_DIR)/$(PHREEQC_LINUX_FILE)
#	cp -f $(PHREEQC_SOURCE_DIR)/$(PHREEQC_SOURCE_FILE) $(PHREEQC_SOURCE_TARGET_DIR)/$(PHREEQC_SOURCE_FILE)
#	cp -f $(PHREEQC_WIN_PUB_DIR)/README.Win.$(PHREEQC_WIN_SVN).TXT $(HTML_PHREEQC)
# uncomment next line
# ERROR with this copy because .exe not available
	cp -f $(PHREEQC_WIN_DIR)/README.Win.$(PHREEQC_WIN_SVN).TXT $(HTML_PHREEQC)
	cp -f $(PHREEQC_LINUX_DIR)/README.Linux.$(PHREEQC_LINUX_SVN).TXT $(HTML_PHREEQC)
	cp -f $(PHREEQC_MAC_DIR)/README.TXT $(HTML_PHREEQC)/README.Mac.TXT
	cp -f $(PHREEQC_WIN_DIR)/RELEASE.$(PHREEQC_LINUX_SVN).TXT $(HTML_PHREEQC)

# PHREEQCI
	sed $(SED_ARGS) < template_phreeqci_index.html > web_phreeqci_index.html
	cp -f web_phreeqci_index.html $(HTML_PHREEQCI)/index.html 
	cd $(HTML_PHREEQCI); rcs -l index.html; ci -l -m"Latest" index.html
#	cp -f $(PHREEQCI_DIR)/$(PHREEQCI_FILE) $(PHREEQCI_TARGET_DIR)/$(PHREEQCI_FILE)

# PHAST
	sed $(SED_ARGS) < template_phast_index.html > web_phast_index.html
	cp -f web_phast_index.html $(HTML_PHAST)/index.html 
	cd $(HTML_PHAST); rcs -l index.html; ci -l -m"Latest" index.html
	cp -f $(PHAST_LINUX_DIR)/README.$(PHAST_LINUX_SVN) $(HTML_PHAST)
	cp -f $(PHAST_LINUX_DIR)/RELEASE.$(PHAST_LINUX_SVN) $(HTML_PHAST)

# NETPATH
	sed $(SED_ARGS) < template_netpath_index.html > web_netpath_index.html
	cp -f web_netpath_index.html $(HTML_NETPATH)/index.html 
	cp -f $(NETPATHXL_DIR)/README.$(NETPATHXL_SVN).TXT $(HTML_NETPATH)

#	mkdir -p temp
#	cp -f $(PHREEQC_LINUX_DIR)/$(PHREEQC_LINUX_FILE) temp
#	cd temp; tar -xzf $(PHREEQC_LINUX_FILE); cp -f ./*/doc/RELEASE.TXT $(HTML_PHREEQC); cd ..; rm -rf temp


clean:
	rm -f web_netpath_index.html web_phast_index.html web_phreeqci_index.html web_phreeqc_index.html