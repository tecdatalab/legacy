#include "Grep.h"

int npats =  178 ;
char* respat[] = {"^.*$","^.*$","^WAT|HOH|H2O|DOD|DIS$","^CA$","^CD$","^.*$","^ACE$", "^.*$", "^.*$", "^.*$", "^.*$", "^.*$", "^ALA$", "^ILE|THR|VAL$", "^.*$", "^ASN|ASP|ASX|HIS|HIP|HIE|HID|HISN|HISL|LEU|PHE|TRP|TYR$", "^ARG|GLU|GLN|GLX|MET$", "^LEU$", "^.*$", "^GLN$", "^GLN$", "^ACE$", "^ARG$", "^ARG$", "^ARG$", "^ARG$", "^ARG$", "^ARG$", "^ASN$", "^ASN$", "^ASN$", "^ASN$", "^ASP$", "^ASP$", "^ASX$", "^ASX$", "^ASX$", "^ASX$", "^ASX$", "^CYS|MET$", "^CY[SXM]$", "^CYH$", "^GLU$", "^GLU$", "^GLU|GLN|GLX$", "^GLN$", "^GLN$", "^GLN|GLX$", "^HIS|HID|HIE|HIP|HISL$", "^HIS|HIE|HISL$", "^HID|HIP$", "^HID|HIP$", "^HIS|HIE|HIP$", "^HIS|HIE|HIP$", "^HID|HISL$", "^HID|HISL$", "^HIS|HID|HIP|HISD$", "^ILE$", "^ILE$", "^ILE$", "^LEU$", "^LEU$", "^LYS$", "^LYS$", "^LYS$", "^MET$", "^MET$", "^PHE$", "^PHE$", "^PRO|CPR$", "^CSO$", "^CSO$", "^CSO$", "^CSO$", "^SER$", "^THR$", "^THR$", "^TRP$", "^TRP$", "^TRP$", "^TRP$", "^TRP$", "^TRP$", "^TRP$", "^TRP$", "^TYR$", "^TYR$", "^TYR$", "^VAL$", "^VAL$", "^.*$", "^.*$", "^FS[34]$", "^FS[34]$", "^FS3$", "^FEO$", "^FEO$", "^HEM$", "^HEM$", "^HEM$", "^HEM$", "^HEM$", "^HEM$", "^HEM$", "^HEM$", "^HEM$", "^HEM$", "^HEM$", "^HEM$", "^HEM$", "^AZI$", "^MPD$", "^MPD$", "^MPD$", "^MPD$", "^MPD$", "^MPD$", "^MPD$", "^MPD$", "^SO4|SUL$", "^SO4|SUL$", "^PO4|PHO$", "^PC$", "^PC$", "^PC$", "^PC$", "^PC$", "^PC$", "^BIG$", "^POI$", "^DOT$", "^.*$", "^.*$", "^.*$", "^.*$", "^.*$", "^.*$", "^.*$", "^.*$", "^.*$", "^FMN$", "^FMN$", "^FMN$", "^FMN$", "^FMN$", "^FMN$", "^FMN$", "^FMN$", "^FMN$", "^FMN$", "^FMN$", "^FMN$", "^FMN$", "^FMN$", "^FMN$", "^ALK|MYR$", "^ALK|MYR$", "^ALK$", "^MYR$", "^ALK|MYR$", "^.*$", "^.*$", "^.*$", "^.*$", "^.*$", "^.*$", "^.*$", "^.*$", "^.*$", "^.*$", "^.*$", "^.*$", "^.*$", "^FAD|NAD|AMX|APU$", "^FAD|NAD|AMX|APU$", "^FAD|NAD|AMX|APU$", "^FAD|NAD|AMX|APU$", "^FAD|NAD|AMX|APU$" }; 
float get_radius(char* resname, char* aname) {

	int pat, h_select=5;	// 4 for water
	CGrep grep;

	if (grep.match("^[ 0-9][HhDd]$", aname))
		aname = "H";

	if (grep.match("^[Hh][^Gg]$", aname))
		aname = "H";

	for(pat=0;pat<npats;pat++) {
		if (grep.match(atmpat[pat], aname) &&
			grep.match(respat[pat], resname))
		break;
	}
	if(pat==npats) {
		for(pat=0;pat<npats;pat++) {
			if (grep.match(atmpat[pat], aname)) break;
		} if(pat==npats) return -1;
	}
	return united_rad[atmnum[pat]-1];
}