#include "GSam.h"
#include "iostream"
#include <gclib/GArgs.h>
#include <gclib/GStr.h>
#include "commons.h"
#include "tmerge.h"

#define VERSION "0.0.1"

using namespace std;

const char* USAGE = "Junction_extractor v" VERSION "\n"
                              "==========================================================================================\n"
                              "This is a simple script to extract all splice junctions from a BAM file or a list of BAM files.\n"
                              "==========================================================================================\n\n"
                              "  usage: junction_extractor [-hv] [-o junctions.bed] input\n\n";

GStr outfname;
FILE* outf=NULL;
GSamWriter* outfile=NULL;
bool verbose=false;
int juncCount=0;

TInputFiles inRecords;

struct CJunc {
	int start, end;
	char strand;
	uint64_t dupcount;
	CJunc(int vs=0, int ve=0, char vstrand='+', uint64_t dcount=1):
	  start(vs), end(ve), strand(vstrand), dupcount(dcount) { }

	bool operator==(const CJunc& a) {
		return (strand==a.strand && start==a.start && end==a.end);
	}

    bool operator<(const CJunc& a) { // sort by strand last
        if (start==a.start){
            if(end==a.end){
                return strand<a.strand;
            }
            else{
                return (end<a.end);
            }
        }
        else{
            return (start<a.start);
        }
    }

	void add(CJunc& j) {
       dupcount+=j.dupcount;
	}

	void write(FILE* f, const char* chr) {
		juncCount++;
		fprintf(f, "%s\t%d\t%d\tJUNC%08d\t%ld\t%c\n",
				chr, start-1, end, juncCount, (long)dupcount, strand);
	}
};

GArray<CJunc> junctions(64, true);

void addJunction(GSamRecord& r, int dupcount) {
	char strand = r.spliceStrand();
//	if (strand!='+' && strand!='-') return; // TODO: should we output .?
	for (int i=1;i<r.exons.Count();i++) {
		CJunc j(r.exons[i-1].end+1, r.exons[i].start-1, strand,
				dupcount);
		int ei;
		int r=junctions.AddIfNew(j, &ei);
		if (r==-1) {//existing junction, update
			junctions[ei].add(j);
		}
	}
}

void flushJuncs(FILE* f, const char* chr) {
    for (int i=0;i<junctions.Count();i++) {
    	junctions[i].write(f,chr);
    }
    junctions.Clear();
    junctions.setCapacity(128);
}

void processOptions(int argc, char* argv[]);

int main(int argc, char*argv[]) {

    // ANSI Shadow
	const char *banner = R"""(
     ██╗██╗   ██╗███╗   ██╗ ██████╗████████╗██╗ ██████╗ ███╗   ██╗    ███████╗██╗  ██╗████████╗██████╗  █████╗  ██████╗████████╗ ██████╗ ██████╗ 
     ██║██║   ██║████╗  ██║██╔════╝╚══██╔══╝██║██╔═══██╗████╗  ██║    ██╔════╝╚██╗██╔╝╚══██╔══╝██╔══██╗██╔══██╗██╔════╝╚══██╔══╝██╔═══██╗██╔══██╗
     ██║██║   ██║██╔██╗ ██║██║        ██║   ██║██║   ██║██╔██╗ ██║    █████╗   ╚███╔╝    ██║   ██████╔╝███████║██║        ██║   ██║   ██║██████╔╝
██   ██║██║   ██║██║╚██╗██║██║        ██║   ██║██║   ██║██║╚██╗██║    ██╔══╝   ██╔██╗    ██║   ██╔══██╗██╔══██║██║        ██║   ██║   ██║██╔══██╗
╚█████╔╝╚██████╔╝██║ ╚████║╚██████╗   ██║   ██║╚██████╔╝██║ ╚████║    ███████╗██╔╝ ██╗   ██║   ██║  ██║██║  ██║╚██████╗   ██║   ╚██████╔╝██║  ██║
 ╚════╝  ╚═════╝ ╚═╝  ╚═══╝ ╚═════╝   ╚═╝   ╚═╝ ╚═════╝ ╚═╝  ╚═══╝    ╚══════╝╚═╝  ╚═╝   ╚═╝   ╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝   ╚═╝    ╚═════╝ ╚═╝  ╚═╝
	)""";
    cout << banner << endl;

    inRecords.setup(VERSION, argc, argv);
	processOptions(argc, argv);
	int numSamples=inRecords.start();
	cout << "numSamples: " << numSamples << endl;
	cout << "outfname  : " << outfname << endl;
	outfile = new GSamWriter(outfname, inRecords.header(), GSamFile_BAM);
	TInputRecord* irec=NULL;
	GSamRecord* brec=NULL;

	/***************************
	 * Creating the output junction bed file
	 ***************************/    
    if (!outfname.is_empty()) {
        if (strcmp(outfname.substr(outfname.length()-4, 4).chars(), ".bed")!=0) {
            outfname.append(".bed");
        }
        outf = fopen(outfname.chars(), "w");
        if (outf==NULL) GError("Error creating file %s\n", outfname.chars());
        fprintf(outf, "track name=junctions\n");
    }

	/***************************
	 * Reading BAM file.
	 ***************************/
	int counter = 0;
    int prev_tid=-1;
    char* prev_refname;
    GVec<uint64_t> bcov(2048*1024);
    std::vector<std::pair<float,uint64_t>> bsam(2048*1024,{0,1}); // number of samples. 1st - current average; 2nd - total number of values
//    std::vector<std::set<int>> bsam_idx(2048*1024,std::set<int>{}); // for indexed runs
    int b_end=0; //bundle start, end (1-based)
    int b_start=0; //1 based

    // cout << "b_end: " << b_end << "  b_start: " << b_start << endl;
	while ((irec=inRecords.next())!=NULL) {
		brec=irec->brec;
        // cout << brec->refId() << endl;
        uint32_t dupcount=0;
        std::vector<int> cur_samples;
        int endpos=brec->end;

        if (brec->refId()!=prev_tid || (int)brec->start>b_end) {
            if (outf) {
                flushJuncs(outf, prev_refname);
            } // TODO: write the last column to 3 dec places
            b_start=brec->start;
            b_end=endpos;
            prev_tid=brec->refId();

            prev_refname=(char*)brec->refName();
        } else { //extending current bundle
            if (b_end<endpos) {
                b_end=endpos;
                bcov.setCount(b_end-b_start+1, (int)0);
            }
        }
        int accYC = 0;
        accYC = brec->tag_int("YC", 1);
        // cout << "accYC: " << accYC << endl;
        if (outf && brec->exons.Count()>1) {
            addJunction(*brec, accYC);
        }
	}
    flushJuncs(outf, prev_refname);
    fclose(outf);
}

void processOptions(int argc, char* argv[]) {
    GArgs args(argc, argv, "help;debug;verbose;version;full;clip;exon;keep-supp;keep-unmap;SMLPEDVho:N:Q:F:G:");
    args.printError(USAGE, true);

    if (args.getOpt('h') || args.getOpt("help")) {
        fprintf(stdout,"%s", USAGE);
        exit(0);
    }

    if (args.getOpt('v') || args.getOpt("version")) {
        fprintf(stdout,"%s\n", VERSION);
        exit(0);
    }

    if (args.startNonOpt()==0) {
        GMessage(USAGE);
        GMessage("\nError: no input provided!\n");
        exit(1);
    }

    outfname=args.getOpt('o');
    if (outfname.is_empty()) {
        GMessage(USAGE);
        GMessage("\nError: output filename must be provided (-o)!\n");
        exit(1);
    }

    verbose=(args.getOpt("verbose")!=NULL || args.getOpt('V')!=NULL);
    if (verbose) {
        fprintf(stderr, "Running intron_matcher " VERSION ". Command line:\n");
        args.printCmdLine(stderr);
    }
    const char* ifn=NULL;
    while ( (ifn=args.nextNonOpt())!=NULL) {
        //input alignment files
        std::string absolute_ifn = get_full_path(ifn);
        cout << "absolute_ifn: " << absolute_ifn << endl;
        inRecords.addFile(absolute_ifn.c_str());
    }


	// if (args.getOpt('G')) {
	// 	guidegff=args.getOpt('G');
	// 	if (fileExists(guidegff.chars())>1) {
	// 		// guided=true;
	// 	} else {
	// 		GError("Error: reference annotation file (%s) not found.\n",
	// 				guidegff.chars());
	// 	}
	// } else {
    //     GMessage(USAGE);
    //     GMessage("\nError: gff reference file must be provided (-G)!\n");
    //     exit(1);
	// }

    // verbose=(args.getOpt("verbose")!=NULL || args.getOpt('V')!=NULL);
    // if (verbose) {
    //     fprintf(stderr, "Running intron_matcher " VERSION ". Command line:\n");
    //     args.printCmdLine(stderr);
    // }
    // const char* ifn=NULL;
    // while ( (ifn=args.nextNonOpt())!=NULL) {
    //     //input alignment files
    //     std::string absolute_ifn = get_full_path(ifn);
    //     cout << "absolute_ifn: " << absolute_ifn << endl;
    //     inRecords.addFile(absolute_ifn.c_str());
    // }

}