// Single-sample, exact raw-R1 bead-barcode compatibility diagnostic.
#include "whitelist.hpp"
#include <zlib.h>
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>

namespace fs = std::filesystem;
constexpr size_t WIDTH = 1536, TOP = 50000;
void require(bool ok, const std::string& message) {
    if (!ok) throw std::runtime_error(message);
}
std::string quote(const std::string& s) {
    std::ostringstream out; out << '"';
    for (unsigned char c : s) {
        switch(c) {
        case '"': out << "\\\""; break;
        case '\\': out << "\\\\"; break;
        case '\n': out << "\\n"; break;
        case '\r': out << "\\r"; break;
        case '\t': out << "\\t"; break;
        default:
            if (c < 32) out << "\\u" << std::hex << std::setw(4) << std::setfill('0') << int(c) << std::dec;
            else out << char(c);
        }
    }
    return out.str() + '"';
}
struct Signature {
    dev_t device; ino_t inode; off_t bytes; timespec modified;
    explicit Signature(const std::string& path) {
        struct stat st{};
        require(stat(path.c_str(), &st) == 0 && S_ISREG(st.st_mode), "Cannot stat regular input: " + path);
        device=st.st_dev; inode=st.st_ino; bytes=st.st_size; modified=st.st_mtim;
    }
    bool operator==(const Signature& b) const {
        return device==b.device && inode==b.inode && bytes==b.bytes &&
            modified.tv_sec==b.modified.tv_sec && modified.tv_nsec==b.modified.tv_nsec;
    }
};
// zlib inflate directly: validate every gzip member, including trailer/CRC,
// and reject trailing garbage. Plain FASTQ is supported as well.
class Reader {
    std::ifstream file;
    std::array<unsigned char, 1<<16> input{}, output{};
    z_stream z{};
    bool gzip=false, initialized=false, member_end=false, done=false;
    size_t pos=0, available=0;
    std::string path;
    bool refill() {
        if(done) return false;
        if(!gzip) {
            file.read(reinterpret_cast<char*>(output.data()), output.size());
            require(!file.bad(), "Input read error: " + path);
            available=file.gcount(); pos=0; done=available==0;
            return available!=0;
        }
        while(true) {
            if(!z.avail_in) {
                file.read(reinterpret_cast<char*>(input.data()), input.size());
                require(!file.bad(), "Input read error: " + path);
                z.avail_in=file.gcount(); z.next_in=input.data();
            }
            if(member_end) {
                if(!z.avail_in) { done=true; return false; }
                auto next=z.next_in; auto count=z.avail_in;
                require(inflateReset2(&z,31)==Z_OK,"gzip reset failed");
                z.next_in=next; z.avail_in=count; member_end=false;
            }
            require(z.avail_in!=0,"Truncated gzip member: " + path);
            z.next_out=output.data(); z.avail_out=output.size();
            auto before=z.avail_in;
            int rc=inflate(&z,Z_NO_FLUSH);
            require(rc==Z_OK || rc==Z_STREAM_END,
                "gzip inflate/CRC error: " + path + (z.msg ? std::string(": ")+z.msg : ""));
            member_end=rc==Z_STREAM_END; available=output.size()-z.avail_out; pos=0;
            if(available) return true;
            require(member_end || z.avail_in!=before,"gzip made no progress: " + path);
        }
    }
public:
    explicit Reader(const std::string& name):file(name,std::ios::binary),path(name) {
        require(file.good(),"Cannot open input: " + path);
        unsigned char magic[2]{}; file.read(reinterpret_cast<char*>(magic),2);
        gzip=file.gcount()==2 && magic[0]==0x1f && magic[1]==0x8b;
        require(gzip || fs::path(name).extension()!=".gz", "Invalid gzip header: " + name);
        file.clear(); file.seekg(0);
        if(gzip) { require(inflateInit2(&z,31)==Z_OK,"gzip init failed"); initialized=true; }
    }
    ~Reader() { if(initialized) inflateEnd(&z); }
    bool line(std::string& s) {
        s.clear();
        while(true) {
            if(pos==available && !refill()) return !s.empty();
            const auto* begin=output.data()+pos;
            const auto* end=static_cast<const unsigned char*>(std::memchr(begin,'\n',available-pos));
            size_t n=end ? size_t(end-begin)+1 : available-pos;
            s.append(reinterpret_cast<const char*>(begin),n); pos+=n;
            if(s.size()>1<<20) throw std::runtime_error("FASTQ line too long: " + path);
            if(end) {
                s.pop_back(); if(!s.empty() && s.back()=='\r') s.pop_back(); return true;
            }
        }
    }
    bool record(std::string& sequence) {
        std::string header, plus, quality;
        if(!line(header)) return false;
        if(!line(sequence) || !line(plus) || !line(quality))
            throw std::runtime_error("Truncated FASTQ record: " + path);
        if(header.size()<=1 || header[0]!='@' || header[1]==' ' || header[1]=='\t')
            throw std::runtime_error("Invalid FASTQ header: " + path);
        if(plus.empty() || plus[0]!='+') throw std::runtime_error("Invalid FASTQ plus line: " + path);
        if(sequence.empty() || sequence.size()!=quality.size())
            throw std::runtime_error("Invalid FASTQ sequence/quality lengths: " + path);
        for(char c:sequence)
            if(c!='A' && c!='C' && c!='G' && c!='T' && c!='N')
                throw std::runtime_error("Non-ACGTN FASTQ base: " + path);
        for(unsigned char c:quality)
            if(c<33 || c>126) throw std::runtime_error("Invalid FASTQ quality: " + path);
        return true;
    }
};
struct Stats {
    uint64_t reads=0, eligible=0, short_read=0, linker=0, barcode=0;
    bool eof=false;
};
class Extractor {
    std::array<std::vector<int16_t>,2> lookup;
    static int code(const char* s) {
        int value=0;
        for(int i=0;i<10;++i) {
            int b=s[i]=='A'?0:s[i]=='C'?1:s[i]=='G'?2:s[i]=='T'?3:-1;
            if(b<0) return -1;
            value=(value<<2)|b;
        }
        return value;
    }
public:
    Extractor() {
        for(auto& v:lookup) v.assign(1<<20,-1);
        for(size_t i=0;i<WIDTH;++i) {
            lookup[0][whitelist::part0[i]]=i; lookup[1][whitelist::part1[i]]=i;
        }
    }
    int extract(const std::string& s, Stats& st) const {
        ++st.reads;
        if(s.size()<32) { ++st.short_read; return -1; }
        constexpr char anchor[]="TCTGCGCCTTCC";
        int observed=0;
        for(int i=0;i<12;++i) {
            char x=s[i<6?i:i+10];
            if(x!='N') { ++observed; if(x!=anchor[i]) { ++st.linker; return -1; } }
        }
        if(observed<3) { ++st.linker; return -1; }
        int a=code(s.data()+6), b=code(s.data()+22);
        if(a<0 || b<0 || lookup[0][a]<0 || lookup[1][b]<0) { ++st.barcode; return -1; }
        ++st.eligible; return lookup[0][a]*WIDTH+lookup[1][b];
    }
};
std::string stats_json(const Stats& st) {
    return "{\"reads\":"+std::to_string(st.reads)+",\"eligible_reads\":"+std::to_string(st.eligible)+
        ",\"short_reads\":"+std::to_string(st.short_read)+",\"invalid_linker_reads\":"+std::to_string(st.linker)+
        ",\"invalid_barcode_reads\":"+std::to_string(st.barcode)+",\"full_input_validated\":"+(st.eof?"true":"false")+"}";
}
// Exclusive output creation; never overwrite inputs or existing results.
void save_exclusive(const std::string& path, const std::string& data) {
    int fd=open(path.c_str(),O_WRONLY|O_CREAT|O_EXCL,0644);
    require(fd>=0,"Cannot create output (possibly exists): " + path);
    size_t pos=0;
    while(pos<data.size()) {
        ssize_t n=write(fd,data.data()+pos,data.size()-pos);
        if(n<=0) { close(fd); throw std::runtime_error("Output write failed: " + path); }
        pos+=n;
    }
    require(close(fd)==0,"Output close failed: " + path);
}
int main(int argc,char** argv) {
    try {
        const auto started=std::chrono::steady_clock::now();
        const std::string help=
            "oligo_cdna_check 1.0.0\nUsage: oligo_cdna_check --sample NAME --cdna-r1 FILE --oligo-r1 FILE\n"
            "                        --out-prefix PATH [--reads 1M|10M|all]\n"
            "One premerged R1 per library. Full cDNA scan; oligo prefix default 1M.\n"
            "Exact whitelist/linkers; matches / eligible oligo reads; warning below 10%.\n"
            "Outputs PATH.json and PATH.tsv. Exit 0: completed (including low-match warning);\n"
            "1: input/runtime/output error; 2: no usable reference or eligible oligo reads.\n";
        if(argc==2 && std::string(argv[1])=="--help") { std::cout<<help; return 0; }
        if(argc==2 && std::string(argv[1])=="--version") { std::cout<<"1.0.0\n"; return 0; }
        std::map<std::string,std::string> args;
        require(argc%2==1,"Expected --key value pairs.\n"+help);
        for(int i=1;i<argc;i+=2) {
            std::string k=argv[i];
            require(k=="--sample"||k=="--cdna-r1"||k=="--oligo-r1"||k=="--out-prefix"||k=="--reads","Unknown option: "+k);
            require(args.emplace(k,argv[i+1]).second,"Duplicate option: "+k);
        }
        for(const char* k:{"--sample","--cdna-r1","--oligo-r1","--out-prefix"})
            require(args.count(k) && !args.at(k).empty(),std::string("Missing ")+k+"\n"+help);
        std::string sample=args.at("--sample"), cdna=args.at("--cdna-r1"), oligo=args.at("--oligo-r1"), prefix=args.at("--out-prefix");
        require(sample.find_first_of("\t\r\n")==std::string::npos,"Sample name contains TSV control characters");
        std::string mode=args.count("--reads")?args.at("--reads"):"1M";
        require(mode=="1M"||mode=="10M"||mode=="all","--reads must be 1M, 10M or all");
        uint64_t limit=mode=="all"?std::numeric_limits<uint64_t>::max():mode=="1M"?1000000:10000000;
        for(const char* ext:{".json",".tsv"}) require(!fs::exists(prefix+ext) && !fs::is_symlink(prefix+ext),"Refuse existing output: "+prefix+ext);
        fs::path parent=fs::path(prefix).parent_path();
        require(parent.empty() || fs::is_directory(parent),"Output parent directory must already exist");
        Signature c_before(cdna), o_before(oligo);
        Extractor extractor; Stats c, o; uint64_t hits=0;
        std::vector<uint64_t> counts(WIDTH*WIDTH);
        std::string s;
        {
            Reader reader(cdna);
            while(reader.record(s)) {
                int k=extractor.extract(s,c); if(k>=0) ++counts[k];
                if(c.reads%10000000==0) std::cerr<<"cDNA: "<<c.reads<<" reads\n";
            }
            c.eof=true;
        }
        std::vector<uint32_t> order;
        for(size_t i=0;i<counts.size();++i) if(counts[i]) order.push_back(i);
        const size_t unique=order.size(), selected=std::min(TOP,unique);
        auto better=[&](uint32_t a,uint32_t b) { return counts[a]!=counts[b]?counts[a]>counts[b]:a>b; };
        std::partial_sort(order.begin(),order.begin()+selected,order.end(),better);
        std::vector<uint8_t> reference(counts.size());
        for(size_t i=0;i<selected;++i) reference[order[i]]=1;
        {
            Reader reader(oligo);
            while(o.reads<limit) {
                if(!reader.record(s)) { o.eof=true; break; }
                int k=extractor.extract(s,o); if(k>=0 && reference[k]) ++hits;
                if(o.reads%10000000==0) std::cerr<<"oligo: "<<o.reads<<" reads\n";
            }
        }
        require(c_before==Signature(cdna) && o_before==Signature(oligo),"Input changed during scan");
        const bool usable=selected>0 && o.eligible>0;
        const double rate=o.eligible?100.0*hits/o.eligible:0;
        // Integer boundary comparison, with no rounding-dependent decision.
        const bool low=usable && hits <= (o.eligible-1)/10;
        const std::string status=!selected?"NO_CDNA_REFERENCE":!o.eligible?"NO_ELIGIBLE_READS":low?"LOW_MATCH_RATE":"OK";
        const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
        std::ostringstream pct; pct<<std::setprecision(12)<<rate;
        std::ostringstream json;
        json<<"{\n  \"version\":\"1.0.0\",\n  \"metric_version\":\"raw_exact_top50000_v1\",\n"
            <<"  \"sample\":"<<quote(sample)<<",\n  \"status\":"<<quote(status)<<",\n"
            <<"  \"warning_low_match\":"<<(low?"true":"false")<<",\n  \"warning_threshold_percent\":10,\n"
            <<"  \"sampling\":"<<quote(mode)<<",\n  \"sampling_method\":\"oligo_prefix_not_random\",\n"
            <<"  \"cdna_reference_source\":\"full_raw_cDNA_R1_exact_no_correction\",\n"
            <<"  \"cdna_r1\":"<<quote(fs::absolute(cdna).string())<<",\n  \"oligo_r1\":"<<quote(fs::absolute(oligo).string())<<",\n"
            <<"  \"cdna_input_bytes\":"<<c_before.bytes<<",\n  \"oligo_input_bytes\":"<<o_before.bytes<<",\n"
            <<"  \"input_stat_signatures_unchanged\":true,\n  \"whitelist_sha256\":"<<quote(whitelist::sha256)<<",\n"
            <<"  \"cdna\":"<<stats_json(c)<<",\n  \"oligo\":"<<stats_json(o)<<",\n"
            <<"  \"cdna_unique_eligible_barcodes\":"<<unique<<",\n  \"cdna_reference_barcodes\":"<<selected<<",\n"
            <<"  \"top_barcode_tie_rule\":\"count_descending_then_barcode_lexicographic_descending\",\n"
            <<"  \"matching_reads\":"<<hits<<",\n  \"matching_percent\":"<<(usable?pct.str():"null")<<",\n"
            <<"  \"elapsed_seconds\":"<<seconds<<"\n}\n";
        std::ostringstream tsv;
        tsv<<"sample\tstatus\twarning_low_match\tsampling\tcdna_reads\tcdna_eligible_reads\tcdna_reference_barcodes\toligo_reads\toligo_eligible_reads\tmatching_reads\tmatching_percent\n"
           <<sample<<'\t'<<status<<'\t'<<(low?"true":"false")<<'\t'<<mode<<'\t'<<c.reads<<'\t'<<c.eligible<<'\t'<<selected<<'\t'<<o.reads<<'\t'<<o.eligible<<'\t'<<hits<<'\t'<<(usable?pct.str():"NA")<<'\n';
        save_exclusive(prefix+".json",json.str()); save_exclusive(prefix+".tsv",tsv.str());
        std::cout<<sample<<"\t"<<status<<"\t"<<(usable?pct.str():"NA")<<"%\n";
        if(low) std::cerr<<"WARNING LOW_MATCH_RATE: "<<pct.str()<<"% < 10% (matched / eligible oligo reads). Review library pairing; this is not proof of a mismatch.\n";
        if(!usable) std::cerr<<"ERROR "<<status<<": cannot compute an informative matching rate.\n";
        return usable?0:2;
    } catch(const std::exception& e) { std::cerr<<"ERROR: "<<e.what()<<'\n'; return 1; }
}
