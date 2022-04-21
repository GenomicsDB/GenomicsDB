#include "tiledb_generic_loader.h"

void read_sam_file(std::string filename) {
  std::cerr << "SAM file is " << filename << std::endl;

  samFile *fp_in = hts_open(filename.c_str(),"r"); //open bam file
  bam_hdr_t *bamHdr = sam_hdr_read(fp_in); //read header
  bam1_t *aln = bam_init1(); //initialize an alignment

  if(!bamHdr) {
    std::cerr << "header is null" << std::endl;
  } else {
    std::cerr << "header is NOT null" << std::endl;
  }
  
  // header parse
  // uint32_t *tar = bamHdr->text ;
  // uint32_t *tarlen = bamHdr->target_len ;

  // printf("%d\n",tar);
  
  int rc;
  std::cerr << "before while" << std::endl;
  while(!(rc = sam_read1(fp_in,bamHdr,aln))){
          
    int32_t pos = aln->core.pos +1; //left most position of alignment in zero based coordinate (+1)
    char *chr = bamHdr->target_name[aln->core.tid] ; //contig name (chromosome)
    uint32_t len = aln->core.l_qseq; //length of the read.
    
    uint8_t *q = bam_get_seq(aln); //quality string
    uint32_t q2 = aln->core.qual ; //mapping quality
    
    char* qname = bam_get_qname(aln);    
    uint16_t flag = aln->core.flag;
    uint32_t* cigar = bam_get_cigar(aln);
    uint32_t n_cigar = aln->core.n_cigar;
    char cigar_codes[] = {'M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X'};
    //uint8_t* qual = bam_get_qual(aln);
    char* qual = (char*)bam_get_qual(aln);
    uint8_t mapq = aln->core.qual;
    //char* seq = (char*)bam_get_seq(aln);
    int32_t rnext = aln->core.mtid;
    int32_t pnext = aln->core.mpos;
    int32_t tlen = aln->core.isize;

    char *qseq = (char *)malloc(len);
    
    for(int i=0; i< len ; i++){
      qseq[i] = seq_nt16_str[bam_seqi(q,i)]; //gets nucleotide id and converts them into IUPAC id.
    }
    
    //printf("chr=%s\tpos=%d\tlen=%d\tqseq=%s\tq=%s\tq2=%d\n",chr,pos,len,qseq,q,q2);
    printf("qname=%s\tflag=%d\tchr=%s\tpos=%d\tmapq=%d\tlen=%d\tqseq=%s\tq2=%d\n",qname,flag,chr,pos,mapq,len,qseq,q2);
    std::cout << "cigar=";
    for(uint32_t i = 0; i < n_cigar; i++) {
      auto op_len = cigar[i] >> 4;
      auto op = cigar_codes[cigar[i] & 0b1111];

      std::cout << op_len << op << ", ";
    }
    std::cout << std::endl;
    std::cout << "rnext=" << rnext << "\tpnext=" << pnext << "\ttlen=" << tlen << std::endl;
    std::cerr << "after while rc is " << rc << std::endl;

    std::cout << "qual=" << std::endl;
    for(uint32_t i = 0; i < len; i++) {
      double q = qual[i];
      std::cout << pow(10, (q/-10)) << ", ";
    }
    std::cout << std::endl << std::endl;

    free(qseq);
  }
  bam_destroy1(aln);
  sam_close(fp_in);
}

std::vector<uint8_t> OmicsMultiCell::as_cell(const VariantArraySchema& schema) {
  return {};
  //construct cell
  /*std::vector<uint8_t> cell(16);
  *(reinterpret_cast<uint64_t*>(cell.data())) = coords[0]; // write row in cell
  *(reinterpret_cast<uint64_t*>(cell.data()) + 1) = coords[1]; // write position in cell

  // reserve space for cell size
  for(int i = 0; i < sizeof(size_t); i++) {
    cell.push_back(0);
  }

  // attributes
  for(int i = 0; i < schema.attribute_num(); i++) {
    std::string attribute_name = schema.attribute_name(i);

    for(auto& sc : subcells) {
      
      cell.push_back();

      if(attribute_name == "NAME") {
        cell.insert(cell.end(), {0, 0, 0, 0});
        *(reinterpret_cast<uint32_t*>(cell.data() + cell.size() - 4)) = name.length();
        for(auto c : name) {
          cell.push_back(c);
        }
      }
    }
  }
  // fill in cell size
  *(reinterpret_cast<size_t*>(cell.data() + 2*sizeof(int64_t))) = cell.size();*/
}

GenericTileDBLoader::GenericTileDBLoader(
                                         const std::string& config_filename,
                                         const int idx,
                                         const bool superimpose
                                        ): m_config_filename(config_filename), m_idx(idx), m_superimpose(superimpose) {}
