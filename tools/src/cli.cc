#include "cli.h"
#include <iostream>
#include <sys/wait.h>
#include <functional>
#include <string>
#include <mpi.h>
#include <cstdlib>
#include <cstring>
#include <vector>

void print_opts(){
    std::cout << "Options: " << std::endl;
    std::cout << "\tcreate_genomicsdb_(w)orkspace" << std::endl;
    std::cout << "\t(g)t_mpi_gather" << std::endl;
    std::cout << "\tvcf(2)genomicsdb" << std::endl;
    std::cout << "\tvcf(d)iff" << std::endl;
    std::cout << "\t(c)onsolidate_genomicsdb_array" << std::endl;
    std::cout << "\tvcf2genomicsdb_(i)nit" << std::endl;
    std::cout << "\tvcf_(h)istogram" << std::endl;
    std::cout << "\t(e)xit" << std::endl;
}

void initialize_MPI(int* procs, int* rank){
    auto rc = MPI_Init(0, 0);
    if (rc != MPI_SUCCESS) {
        std::cerr << "Error starting MPI program. Terminating." << std::endl;
        MPI_Abort(MPI_COMM_WORLD, rc);
    }
    MPI_Comm_size(MPI_COMM_WORLD, procs);
    MPI_Comm_rank(MPI_COMM_WORLD, rank);
}

int main(int argc, char* argv[]){
    int procs, rank;
    initialize_MPI(&procs, &rank);

    if(argc != 1){
        int rv = 0;
        if(!std::strcmp(argv[1], "create_genomicsdb_workspace")){
            if(!rank){
                rv = create_genomicsdb_workspace_main(argc - 1, argv + 1);
            }
            exit(rv);
        }
        if(!std::strcmp(argv[1], "gt_mpi_gather")){
            rv = gt_mpi_gather_main(argc - 1, argv + 1);
            exit(rv);
        }
        if(!std::strcmp(argv[1], "vcf2genomicsdb")){
            rv = vcf2genomicsdb_main(argc - 1, argv + 1);
            exit(rv);
        }
        if(!std::strcmp(argv[1], "vcfdiff")){
            rv = vcfdiff_main(argc - 1, argv + 1);
            exit(rv);
        }
        if(!std::strcmp(argv[1], "consolidate_genomicsdb_array")){
            if(!rank){
                rv = consolidate_genomicsdb_array_main(argc - 1, argv + 1);
            }
            exit(rv);
        }
        if(!std::strcmp(argv[1], "vcf2genomicsdb_init")){
            if(!rank){
                rv = vcf2genomicsdb_init_main(argc - 1, argv + 1);
            }
            exit(rv);
        }
        if(!std::strcmp(argv[1], "vcf_histogram")){
            if(!rank){
                vcf_histogram_main(argc - 1, argv + 1);
            }
            exit(rv);
        }
        
        if(!rank){
        std::cout << "Invalid tool - supported tools: " << std::endl;
            std::cout << "\tcreate_genomicsdb_workspace" << std::endl;
            std::cout << "\tgt_mpi_gather" << std::endl;
            std::cout << "\tvcf2genomicsdb" << std::endl;
            std::cout << "\tvcfdiff" << std::endl;
            std::cout << "\tconsolidate_genomicsdb_array" << std::endl;
            std::cout << "\tvcf2genomicsdb_init" << std::endl;
            std::cout << "\tvcf_histogram" << std::endl;
        }
        exit(1);
    }

    while(1){
        char rec_buf;
        std::string args;
        int arg_len;
        char* args_ptr;

        if(!rank){
            std::string line;

            print_opts();        
            std::cout << "Choice: ";
            std::getline(std::cin, line);

            if(!line.length()){
                std::cout << "Invalid choice" << std::endl;
                continue;
            }

            rec_buf = line[0];

            // usage
            switch(line[0]){
                case 'w': create_genomicsdb_workspace_main(argc, argv); break;
                case 'g': gt_mpi_gather_main(argc, argv); break;
                case '2': vcf2genomicsdb_main(argc, argv); break;
                case 'd': vcfdiff_main(argc, argv); break;
                case 'c': consolidate_genomicsdb_array_main(argc, argv); break;
                case 'i': vcf2genomicsdb_init_main(argc, argv); break;
                case 'h': vcf_histogram_main(argc, argv); break;
                case 'e': exit(0);
                default: std::cout << "Invalid choice" << std::endl; continue;
            }

            // get args
            std::cout << "Enter arguments: ";
            std::getline(std::cin, args);
            arg_len = args.length() + 1;
        }


        // transmit tool
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(&rec_buf, 1, MPI_CHAR, 0, MPI_COMM_WORLD);

        // transmit args (size and str)
        MPI_Bcast(&arg_len, 1, MPI_INT, 0, MPI_COMM_WORLD);
        args_ptr = new char [arg_len];

        if(!rank){
            std::strcpy(args_ptr, args.c_str());
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(args_ptr, arg_len, MPI_CHAR, 0, MPI_COMM_WORLD);

        // split by space
        std::vector<std::string> argv_vector;

        int ind = 0, newind;
        while((size_t)(newind = args.find(" ", ind)) != std::string::npos){
            argv_vector.push_back(args.substr(ind, newind - ind));
            ind = newind + 1;
        }
        argv_vector.push_back(args.substr(ind));

        // construct argv
        auto argv_ptr = new char* [argv_vector.size() + 1];
        auto temp = new char [5];
        std::strcpy(temp, "tool");
        argv_ptr[0] = temp;
        for(int i = 0; i < (int)argv_vector.size(); i++){
            temp = new char [argv_vector[i].length() + 1];
            std::strcpy(temp, argv_vector[i].c_str());
            argv_ptr[i + 1] = temp;
        }

        int argc_int = argv_vector.size() + 1;

        // parallel tools
        int code;
        bool val = 0;
        switch(rec_buf){
            case 'g': val = true; code = gt_mpi_gather_main(argc_int, argv_ptr); break;
            case '2': val = true; code = vcf2genomicsdb_main(argc_int, argv_ptr); break;
            case 'd': val = true; code = vcfdiff_main(argc_int, argv_ptr); break;
            case 'e': MPI_Finalize(); exit(0);
        }

        // serial tools
        if(!rank){
            switch(rec_buf){
                case 'w': code = create_genomicsdb_workspace_main(argc_int, argv_ptr); break;
                case 'c': code = consolidate_genomicsdb_array_main(argc_int, argv_ptr); break;
                case 'i': code = vcf2genomicsdb_init_main(argc_int, argv_ptr); break;
                case 'h': code = vcf_histogram_main(argc_int, argv_ptr); break;
                default: if(!val) {std::cout << "Invalid choice" << std::endl; }
            }
        }

        delete [] args_ptr;

        // clean up argv
        for(int i = 0; i < argc_int; i++){
            delete [] argv_ptr[i];
        }
        delete [] argv_ptr;
    }   
}
