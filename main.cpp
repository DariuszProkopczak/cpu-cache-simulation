#include <iostream>
#include <queue>
#include <string>
#include "cpu.cpp"
#define DEBUG 0


// create ram object from user arguments
Ram Initialize_Ram(int data_block_size_bytes, int matrix_dim, int double_size_bytes, int cache_size_bytes, int n_way_associativity,  int algorithm)
{
    int block_size_doubles = (data_block_size_bytes / double_size_bytes);
    int cache_size_blocks = cache_size_bytes / data_block_size_bytes;
    int num_cache_sets = cache_size_blocks / n_way_associativity;
    int blocks_in_ram;

    if (algorithm == 0) // daxpy: size of memory needed is 3 * <dimension size>
    {
        blocks_in_ram = (int) (3 * ceil(((double) matrix_dim) / block_size_doubles));
        return Ram(blocks_in_ram, num_cache_sets, block_size_doubles);
    }
    else // matrix multiply: size of memory needed is 3 * <dimension size>^2
    {
        blocks_in_ram = (int) (3 * ceil( ((double) (matrix_dim * matrix_dim)) / block_size_doubles));
        return Ram(blocks_in_ram, num_cache_sets, block_size_doubles);
    }
}

// create ram object from user arguments
Cache Initialize_Cache(Ram * ram, int data_block_size_bytes, int double_size_bytes, int cache_size_bytes, int n_way_associativity, int block_replacement_strategy)
{
    int block_size_doubles = (data_block_size_bytes / double_size_bytes);
    int cache_size_blocks = cache_size_bytes / data_block_size_bytes;

    int block_placement_strategy;
    if (n_way_associativity == 1)
        block_placement_strategy = 0;
    else if (n_way_associativity == cache_size_blocks)
        block_placement_strategy = 1;
    else
        block_placement_strategy = 2;

    return Cache(ram, n_way_associativity, cache_size_blocks, block_size_doubles, block_placement_strategy, block_replacement_strategy);

}

int main(int argc, char* argv[]) {

    int cache_size_bytes = 65536;
    int data_block_size_bytes = 64;
    int n_way_associativity = 2;
    int replacement_policy = 2; // 0: random, 1: FIFO, 2: LRU
    int algorithm = 2; // 0: daxpy, 1: mxm, 2: mxm_block
    int matrix_dim = 480;
    int blocking_factor = 32;
    int print_solution_matrix = 0;

    int double_size_bytes = 8; // hard-code size of double to 8 bytes

    // parse arguments
    for (int i = 1; i < argc; i++) {
        if (argv[i][1] == 'c') {
            cache_size_bytes = atoi(argv[i + 1]);
        } else if (argv[i][1] == 'b') {
            data_block_size_bytes = atoi(argv[i + 1]);
        } else if (argv[i][1] == 'n') {
            n_way_associativity = atoi(argv[i + 1]);
        } else if (argv[i][1] == 'r') {
            if (strcmp(argv[i + 1], "random") == 0) {
                replacement_policy = 0;
            } else if (strcmp(argv[i + 1], "FIFO") == 0) {
                replacement_policy = 1;
            } else if (strcmp(argv[i + 1], "LRU") == 0) {
                replacement_policy = 2;
            }
        } else if (argv[i][1] == 'a') {
            if (strcmp(argv[i + 1], "daxpy") == 0) {
                algorithm = 0;
            } else if (strcmp(argv[i + 1], "mxm") == 0) {
                algorithm = 1;
            } else if (strcmp(argv[i + 1], "mxm_block") == 0) {
                algorithm = 2;
            }
        } else if (argv[i][1] == 'd') {
            matrix_dim = atoi(argv[i + 1]);
        } else if (argv[i][1] == 'f') {
            blocking_factor = atoi(argv[i + 1]);
        } else if (argv[i][1] == 'p') {
            print_solution_matrix = 1;
            continue; // continue to skip i++ below
        }

        i++; // skip an argument
    }

    Ram ram = Initialize_Ram(data_block_size_bytes, matrix_dim, double_size_bytes, cache_size_bytes,
                             n_way_associativity, algorithm);
    int ram_size_blocks = ram.getRamSize();
    Cache cache = Initialize_Cache(&ram, data_block_size_bytes, double_size_bytes, cache_size_bytes,
                                   n_way_associativity, replacement_policy);
    int num_cache_sets = (cache_size_bytes / data_block_size_bytes) / n_way_associativity;
    int block_size_doubles = data_block_size_bytes / double_size_bytes;

    Cpu cpu = Cpu(&cache, num_cache_sets, block_size_doubles, double_size_bytes);

    int read_hits;
    int read_misses;
    int write_hits;
    int write_misses;
    int instruction_count;

    // execute daxpy
    if (algorithm == 0) {

        std::vector<Address *> a;
        std::vector<Address *> b;
        std::vector<Address *> c;

        double register0 = 3;
        double register1;
        double register2;
        Address * temp_address;

         // addresses
        for (int i = 0; i < matrix_dim; i++) {
            temp_address = new Address(i * double_size_bytes, num_cache_sets, data_block_size_bytes / double_size_bytes, double_size_bytes);
            a.emplace_back(temp_address);
            temp_address = new Address( (matrix_dim * double_size_bytes) + (i * double_size_bytes), num_cache_sets, data_block_size_bytes / double_size_bytes, double_size_bytes);
            b.emplace_back(temp_address);
            temp_address = new Address((2 * matrix_dim * double_size_bytes) + (i * double_size_bytes), num_cache_sets, data_block_size_bytes / double_size_bytes, double_size_bytes);
            c.emplace_back(temp_address);
        }

        for (int i = 0; i < matrix_dim; i++) {
            cpu.storeDouble(a[i], (double) i);
            cpu.storeDouble(b[i], (double) (2 * i));
            cpu.storeDouble(c[i], (double) 0);
        }

        #if DEBUG
        for (int i = 0; i < matrix_dim; i++) {
            double test1 = cpu.loadDouble(a[i]);
            double test2 = cpu.loadDouble(b[i]);
            printf("Val %d: %f\n", i, test1);
            printf("Val %d: %f\n", i, test2);
        }
        #endif

        cpu.reset_logging();

        for (int i = 0; i < matrix_dim; i++) {
            #if DEBUG
            if ((i%50) == 0)
                printf("i: %d\n", i);
            #endif
            register1 = cpu.loadDouble(a[i]);
            register1 = cpu.multDouble(register1, register0);
            register2 = cpu.loadDouble(b[i]);
            register2 = cpu.multDouble(register1, register2);
            cpu.storeDouble(c[i], register2);
        }

        read_hits = cpu.get_read_hits();
        read_misses = cpu.get_read_misses();
        write_hits = cpu.get_write_hits();
        write_misses = cpu.get_write_misses();
        instruction_count = cpu.get_instruction_count();

        if (print_solution_matrix) {
            printf("RESULT ARRAY:\n");
            for (int i = 0; i < matrix_dim; i++) {
                register1 = cpu.loadDouble(c[i]);
                printf("%f\n", register1);
            }
            printf("\n");
        }

    }
    else if (algorithm == 1) // mxm
    {

        std::vector<Address *> a;
        std::vector<Address *> b;
        std::vector<Address *> c;

        double register0;
        double register1;
        double register2;
        double register3;
        Address * temp_address;
        // addresses
        for (int i = 0; i < (matrix_dim * matrix_dim); i++) {
            temp_address = new Address(i * double_size_bytes, num_cache_sets, data_block_size_bytes / double_size_bytes, double_size_bytes);
            a.emplace_back(temp_address);
            temp_address = new Address(((matrix_dim * matrix_dim) * double_size_bytes) + (i * double_size_bytes), num_cache_sets, data_block_size_bytes / double_size_bytes, double_size_bytes);
            b.emplace_back(temp_address);
            temp_address = new Address((2 * (matrix_dim * matrix_dim) * double_size_bytes) + (i * double_size_bytes), num_cache_sets, data_block_size_bytes / double_size_bytes, double_size_bytes);
            c.emplace_back(temp_address);
        }

        for (int i = 0; i < (matrix_dim * matrix_dim); i++)
        {
                cpu.storeDouble(a[i], (double) i);
                cpu.storeDouble(b[i], (double) (2 * i));
                cpu.storeDouble(c[i], (double) 0);
        }

        #if DEBUG
        for (int i = 0; i < (matrix_dim * matrix_dim); i++) {
            double test1 = cpu.loadDouble(a[i]);
            double test2 = cpu.loadDouble(b[i]);
            printf("Val %d: %f\n", i, test1);
            printf("Val %d: %f\n", i, test2);
        }
        #endif

        cpu.reset_logging();

        for (int i = 0; i < matrix_dim; i++) {
            #if DEBUG
            printf("i: %d\n", i);
            #endif
            for (int j = 0; j < matrix_dim; j++) {
                register0 = 0; // sum
                for (int k = 0; k < matrix_dim; k++) {
                    register1 = cpu.loadDouble(a[(matrix_dim * i) + k]);
                    register2 = cpu.loadDouble(b[(matrix_dim * k) + j]);
                    register3 = cpu.multDouble(register1, register2);
                    register0 = cpu.addDouble(register0, register3);
                }
                cpu.storeDouble(c[(matrix_dim * i) + j], register0);
            }
        }

        read_hits = cpu.get_read_hits();
        read_misses = cpu.get_read_misses();
        write_hits = cpu.get_write_hits();
        write_misses = cpu.get_write_misses();
        instruction_count = cpu.get_instruction_count();

        if (print_solution_matrix) {
            printf("RESULT MATRIX:\n");
            for (int i = 0; i < matrix_dim; i++) {
                for (int j = 0; j < matrix_dim; j++) {
                    register1 = cpu.loadDouble(c[(matrix_dim * i) + j]);
                    printf("%f ", register1);
                }
                printf("\n");
            }
            printf("\n");
        }

    }
    else // mxm block
    {

        std::vector<Address *> a;
        std::vector<Address *> b;
        std::vector<Address *> c;

        double register0;
        double register1;
        double register2;
        double register3;
        Address * temp_address;
        // addresses
        for (int i = 0; i < (matrix_dim * matrix_dim); i++) {
            temp_address = new Address(i * double_size_bytes, num_cache_sets, data_block_size_bytes / double_size_bytes, double_size_bytes);
            a.emplace_back(temp_address);
            temp_address = new Address(((matrix_dim * matrix_dim) * double_size_bytes) + (i * double_size_bytes), num_cache_sets, data_block_size_bytes / double_size_bytes, double_size_bytes);
            b.emplace_back(temp_address);
            temp_address = new Address((2 * (matrix_dim * matrix_dim) * double_size_bytes) + (i * double_size_bytes), num_cache_sets, data_block_size_bytes / double_size_bytes, double_size_bytes);
            c.emplace_back(temp_address);
        }

        for (int i = 0; i < (matrix_dim * matrix_dim); i++)
        {
            cpu.storeDouble(a[i], (double) i);
            cpu.storeDouble(b[i], (double) (2 * i));
            cpu.storeDouble(c[i], (double) 0);
        }

        #if DEBUG
        for (int i = 0; i < (matrix_dim * matrix_dim); i++) {
            double test1 = cpu.loadDouble(a[i]);
            double test2 = cpu.loadDouble(b[i]);
            printf("Val %d: %f\n", i, test1);
            printf("Val %d: %f\n", i, test2);
        }
        #endif

        cpu.reset_logging();

        for (int sub_i=0; sub_i<(matrix_dim/blocking_factor); sub_i++) {
            for (int sub_j=0; sub_j<(matrix_dim/blocking_factor); sub_j++) {
                for (int i = (sub_i*blocking_factor); i < ((sub_i*blocking_factor) + blocking_factor); i++) {
                    for (int j = (sub_j*blocking_factor); j < ((sub_j*blocking_factor) + blocking_factor); j++) {
                        register0 = 0; // sum
                        for (int k = 0; k < matrix_dim; k++) {
                            #if DEBUG
                            printf("i: %d\n", i);
                            #endif
                            register1 = cpu.loadDouble(a[(matrix_dim * i) + k]);
                            register2 = cpu.loadDouble(b[(matrix_dim * k) + j]);
                            register3 = cpu.multDouble(register1, register2);
                            register0 = cpu.addDouble(register0, register3);
                        }
                        cpu.storeDouble(c[(matrix_dim * i) + j], register0);
                    }
                }
            }
        }

        read_hits = cpu.get_read_hits();
        read_misses = cpu.get_read_misses();
        write_hits = cpu.get_write_hits();
        write_misses = cpu.get_write_misses();
        instruction_count = cpu.get_instruction_count();

        if (print_solution_matrix)
        {
            printf("RESULT MATRIX:\n");
            for (int i=0; i<matrix_dim; i++)
            {
                for (int j=0; j<matrix_dim; j++)
                {
                    register1 = cpu.loadDouble(c[(matrix_dim * i) + j]);
                    printf("%f ", register1);
                }
                printf("\n");
            }
            printf("\n");
        }
    }

    printf("INPUTS=====================================\n");
    printf("Ram Size =                   %d bytes\n", ram_size_blocks * data_block_size_bytes);
    printf("Cache Size =                 %d bytes\n", cache_size_bytes);
    printf("Block Size =                 %d bytes\n", data_block_size_bytes);
    printf("Total Blocks in Cache =      %d\n", cache_size_bytes/data_block_size_bytes);
    printf("Associativity =              %d\n", n_way_associativity);
    printf("Number of Sets =             %d\n", (cache_size_bytes/data_block_size_bytes) / n_way_associativity);

    if (replacement_policy == 0)
        printf("Replacement Policy =         random\n");
    else if (replacement_policy == 1)
        printf("Replacement Policy =         FIFO\n");
    else
        printf("Replacement Policy =         LRU\n");

    if (algorithm == 0)
        printf("Algorithm =                  daxpy\n");
    else if (algorithm == 1)
        printf("Algorithm =                  mxm\n");
    else
        printf("Algorithm =                  blocked mxm\n");

    printf("MXM Blocking Factor =        %d\n", blocking_factor);
    printf("Matrix or Vector Dimension = %d\n", matrix_dim);

    printf("RESULTS=====================================\n");
    printf("Instruction Count =          %d\n", instruction_count);
    printf("Read Hits =                  %d\n", read_hits);
    printf("Read Misses =                %d\n", read_misses);
    printf("Read Miss Rate =             %f %%\n", (read_misses / ((double) read_hits + read_misses)) * 100);
    // TODO: Getting strange results with writes. May need to debug.
    // printf("Write Hits =                 %d\n", write_hits);
    // printf("Write Misses =               %d\n", write_misses);
    // printf("Write Miss Rate =            %f %%\n", (write_misses / ((double) write_hits + write_misses)) * 100);

    return 0;
}