#include<iostream>
#include<queue>
#include<cmath>
#include<string.h>
#define DEBUG 0

class Address {

private:

    int tag;
    int index;
    int offset;

public:

    int address;

    Address(int address, int num_cache_sets, int block_size, int bytes_per_value) {

        this->address = address;

        int num_bits_after_tag = (int) (log2(num_cache_sets) + log2(block_size) + log2(bytes_per_value));
        this->tag = address >> num_bits_after_tag;

        int num_bits_after_index = (int) (log2(block_size) + log2(bytes_per_value));
        int temp_address = address >> num_bits_after_index;
        this->index = temp_address % num_cache_sets;

        temp_address = address >> ((int) log2(bytes_per_value));
        int num_bits_offset = (int) log2(block_size);
        int offset_mask = (int) (pow(2, num_bits_offset) - 1);
        this->offset = temp_address & offset_mask;

    }

    int getTag()
    {
        return this->tag;
    }

    int getIndex()
    {
        return this->index;
    }

    int getOffset()
    {
        return this->offset;

    }
};

class DataBlock {

public:

    int size;
    double * data;

    DataBlock(int size)
    {
        this->size = size;
        this->data = new double[this->size];
    }

};

class Ram {

private:

    int num_blocks;
    int num_cache_sets;
    // TODO: free vector
    std::vector<DataBlock> data;

    int get_block_position(Address * address)
    {
        int tag = address->getTag();
        int index = address->getIndex();

        int num_bits_index = (int) log2(this->num_cache_sets);
        return (tag << num_bits_index) | index;
    }

public:

    Ram(int num_blocks, int num_cache_sets, int block_size) {
        this->num_blocks = num_blocks;
        this->num_cache_sets = num_cache_sets;

        for (int i=0; i<this->num_blocks; i++)
            this->data.emplace_back(DataBlock(block_size));
    }

    ~Ram()
    {
        for (int i=0; i<this->num_blocks; i++)
            free(this->data[i].data);
    }

    int getRamSize()
    {
        return this->num_blocks;
    }

    DataBlock getBlock(Address * address)
    {
        return this->data[this->get_block_position(address)];
    }

    //TODO: remove
//    void setBlock(Address address, DataBlock dataBlock)
//    {
//        this->data[this->get_block_position(address)] = dataBlock;
//    }
};

class Cache {

private:

    // TODO: free arrays
    std::vector<std::vector<DataBlock>> blocks;
    std::vector<std::vector<int>> last_used;   // for least recently used replacement strategy. lower value signifies
                                               //     more recently used block
    std::vector<int> blocks_to_replace;        // used by the "first in first out" block replacement strategy
                                               //     to keep track of which block in the set was placed last
    std::vector<std::vector<int>> tags;
    std::vector<std::vector<int>> valid_bits;
    int block_replacement_strategy;
    int set_size;
    int num_blocks;
    int read_hits;         // user can reset to 0 with "reset_cache_logging()" and get value by calling "get_read_hits()"
    int read_misses;       // user can reset to 0 with "reset_cache_logging()" and get value by calling "get_read_misses()"
    int write_hits;        // user can reset to 0 with "reset_cache_logging()" and get value by calling "get_write_hits()"
    int write_misses;      // user can reset to 0 with "reset_cache_logging()" and get value by calling "get_write_misses()"
    Ram * ram; // reference to ram object for getting and setting blocks from ram

    // initializes vector used to keep track of which data block in a set was last used
    //    the vector is only initialized if the "least recently used" block replacement strategy
    //    is chosen
    void initialize_last_used_vector()
    {
        int num_sets = this->num_blocks/this->set_size;
        for (int set_idx = 0; set_idx<num_sets; set_idx++)
        {
            std::vector<int> row;
            for (int block_idx = 0; block_idx<this->set_size; block_idx++)
            {
                row.emplace_back(this->set_size - 1); // set to max value initially
            }
            this->last_used.emplace_back(row);
        }
    }

    // used by "first in, first out" block replacement strategy to keep track block order
    void initialize_last_placed_block_vector()
    {
        int num_sets = this->num_blocks/this->set_size;
        for (int set_idx = 0; set_idx<num_sets; set_idx++)
            this->blocks_to_replace.push_back(0);
    }

    // returns the block to use by the fifo replacement strategy, and updates block order in set
    int get_and_update_fifo_block(int index)
    {
        int block_to_replace = this->blocks_to_replace[index];
        this->blocks_to_replace[index] = (block_to_replace + 1) % this->set_size;
        return block_to_replace;
    }

    // used by "least recently used" block replacement strategy to update block order
    void update_last_used_block_order(int index, int block)
    {
        int curr_block_precedence = this->last_used[index][block];
        for (int i=0; i<this->set_size; i++)
        {
            // decrease the precedence of other blocks that originally had a higher precedence than the selected block
            if ((this->last_used[index][i] < curr_block_precedence) && (i != block))
            {
                this->last_used[index][i] = this->last_used[index][i] + 1;
            }
        }
        this->last_used[index][block] = 0; // the block now has the highest precedence
    }

    // used by "least recently used" block replacement strategy
    int get_least_recently_used_block(int set_idx)
    {
        for (int i=0; i<this->set_size; i++)
        {
            if (this->last_used[set_idx][i] == (this->set_size-1))
                return i;
        }
        return -1;  // something went wrong
    }

    DataBlock getBlock(Address * address)
    {

        int index = address->getIndex();
        int tag = address->getTag();

        // cycle through each block of the set and check valid bit and tag
        for (int block_idx=0; block_idx<this->set_size; block_idx++)
        {
            if ((this->valid_bits[index][block_idx] != 0) & (tag == this->tags[index][block_idx]))
            {
                if (this->block_replacement_strategy == 2)
                    update_last_used_block_order(index, block_idx);

                this->read_hits += 1;

                return this->blocks[index][block_idx];
            }
        }

        // if this point was reached, block was not in cache. Will get block from RAM and place in cache
        this->read_misses += 1;
        DataBlock block = this->ram->getBlock(address);
        this->setBlock(address, block);

        return block;
    }

    void setBlock(Address * address, DataBlock block)
    {
        int block_idx;

        int tag = address->getTag();
        int set_idx = address->getIndex();
        if (this->block_replacement_strategy == 0)    // random
        {
            block_idx = rand() % this->set_size;
        }
        else if (this->block_replacement_strategy == 1)// first in, first out
        {
            block_idx = this->get_and_update_fifo_block(set_idx);
        }
        else      // least recently used
        {
            block_idx = this->get_least_recently_used_block(set_idx);
            this->update_last_used_block_order(set_idx, block_idx);
        }
        this->tags[set_idx][block_idx] = tag;
        this->valid_bits[set_idx][block_idx] = 1;
        this->blocks[set_idx][block_idx] = block;
    }

public:

    Cache(Ram * ram, int set_size, int num_blocks, int block_size, int block_placement_strategy, int block_replacement_strategy)
    {
        /*
         * block_placement_strategy: 0: direct-mapped, 1: fully associate, 2: n-way set associative
         * block_replacement_strategy: 0: random, 1: first in, first out, 2: least recently used
         */

        this->ram = ram; // set reference to ram object for getting and setting blocks from ram

        this->read_hits = 0;
        this->read_misses = 0;
        this->write_hits = 0;
        this->write_misses = 0;

        if (block_placement_strategy == 0)
            this->set_size = 1; // direct-mapped
        else if (block_placement_strategy == 1)
            this->set_size = num_blocks; // fully associative
        else
            this->set_size = set_size; // n-way set associative

        this->num_blocks = num_blocks;
        int num_sets = num_blocks / set_size;

        if (block_replacement_strategy == 0)
            this->block_replacement_strategy = 0; // random
        else if (block_replacement_strategy == 1)
        {
            this->block_replacement_strategy = 1; // First in, first out
            initialize_last_placed_block_vector();
        }
        else
        {
            this->block_replacement_strategy = 2; // Least recently used. Will use array to keep track of least recently
            //     used block
            initialize_last_used_vector();
        }

        // initialize cache with data blocks
        for (int set_idx=0; set_idx<num_sets; set_idx++)
        {
            std::vector<DataBlock> blocks_row;
            std::vector<int> tags_row;
            std::vector<int> valid_bits_row;
            for (int block_idx=0; block_idx<set_size; block_idx++)
            {
                blocks_row.emplace_back(DataBlock(block_size));
                tags_row.emplace_back(0);
                valid_bits_row.emplace_back(0);
            }
            this->blocks.emplace_back(blocks_row);
            this->tags.emplace_back(tags_row);
            this->valid_bits.emplace_back(valid_bits_row);
        }
    }

    void reset_cache_logging()
    {
        this->read_hits = 0;
        this->read_misses = 0;
        this->write_hits = 0;
        this->write_misses = 0;
    }

    int get_read_hits()
    {
        return this->read_hits;
    }

    int get_read_misses()
    {
        return this->read_misses;
    }

    int get_write_hits()
    {
        return this->write_hits;
    }

    int get_write_misses()
    {
        return this->write_misses;
    }

    double getDouble(Address * address)
    {
        DataBlock block = this->getBlock(address);
        return block.data[address->getOffset()];
    }

    void setDouble(Address * address, double value)
    {
        int index = address->getIndex();
        int tag = address->getTag();
        int offset = address->getOffset();

        DataBlock block = this->ram->getBlock(address);
        block.data[offset] = value;     // updates value in ram

        // cycle through each block of the set and check valid bit and tag
        for (int block_idx=0; block_idx<this->set_size; block_idx++)
        {
            if ((this->valid_bits[index][block_idx] != 0) & (tag == this->tags[index][block_idx]))
            {
                if (this->block_replacement_strategy == 2)
                    update_last_used_block_order(index, block_idx);

                this->write_hits += 1;

                this->blocks[index][block_idx].data[offset] = value;
                return;
            }
        }

        // if this point is reached, block is not in cache, set using block taken from ram
        this->write_misses += 1;
        this->setBlock(address, block);

    }

};

class Cpu {

private:

    Cache * cache;
    int num_cache_sets;
    int block_size_doubles;
    int double_size_bytes;
    int instruction_count;

public:

    Cpu(Cache * cache, int num_cache_sets, int block_size_doubles, int double_size_bytes)
    {
        this->cache = cache;
        this->num_cache_sets = num_cache_sets;
        this->block_size_doubles = block_size_doubles;
        this->double_size_bytes = double_size_bytes;
        this->instruction_count = 0;
    }

    void reset_logging()
    {
        this->cache->reset_cache_logging();
        this->instruction_count = 0;
    }

    int get_read_hits()
    {
        return this->cache->get_read_hits();
    }

    int get_read_misses()
    {
        return this->cache->get_read_misses();
    }

    int get_write_hits()
    {
        return this->cache->get_write_hits();
    }

    int get_write_misses()
    {
        return this->cache->get_write_misses();
    }

    int get_instruction_count()
    {
        return this->instruction_count;
    }

    double loadDouble(Address * address)
    {
        this->instruction_count += 1;
        return this->cache->getDouble(address);
    }

    void storeDouble(Address * address, double value)
    {
        this->instruction_count += 1;
        this->cache->setDouble(address, value);
    }

    double addDouble(double value1, double value2)
    {
        this->instruction_count += 1;
        return value1 + value2;
    }

    double multDouble(double value1, double value2)
    {
        this->instruction_count += 1;
        return value1 * value2;
    }
};

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