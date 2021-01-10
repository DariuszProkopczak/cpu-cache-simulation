#include "ram.cpp"

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
