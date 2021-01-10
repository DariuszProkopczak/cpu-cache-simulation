#include "address.cpp"
#include "datablock.cpp"
#include <vector>

class Ram {

private:

    int num_blocks;
    int num_cache_sets;
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
